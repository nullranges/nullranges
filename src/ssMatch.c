#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

/* ssMatch's remaining focal and pool, kept in one list that only ssMatch
 * holds and that ss_pass updates in place. */
enum {
	S_F, S_FIDX, S_P, S_PIDX, S_CNT, S_CODE, S_LENS, S_OUTF, S_OUTP,
	S_FBIN, S_PBIN, S_FN, S_PN, S_RANK, S_BINOF, S_POFF, S_FOFF, S_PGRP,
	S_FGRP, S_BR, S_N
};
/* S_LENS: remaining focal, remaining pool, nLeft, matched so far */

SEXP ss_init(SEXP fps, SEXP pps, SEXP cnt, SEXP code, SEXP nLeft)
{
	if (!isReal(fps) || !isReal(pps) || !isInteger(cnt) ||
	    !isInteger(code) || !isInteger(nLeft) || XLENGTH(nLeft) != 1)
		error("ss_init: invalid arguments");
	R_xlen_t nf = XLENGTH(fps), np = XLENGTH(pps);
	if (nf + np > INT_MAX / 2 || XLENGTH(code) != nf + np ||
	    XLENGTH(cnt) != nf + np)
		error("ss_init: invalid lengths");
	R_xlen_t cap = nf + np + 2, capr = nf + 2;

	SEXP st = PROTECT(allocVector(VECSXP, S_N));
	SET_VECTOR_ELT(st, S_F, allocVector(REALSXP, nf));
	SET_VECTOR_ELT(st, S_FIDX, allocVector(INTSXP, nf));
	SET_VECTOR_ELT(st, S_P, allocVector(REALSXP, np));
	SET_VECTOR_ELT(st, S_PIDX, allocVector(INTSXP, np));
	SET_VECTOR_ELT(st, S_CNT, allocVector(INTSXP, nf + np));
	SET_VECTOR_ELT(st, S_CODE, code);
	SET_VECTOR_ELT(st, S_LENS, allocVector(INTSXP, 4));
	SET_VECTOR_ELT(st, S_OUTF, allocVector(INTSXP, nf));
	SET_VECTOR_ELT(st, S_OUTP, allocVector(INTSXP, nf));
	SET_VECTOR_ELT(st, S_FBIN, allocVector(INTSXP, nf));
	SET_VECTOR_ELT(st, S_PBIN, allocVector(INTSXP, np));
	SET_VECTOR_ELT(st, S_FN, allocVector(INTSXP, cap));
	SET_VECTOR_ELT(st, S_PN, allocVector(INTSXP, cap));
	SET_VECTOR_ELT(st, S_RANK, allocVector(INTSXP, cap));
	SET_VECTOR_ELT(st, S_BINOF, allocVector(INTSXP, capr));
	SET_VECTOR_ELT(st, S_POFF, allocVector(INTSXP, capr));
	SET_VECTOR_ELT(st, S_FOFF, allocVector(INTSXP, capr));
	SET_VECTOR_ELT(st, S_PGRP, allocVector(INTSXP, np));
	SET_VECTOR_ELT(st, S_FGRP, allocVector(INTSXP, nf));
	SET_VECTOR_ELT(st, S_BR, allocVector(REALSXP, cap));

	if (nf)
		memcpy(REAL(VECTOR_ELT(st, S_F)), REAL(fps), nf * sizeof(double));
	if (np)
		memcpy(REAL(VECTOR_ELT(st, S_P)), REAL(pps), np * sizeof(double));
	if (nf + np)
		memcpy(INTEGER(VECTOR_ELT(st, S_CNT)), INTEGER(cnt),
			(nf + np) * sizeof(int));
	int *fi = INTEGER(VECTOR_ELT(st, S_FIDX));
	for (R_xlen_t j = 0; j < nf; j++)
		fi[j] = (int) j + 1;
	int *pi = INTEGER(VECTOR_ELT(st, S_PIDX));
	for (R_xlen_t j = 0; j < np; j++)
		pi[j] = (int) j + 1;
	int *lens = INTEGER(VECTOR_ELT(st, S_LENS));
	lens[0] = (int) nf;
	lens[1] = (int) np;
	lens[2] = INTEGER(nLeft)[0];
	lens[3] = 0;
	UNPROTECT(1);
	return st;
}

/* min and max of c(f, p) as R's min() and max() give them for doubles
 * (the first of equal values, so the first of 0 and -0); NaN when
 * either holds NA or NaN, which seq() then rejects either way. */
SEXP ss_range(SEXP st)
{
	const int *lens = INTEGER(VECTOR_ELT(st, S_LENS));
	const double *f = REAL(VECTOR_ELT(st, S_F)), *p = REAL(VECTOR_ELT(st, S_P));
	int lf = lens[0], lp = lens[1], any = 0;
	double mn = R_PosInf, mx = R_NegInf;
	for (int s = 0; s < 2; s++) {
		const double *x = s ? p : f;
		int n = s ? lp : lf;
		for (int j = 0; j < n; j++) {
			double v = x[j];
			if (ISNAN(v)) {
				mn = mx = R_NaN;
				goto done;
			}
			if (!any) {
				mn = mx = v;
				any = 1;
			} else {
				if (v < mn)
					mn = v;
				if (v > mx)
					mx = v;
			}
		}
	}
done:;
	SEXP out = PROTECT(allocVector(REALSXP, 2));
	REAL(out)[0] = mn;
	REAL(out)[1] = mx;
	UNPROTECT(1);
	return out;
}

/* .bincode(x, br, right = TRUE, include.lowest = TRUE) for sorted breaks
 * without NaN, NA given as nb: a guess from the spacing, checked, else
 * R's binary search. */
static inline int bin_of(double x, const double *br, int nb, double mn,
	double by)
{
	if (ISNAN(x) || x < br[0] || br[nb - 1] < x)
		return nb;
	double t = ceil((x - mn) / by);
	if (!(t >= 1))
		t = 1;
	if (t > nb - 1)
		t = nb - 1;
	int g = (int) t;
	if (br[g - 1] < x && x <= br[g])
		return g;
	int lo = 0, hi = nb - 1;
	while (hi - lo >= 2) {
		int mid = (hi + lo) / 2;
		if (x > br[mid])
			lo = mid;
		else
			hi = mid;
	}
	return lo + 1;
}

static inline void count_out(int *cnt, const int *code, int row, int *nLeft)
{
	int c = code[row] - 1;
	if (--cnt[c] == 0)
		(*nLeft)--;
}

static SEXP run_pass(SEXP st, const double *br, int nb, double mn,
	double by, SEXP replace_, SEXP count_);

/* One pass of ssMatch after seq() made the breaks: bin the remaining
 * focal and pool as stratify does, draw each qualifying bin's pool rows
 * (bins in order of first appearance in the pool) exactly as
 * sample(length(x), size, replace) does, record the pairs, count the
 * matched rows out of nLeft when `count`, and drop them from the
 * remaining focal and pool. Returns c(matched, nLeft): matched is 0 when
 * no bin qualifies and NA when the breaks fail binCode's guard; neither
 * changes the state or draws a number. */
SEXP ss_pass(SEXP st, SEXP br_, SEXP mn_, SEXP by_, SEXP replace_,
	SEXP count_)
{
	/* binCode's guard, after which .bincode runs on the whole vector */
	double mn = asReal(mn_), by = asReal(by_);
	int ok = isReal(br_) && XLENGTH(br_) >= 2 &&
		XLENGTH(br_) <= INT_MAX - 2 && R_FINITE(mn) && R_FINITE(by) &&
		by > 0;
	const double *br = ok ? REAL(br_) : NULL;
	int nb = ok ? (int) XLENGTH(br_) : 0;
	for (int b = 0; ok && b < nb; b++)
		if (ISNAN(br[b]) || (b && br[b - 1] > br[b]))
			ok = 0;
	if (!ok) {
		SEXP out = PROTECT(allocVector(INTSXP, 2));
		INTEGER(out)[0] = NA_INTEGER;
		INTEGER(out)[1] = INTEGER(VECTOR_ELT(st, S_LENS))[2];
		UNPROTECT(1);
		return out;
	}
	return run_pass(st, br, nb, mn, by, replace_, count_);
}

/* ss_pass with the breaks seq(from = mn, to = mx, by = by) made here, in
 * the case where seq.default computes pmin(from + (0L:n) * by, to) with
 * n = as.integer((to - from)/by + 1e-10) and returns at least two
 * breaks; those are sorted and finite, so binCode's guard holds. In any
 * other case it returns NULL, changes nothing, and the caller runs seq()
 * and ss_pass. */
SEXP ss_pass_by(SEXP st, SEXP mn_, SEXP mx_, SEXP by_, SEXP replace_,
	SEXP count_)
{
	if (!isReal(mn_) || !isReal(mx_) || !isReal(by_) || XLENGTH(mn_) != 1 ||
	    XLENGTH(mx_) != 1 || XLENGTH(by_) != 1)
		return R_NilValue;
	double from = REAL(mn_)[0], to = REAL(mx_)[0], by = REAL(by_)[0];
	if (!R_FINITE(from) || !R_FINITE(to) || !R_FINITE(by) || !(by > 0))
		return R_NilValue;
	double del = to - from;
	if (!R_FINITE(del) || !(del > 0))
		return R_NilValue;
	double n = del / by;
	if (!R_FINITE(n) || n > INT_MAX)
		return R_NilValue;
	if (fabs(del) / fmax(fabs(to), fabs(from)) < 100 * DBL_EPSILON)
		return R_NilValue;
	SEXP buf = VECTOR_ELT(st, S_BR);
	double t = n + 1e-10;
	if (!(t >= 1) || !(t < (double) XLENGTH(buf) - 1))
		return R_NilValue;
	int ni = (int) t;
	double *br = REAL(buf);
	for (int k = 0; k <= ni; k++) {
		/* R rounds (0L:n) * by and from + that separately */
		volatile double kb = (double) k * by;
		double v = from + kb;
		br[k] = to < v ? to : v;
	}
	return run_pass(st, br, ni + 1, from, by, replace_, count_);
}

static SEXP run_pass(SEXP st, const double *br, int nb, double mn,
	double by, SEXP replace_, SEXP count_)
{
	int *lens = INTEGER(VECTOR_ELT(st, S_LENS));
	int lf = lens[0], lp = lens[1], nLeft = lens[2], nDone = lens[3];
	SEXP out = PROTECT(allocVector(INTSXP, 2));
	INTEGER(out)[0] = NA_INTEGER;
	INTEGER(out)[1] = nLeft;

	const double *f = REAL(VECTOR_ELT(st, S_F));
	const double *p = REAL(VECTOR_ELT(st, S_P));
	int *fbin = INTEGER(VECTOR_ELT(st, S_FBIN));
	int *pbin = INTEGER(VECTOR_ELT(st, S_PBIN));
	int *fN = INTEGER(VECTOR_ELT(st, S_FN));
	int *pN = INTEGER(VECTOR_ELT(st, S_PN));
	int *rank = INTEGER(VECTOR_ELT(st, S_RANK));
	if ((R_xlen_t) nb + 1 > XLENGTH(VECTOR_ELT(st, S_FN))) {
		fN = (int *) R_alloc(nb + 1, sizeof(int));
		pN = (int *) R_alloc(nb + 1, sizeof(int));
		rank = (int *) R_alloc(nb + 1, sizeof(int));
	}
	memset(fN, 0, (nb + 1) * sizeof(int));
	memset(pN, 0, (nb + 1) * sizeof(int));

	/* Bin codes and counts; NA (outside the breaks) is bin nb */
	for (int j = 0; j < lf; j++)
		fN[fbin[j] = bin_of(f[j], br, nb, mn, by)]++;
	for (int j = 0; j < lp; j++)
		pN[pbin[j] = bin_of(p[j], br, nb, mn, by)]++;

	/* A bin qualifies with at least one focal and no fewer pool */
	int any = 0;
	for (int j = 0; j < lf; j++)
		if (fN[fbin[j]] <= pN[fbin[j]]) {
			any = 1;
			break;
		}
	if (!any) {
		INTEGER(out)[0] = 0;
		UNPROTECT(1);
		return out;
	}

	/* Qualifying bins ranked by first appearance in the pool */
	int *binOf = INTEGER(VECTOR_ELT(st, S_BINOF));
	int *poff = INTEGER(VECTOR_ELT(st, S_POFF));
	int *foff = INTEGER(VECTOR_ELT(st, S_FOFF));
	memset(rank, 0, (nb + 1) * sizeof(int));
	int K = 0;
	for (int j = 0; j < lp; j++) {
		int b = pbin[j];
		if (!rank[b] && fN[b] > 0 && fN[b] <= pN[b]) {
			rank[b] = ++K;
			binOf[K] = b;
		}
	}

	/* Positions of each qualifying bin, in row order, grouped by rank;
	 * fN and pN become each bin's write cursor */
	poff[1] = foff[1] = 0;
	for (int k = 1; k <= K; k++) {
		int b = binOf[k];
		poff[k + 1] = poff[k] + pN[b];
		foff[k + 1] = foff[k] + fN[b];
		pN[b] = poff[k];
		fN[b] = foff[k];
	}
	int *pgrp = INTEGER(VECTOR_ELT(st, S_PGRP));
	int *fgrp = INTEGER(VECTOR_ELT(st, S_FGRP));
	for (int j = 0; j < lp; j++)
		if (rank[pbin[j]])
			pgrp[pN[pbin[j]]++] = j;
	for (int j = 0; j < lf; j++)
		if (rank[fbin[j]])
			fgrp[fN[fbin[j]]++] = j;

	/* Each bin's draw as R's do_sample makes it: one R_unif_index per
	 * element with replacement or for a single element, otherwise its
	 * partial shuffle; sampled pool rows are marked with bin 0 */
	int rep = asLogical(replace_);
	int *fIdx = INTEGER(VECTOR_ELT(st, S_FIDX));
	int *pIdx = INTEGER(VECTOR_ELT(st, S_PIDX));
	int *outF = INTEGER(VECTOR_ELT(st, S_OUTF)) + nDone;
	int *outP = INTEGER(VECTOR_ELT(st, S_OUTP)) + nDone;
	GetRNGstate();
	for (int k = 1; k <= K; k++) {
		int m = poff[k + 1] - poff[k], s = foff[k + 1] - foff[k];
		int *seg = pgrp + poff[k];
		const int *fseg = fgrp + foff[k];
		if (rep || s < 2) {
			double dn = (double) m;
			for (int t = 0; t < s; t++) {
				int pos = seg[(int) R_unif_index(dn)];
				*outF++ = fIdx[fseg[t]];
				*outP++ = pIdx[pos];
				pbin[pos] = 0;
			}
		} else {
			int mm = m;
			for (int t = 0; t < s; t++) {
				int j = (int) R_unif_index((double) mm);
				int pos = seg[j];
				seg[j] = seg[--mm];
				*outF++ = fIdx[fseg[t]];
				*outP++ = pIdx[pos];
				pbin[pos] = 0;
			}
		}
	}
	PutRNGstate();

	/* Count the matched rows out and drop them, keeping row order */
	int count = asLogical(count_) == 1;
	double *fw = REAL(VECTOR_ELT(st, S_F)), *pw = REAL(VECTOR_ELT(st, S_P));
	int *cnt = INTEGER(VECTOR_ELT(st, S_CNT));
	const int *code = INTEGER(VECTOR_ELT(st, S_CODE));
	int nf = (int) XLENGTH(VECTOR_ELT(st, S_F));
	int w = 0;
	for (int j = 0; j < lf; j++) {
		if (rank[fbin[j]]) {
			if (count)
				count_out(cnt, code, fIdx[j] - 1, &nLeft);
		} else {
			fw[w] = fw[j];
			fIdx[w++] = fIdx[j];
		}
	}
	int matched = lf - w;
	lens[0] = w;
	w = 0;
	for (int j = 0; j < lp; j++) {
		if (!pbin[j]) {
			if (count)
				count_out(cnt, code, nf + pIdx[j] - 1, &nLeft);
		} else {
			pw[w] = pw[j];
			pIdx[w++] = pIdx[j];
		}
	}
	lens[1] = w;
	lens[2] = nLeft;
	lens[3] = nDone + matched;
	INTEGER(out)[0] = matched;
	INTEGER(out)[1] = nLeft;
	UNPROTECT(1);
	return out;
}

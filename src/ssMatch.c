#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <stdint.h>

/* ssMatch's remaining focal and pool, kept in one list that only ssMatch
 * holds and that ss_pass updates in place. */
enum {
	S_F, S_FIDX, S_P, S_PIDX, S_CNT, S_CODE, S_LENS, S_OUTF, S_OUTP,
	S_FBIN, S_PBIN, S_FN, S_PN, S_RANK, S_BINOF, S_POFF, S_FOFF, S_PGRP,
	S_FGRP, S_BR, S_N
};
/* S_LENS: remaining focal, remaining pool, nLeft, matched so far */

/* The state list with the scores copied and the row indices filled;
 * S_CODE is code, S_CNT is allocated but not filled */
static SEXP state_new(SEXP fps, SEXP pps, SEXP code)
{
	R_xlen_t nf = XLENGTH(fps), np = XLENGTH(pps);
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
	int *fi = INTEGER(VECTOR_ELT(st, S_FIDX));
	for (R_xlen_t j = 0; j < nf; j++)
		fi[j] = (int) j + 1;
	int *pi = INTEGER(VECTOR_ELT(st, S_PIDX));
	for (R_xlen_t j = 0; j < np; j++)
		pi[j] = (int) j + 1;
	int *lens = INTEGER(VECTOR_ELT(st, S_LENS));
	lens[0] = (int) nf;
	lens[1] = (int) np;
	lens[3] = 0;
	UNPROTECT(1);
	return st;
}

SEXP ss_init(SEXP fps, SEXP pps, SEXP cnt, SEXP code, SEXP nLeft)
{
	if (!isReal(fps) || !isReal(pps) || !isInteger(cnt) ||
	    !isInteger(code) || !isInteger(nLeft) || XLENGTH(nLeft) != 1)
		error("ss_init: invalid arguments");
	R_xlen_t nf = XLENGTH(fps), np = XLENGTH(pps);
	if (nf + np > INT_MAX / 2 || XLENGTH(code) != nf + np ||
	    XLENGTH(cnt) != nf + np)
		error("ss_init: invalid lengths");

	SEXP st = PROTECT(state_new(fps, pps, code));
	if (nf + np)
		memcpy(INTEGER(VECTOR_ELT(st, S_CNT)), INTEGER(cnt),
			(nf + np) * sizeof(int));
	INTEGER(VECTOR_ELT(st, S_LENS))[2] = INTEGER(nLeft)[0];
	UNPROTECT(1);
	return st;
}

/* Sorted mode, held in the state's "ss_sorted" attribute (the state list
 * keeps its 20 elements): the focal and the pool ordered by score, with
 * alive flags by row. S_F, S_FIDX, S_P and S_PIDX then stay the full
 * input, S_LENS[0:1] count the alive rows, and ss_compact turns the state
 * back into the compacted form the scan and the R fallback read. */
enum {
	Q_FORD, Q_FSV, Q_PORD, Q_PSV, Q_FAL, Q_PAL, Q_PTR, Q_RUN, Q_PTAG,
	Q_KEYS, Q_N
};
/* Q_PTR: first and one past the last possibly alive sorted position,
 * focal then pool; Q_RUN: R_W ints per focal bin of a pass */
enum { R_B, R_FS, R_FE, R_FN, R_LO, R_HI, R_PN, R_W };

static SEXP sorted_sym(void)
{
	static SEXP s = NULL;
	if (!s)
		s = install("ss_sorted");
	return s;
}

static SEXP sorted_of(SEXP st)
{
	SEXP q = getAttrib(st, sorted_sym());
	return TYPEOF(q) == VECSXP ? q : NULL;
}

/* Write the alive rows of S_F/S_FIDX and S_P/S_PIDX to their fronts in
 * row order, as run_pass's compaction leaves them, and leave sorted mode */
SEXP ss_compact(SEXP st)
{
	SEXP q = sorted_of(st);
	if (!q)
		return R_NilValue;
	int *lens = INTEGER(VECTOR_ELT(st, S_LENS));
	for (int s = 0; s < 2; s++) {
		double *x = REAL(VECTOR_ELT(st, s ? S_P : S_F));
		int *idx = INTEGER(VECTOR_ELT(st, s ? S_PIDX : S_FIDX));
		const Rbyte *al = RAW(VECTOR_ELT(q, s ? Q_PAL : Q_FAL));
		int n = (int) XLENGTH(VECTOR_ELT(st, s ? S_P : S_F)), w = 0;
		for (int j = 0; j < n; j++)
			if (al[j]) {
				x[w] = x[j];
				idx[w++] = idx[j];
			}
		lens[s] = w;
	}
	setAttrib(st, sorted_sym(), R_NilValue);
	return R_NilValue;
}

typedef struct {
	uint64_t k;
	int i;
} kv;

static int cmp_kv(const void *a, const void *b)
{
	uint64_t x = ((const kv *) a)->k, y = ((const kv *) b)->k;
	return (x > y) - (x < y);
}

/* An unsigned key in the order of the double x (not NaN); 0 and -0 share
 * one */
static inline uint64_t dkey(double x)
{
	uint64_t u;
	if (x == 0)
		x = 0.0;
	memcpy(&u, &x, sizeof u);
	return (u >> 63) ? ~u : u | ((uint64_t) 1 << 63);
}

/* ord (0-based rows) and sv: x in non-decreasing order, ties in any
 * order, by buckets on the top bits of the key then a sort per bucket;
 * 0 when x holds NaN or the result is not in order */
#ifndef SS_SORT_BITS
#define SS_SORT_BITS 20
#endif
static int sort_scores(const double *x, int n, int *ord, double *sv)
{
	if (n == 0)
		return 1;
	kv *a = (kv *) R_alloc(n, sizeof(kv)), *b = (kv *) R_alloc(n, sizeof(kv));
	uint64_t lo = UINT64_MAX, hi = 0;
	for (int i = 0; i < n; i++) {
		if (ISNAN(x[i]))
			return 0;
		uint64_t k = dkey(x[i]);
		a[i].k = k;
		a[i].i = i;
		if (k < lo)
			lo = k;
		if (k > hi)
			hi = k;
	}
	int bits = 1;
	while ((1 << bits) < n && bits < SS_SORT_BITS)
		bits++;
	uint64_t span = hi - lo;
	int top = 0;
	while (top < 64 && (span >> top) != 0)
		top++;
	int shift = top > bits ? top - bits : 0;
	int nb = (int) (span >> shift) + 1;
	int *cnt = (int *) R_alloc(nb + 1, sizeof(int));
	memset(cnt, 0, (nb + 1) * sizeof(int));
	for (int i = 0; i < n; i++)
		cnt[((a[i].k - lo) >> shift) + 1]++;
	for (int q = 1; q <= nb; q++)
		cnt[q] += cnt[q - 1];
	for (int i = 0; i < n; i++)
		b[cnt[(a[i].k - lo) >> shift]++] = a[i];
	/* bucket q now ends at cnt[q] */
	for (int q = 0, s = 0; q < nb; s = cnt[q++]) {
		int m = cnt[q] - s;
		kv *v = b + s;
		if (m <= 1)
			continue;
		if (m <= 16) {
			for (int t = 1; t < m; t++) {
				kv w = v[t];
				int j = t;
				while (j > 0 && v[j - 1].k > w.k) {
					v[j] = v[j - 1];
					j--;
				}
				v[j] = w;
			}
		} else {
			int same = 1;
			for (int t = 1; t < m && same; t++)
				same = v[t].k == v[0].k;
			if (!same)
				qsort(v, m, sizeof(kv), cmp_kv);
		}
	}
	for (int t = 0; t < n; t++) {
		ord[t] = b[t].i;
		sv[t] = x[b[t].i];
		if (t && !(sv[t - 1] <= sv[t]))
			return 0;
	}
	return 1;
}

/* The "ss_sorted" list, its vectors allocated */
static SEXP sorted_new(int nf, int np)
{
	SEXP q = PROTECT(allocVector(VECSXP, Q_N));
	SET_VECTOR_ELT(q, Q_FORD, allocVector(INTSXP, nf));
	SET_VECTOR_ELT(q, Q_FSV, allocVector(REALSXP, nf));
	SET_VECTOR_ELT(q, Q_PORD, allocVector(INTSXP, np));
	SET_VECTOR_ELT(q, Q_PSV, allocVector(REALSXP, np));
	SET_VECTOR_ELT(q, Q_FAL, allocVector(RAWSXP, nf));
	SET_VECTOR_ELT(q, Q_PAL, allocVector(RAWSXP, np));
	SET_VECTOR_ELT(q, Q_PTR, allocVector(INTSXP, 4));
	SET_VECTOR_ELT(q, Q_RUN, allocVector(INTSXP, (R_xlen_t) R_W * (nf + 1)));
	SET_VECTOR_ELT(q, Q_PTAG, allocVector(INTSXP, np));
	SET_VECTOR_ELT(q, Q_KEYS, allocVector(INTSXP, 2 * ((R_xlen_t) nf + 1)));
	UNPROTECT(1);
	return q;
}

/* Sort f and p into q; FALSE (q partly written) when either holds NaN */
static int sorted_fill(SEXP q, const double *f, int nf, const double *p,
	int np)
{
	for (int s = 0; s < 2; s++) {
		int n = s ? np : nf;
		if (!sort_scores(s ? p : f, n,
			INTEGER(VECTOR_ELT(q, s ? Q_PORD : Q_FORD)),
			REAL(VECTOR_ELT(q, s ? Q_PSV : Q_FSV))))
			return 0;
		if (n)
			memset(RAW(VECTOR_ELT(q, s ? Q_PAL : Q_FAL)), 1, n);
	}
	int *ptr = INTEGER(VECTOR_ELT(q, Q_PTR));
	ptr[0] = 0;
	ptr[1] = nf;
	ptr[2] = 0;
	ptr[3] = np;
	if (np)
		memset(INTEGER(VECTOR_ELT(q, Q_PTAG)), 0, np * sizeof(int));
	return 1;
}

/* Enter sorted mode on a state no pass has changed, when no score is NaN;
 * TRUE when it did */
SEXP ss_sort(SEXP st)
{
	const int *lens = INTEGER(VECTOR_ELT(st, S_LENS));
	int nf = (int) XLENGTH(VECTOR_ELT(st, S_F));
	int np = (int) XLENGTH(VECTOR_ELT(st, S_P));
	if (sorted_of(st) || lens[0] != nf || lens[1] != np || lens[3] != 0)
		return ScalarLogical(FALSE);
	SEXP q = PROTECT(sorted_new(nf, np));
	if (!sorted_fill(q, REAL(VECTOR_ELT(st, S_F)), nf,
		REAL(VECTOR_ELT(st, S_P)), np)) {
		UNPROTECT(1);
		return ScalarLogical(FALSE);
	}
	setAttrib(st, sorted_sym(), q);
	UNPROTECT(1);
	return ScalarLogical(TRUE);
}

/* ss_init and ss_sort in one, with code, cnt and nLeft taken from the
 * sorted scores rather than from match(x, x) on x = c(fps, pps): each
 * maximal run of equal scores (==, so 0 and -0 together) across the
 * sorted focal and pool is one score, whose code is its first position
 * in x, 1 + its smallest focal row or else nf + 1 + its smallest pool row
 * (match() on doubles also ties 0 and -0). NULL when either side holds
 * NaN, with nothing changed: the caller then calls match(), ss_init and
 * ss_sort. */
SEXP ss_init_sorted(SEXP fps, SEXP pps)
{
	if (!isReal(fps) || !isReal(pps))
		error("ss_init_sorted: invalid arguments");
	R_xlen_t nf_ = XLENGTH(fps), np_ = XLENGTH(pps);
	if (nf_ + np_ > INT_MAX / 2)
		error("ss_init_sorted: invalid lengths");
	int nf = (int) nf_, np = (int) np_;
	SEXP q = PROTECT(sorted_new(nf, np));
	if (!sorted_fill(q, REAL(fps), nf, REAL(pps), np)) {
		UNPROTECT(1);
		return R_NilValue;
	}
	SEXP code_ = PROTECT(allocVector(INTSXP, nf_ + np_));
	SEXP st = PROTECT(state_new(fps, pps, code_));
	int *code = INTEGER(code_);
	int *cnt = INTEGER(VECTOR_ELT(st, S_CNT));
	if (nf + np)
		memset(cnt, 0, (nf_ + np_) * sizeof(int));
	const int *ford = INTEGER(VECTOR_ELT(q, Q_FORD));
	const int *pord = INTEGER(VECTOR_ELT(q, Q_PORD));
	const double *fsv = REAL(VECTOR_ELT(q, Q_FSV));
	const double *psv = REAL(VECTOR_ELT(q, Q_PSV));
	int i = 0, j = 0, groups = 0;
	while (i < nf || j < np) {
		double v = i < nf && (j >= np || fsv[i] <= psv[j]) ? fsv[i] : psv[j];
		int i0 = i, j0 = j, fmin = INT_MAX, pmin = INT_MAX;
		/* sort_scores leaves ties in any order */
		for (; i < nf && fsv[i] == v; i++)
			if (ford[i] < fmin)
				fmin = ford[i];
		for (; j < np && psv[j] == v; j++)
			if (pord[j] < pmin)
				pmin = pord[j];
		int c = i > i0 ? fmin + 1 : nf + pmin + 1;
		for (int t = i0; t < i; t++)
			code[ford[t]] = c;
		for (int t = j0; t < j; t++)
			code[nf + pord[t]] = c;
		cnt[c - 1] = (i - i0) + (j - j0);
		groups++;
	}
	INTEGER(VECTOR_ELT(st, S_LENS))[2] = groups;
	setAttrib(st, sorted_sym(), q);
	UNPROTECT(3);
	return st;
}

/* Move the alive bounds of the sorted focal (s = 0) or pool (s = 1) past
 * rows that have been matched */
static void advance(SEXP q, int s)
{
	int *ptr = INTEGER(VECTOR_ELT(q, Q_PTR)) + 2 * s;
	const int *ord = INTEGER(VECTOR_ELT(q, s ? Q_PORD : Q_FORD));
	const Rbyte *al = RAW(VECTOR_ELT(q, s ? Q_PAL : Q_FAL));
	int a = ptr[0], e = ptr[1];
	while (a < e && !al[ord[a]])
		a++;
	while (e > a && !al[ord[e - 1]])
		e--;
	ptr[0] = a;
	ptr[1] = e;
}

/* min and max of c(f, p) as R's min() and max() give them for doubles
 * (the first of equal values, so the first of 0 and -0); NaN when
 * either holds NA or NaN, which seq() then rejects either way. */
SEXP ss_range(SEXP st)
{
	SEXP q = sorted_of(st);
	if (q) {
		/* the ends of the sorted alive rows; a zero may be 0 or -0, and
		 * which comes first in row order decides, so compact for that */
		advance(q, 0);
		advance(q, 1);
		const int *ptr = INTEGER(VECTOR_ELT(q, Q_PTR));
		const double *fsv = REAL(VECTOR_ELT(q, Q_FSV));
		const double *psv = REAL(VECTOR_ELT(q, Q_PSV));
		int any = 0;
		double mn = R_PosInf, mx = R_NegInf;
		if (ptr[0] < ptr[1]) {
			mn = fsv[ptr[0]];
			mx = fsv[ptr[1] - 1];
			any = 1;
		}
		if (ptr[2] < ptr[3]) {
			double a = psv[ptr[2]], b = psv[ptr[3] - 1];
			if (!any || a < mn)
				mn = a;
			if (!any || b > mx)
				mx = b;
			any = 1;
		}
		if (!any || (mn != 0 && mx != 0)) {
			SEXP out = PROTECT(allocVector(REALSXP, 2));
			REAL(out)[0] = mn;
			REAL(out)[1] = mx;
			UNPROTECT(1);
			return out;
		}
		ss_compact(st);
	}
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
static SEXP run_pass_sorted(SEXP st, SEXP q, const double *br, int nb,
	double mn, double by, SEXP replace_, SEXP count_);

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
	/* this scans the compacted rows, as does the R fallback after NA */
	ss_compact(st);
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
	SEXP q = sorted_of(st);
	if (q)
		return run_pass_sorted(st, q, br, ni + 1, from, by, replace_, count_);
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

/* The first position in [s, e) of the non-decreasing v whose value is
 * above t (or at least t when !strict), galloping forward from s; e when
 * there is none */
static inline int gallop(const double *v, int s, int e, double t, int strict)
{
	if (s >= e || (strict ? v[s] > t : v[s] >= t))
		return s;
	long long lo = s, hi, step = 1;
	for (;;) {
		hi = lo + step;
		if (hi >= e) {
			hi = e;
			break;
		}
		if (strict ? v[hi] > t : v[hi] >= t)
			break;
		lo = hi;
		step *= 2;
	}
	while (hi - lo > 1) {
		long long mid = lo + (hi - lo) / 2;
		if (strict ? v[mid] > t : v[mid] >= t)
			hi = mid;
		else
			lo = mid;
	}
	return (int) hi;
}

static int cmp_int(const void *a, const void *b)
{
	int x = *(const int *) a, y = *(const int *) b;
	return (x > y) - (x < y);
}

static void sort_ints(int *x, int n)
{
	if (n > 32) {
		qsort(x, n, sizeof(int), cmp_int);
		return;
	}
	for (int i = 1; i < n; i++) {
		int v = x[i], j = i;
		while (j > 0 && x[j - 1] > v) {
			x[j] = x[j - 1];
			j--;
		}
		x[j] = v;
	}
}

/* Pool rows of the qualifying bins over this share of the alive pool are
 * grouped by one sweep over the row indices rather than sorted per bin */
#ifndef SS_SWEEP_SHARE
#define SS_SWEEP_SHARE 16
#endif

/* run_pass on the sorted state: the alive focal rows in score order fall
 * in non-decreasing bins, and the pool rows of bin b are one range of the
 * sorted pool, br[b-1] < x <= br[b] (br[0] <= x for b = 1, x > br[nb-1]
 * for the NA bin nb), the rows bin_of puts there. Bins are ranked by
 * their smallest alive pool row, the first appearance in the compacted
 * pool, and each lists its rows in row order, so the draws and pairs are
 * run_pass's; matched rows are marked dead rather than compacted. */
static SEXP run_pass_sorted(SEXP st, SEXP q, const double *br, int nb,
	double mn, double by, SEXP replace_, SEXP count_)
{
	int *lens = INTEGER(VECTOR_ELT(st, S_LENS));
	int nLeft = lens[2], nDone = lens[3];
	advance(q, 0);
	advance(q, 1);
	const int *ptr = INTEGER(VECTOR_ELT(q, Q_PTR));
	int fa = ptr[0], fz = ptr[1], pa = ptr[2], pz = ptr[3];
	const int *ford = INTEGER(VECTOR_ELT(q, Q_FORD));
	const int *pord = INTEGER(VECTOR_ELT(q, Q_PORD));
	const double *fsv = REAL(VECTOR_ELT(q, Q_FSV));
	const double *psv = REAL(VECTOR_ELT(q, Q_PSV));
	Rbyte *fal = RAW(VECTOR_ELT(q, Q_FAL));
	Rbyte *pal = RAW(VECTOR_ELT(q, Q_PAL));

	SEXP out = PROTECT(allocVector(INTSXP, 2));
	INTEGER(out)[0] = 0;
	INTEGER(out)[1] = nLeft;
	if (fa >= fz) {
		UNPROTECT(1);
		return out;
	}
	/* bins rise with the score only above br[0], which ss_range's min
	 * guarantees; otherwise scan */
	if (fsv[fa] < br[0] || (pa < pz && psv[pa] < br[0])) {
		UNPROTECT(1);
		ss_compact(st);
		return run_pass(st, br, nb, mn, by, replace_, count_);
	}

	/* Each focal bin: its sorted focal span, its pool range and counts */
	int *run = INTEGER(VECTOR_ELT(q, Q_RUN));
	int nr = 0, cur = pa;
	for (int s = fa; s < fz;) {
		if (!fal[ford[s]]) {
			s++;
			continue;
		}
		int b = bin_of(fsv[s], br, nb, mn, by);
		double top = b < nb ? br[b] : R_PosInf;
		int fs = s, fN = 1;
		for (s++; s < fz; s++)
			if (fal[ford[s]]) {
				if (!(fsv[s] <= top))
					break;
				fN++;
			}
		int lo, hi;
		if (b == nb) {
			lo = gallop(psv, cur, pz, br[nb - 1], 1);
			hi = pz;
		} else {
			lo = b == 1 ? gallop(psv, cur, pz, br[0], 0) :
				gallop(psv, cur, pz, br[b - 1], 1);
			hi = gallop(psv, lo, pz, br[b], 1);
		}
		cur = hi;
		int pN = 0;
		for (int t = lo; t < hi; t++)
			pN += pal[pord[t]];
		int *r = run + (R_xlen_t) R_W * nr++;
		r[R_B] = b;
		r[R_FS] = fs;
		r[R_FE] = s;
		r[R_FN] = fN;
		r[R_LO] = lo;
		r[R_HI] = hi;
		r[R_PN] = pN;
	}

	/* A bin qualifies with at least one focal and no fewer pool; keep the
	 * qualifying ones at the front, in score order */
	int K = 0, Fq = 0;
	long long Wq = 0;
	for (int i = 0; i < nr; i++) {
		int *r = run + (R_xlen_t) R_W * i;
		if (r[R_FN] <= r[R_PN]) {
			if (K != i)
				memcpy(run + (R_xlen_t) R_W * K, r, R_W * sizeof(int));
			K++;
			Fq += r[R_FN];
			Wq += r[R_PN];
		}
	}
	if (!K) {
		UNPROTECT(1);
		return out;
	}

	/* Each qualifying bin's pool and focal rows in row order, and the bins
	 * ranked by their first pool row (ord) */
	int *poff = INTEGER(VECTOR_ELT(st, S_POFF));
	int *foff = INTEGER(VECTOR_ELT(st, S_FOFF));
	int *ord = INTEGER(VECTOR_ELT(st, S_BINOF));
	int *pg = INTEGER(VECTOR_ELT(st, S_PGRP));
	int *fg = INTEGER(VECTOR_ELT(st, S_FGRP));
	int *keys = INTEGER(VECTOR_ELT(q, Q_KEYS));
	poff[0] = foff[0] = 0;
	for (int k = 0; k < K; k++) {
		const int *r = run + (R_xlen_t) R_W * k;
		poff[k + 1] = poff[k] + r[R_PN];
		foff[k + 1] = foff[k] + r[R_FN];
		int w = foff[k];
		for (int t = r[R_FS]; t < r[R_FE]; t++)
			if (fal[ford[t]])
				fg[w++] = ford[t];
		sort_ints(fg + foff[k], r[R_FN]);
	}
	if (Wq * SS_SWEEP_SHARE > lens[1]) {
		int *ptag = INTEGER(VECTOR_ELT(q, Q_PTAG));
		int jmin = INT_MAX, jmax = -1;
		for (int k = 0; k < K; k++) {
			const int *r = run + (R_xlen_t) R_W * k;
			for (int t = r[R_LO]; t < r[R_HI]; t++) {
				int j = pord[t];
				if (pal[j]) {
					ptag[j] = k + 1;
					if (j < jmin)
						jmin = j;
					if (j > jmax)
						jmax = j;
				}
			}
			keys[k] = poff[k];
		}
		int nk = 0;
		for (int j = jmin; j <= jmax; j++) {
			int g = ptag[j];
			if (g) {
				g--;
				if (keys[g] == poff[g])
					ord[nk++] = g;
				pg[keys[g]++] = j;
				ptag[j] = 0;
			}
		}
	} else {
		for (int k = 0; k < K; k++) {
			const int *r = run + (R_xlen_t) R_W * k;
			int w = poff[k];
			for (int t = r[R_LO]; t < r[R_HI]; t++)
				if (pal[pord[t]])
					pg[w++] = pord[t];
			sort_ints(pg + poff[k], r[R_PN]);
			keys[2 * k] = pg[poff[k]];
			keys[2 * k + 1] = k;
		}
		qsort(keys, K, 2 * sizeof(int), cmp_int);
		for (int k = 0; k < K; k++)
			ord[k] = keys[2 * k + 1];
	}

	/* Each bin's draw as R's do_sample makes it, as in run_pass; drawn
	 * pool rows and the focal rows are counted out once and marked dead */
	int rep = asLogical(replace_);
	int count = asLogical(count_) == 1;
	int nf = (int) XLENGTH(VECTOR_ELT(st, S_F));
	int *cnt = INTEGER(VECTOR_ELT(st, S_CNT));
	const int *code = INTEGER(VECTOR_ELT(st, S_CODE));
	int *outF = INTEGER(VECTOR_ELT(st, S_OUTF)) + nDone;
	int *outP = INTEGER(VECTOR_ELT(st, S_OUTP)) + nDone;
	int pdead = 0;
	GetRNGstate();
	for (int i = 0; i < K; i++) {
		int k = ord[i];
		int m = poff[k + 1] - poff[k], s = foff[k + 1] - foff[k];
		int *seg = pg + poff[k];
		const int *fseg = fg + foff[k];
		if (rep || s < 2) {
			double dn = (double) m;
			for (int t = 0; t < s; t++) {
				int j = seg[(int) R_unif_index(dn)];
				*outF++ = fseg[t] + 1;
				*outP++ = j + 1;
				if (pal[j]) {
					pal[j] = 0;
					pdead++;
					if (count)
						count_out(cnt, code, nf + j, &nLeft);
				}
			}
		} else {
			int mm = m;
			for (int t = 0; t < s; t++) {
				int jj = (int) R_unif_index((double) mm);
				int j = seg[jj];
				seg[jj] = seg[--mm];
				*outF++ = fseg[t] + 1;
				*outP++ = j + 1;
				pal[j] = 0;
				pdead++;
				if (count)
					count_out(cnt, code, nf + j, &nLeft);
			}
		}
	}
	PutRNGstate();
	for (int t = 0; t < Fq; t++) {
		fal[fg[t]] = 0;
		if (count)
			count_out(cnt, code, fg[t], &nLeft);
	}
	lens[0] -= Fq;
	lens[1] -= pdead;
	lens[2] = nLeft;
	lens[3] = nDone + Fq;
	INTEGER(out)[0] = Fq;
	INTEGER(out)[1] = nLeft;
	UNPROTECT(1);
	return out;
}

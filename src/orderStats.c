#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <math.h>

/* Order statistics of a double vector without sorting a copy of it.
 *
 * order_stats(x, ranks) returns, for each rank r in ranks (1-based, strictly
 * increasing, within 1..length(x)), the value sort(x)[r], or NULL when x
 * holds a NaN or NA. The values are exact: the result equals what
 * sort(x, partial = ranks)[ranks] gives, up to the sign of a zero where -0
 * and 0 tie.
 *
 * Large x: a fixed-stride sample of x, sorted, brackets each rank by two
 * sample values [L, U]; a pass per bracket counts the elements below L and
 * copies those in [L, U] into a buffer; a rank r with
 * less < r <= less + count is then the (r - less)-th smallest of its buffer.
 * A rank the brackets miss (an unrepresentative sample) or a full buffer
 * (many ties) sends every rank to the exact path below. Small x, or that
 * fallback: select on a copy. */

#define OS_SAMPLE 4096
#define OS_SMALL (8 * OS_SAMPLE)

static void os_swap(double *a, R_xlen_t i, R_xlen_t j)
{
	double t = a[i];
	a[i] = a[j];
	a[j] = t;
}

/* Rearrange a[lo..hi] so that a[k] is the k-th smallest for every k in
 * ks[0..nk-1] (0-based, increasing, within lo..hi): median-of-three Hoare
 * partitions, a sort of the range when it is small or after depth levels.
 * No element may be NaN. */
static void os_select(double *a, R_xlen_t lo, R_xlen_t hi, const R_xlen_t *ks,
	int nk, int depth)
{
	while (nk > 0) {
		if (hi - lo < 16 || depth-- <= 0) {
			R_qsort(a + lo, 1, (size_t) (hi - lo + 1));
			return;
		}
		R_xlen_t mid = lo + (hi - lo) / 2;
		if (a[mid] < a[lo]) os_swap(a, mid, lo);
		if (a[hi] < a[lo]) os_swap(a, hi, lo);
		if (a[hi] < a[mid]) os_swap(a, hi, mid);
		double p = a[mid];
		R_xlen_t i = lo - 1, j = hi + 1;
		for (;;) {
			do i++; while (a[i] < p);
			do j--; while (a[j] > p);
			if (i >= j) break;
			os_swap(a, i, j);
		}
		/* a[lo..j] <= p <= a[j+1..hi], lo <= j < hi */
		int nl = 0;
		while (nl < nk && ks[nl] <= j) nl++;
		if (nl < nk - nl) {
			os_select(a, lo, j, ks, nl, depth);
			ks += nl;
			nk -= nl;
			lo = j + 1;
		} else {
			os_select(a, j + 1, hi, ks + nl, nk - nl, depth);
			nk = nl;
			hi = j;
		}
	}
}

static int os_depth(R_xlen_t n)
{
	int d = 0;
	while (n > 1) {
		n >>= 1;
		d++;
	}
	return 2 * d + 4;
}

/* min(x) and max(x) as R's rmin and rmax take them from a vector without
 * NaN: from x[0], replaced only on a strict compare, so of tied values (-0
 * and 0) the first in row order */
static void os_minmax(const double *x, R_xlen_t n, double *mm)
{
	double mn = x[0], mx = x[0];
	for (R_xlen_t i = 1; i < n; i++) {
		double v = x[i];
		if (v < mn) mn = v;
		if (v > mx) mx = v;
	}
	mm[0] = mn;
	mm[1] = mx;
}

/* the exact path: a copy of x, selected; 0 when x holds a NaN. With mm,
 * also x's min and max into mm[0] and mm[1]. */
static int os_full(const double *x, R_xlen_t n, const int *r, int nr,
	double *out, double *mm)
{
	double *a = (double *) R_alloc((size_t) n, sizeof(double));
	R_xlen_t *ks = (R_xlen_t *) R_alloc((size_t) nr, sizeof(R_xlen_t));
	for (R_xlen_t i = 0; i < n; i++) {
		if (ISNAN(x[i])) return 0;
		a[i] = x[i];
	}
	if (mm) os_minmax(x, n, mm);
	for (int k = 0; k < nr; k++)
		ks[k] = (R_xlen_t) r[k] - 1;
	os_select(a, 0, n - 1, ks, nr, os_depth(n));
	for (int k = 0; k < nr; k++)
		out[k] = a[ks[k]];
	return 1;
}

/* with_range: the result also holds min(x) and max(x) after the order
 * statistics, from the pass that reads x first (x of length >= 1) */
static SEXP os_impl(SEXP x, SEXP ranks, int with_range)
{
	R_xlen_t n = XLENGTH(x);
	int nr = LENGTH(ranks);
	const double *px = REAL(x);
	const int *r = INTEGER(ranks);
	if (with_range && n < 1)
		return R_NilValue;
	SEXP out = PROTECT(allocVector(REALSXP, nr + (with_range ? 2 : 0)));
	double *po = REAL(out);
	double *mm = with_range ? po + nr : NULL;
	if (nr == 0) {
		if (mm) {
			for (R_xlen_t i = 0; i < n; i++)
				if (ISNAN(px[i])) {
					UNPROTECT(1);
					return R_NilValue;
				}
			os_minmax(px, n, mm);
		}
		UNPROTECT(1);
		return out;
	}
	if (n < OS_SMALL) {
		int ok = os_full(px, n, r, nr, po, mm);
		UNPROTECT(1);
		return ok ? out : R_NilValue;
	}

	/* the sample, sorted */
	int m = OS_SAMPLE;
	R_xlen_t stride = n / m;
	double *s = (double *) R_alloc((size_t) m, sizeof(double));
	for (int j = 0; j < m; j++) {
		s[j] = px[(R_xlen_t) j * stride + stride / 2];
		if (ISNAN(s[j])) {
			UNPROTECT(1);
			return R_NilValue;
		}
	}
	R_qsort(s, 1, (size_t) m);

	/* one bracket per run of ranks whose sample windows overlap: about
	 * four standard deviations of a sample rank either side */
	int *blo = (int *) R_alloc((size_t) nr, sizeof(int));
	int *bhi = (int *) R_alloc((size_t) nr, sizeof(int));
	int *first = (int *) R_alloc((size_t) nr + 1, sizeof(int));
	int nb = 0;
	for (int k = 0; k < nr; k++) {
		double q = (double) r[k] / (double) n;
		double pos = q * m - 0.5;
		double d = 4.0 * sqrt(m * q * (1.0 - q)) + 8.0;
		double flo = floor(pos - d), fhi = ceil(pos + d);
		int jlo = flo < 0 ? -1 : (int) flo;
		int jhi = fhi > m - 1 ? m : (int) fhi;
		if (nb > 0 && jlo <= bhi[nb - 1]) {
			if (jhi > bhi[nb - 1]) bhi[nb - 1] = jhi;
		} else {
			first[nb] = k;
			blo[nb] = jlo;
			bhi[nb] = jhi;
			nb++;
		}
	}
	first[nb] = nr;

	/* L = -Inf and U = Inf past the sample's ends; each buffer holds one and
	 * a half times the expected count, plus slack, and one spare slot */
	double *L = (double *) R_alloc((size_t) nb, sizeof(double));
	double *U = (double *) R_alloc((size_t) nb, sizeof(double));
	R_xlen_t *less = (R_xlen_t *) R_alloc((size_t) nb, sizeof(R_xlen_t));
	R_xlen_t *cnt = (R_xlen_t *) R_alloc((size_t) nb, sizeof(R_xlen_t));
	R_xlen_t *cap = (R_xlen_t *) R_alloc((size_t) nb, sizeof(R_xlen_t));
	double **buf = (double **) R_alloc((size_t) nb, sizeof(double *));
	R_xlen_t total = 0;
	for (int b = 0; b < nb; b++) {
		L[b] = blo[b] < 0 ? R_NegInf : s[blo[b]];
		U[b] = bhi[b] >= m ? R_PosInf : s[bhi[b]];
		double c = 1.5 * (double) (bhi[b] - blo[b] + 1) * (double) stride + 1024.0;
		cap[b] = c > (double) n ? n : (R_xlen_t) c;
		total += cap[b];
		less[b] = 0;
		cnt[b] = 0;
	}
	if (total > n / 2) {
		int ok = os_full(px, n, r, nr, po, mm);
		UNPROTECT(1);
		return ok ? out : R_NilValue;
	}
	for (int b = 0; b < nb; b++)
		buf[b] = (double *) R_alloc((size_t) cap[b] + 1, sizeof(double));

	/* one pass per bracket, without branches on the data: every element is
	 * stored at the buffer's next free slot (or its spare last slot once
	 * full) and the slot advances only for an element inside [L, U] */
	int nan = 0;
	for (int b = 0; b < nb; b++) {
		const double Lb = L[b], Ub = U[b];
		double *bb = buf[b];
		R_xlen_t lt = 0, c = 0, cp = cap[b];
		int nn = 0;
		if (b == 0 && mm) {
			double mn = px[0], mx = px[0];
			for (R_xlen_t i = 0; i < n; i++) {
				double v = px[i];
				nn |= (v != v);
				lt += (v < Lb);
				bb[c < cp ? c : cp] = v;
				c += (v >= Lb) & (v <= Ub);
				mn = (v < mn) ? v : mn;
				mx = (v > mx) ? v : mx;
			}
			mm[0] = mn;
			mm[1] = mx;
		} else {
			for (R_xlen_t i = 0; i < n; i++) {
				double v = px[i];
				nn |= (v != v);
				lt += (v < Lb);
				bb[c < cp ? c : cp] = v;
				c += (v >= Lb) & (v <= Ub);
			}
		}
		less[b] = lt;
		cnt[b] = c;
		nan |= nn;
	}
	if (nan) {
		UNPROTECT(1);
		return R_NilValue;
	}

	int hit = 1;
	for (int b = 0; b < nb && hit; b++) {
		if (cnt[b] > cap[b]) hit = 0;
		for (int k = first[b]; k < first[b + 1] && hit; k++)
			if (!(less[b] < r[k] && r[k] <= less[b] + cnt[b])) hit = 0;
	}
	if (!hit) {
		int ok = os_full(px, n, r, nr, po, mm);
		UNPROTECT(1);
		return ok ? out : R_NilValue;
	}

	R_xlen_t *ks = (R_xlen_t *) R_alloc((size_t) nr, sizeof(R_xlen_t));
	for (int b = 0; b < nb; b++) {
		int k0 = first[b], nk = first[b + 1] - k0;
		for (int k = 0; k < nk; k++)
			ks[k0 + k] = (R_xlen_t) r[k0 + k] - 1 - less[b];
		os_select(buf[b], 0, cnt[b] - 1, ks + k0, nk, os_depth(cnt[b]));
		for (int k = 0; k < nk; k++)
			po[k0 + k] = buf[b][ks[k0 + k]];
	}
	UNPROTECT(1);
	return out;
}

SEXP order_stats(SEXP x, SEXP ranks)
{
	return os_impl(x, ranks, 0);
}

/* order_stats(x, ranks) followed by min(x) and max(x) */
SEXP order_stats_range(SEXP x, SEXP ranks)
{
	return os_impl(x, ranks, 1);
}

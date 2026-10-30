#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <string.h>

/* nnMatch's candidate pool rows. With u the distinct focal scores sorted,
 * pool score p lies in interval k, u[k - 1] < p <= u[k] (k = 0 at or below
 * u[0], k = m above u[m - 1]). The rolling nearest join gives each u the
 * rows of the largest pool score <= u or of the smallest >= u, and each of
 * those is its interval's max or min. Returns the 1-based pool rows whose
 * score equals its interval's min or max, in increasing order. NULL when
 * fps or pps is not a double vector, either is empty or longer than
 * INT_MAX, or either holds a value that is not finite. */

/* buckets per distinct focal score in the index over u */
#define NN_BUCKETS 16

/* finite: not NaN, NA or infinite, without a call into R */
#define NN_FINITE(x) (fabs(x) <= DBL_MAX)

SEXP nn_keep(SEXP fps, SEXP pps)
{
	if (TYPEOF(fps) != REALSXP || TYPEOF(pps) != REALSXP)
		return R_NilValue;
	R_xlen_t nf_ = XLENGTH(fps), np_ = XLENGTH(pps);
	if (nf_ < 1 || np_ < 1 || nf_ > INT_MAX / (2 * NN_BUCKETS) ||
	    np_ > INT_MAX)
		return R_NilValue;
	int nf = (int) nf_, np = (int) np_;
	const double *f = REAL_RO(fps), *p = REAL_RO(pps);

	double *u = (double *) R_alloc(nf, sizeof(double));
	for (int i = 0; i < nf; i++) {
		if (!NN_FINITE(f[i]))
			return R_NilValue;
		u[i] = f[i];
	}
	R_qsort(u, 1, (size_t) nf);
	int m = 1;
	for (int i = 1; i < nf; i++)
		if (u[i] != u[m - 1])
			u[m++] = u[i];

	/* start[b]: the number of u whose bucket is below b, so the interval of
	 * a p in bucket b is at least start[b] (the bucket is monotone in p) */
	double lo = u[0], hi = u[m - 1], span = hi - lo;
	int B = m > 1 ? NN_BUCKETS * m : 1;
	double scale = m > 1 ? (double) B / span : 0;
	int use_b = m > 1 && R_FINITE(span) && span > 0 && R_FINITE(scale);
	int *start = NULL;
	if (use_b) {
		start = (int *) R_alloc((size_t) B + 1, sizeof(int));
		memset(start, 0, ((size_t) B + 1) * sizeof(int));
		for (int j = 0; j < m; j++) {
			int b = (int) ((u[j] - lo) * scale);
			if (b > B - 1)
				b = B - 1;
			start[b + 1]++;
		}
		for (int b = 1; b <= B; b++)
			start[b] += start[b - 1];
	}

	int *kk = (int *) R_alloc(np, sizeof(int));
	double *mn = (double *) R_alloc((size_t) m + 1, sizeof(double));
	double *mx = (double *) R_alloc((size_t) m + 1, sizeof(double));
	for (int k = 0; k <= m; k++) {
		mn[k] = R_PosInf;
		mx[k] = R_NegInf;
	}
	for (int i = 0; i < np; i++) {
		double x = p[i];
		if (!NN_FINITE(x))
			return R_NilValue;
		int k;
		if (x <= lo)
			k = 0;
		else if (x > hi)
			k = m;
		else if (use_b) {
			int b = (int) ((x - lo) * scale);
			if (b > B - 1)
				b = B - 1;
			/* two steps without a branch, as most scans are short; k stays
			 * below m, since u[m - 1] = hi >= x */
			k = start[b];
			k += u[k] < x;
			k += u[k] < x;
			while (u[k] < x)
				k++;
		} else {
			/* the first k with u[k] >= x; u[m - 1] = hi >= x */
			int a = 0, z = m - 1;
			while (a < z) {
				int mid = a + (z - a) / 2;
				if (u[mid] < x)
					a = mid + 1;
				else
					z = mid;
			}
			k = a;
		}
		kk[i] = k;
		double a = mn[k], z = mx[k];
		mn[k] = x < a ? x : a;
		mx[k] = x > z ? x : z;
	}

	/* kept rows written over kk, which is read ahead of the writes */
	int w = 0;
	for (int i = 0; i < np; i++) {
		int k = kk[i];
		if (p[i] == mn[k] || p[i] == mx[k])
			kk[w++] = i + 1;
	}
	SEXP out = PROTECT(allocVector(INTSXP, w));
	if (w)
		memcpy(INTEGER(out), kk, (size_t) w * sizeof(int));
	UNPROTECT(1);
	return out;
}

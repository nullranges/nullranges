#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <stdint.h>

/* The 1-based indices, in order, of the ranges (code, start, width > 0)
 * that overlap none of the intervals of their sequence: sequence c holds
 * intervals seg[c - 1] to seg[c] - 1 of i_start and i_end, sorted and
 * disjoint. NULL for a code outside 1..length(seg) - 1 or more ranges
 * than an integer index holds. */
SEXP outside_intervals(SEXP code, SEXP start, SEXP width, SEXP seg,
	SEXP i_start, SEXP i_end)
{
	R_xlen_t n = XLENGTH(code);
	int n_lev = LENGTH(seg) - 1;
	if (n > INT_MAX || XLENGTH(start) != n || XLENGTH(width) != n ||
		n_lev < 0)
		return R_NilValue;
	const int *c = INTEGER(code), *s = INTEGER(start), *w = INTEGER(width);
	const int *sg = INTEGER(seg), *is = INTEGER(i_start), *ie = INTEGER(i_end);
	unsigned char *hit = (unsigned char *) R_alloc((size_t) n + 1, 1);
	R_xlen_t n_keep = 0;
	int prev = 0, k = 0;
	for (R_xlen_t r = 0; r < n; r++) {
		int cc = c[r];
		if (cc < 1 || cc > n_lev)
			return R_NilValue;
		int lo = sg[cc - 1], hi = sg[cc], x = s[r];
		/* k: the first interval of the sequence ending at or after x; the
		 * previous range's k when it still is */
		if (cc != prev || (k > lo && ie[k - 1] >= x) ||
			(k < hi && ie[k] < x)) {
			int a = lo, b = hi;
			while (a < b) {
				int m = a + (b - a) / 2;
				if (ie[m] < x)
					a = m + 1;
				else
					b = m;
			}
			k = a;
			prev = cc;
		}
		hit[r] = k < hi && (int64_t) is[k] <= (int64_t) x + w[r] - 1;
		n_keep += !hit[r];
	}
	SEXP ans = PROTECT(allocVector(INTSXP, n_keep));
	int *a = INTEGER(ans);
	R_xlen_t j = 0;
	for (R_xlen_t r = 0; r < n; r++)
		if (!hit[r])
			a[j++] = (int) (r + 1);
	UNPROTECT(1);
	return ans;
}

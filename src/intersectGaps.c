#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <stdint.h>

/* The gaps of [lo, hi) with an end at or after x start at the returned
 * index; k, the previous range's answer, is kept when it still is. */
static inline int first_gap(const int *ge, int lo, int hi, int k, int same,
	int x)
{
	if (same && (k == lo || ge[k - 1] < x) && (k == hi || ge[k] >= x))
		return k;
	int a = lo, b = hi;
	while (a < b) {
		int m = a + (b - a) / 2;
		if (ge[m] < x)
			a = m + 1;
		else
			b = m;
	}
	return a;
}

/* Each range (code, start, end) clipped to every gap of its sequence it
 * overlaps, by range and then by gap: sequence c holds gaps seg[c - 1] to
 * seg[c] - 1 of g_start and g_end, sorted and disjoint, and a code
 * outside 1..length(seg) - 1 has none. Returns list(q, start, width): the
 * 1-based index of each piece's range, its start and its width. NULL for
 * an argument not integer, lengths that disagree, an NA code, start or
 * end, a range whose width an integer cannot hold, or more pieces than an
 * integer vector indexes. */
SEXP intersect_gaps(SEXP code, SEXP start, SEXP end, SEXP seg,
	SEXP g_start, SEXP g_end)
{
	if (TYPEOF(code) != INTSXP || TYPEOF(start) != INTSXP ||
		TYPEOF(end) != INTSXP || TYPEOF(seg) != INTSXP ||
		TYPEOF(g_start) != INTSXP || TYPEOF(g_end) != INTSXP)
		return R_NilValue;
	R_xlen_t n = XLENGTH(code);
	int n_lev = LENGTH(seg) - 1;
	if (n > INT_MAX || XLENGTH(start) != n || XLENGTH(end) != n ||
		n_lev < 0 || XLENGTH(g_end) != XLENGTH(g_start) ||
		INTEGER(seg)[n_lev] != XLENGTH(g_start))
		return R_NilValue;
	const int *c = INTEGER(code), *s = INTEGER(start), *e = INTEGER(end);
	const int *sg = INTEGER(seg), *gs = INTEGER(g_start), *ge = INTEGER(g_end);

	/* pass 1: count the pieces */
	R_xlen_t total = 0;
	int prev = 0, k = 0;
	for (R_xlen_t r = 0; r < n; r++) {
		int cc = c[r], x = s[r], y = e[r];
		if (cc == NA_INTEGER || x == NA_INTEGER || y == NA_INTEGER ||
			(int64_t) y - x >= INT_MAX)
			return R_NilValue;
		if (cc < 1 || cc > n_lev)
			continue;
		int hi = sg[cc];
		k = first_gap(ge, sg[cc - 1], hi, k, cc == prev, x);
		prev = cc;
		int j = k;
		while (j < hi && gs[j] <= y)
			j++;
		total += j - k;
	}
	if (total > INT_MAX)
		return R_NilValue;

	/* pass 2: fill them */
	SEXP ans = PROTECT(allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ans, 0, allocVector(INTSXP, total));
	SET_VECTOR_ELT(ans, 1, allocVector(INTSXP, total));
	SET_VECTOR_ELT(ans, 2, allocVector(INTSXP, total));
	int *q = INTEGER(VECTOR_ELT(ans, 0)), *ns = INTEGER(VECTOR_ELT(ans, 1)),
		*w = INTEGER(VECTOR_ELT(ans, 2));
	R_xlen_t o = 0;
	prev = 0;
	k = 0;
	for (R_xlen_t r = 0; r < n; r++) {
		int cc = c[r], x = s[r], y = e[r];
		if (cc < 1 || cc > n_lev)
			continue;
		int hi = sg[cc];
		k = first_gap(ge, sg[cc - 1], hi, k, cc == prev, x);
		prev = cc;
		for (int j = k; j < hi && gs[j] <= y; j++) {
			int a = x > gs[j] ? x : gs[j];
			int b = y < ge[j] ? y : ge[j];
			q[o] = (int) (r + 1);
			ns[o] = a;
			w[o] = b - a + 1;
			o++;
		}
	}
	UNPROTECT(1);
	return ans;
}

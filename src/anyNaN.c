#include <R.h>
#include <Rinternals.h>

/* any_nan(x): anyNA(x) for a double vector with a data pointer, else NULL
 * (another type, or an ALTREP vector without a pointer), and the caller
 * uses anyNA. NA_real_ is a NaN, and v != v holds exactly for NaN, so an OR
 * of it per block, with no branch per element, gives the same TRUE/FALSE. */

#define NAN_BLOCK 512

SEXP any_nan(SEXP x)
{
	if (TYPEOF(x) != REALSXP)
		return R_NilValue;
	const double *p = REAL_OR_NULL(x);
	if (p == NULL)
		return R_NilValue;
	R_xlen_t n = XLENGTH(x), i = 0;
	for (; i + NAN_BLOCK <= n; i += NAN_BLOCK) {
		const double *q = p + i;
		int acc = 0;
		for (int k = 0; k < NAN_BLOCK; k++)
			acc |= (q[k] != q[k]);
		if (acc)
			return ScalarLogical(TRUE);
	}
	int acc = 0;
	for (; i < n; i++)
		acc |= (p[i] != p[i]);
	return ScalarLogical(acc);
}

#include <R.h>
#include <Rinternals.h>
#include <string.h>
#include <limits.h>

/* propensityMatch's row index into covarData for matchedData:
 * c(seq_len(nf), nf + index, nf + seq_len(np), nf + which(!drawn)), where
 * drawn marks the pool rows index holds. Returns NULL unless nf and np are
 * non-negative integer scalars with nf + np <= INT_MAX and index is an
 * integer vector of values in 1..np. */
SEXP matched_rows(SEXP s_nf, SEXP index, SEXP s_np)
{
	if (TYPEOF(s_nf) != INTSXP || XLENGTH(s_nf) != 1 ||
		TYPEOF(s_np) != INTSXP || XLENGTH(s_np) != 1 ||
		TYPEOF(index) != INTSXP)
		return R_NilValue;
	int nf = INTEGER_RO(s_nf)[0], np = INTEGER_RO(s_np)[0];
	if (nf == NA_INTEGER || np == NA_INTEGER || nf < 0 || np < 0 ||
		(double) nf + (double) np > (double) INT_MAX)
		return R_NilValue;
	R_xlen_t m = XLENGTH(index);
	const int *ix = INTEGER_RO(index);
	unsigned char *drawn = (unsigned char *) R_alloc((size_t) np + 1, 1);
	memset(drawn, 0, (size_t) np + 1);
	for (R_xlen_t k = 0; k < m; k++) {
		int v = ix[k];
		if (v == NA_INTEGER || v < 1 || v > np)
			return R_NilValue;
		drawn[v - 1] = 1;
	}
	R_xlen_t nu = 0;
	for (int j = 0; j < np; j++)
		nu += !drawn[j];
	SEXP out = PROTECT(allocVector(INTSXP, (R_xlen_t) nf + m + np + nu));
	int *o = INTEGER(out);
	R_xlen_t at = 0;
	for (int i = 1; i <= nf; i++)
		o[at++] = i;
	for (R_xlen_t k = 0; k < m; k++)
		o[at++] = nf + ix[k];
	for (int j = 1; j <= np; j++)
		o[at++] = nf + j;
	for (int j = 0; j < np; j++)
		if (!drawn[j])
			o[at++] = nf + j + 1;
	UNPROTECT(1);
	return out;
}

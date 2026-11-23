#include <R.h>
#include <Rinternals.h>
#include <string.h>

/* propensityMatch's focal and pool scores, list(ps[seq_len(nf)],
 * ps[nf + seq_len(n - nf)]), as two block copies, for an attribute-free
 * double ps of length n (the caller checks the attributes), an integer
 * scalar nf in 0..n and an integer scalar n. Returns NULL otherwise. */
SEXP split_scores(SEXP ps, SEXP s_nf, SEXP s_n)
{
	if (TYPEOF(ps) != REALSXP || isObject(ps) ||
		TYPEOF(s_nf) != INTSXP || XLENGTH(s_nf) != 1 ||
		TYPEOF(s_n) != INTSXP || XLENGTH(s_n) != 1)
		return R_NilValue;
	int nf = INTEGER_RO(s_nf)[0], n = INTEGER_RO(s_n)[0];
	if (nf == NA_INTEGER || n == NA_INTEGER || nf < 0 || nf > n ||
		XLENGTH(ps) != (R_xlen_t) n)
		return R_NilValue;
	/* a wrapper's (names(ps) <- NULL) read pointer is its payload's */
	const double *x = REAL_RO(ps);
	SEXP out = PROTECT(allocVector(VECSXP, 2));
	SEXP fps = allocVector(REALSXP, nf);
	SET_VECTOR_ELT(out, 0, fps);
	SEXP pps = allocVector(REALSXP, n - nf);
	SET_VECTOR_ELT(out, 1, pps);
	if (nf > 0)
		memcpy(REAL(fps), x, (size_t) nf * sizeof(double));
	if (n - nf > 0)
		memcpy(REAL(pps), x + nf, (size_t) (n - nf) * sizeof(double));
	UNPROTECT(1);
	return out;
}

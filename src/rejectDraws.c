#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <float.h>

/* rejectSample's acceptance draws: for each pool score, the focal density fx
 * and pool density px at it, the acceptance probability
 * pmin(ifelse(fx > 1e-3, fx, 0) / (scale * px), 1) and one
 * rbinom(1, size = 1, prob = p), as stats' rbinom makes them (Rmath's
 * rbinom per element between one GetRNGstate and one PutRNGstate, returned
 * as integer).
 *
 * Returns NULL, before touching the RNG, unless fx and px are doubles of one
 * length n >= 1, scale is one double, no fx is NA or NaN, and every
 * scale * px is finite and > 0: then no probability is NaN and no draw can
 * be NA, and R's route computes the same values. */
SEXP reject_draws(SEXP fx, SEXP px, SEXP s_scale)
{
	if (TYPEOF(fx) != REALSXP || TYPEOF(px) != REALSXP ||
		TYPEOF(s_scale) != REALSXP || XLENGTH(s_scale) != 1)
		return R_NilValue;
	R_xlen_t n = XLENGTH(fx);
	if (n < 1 || XLENGTH(px) != n)
		return R_NilValue;
	const double *f = REAL_RO(fx), *p = REAL_RO(px);
	const double scale = REAL_RO(s_scale)[0];
	for (R_xlen_t i = 0; i < n; i++) {
		double d = scale * p[i];
		if (ISNAN(f[i]) || !(d > 0 && d <= DBL_MAX))
			return R_NilValue;
	}
	SEXP out = PROTECT(allocVector(INTSXP, n));
	int *o = INTEGER(out);
	GetRNGstate();
	for (R_xlen_t i = 0; i < n; i++) {
		double t = f[i] > 1e-3 ? f[i] : 0.0;
		double q = t / (scale * p[i]);
		double pr = (1.0 < q) ? 1.0 : q;
		o[i] = (int) rbinom(1.0, pr);
	}
	PutRNGstate();
	UNPROTECT(1);
	return out;
}

#define USE_FC_LEN_T
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Applic.h>
#include <R_ext/BLAS.h>
#include <float.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#ifndef FCONE
# define FCONE
#endif

/* the elementwise loops run on at most this many threads, from this many
 * rows */
#define LOGIT_MAX_THREADS 4
#define LOGIT_PAR_MIN 16384
/* rows per unit of work in the z, w and x * w loop, whose w stays in cache */
#define LOGIT_BLOCK 2048

/* stats::glm.fit's IRLS loop for binomial("logit") with prior weights of 1
 * and no offset. Every value is computed with the operations, in the order,
 * that glm.fit's R code and the C functions it calls use, so the fitted
 * values are bitwise those of glm.fit. The helpers below copy stats'
 * family.c (logit_linkinv, logit_mu_eta, binomial_dev_resids) expression for
 * expression, so a compiler that contracts multiply-adds does so in the same
 * places. */

static const double THRESH = 30.;
static const double MTHRESH = -30.;
static const double INVEPS = 1/DBL_EPSILON;

/* mustart is 0.25 or 0.75 here, so family.c's range error cannot arise */
static R_INLINE double x_d_omx(double x) {
	return x/(1 - x);
}

static R_INLINE double x_d_opx(double x) {return x/(1 + x);}

/* linkinv1 also stores the exp(etai) it computes in *expE, and mu_eta1 reads
 * that value back in the next IRLS step instead of calling exp(etai) again.
 * mu_eta1 uses expE only when MTHRESH <= etai <= THRESH, exactly where
 * linkinv1 called exp(etai). */
static R_INLINE double linkinv1(double etai, double *expE)
{
	double tmp;
	tmp = (etai < MTHRESH) ? DBL_EPSILON :
		((etai > THRESH) ? INVEPS : exp(etai));
	*expE = tmp;
	return x_d_opx(tmp);
}

static R_INLINE double mu_eta1(double etai, double expE)
{
	double opexp = 1 + expE;
	return (etai > THRESH || etai < MTHRESH) ? DBL_EPSILON : expE/(opexp * opexp);
}

static R_INLINE double y_log_y(double y, double mu)
{
	return (y != 0.) ? (y * log(y/mu)) : 0;
}

static R_INLINE double dev_resid1(double yi, double mui, double wt)
{
	return 2 * wt *
		(y_log_y(yi, mui) + y_log_y(1 - yi, 1 - mui));
}

/* sum(): summary.c's rsum, a long double accumulator */
static R_INLINE double rsum_value(long double s)
{
	if (s > DBL_MAX) return R_PosInf;
	else if (s < -DBL_MAX) return R_NegInf;
	else return (double) s;
}

/* array.c's test for whether %*% uses the BLAS (false) or its own loop */
static int may_have_nan_or_inf(const double *x, R_xlen_t n)
{
	if ((n&1) != 0 && !R_FINITE(x[0]))
		return 1;
	for (R_xlen_t i = n&1; i < n; i += 2)
		if (!R_FINITE(x[i]+x[i+1]))
			return 1;
	return 0;
}

/* x: glm's model matrix (double, n x p); y: the response factor's codes.
 * Returns list(fitted values, converged, any fitted value within 10 eps of
 * 0 or 1), or NULL wherever glm.fit would stop, warn inside its loop, halve
 * a step or drop a row, and wherever %*% would not call dgemv.
 * The elementwise loops run on up to n_threads threads; each element is
 * computed by the same expression on any thread, and the deviance is summed
 * afterwards in one pass in row order, so the result does not depend on the
 * thread count. */
static SEXP logit_impl(SEXP x, SEXP ycode, SEXP s_eps, SEXP s_maxit,
	SEXP s_tol, int n_threads)
{
	SEXP dim = getAttrib(x, R_DimSymbol);
	if (TYPEOF(x) != REALSXP || TYPEOF(ycode) != INTSXP || LENGTH(dim) != 2)
		return R_NilValue;
	int n = INTEGER(dim)[0], p = INTEGER(dim)[1];
	R_xlen_t N = n, NP = (R_xlen_t) n * p;
	if (n < 1 || p < 1 || XLENGTH(ycode) != N)
		return R_NilValue;
	const double *X = REAL(x);
	const int *yc = INTEGER(ycode);
	double epsilon = asReal(s_eps), tol = asReal(s_tol);
	int maxit = asInteger(s_maxit);
	if (may_have_nan_or_inf(X, NP))
		return R_NilValue;

	/* the n-long buffers as glm.fit's .Call(C_Cdqrls) allocates its own */
	SEXP s_mu = PROTECT(allocVector(REALSXP, N));
	SEXP s_qr = PROTECT(allocVector(REALSXP, NP));
	SEXP s_zw = PROTECT(allocVector(REALSXP, N));
	SEXP s_rsd = PROTECT(allocVector(REALSXP, N));
	SEXP s_qty = PROTECT(allocVector(REALSXP, N));
	SEXP s_eta = PROTECT(allocVector(REALSXP, N));
	SEXP s_y = PROTECT(allocVector(REALSXP, N));
	SEXP s_w = PROTECT(allocVector(REALSXP, N));
	double *mu = REAL(s_mu), *qr = REAL(s_qr), *zw = REAL(s_zw),
		*rsd = REAL(s_rsd), *qty = REAL(s_qty), *eta = REAL(s_eta),
		*y = REAL(s_y), *w = REAL(s_w);
	/* w[i] holds linkinv1's exp(eta[i]) from the end of one step until the
	 * next step reads it for mu_eta1 and overwrites it with the weight */
	double *coef = (double *) R_alloc(p, sizeof(double));
	double *start = (double *) R_alloc(p, sizeof(double));
	double *qraux = (double *) R_alloc(p, sizeof(double));
	double *work = (double *) R_alloc(2 * (size_t) p, sizeof(double));
	int *pivot = (int *) R_alloc(p, sizeof(int));
	int conv = 0, bound = 0, bad;
	long double s;
	double dev, devold;
	/* the deviance residuals, summed after the linkinv loop; dqrls
	 * overwrites rsd before reading it */
	double *dr = rsd;

	int n_thr = 1;
#ifdef _OPENMP
	if (N >= LOGIT_PAR_MIN && n_threads > 1) {
		n_thr = omp_get_max_threads();
		if (n_thr > n_threads)
			n_thr = n_threads;
		if (n_thr > LOGIT_MAX_THREADS)
			n_thr = LOGIT_MAX_THREADS;
		if (n_thr < 1)
			n_thr = 1;
	}
#else
	(void) n_threads;
#endif

	/* binomial()$initialize for a factor: y != levels(y)[1L], then
	 * mustart <- (weights * y + 0.5)/(weights + 1); eta <- linkfun(mustart).
	 * One thread: mustart takes two values, and the compiler folds each
	 * branch to constants */
	s = 0.0;
	for (R_xlen_t i = 0; i < N; i++) {
		if (yc[i] == NA_INTEGER)
			goto fallback;
		double yi = (yc[i] != 1) ? 1.0 : 0.0;
		y[i] = yi;
		double mustart = (1.0 * yi + 0.5)/(1.0 + 1);
		double etai = log(x_d_omx(mustart));
		eta[i] = etai;
		double mui = linkinv1(etai, w + i);
		mu[i] = mui;
		if (!(mui > 0 && mui < 1))
			goto fallback;
		s += dev_resid1(yi, mui, 1.0);
	}
	devold = rsum_value(s);

	for (int iter = 1; iter <= maxit; iter++) {
		R_CheckUserInterrupt();
		/* varmu, mu.eta.val, z and w; every row is good while weights are 1
		 * and mu.eta.val is not 0. Then the row of x[good, , drop = FALSE]
		 * * w; Cdqrls stops on a non-finite value. In package code R_FINITE
		 * is a call into libR's R_finite, so that test uses C99 isfinite,
		 * which agrees on every double */
		bad = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_thr) schedule(static) if(n_thr > 1) \
	reduction(|:bad)
#endif
		for (R_xlen_t b = 0; b < N; b += LOGIT_BLOCK) {
			R_xlen_t e = (N - b > LOGIT_BLOCK) ? b + LOGIT_BLOCK : N;
			for (R_xlen_t i = b; i < e; i++) {
				double mui = mu[i], etai = eta[i];
				double varmu = mui * (1 - mui);
				if (ISNAN(varmu) || varmu == 0) {
					bad = 1;
					break;
				}
				double mev = mu_eta1(etai, w[i]);
				if (ISNAN(mev) || mev == 0) {
					bad = 1;
					break;
				}
				double z = (etai - 0.0) + (y[i] - mui)/mev;
				double wi = sqrt((1.0 * (mev * mev))/(mui * (1 - mui)));
				w[i] = wi;
				double zwi = z * wi;
				if (!R_FINITE(zwi)) {
					bad = 1;
					break;
				}
				zw[i] = zwi;
			}
			if (bad)
				continue;
			for (int j = 0; j < p; j++) {
				const double *xj = X + (R_xlen_t) j * N;
				double *qj = qr + (R_xlen_t) j * N;
				int nf = 0;
				for (R_xlen_t i = b; i < e; i++) {
					double v = xj[i] * w[i];
					nf |= !isfinite(v);
					qj[i] = v;
				}
				bad |= nf;
			}
		}
		if (bad)
			goto fallback;
		/* Cdqrls */
		int ny = 1, rank;
		for (int j = 0; j < p; j++)
			pivot[j] = j + 1;
		F77_CALL(dqrls)(qr, &n, &p, zw, &ny, &tol, coef, rsd, qty, &rank,
			pivot, qraux, work);
		for (int j = 0; j < p; j++)
			if (!R_FINITE(coef[j]))
				goto fallback;
		/* start[fit$pivot] <- fit$coefficients */
		for (int j = 0; j < p; j++)
			start[pivot[j] - 1] = coef[j];
		/* eta <- drop(x %*% start): matprod's dgemv for a vector */
		if (may_have_nan_or_inf(start, p))
			goto fallback;
		{
			const char *transN = "N";
			double one = 1.0, zero = 0.0;
			int ione = 1;
			F77_CALL(dgemv)(transN, &n, &p, &one, X, &n, start, &ione, &zero,
				eta, &ione FCONE);
		}
		/* mu <- linkinv(eta <- eta + offset); dev <- sum(dev.resids(...)) */
		bad = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_thr) schedule(static) if(n_thr > 1) \
	reduction(|:bad)
#endif
		for (R_xlen_t i = 0; i < N; i++) {
			double etai = eta[i] + 0.0;
			eta[i] = etai;
			double mui = linkinv1(etai, w + i);
			mu[i] = mui;
			dr[i] = dev_resid1(y[i], mui, 1.0);
			if (!(mui > 0 && mui < 1))
				bad = 1;
		}
		s = 0.0;
		for (R_xlen_t i = 0; i < N; i++)
			s += dr[i];
		dev = rsum_value(s);
		if (!R_FINITE(dev) || bad)
			goto fallback;
		if (fabs(dev - devold)/(0.1 + fabs(dev)) < epsilon) {
			conv = 1;
			break;
		}
		devold = dev;
	}

	/* glm.fit's "fitted probabilities numerically 0 or 1" test */
	double eps = 10 * DBL_EPSILON;
	for (R_xlen_t i = 0; i < N; i++)
		if (mu[i] > 1 - eps || mu[i] < eps) {
			bound = 1;
			break;
		}

	SEXP ans = PROTECT(allocVector(VECSXP, 3));
	SET_VECTOR_ELT(ans, 0, s_mu);
	SET_VECTOR_ELT(ans, 1, ScalarLogical(conv));
	SET_VECTOR_ELT(ans, 2, ScalarLogical(bound));
	UNPROTECT(9);
	return ans;

fallback:
	UNPROTECT(8);
	return R_NilValue;
}

/* logit_impl() on one thread */
SEXP logit_irls(SEXP x, SEXP ycode, SEXP s_eps, SEXP s_maxit, SEXP s_tol)
{
	return logit_impl(x, ycode, s_eps, s_maxit, s_tol, 1);
}

/* logit_impl() on up to n_threads_ threads (at most 4, and only with
 * OpenMP); the caller passes data.table's thread count, which is 1 in a
 * forked child, where a parallel region would wait for the parent's OpenMP
 * threads */
SEXP logit_irls_par(SEXP x, SEXP ycode, SEXP s_eps, SEXP s_maxit, SEXP s_tol,
	SEXP n_threads_)
{
	int n_threads = 1;
	if (TYPEOF(n_threads_) == INTSXP && XLENGTH(n_threads_) == 1 &&
		INTEGER(n_threads_)[0] != NA_INTEGER)
		n_threads = INTEGER(n_threads_)[0];
	return logit_impl(x, ycode, s_eps, s_maxit, s_tol, n_threads);
}

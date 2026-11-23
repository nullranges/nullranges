#include <R.h>
#include <Rinternals.h>

/* KernSmooth's linbin(sx, gpoints, truncate = TRUE) on sx = (x - m) / s,
 * without storing sx; the same double operations in the same order as its
 * src/linbin.f, whose notice follows:
 *
 *   Part of R package KernSmooth
 *   Copyright (C) 1995  M. P. Wand
 *
 *   Unlimited use and distribution (see LICENCE).
 *
 * KernSmooth's LICENCE.note: "By agreement with Dr Wand (1998-June-22), the
 * R port can be used and distributed freely".
 *
 * The range test on lxi before the cast replaces Fortran's test on
 * li = int(lxi): for finite lxi they agree, and a NaN or infinite lxi,
 * whose cast is undefined in C, is skipped as linbin skips it on x86. */
SEXP scaled_linbin(SEXP x, SEXP s_m, SEXP s_s, SEXP s_a, SEXP s_b, SEXP s_M)
{
	R_xlen_t n = XLENGTH(x);
	const double *xp = REAL(x);
	double m = asReal(s_m), s = asReal(s_s), a = asReal(s_a), b = asReal(s_b);
	int M = asInteger(s_M);
	SEXP out = PROTECT(allocVector(REALSXP, M));
	double *g = REAL(out);
	for (int i = 0; i < M; i++)
		g[i] = 0.0;
	double delta = (b - a) / (M - 1);
	for (R_xlen_t i = 0; i < n; i++) {
		double sx = (xp[i] - m) / s;
		double lxi = ((sx - a) / delta) + 1;
		if (!(lxi >= 1.0 && lxi < (double) M))
			continue;
		int li = (int) lxi;
		double rem = lxi - li;
		g[li - 1] = g[li - 1] + (1 - rem);
		g[li] = g[li] + rem;
	}
	UNPROTECT(1);
	return out;
}

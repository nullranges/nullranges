#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP moved_positions(SEXP b_code, SEXP b_start, SEXP b_width, SEXP n_rep,
	SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width, SEXP circular,
	SEXP seqlen, SEXP bounded, SEXP shift, SEXP t_code);
SEXP moved_rows(SEXP b_code, SEXP b_start, SEXP b_width, SEXP n_rep,
	SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width, SEXP circular,
	SEXP seqlen, SEXP bounded, SEXP shift, SEXP t_code, SEXP mode,
	SEXP xseg, SEXP xs, SEXP xe);
SEXP moved_rows_lens(SEXP b_code, SEXP b_start, SEXP b_width, SEXP n_rep,
	SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width, SEXP circular,
	SEXP seqlen, SEXP bounded, SEXP shift, SEXP t_code, SEXP mode,
	SEXP xseg, SEXP xs, SEXP xe, SEXP cols);
SEXP ss_init(SEXP fps, SEXP pps, SEXP cnt, SEXP code, SEXP nLeft);
SEXP ss_range(SEXP st);
SEXP ss_sort(SEXP st);
SEXP ss_compact(SEXP st);
SEXP ss_pass(SEXP st, SEXP br, SEXP mn, SEXP by, SEXP replace, SEXP count);
SEXP ss_pass_by(SEXP st, SEXP mn, SEXP mx, SEXP by, SEXP replace,
	SEXP count);
SEXP outside_intervals(SEXP code, SEXP start, SEXP width, SEXP seg,
	SEXP i_start, SEXP i_end);
SEXP logit_irls(SEXP x, SEXP ycode, SEXP s_eps, SEXP s_maxit, SEXP s_tol);
SEXP logit_irls_par(SEXP x, SEXP ycode, SEXP s_eps, SEXP s_maxit, SEXP s_tol,
	SEXP n_threads);
SEXP sample_groups(SEXP index, SEXP n_rows, SEXP size);
SEXP intersect_gaps(SEXP code, SEXP start, SEXP end, SEXP seg,
	SEXP g_start, SEXP g_end);
SEXP seg_draws(SEXP s_n_rep, SEXP prob, SEXP rmin, SEXP rmax, SEXP gidx,
	SEXP off, SEXP n_draw, SEXP st, SEXP rearr);

static const R_CallMethodDef callMethods[] = {
	{"moved_positions", (DL_FUNC) &moved_positions, 13},
	{"moved_rows", (DL_FUNC) &moved_rows, 17},
	{"moved_rows_lens", (DL_FUNC) &moved_rows_lens, 18},
	{"ss_init", (DL_FUNC) &ss_init, 5},
	{"ss_range", (DL_FUNC) &ss_range, 1},
	{"ss_sort", (DL_FUNC) &ss_sort, 1},
	{"ss_compact", (DL_FUNC) &ss_compact, 1},
	{"ss_pass", (DL_FUNC) &ss_pass, 6},
	{"ss_pass_by", (DL_FUNC) &ss_pass_by, 6},
	{"outside_intervals", (DL_FUNC) &outside_intervals, 6},
	{"logit_irls", (DL_FUNC) &logit_irls, 5},
	{"logit_irls_par", (DL_FUNC) &logit_irls_par, 6},
	{"sample_groups", (DL_FUNC) &sample_groups, 3},
	{"intersect_gaps", (DL_FUNC) &intersect_gaps, 6},
	{"seg_draws", (DL_FUNC) &seg_draws, 9},
	{NULL, NULL, 0}
};

void R_init_nullranges(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

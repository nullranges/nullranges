#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP moved_positions(SEXP b_code, SEXP b_start, SEXP b_width, SEXP n_rep,
	SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width, SEXP circular,
	SEXP seqlen, SEXP bounded, SEXP shift, SEXP t_code);
SEXP ss_init(SEXP fps, SEXP pps, SEXP cnt, SEXP code, SEXP nLeft);
SEXP ss_range(SEXP st);
SEXP ss_pass(SEXP st, SEXP br, SEXP mn, SEXP by, SEXP replace, SEXP count);
SEXP outside_intervals(SEXP code, SEXP start, SEXP width, SEXP seg,
	SEXP i_start, SEXP i_end);

static const R_CallMethodDef callMethods[] = {
	{"moved_positions", (DL_FUNC) &moved_positions, 13},
	{"ss_init", (DL_FUNC) &ss_init, 5},
	{"ss_range", (DL_FUNC) &ss_range, 1},
	{"ss_pass", (DL_FUNC) &ss_pass, 6},
	{"outside_intervals", (DL_FUNC) &outside_intervals, 6},
	{NULL, NULL, 0}
};

void R_init_nullranges(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

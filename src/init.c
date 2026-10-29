#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP moved_positions(SEXP b_code, SEXP b_start, SEXP b_width, SEXP n_rep,
	SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width, SEXP circular,
	SEXP seqlen, SEXP bounded, SEXP shift, SEXP t_code);

static const R_CallMethodDef callMethods[] = {
	{"moved_positions", (DL_FUNC) &moved_positions, 13},
	{NULL, NULL, 0}
};

void R_init_nullranges(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[] = {
	{NULL, NULL, 0}
};

void R_init_nullranges(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

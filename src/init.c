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
SEXP moved_rows_opt(SEXP b_code, SEXP b_start, SEXP b_width, SEXP n_rep,
	SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width, SEXP circular,
	SEXP seqlen, SEXP bounded, SEXP shift, SEXP t_code, SEXP mode,
	SEXP xseg, SEXP xs, SEXP xe, SEXP cols, SEXP with_code);
SEXP moved_rows_opt_par(SEXP b_code, SEXP b_start, SEXP b_width, SEXP n_rep,
	SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width, SEXP circular,
	SEXP seqlen, SEXP bounded, SEXP shift, SEXP t_code, SEXP mode,
	SEXP xseg, SEXP xs, SEXP xe, SEXP cols, SEXP with_code, SEXP n_threads);
SEXP bait_blocks(SEXP start, SEXP width, SEXP target, SEXP n_rep);
SEXP moved_rows_bait(SEXP b_seg, SEXP seg_code, SEXP b_start, SEXP b_width,
	SEXP n_rep, SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width,
	SEXP circular, SEXP seqlen, SEXP bounded, SEXP t_start, SEXP t_code,
	SEXP mode, SEXP xseg, SEXP xs, SEXP xe, SEXP cols, SEXP with_code);
SEXP moved_rows_bait_par(SEXP b_seg, SEXP seg_code, SEXP b_start,
	SEXP b_width, SEXP n_rep, SEXP y_code, SEXP y_start, SEXP y_end,
	SEXP y_width, SEXP circular, SEXP seqlen, SEXP bounded, SEXP t_start,
	SEXP t_code, SEXP mode, SEXP xseg, SEXP xs, SEXP xe, SEXP cols,
	SEXP with_code, SEXP n_threads);
SEXP outside_intervals(SEXP code, SEXP start, SEXP width, SEXP seg,
	SEXP i_start, SEXP i_end);
SEXP intersect_gaps(SEXP code, SEXP start, SEXP end, SEXP seg,
	SEXP g_start, SEXP g_end);
SEXP seg_draws(SEXP s_n_rep, SEXP prob, SEXP rmin, SEXP rmax, SEXP gidx,
	SEXP off, SEXP n_draw, SEXP st, SEXP rearr);
SEXP seg_draws_moved(SEXP s_n_rep, SEXP prob, SEXP rmin, SEXP rmax,
	SEXP gidx, SEXP off, SEXP n_draw, SEXP st, SEXP rearr, SEXP seg_code,
	SEXP b_width, SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width,
	SEXP circular, SEXP seqlen, SEXP bounded, SEXP t_start, SEXP t_code,
	SEXP mode, SEXP xseg, SEXP xs, SEXP xe, SEXP cols, SEXP with_code);
SEXP seg_draws_moved_par(SEXP s_n_rep, SEXP prob, SEXP rmin, SEXP rmax,
	SEXP gidx, SEXP off, SEXP n_draw, SEXP st, SEXP rearr, SEXP seg_code,
	SEXP b_width, SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width,
	SEXP circular, SEXP seqlen, SEXP bounded, SEXP t_start, SEXP t_code,
	SEXP mode, SEXP xseg, SEXP xs, SEXP xe, SEXP cols, SEXP with_code,
	SEXP n_threads);

static const R_CallMethodDef callMethods[] = {
	{"moved_positions", (DL_FUNC) &moved_positions, 13},
	{"moved_rows", (DL_FUNC) &moved_rows, 17},
	{"moved_rows_lens", (DL_FUNC) &moved_rows_lens, 18},
	{"moved_rows_opt", (DL_FUNC) &moved_rows_opt, 19},
	{"moved_rows_opt_par", (DL_FUNC) &moved_rows_opt_par, 20},
	{"bait_blocks", (DL_FUNC) &bait_blocks, 4},
	{"moved_rows_bait", (DL_FUNC) &moved_rows_bait, 20},
	{"moved_rows_bait_par", (DL_FUNC) &moved_rows_bait_par, 21},
	{"outside_intervals", (DL_FUNC) &outside_intervals, 6},
	{"intersect_gaps", (DL_FUNC) &intersect_gaps, 6},
	{"seg_draws", (DL_FUNC) &seg_draws, 9},
	{"seg_draws_moved", (DL_FUNC) &seg_draws_moved, 26},
	{"seg_draws_moved_par", (DL_FUNC) &seg_draws_moved_par, 27},
	{NULL, NULL, 0}
};

void R_init_nullranges(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

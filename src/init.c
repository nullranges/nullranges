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
SEXP bait_blocks(SEXP start, SEXP width, SEXP target, SEXP n_rep);
SEXP moved_rows_bait(SEXP b_seg, SEXP seg_code, SEXP b_start, SEXP b_width,
	SEXP n_rep, SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width,
	SEXP circular, SEXP seqlen, SEXP bounded, SEXP t_start, SEXP t_code,
	SEXP mode, SEXP xseg, SEXP xs, SEXP xe, SEXP cols, SEXP with_code);
SEXP ss_init(SEXP fps, SEXP pps, SEXP cnt, SEXP code, SEXP nLeft);
SEXP ss_range(SEXP st);
SEXP ss_sort(SEXP st);
SEXP ss_init_sorted(SEXP fps, SEXP pps);
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
SEXP seg_draws_moved(SEXP s_n_rep, SEXP prob, SEXP rmin, SEXP rmax,
	SEXP gidx, SEXP off, SEXP n_draw, SEXP st, SEXP rearr, SEXP seg_code,
	SEXP b_width, SEXP y_code, SEXP y_start, SEXP y_end, SEXP y_width,
	SEXP circular, SEXP seqlen, SEXP bounded, SEXP t_start, SEXP t_code,
	SEXP mode, SEXP xseg, SEXP xs, SEXP xe, SEXP cols, SEXP with_code);
SEXP scaled_linbin(SEXP x, SEXP s_m, SEXP s_s, SEXP s_a, SEXP s_b, SEXP s_M);
SEXP order_stats(SEXP x, SEXP ranks);
SEXP order_stats_range(SEXP x, SEXP ranks);
SEXP reject_draws(SEXP fx, SEXP px, SEXP s_scale);
SEXP matched_rows(SEXP s_nf, SEXP index, SEXP s_np);
SEXP any_nan(SEXP x);
SEXP accept_index(SEXP x);
SEXP nn_keep(SEXP fps, SEXP pps);
SEXP nn_draws(SEXP fps, SEXP pps);
SEXP split_scores(SEXP ps, SEXP s_nf, SEXP s_n);

static const R_CallMethodDef callMethods[] = {
	{"moved_positions", (DL_FUNC) &moved_positions, 13},
	{"moved_rows", (DL_FUNC) &moved_rows, 17},
	{"moved_rows_lens", (DL_FUNC) &moved_rows_lens, 18},
	{"moved_rows_opt", (DL_FUNC) &moved_rows_opt, 19},
	{"bait_blocks", (DL_FUNC) &bait_blocks, 4},
	{"moved_rows_bait", (DL_FUNC) &moved_rows_bait, 20},
	{"ss_init", (DL_FUNC) &ss_init, 5},
	{"ss_range", (DL_FUNC) &ss_range, 1},
	{"ss_sort", (DL_FUNC) &ss_sort, 1},
	{"ss_init_sorted", (DL_FUNC) &ss_init_sorted, 2},
	{"ss_compact", (DL_FUNC) &ss_compact, 1},
	{"ss_pass", (DL_FUNC) &ss_pass, 6},
	{"ss_pass_by", (DL_FUNC) &ss_pass_by, 6},
	{"outside_intervals", (DL_FUNC) &outside_intervals, 6},
	{"logit_irls", (DL_FUNC) &logit_irls, 5},
	{"logit_irls_par", (DL_FUNC) &logit_irls_par, 6},
	{"sample_groups", (DL_FUNC) &sample_groups, 3},
	{"intersect_gaps", (DL_FUNC) &intersect_gaps, 6},
	{"seg_draws", (DL_FUNC) &seg_draws, 9},
	{"seg_draws_moved", (DL_FUNC) &seg_draws_moved, 26},
	{"scaled_linbin", (DL_FUNC) &scaled_linbin, 6},
	{"order_stats", (DL_FUNC) &order_stats, 2},
	{"order_stats_range", (DL_FUNC) &order_stats_range, 2},
	{"reject_draws", (DL_FUNC) &reject_draws, 3},
	{"matched_rows", (DL_FUNC) &matched_rows, 3},
	{"any_nan", (DL_FUNC) &any_nan, 1},
	{"accept_index", (DL_FUNC) &accept_index, 1},
	{"nn_keep", (DL_FUNC) &nn_keep, 2},
	{"nn_draws", (DL_FUNC) &nn_draws, 2},
	{"split_scores", (DL_FUNC) &split_scores, 3},
	{NULL, NULL, 0}
};

void R_init_nullranges(DllInfo *dll)
{
	R_registerRoutines(dll, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(dll, FALSE);
	R_forceSymbols(dll, TRUE);
}

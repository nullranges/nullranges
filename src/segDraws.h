#ifndef NULLRANGES_SEGDRAWS_H
#define NULLRANGES_SEGDRAWS_H

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

/* a table for sample(m, k, replace = TRUE, prob = w) */
typedef struct {
	int n;
	int walker;
	double *p; /* cumulative p in revsort order, or Walker's q */
	int *a;    /* the permutation, or Walker's aliases */
} prob_table;

/* seg_draws()' tables and buffers, and where its draws go: os, og and ol,
 * the starts, segment indices and counts per replicate, w the draws
 * written */
typedef struct {
	int n_rep, n_tab, noprop, n_state;
	R_xlen_t per_rep, total, w;
	const int *o, *nd, *stc, *re, *gi;
	const double *lo, *hi;
	prob_table *tab;
	int *idx, *sbuf, *cnt, *pos, *ord;
	double *ubuf;
	int *os, *og, *ol;
	Sampletype kind;
} seg_plan;

int seg_plan_init(seg_plan *p, SEXP s_n_rep, SEXP prob, SEXP rmin,
	SEXP rmax, SEXP gidx, SEXP off, SEXP n_draw, SEXP st, SEXP rearr);
SEXP seg_plan_out(seg_plan *p, int n_el);
void seg_plan_rep(seg_plan *p, int r);
void seg_plan_trim(const seg_plan *p, SEXP out);

#endif

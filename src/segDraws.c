#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <string.h>
#include "segDraws.h"

/* sample(m, k, replace = TRUE, prob = w) as R 4.6.1's do_sample makes it
 * (src/main/random.c): FixupProb, then walker_ProbSampleReplace when more
 * than 200 of m * p exceed 0.1, else ProbSampleReplace. The tables depend
 * on w alone, so one is built per w and reused by every call. The caller
 * has checked that w is finite, non-negative and has a positive entry. */
static void table_init(prob_table *t, const double *w, int n)
{
	double *p = (double *) R_alloc(n, sizeof(double));
	int *a = (int *) R_alloc(n, sizeof(int));
	double sum = 0.0;
	for (int i = 0; i < n; i++) {
		p[i] = w[i];
		if (p[i] > 0.0)
			sum += p[i];
	}
	for (int i = 0; i < n; i++)
		p[i] /= sum;
	int nc = 0;
	for (int i = 0; i < n; i++)
		if (n * p[i] > 0.1)
			nc++;
	t->n = n;
	t->walker = nc > 200;
	if (t->walker) {
		/* R leaves an alias it never assigns uninitialized; such an
		 * entry is read only when rounding left its q[i] below 1 */
		memset(a, 0, n * sizeof(int));
		int *HL = (int *) R_alloc(n, sizeof(int));
		double *q = (double *) R_alloc(n, sizeof(double));
		int *H = HL, *L = HL + n;
		for (int i = 0; i < n; i++) {
			q[i] = p[i] * n;
			if (q[i] < 1.) *H++ = i; else *--L = i;
		}
		if (H > HL && L < HL + n) {
			for (int k = 0; k < n - 1; k++) {
				int i = HL[k];
				int j = *L;
				a[i] = j;
				q[j] += q[i] - 1.;
				if (q[j] < 1.) L++;
				if (L >= HL + n) break;
			}
		}
		for (int i = 0; i < n; i++)
			q[i] += i;
		t->p = q;
	} else {
		for (int i = 0; i < n; i++)
			a[i] = i + 1;
		revsort(p, a, n);
		for (int i = 1; i < n; i++)
			p[i] += p[i - 1];
		t->p = p;
	}
	t->a = a;
}

/* one draw from a Walker table, 1-based */
static int walker_draw(const prob_table *t, Sampletype kind)
{
	int n = t->n, k;
	double rU;
	if (kind == ROUNDING) {
		rU = unif_rand() * n;
		k = (int) rU;
	} else {
		k = (int) R_unif_index(n);
		rU = k + unif_rand();
	}
	return (rU < t->p[k]) ? k + 1 : t->a[k] + 1;
}

/* R scans for the first j < n - 1 with rU <= p[j] (else n - 1); p is
 * non-decreasing, so that j is the count of p[0..len - 1] below rU, len =
 * n - 1, found by a lower bound written without branches (the comparisons
 * are unpredictable) */
static inline int cum_index(const double *p, int len, double rU)
{
	const double *base = p;
	if (len > 0) {
		while (len > 1) {
			int half = len / 2;
			base = (base[half] < rU) ? base + half : base;
			len -= half;
		}
		base += (*base < rU);
	}
	return (int) (base - p);
}

#define STEP(b, x) b = ((b)[half] < (x)) ? (b) + half : (b)
#define LAST(b, x) b += (*(b) < (x))

/* k draws from a cumulative (non-Walker) table into idx, each off + its
 * 1-based entry. Each draw takes one unif_rand() and no other, so all k are
 * drawn first, in order, into u; then the searches run 8 at a time, each
 * step updating 8 independent lower bounds, whose step sizes depend on n
 * alone, so their chains of loads overlap. */
static void cum_draws(const prob_table *t, int k, double *u, int *idx, int off)
{
	const double *p = t->p;
	const int *a = t->a;
	int len1 = t->n - 1, i = 0;
	for (int j = 0; j < k; j++)
		u[j] = unif_rand();
	if (len1 > 0) {
		for (; i + 8 <= k; i += 8) {
			const double *b0 = p, *b1 = p, *b2 = p, *b3 = p,
				*b4 = p, *b5 = p, *b6 = p, *b7 = p;
			double x0 = u[i], x1 = u[i + 1], x2 = u[i + 2], x3 = u[i + 3],
				x4 = u[i + 4], x5 = u[i + 5], x6 = u[i + 6], x7 = u[i + 7];
			int len = len1;
			while (len > 1) {
				int half = len / 2;
				STEP(b0, x0); STEP(b1, x1); STEP(b2, x2); STEP(b3, x3);
				STEP(b4, x4); STEP(b5, x5); STEP(b6, x6); STEP(b7, x7);
				len -= half;
			}
			LAST(b0, x0); LAST(b1, x1); LAST(b2, x2); LAST(b3, x3);
			LAST(b4, x4); LAST(b5, x5); LAST(b6, x6); LAST(b7, x7);
			idx[i] = off + a[b0 - p];
			idx[i + 1] = off + a[b1 - p];
			idx[i + 2] = off + a[b2 - p];
			idx[i + 3] = off + a[b3 - p];
			idx[i + 4] = off + a[b4 - p];
			idx[i + 5] = off + a[b5 - p];
			idx[i + 6] = off + a[b6 - p];
			idx[i + 7] = off + a[b7 - p];
		}
	}
	for (; i < k; i++)
		idx[i] = off + a[cum_index(p, len1, u[i])];
}

#undef STEP
#undef LAST

/* seg_draws()' checks of its arguments and its tables and buffers, in p;
 * 0, before any draw, for arguments out of shape */
int seg_plan_init(seg_plan *p, SEXP s_n_rep, SEXP prob, SEXP rmin,
	SEXP rmax, SEXP gidx, SEXP off, SEXP n_draw, SEXP st, SEXP rearr)
{
	if (!isReal(prob) || !isReal(rmin) || !isReal(rmax) ||
	    !isInteger(gidx) || !isInteger(off) || !isInteger(n_draw) ||
	    !isInteger(s_n_rep) || XLENGTH(s_n_rep) != 1 || XLENGTH(off) < 2)
		return 0;
	int n_rep = INTEGER(s_n_rep)[0];
	int n_tab = LENGTH(off) - 1;
	R_xlen_t m_all = XLENGTH(prob);
	const int *o = INTEGER(off), *nd = INTEGER(n_draw);
	int noprop = !isNull(st);
	if (n_rep < 1 || XLENGTH(rmin) != m_all || XLENGTH(rmax) != m_all ||
	    XLENGTH(gidx) != m_all || XLENGTH(n_draw) != n_tab ||
	    o[0] != 0 || o[n_tab] != m_all || (noprop && n_tab != 1))
		return 0;
	int max_draw = 0;
	R_xlen_t per_rep = 0;
	for (int s = 0; s < n_tab; s++) {
		if (o[s + 1] - o[s] < 1 || nd[s] < 1)
			return 0;
		if (nd[s] > max_draw)
			max_draw = nd[s];
		per_rep += nd[s];
	}
	int n_state = 0;
	const int *stc = NULL, *re = NULL;
	if (noprop) {
		if (!isInteger(st) || XLENGTH(st) != m_all || !isInteger(rearr) ||
		    XLENGTH(rearr) < 1)
			return 0;
		n_state = LENGTH(rearr);
		stc = INTEGER(st);
		re = INTEGER(rearr);
		for (R_xlen_t j = 0; j < m_all; j++)
			if (stc[j] < 0 || stc[j] > n_state)
				return 0;
		per_rep = 0;
		for (int s = 0; s < n_state; s++) {
			if (re[s] < 0)
				return 0;
			per_rep += re[s];
		}
	}
	if (per_rep > INT_MAX / n_rep)
		return 0;

	p->n_rep = n_rep;
	p->n_tab = n_tab;
	p->noprop = noprop;
	p->n_state = n_state;
	p->per_rep = per_rep;
	p->total = per_rep * n_rep;
	p->w = 0;
	p->o = o;
	p->nd = nd;
	p->stc = stc;
	p->re = re;
	p->tab = (prob_table *) R_alloc(n_tab, sizeof(prob_table));
	const double *pr = REAL(prob);
	p->lo = REAL(rmin);
	p->hi = REAL(rmax);
	p->gi = INTEGER(gidx);
	for (int s = 0; s < n_tab; s++)
		table_init(&p->tab[s], pr + o[s], o[s + 1] - o[s]);
	p->idx = (int *) R_alloc(max_draw, sizeof(int));
	p->ubuf = (double *) R_alloc(max_draw, sizeof(double));
	p->sbuf = p->cnt = p->pos = p->ord = NULL;
	if (noprop) {
		p->sbuf = (int *) R_alloc(max_draw, sizeof(int));
		p->ord = (int *) R_alloc(max_draw, sizeof(int));
		p->cnt = (int *) R_alloc(n_state + 1, sizeof(int));
		p->pos = (int *) R_alloc(n_state + 2, sizeof(int));
	}
	p->os = p->og = p->ol = NULL;
	return 1;
}

/* The value, a list of n_el (at least 3) whose first three elements take
 * the draws: start, seg and len; the caller protects it */
SEXP seg_plan_out(seg_plan *p, int n_el)
{
	SEXP out = PROTECT(allocVector(VECSXP, n_el));
	SEXP start = allocVector(INTSXP, p->total);
	SET_VECTOR_ELT(out, 0, start);
	SEXP seg = allocVector(INTSXP, p->total);
	SET_VECTOR_ELT(out, 1, seg);
	SEXP len = allocVector(INTSXP, p->n_rep);
	SET_VECTOR_ELT(out, 2, len);
	p->os = INTEGER(start);
	p->og = INTEGER(seg);
	p->ol = INTEGER(len);
	UNPROTECT(1);
	return out;
}

/* Replicate r's draws, after those of replicates 0 to r - 1, between
 * GetRNGstate() and PutRNGstate(), p->kind set to R_sample_kind() */
void seg_plan_rep(seg_plan *p, int r)
{
	const int *o = p->o, *nd = p->nd, *stc = p->stc, *re = p->re;
	const int *gi = p->gi;
	const double *lo = p->lo, *hi = p->hi;
	prob_table *tab = p->tab;
	int *idx = p->idx, *sbuf = p->sbuf, *cnt = p->cnt, *pos = p->pos;
	int *ord = p->ord, *os = p->os, *og = p->og;
	double *ubuf = p->ubuf;
	int n_tab = p->n_tab, n_state = p->n_state;
	Sampletype kind = p->kind;
	R_xlen_t w = p->w, w0 = w;
	if (!p->noprop) {
		for (int s = 0; s < n_tab; s++) {
			int k = nd[s], base = o[s];
			if (tab[s].walker) {
				for (int t = 0; t < k; t++)
					idx[t] = base + walker_draw(&tab[s], kind) - 1;
			} else
				cum_draws(&tab[s], k, ubuf, idx, base - 1);
			for (int t = 0; t < k; t++) {
				int j = idx[t];
				os[w] = (int) fround(runif(lo[j], hi[j]), 0.);
				og[w] = gi[j];
				w++;
			}
		}
	} else {
		int k = nd[0];
		if (tab[0].walker) {
			for (int t = 0; t < k; t++)
				idx[t] = walker_draw(&tab[0], kind) - 1;
		} else
			cum_draws(&tab[0], k, ubuf, idx, -1);
		for (int t = 0; t < k; t++) {
			int j = idx[t];
			sbuf[t] = (int) fround(runif(lo[j], hi[j]), 0.);
		}
		/* the draws of each state, in draw order */
		memset(cnt, 0, (n_state + 1) * sizeof(int));
		for (int t = 0; t < k; t++)
			cnt[stc[idx[t]]]++;
		pos[0] = 0;
		for (int c = 0; c <= n_state; c++)
			pos[c + 1] = pos[c] + cnt[c];
		for (int c = 0; c <= n_state; c++)
			cnt[c] = pos[c];
		for (int t = 0; t < k; t++)
			ord[cnt[stc[idx[t]]]++] = t;
		for (int s = 1; s <= n_state; s++) {
			const int *lst = ord + pos[s];
			int rand_n = pos[s + 1] - pos[s], rearr_n = re[s - 1];
			if (rand_n < rearr_n) {
				double dn = rand_n > 0 ? (double) rand_n : 1.0;
				for (int u = 0; u < rearr_n; u++) {
					int i = (int) R_unif_index(dn);
					if (rand_n > 0) {
						int t = lst[i];
						os[w] = sbuf[t];
						og[w] = gi[idx[t]];
						w++;
					}
				}
			} else {
				for (int u = 0; u < rearr_n; u++) {
					int t = lst[u];
					os[w] = sbuf[t];
					og[w] = gi[idx[t]];
					w++;
				}
			}
		}
	}
	p->ol[r] = (int) (w - w0);
	p->w = w;
}

/* The start and seg vectors of out cut to the draws made, after
 * PutRNGstate() */
void seg_plan_trim(const seg_plan *p, SEXP out)
{
	if (p->w < p->total) {
		SET_VECTOR_ELT(out, 0, xlengthgets(VECTOR_ELT(out, 0), p->w));
		SET_VECTOR_ELT(out, 1, xlengthgets(VECTOR_ELT(out, 1), p->w));
	}
}

/* Every replicate's bait blocks for seg_bootstrap_prop() (st = NULL) or
 * seg_bootstrap_noprop(), drawn in the order and with the arguments the R
 * loop uses, between one GetRNGstate and one PutRNGstate. Segment tables:
 * prob, rmin, rmax and gidx (the segment's index in seg), table s holding
 * entries off[s] to off[s + 1] - 1; n_draw[s] draws from table s. Per
 * replicate and table: sample(), then round(runif(n, rmin, rmax)) over the
 * drawn segments. noprop has one table; st is each segment's state (0 for
 * none) and rearr the number of tiling blocks per state, which the drawn
 * blocks of each state are cut or up-sampled to (sample(rand_n, rearr_n,
 * TRUE), with rand_n = 0 drawing from one element and keeping nothing).
 * The caller has checked every argument a draw could fail or warn on and
 * that each start lies in the integer range. Returns list(start, seg, len):
 * the starts and segment indices of all replicates, and the number per
 * replicate; NULL, before any draw, for arguments out of shape. */
SEXP seg_draws(SEXP s_n_rep, SEXP prob, SEXP rmin, SEXP rmax, SEXP gidx,
	SEXP off, SEXP n_draw, SEXP st, SEXP rearr)
{
	seg_plan p;
	if (!seg_plan_init(&p, s_n_rep, prob, rmin, rmax, gidx, off, n_draw, st,
		rearr))
		return R_NilValue;
	SEXP out = PROTECT(seg_plan_out(&p, 3));
	GetRNGstate();
	p.kind = R_sample_kind();
	for (int r = 0; r < p.n_rep; r++)
		seg_plan_rep(&p, r);
	PutRNGstate();
	seg_plan_trim(&p, out);
	UNPROTECT(1);
	return out;
}

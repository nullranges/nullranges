#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Random.h>
#include <R_ext/Utils.h>
#include <limits.h>
#include <string.h>

/* sample(m, k, replace = TRUE, prob = w) as R 4.6.1's do_sample makes it
 * (src/main/random.c): FixupProb, then walker_ProbSampleReplace when more
 * than 200 of m * p exceed 0.1, else ProbSampleReplace. The tables depend
 * on w alone, so one is built per w and reused by every call. The caller
 * has checked that w is finite, non-negative and has a positive entry. */
typedef struct {
	int n;
	int walker;
	double *p; /* cumulative p in revsort order, or Walker's q */
	int *a;    /* the permutation, or Walker's aliases */
} prob_table;

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

/* one draw, 1-based */
static int table_draw(const prob_table *t, Sampletype kind)
{
	if (t->walker) {
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
	double rU = unif_rand();
	/* R scans for the first j < n - 1 with rU <= p[j] (else n - 1); p is
	 * non-decreasing, so that j is the count of p[0..n - 2] below rU, found
	 * by a lower bound written without branches (the comparisons are
	 * unpredictable) */
	const double *base = t->p;
	int len = t->n - 1;
	if (len > 0) {
		while (len > 1) {
			int half = len / 2;
			base = (base[half] < rU) ? base + half : base;
			len -= half;
		}
		base += (*base < rU);
	}
	return t->a[base - t->p];
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
	if (!isReal(prob) || !isReal(rmin) || !isReal(rmax) ||
	    !isInteger(gidx) || !isInteger(off) || !isInteger(n_draw) ||
	    !isInteger(s_n_rep) || XLENGTH(s_n_rep) != 1 || XLENGTH(off) < 2)
		return R_NilValue;
	int n_rep = INTEGER(s_n_rep)[0];
	int n_tab = LENGTH(off) - 1;
	R_xlen_t m_all = XLENGTH(prob);
	const int *o = INTEGER(off), *nd = INTEGER(n_draw);
	int noprop = !isNull(st);
	if (n_rep < 1 || XLENGTH(rmin) != m_all || XLENGTH(rmax) != m_all ||
	    XLENGTH(gidx) != m_all || XLENGTH(n_draw) != n_tab ||
	    o[0] != 0 || o[n_tab] != m_all || (noprop && n_tab != 1))
		return R_NilValue;
	int max_draw = 0;
	R_xlen_t per_rep = 0;
	for (int s = 0; s < n_tab; s++) {
		if (o[s + 1] - o[s] < 1 || nd[s] < 1)
			return R_NilValue;
		if (nd[s] > max_draw)
			max_draw = nd[s];
		per_rep += nd[s];
	}
	int n_state = 0;
	const int *stc = NULL, *re = NULL;
	if (noprop) {
		if (!isInteger(st) || XLENGTH(st) != m_all || !isInteger(rearr) ||
		    XLENGTH(rearr) < 1)
			return R_NilValue;
		n_state = LENGTH(rearr);
		stc = INTEGER(st);
		re = INTEGER(rearr);
		for (R_xlen_t j = 0; j < m_all; j++)
			if (stc[j] < 0 || stc[j] > n_state)
				return R_NilValue;
		per_rep = 0;
		for (int s = 0; s < n_state; s++) {
			if (re[s] < 0)
				return R_NilValue;
			per_rep += re[s];
		}
	}
	if (per_rep > INT_MAX / n_rep)
		return R_NilValue;
	R_xlen_t total = per_rep * n_rep;

	prob_table *tab = (prob_table *) R_alloc(n_tab, sizeof(prob_table));
	const double *pr = REAL(prob), *lo = REAL(rmin), *hi = REAL(rmax);
	const int *gi = INTEGER(gidx);
	for (int s = 0; s < n_tab; s++)
		table_init(&tab[s], pr + o[s], o[s + 1] - o[s]);
	int *idx = (int *) R_alloc(max_draw, sizeof(int));
	int *sbuf = NULL, *cnt = NULL, *pos = NULL, *ord = NULL;
	if (noprop) {
		sbuf = (int *) R_alloc(max_draw, sizeof(int));
		ord = (int *) R_alloc(max_draw, sizeof(int));
		cnt = (int *) R_alloc(n_state + 1, sizeof(int));
		pos = (int *) R_alloc(n_state + 2, sizeof(int));
	}
	SEXP out = PROTECT(allocVector(VECSXP, 3));
	SEXP start = allocVector(INTSXP, total);
	SET_VECTOR_ELT(out, 0, start);
	SEXP seg = allocVector(INTSXP, total);
	SET_VECTOR_ELT(out, 1, seg);
	SEXP len = allocVector(INTSXP, n_rep);
	SET_VECTOR_ELT(out, 2, len);
	int *os = INTEGER(start), *og = INTEGER(seg), *ol = INTEGER(len);

	R_xlen_t w = 0;
	GetRNGstate();
	Sampletype kind = R_sample_kind();
	for (int r = 0; r < n_rep; r++) {
		R_xlen_t w0 = w;
		if (!noprop) {
			for (int s = 0; s < n_tab; s++) {
				int k = nd[s], base = o[s];
				for (int t = 0; t < k; t++)
					idx[t] = base + table_draw(&tab[s], kind) - 1;
				for (int t = 0; t < k; t++) {
					int j = idx[t];
					os[w] = (int) fround(runif(lo[j], hi[j]), 0.);
					og[w] = gi[j];
					w++;
				}
			}
		} else {
			int k = nd[0];
			for (int t = 0; t < k; t++)
				idx[t] = table_draw(&tab[0], kind) - 1;
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
		ol[r] = (int) (w - w0);
	}
	PutRNGstate();
	if (w < total) {
		SET_VECTOR_ELT(out, 0, xlengthgets(start, w));
		SET_VECTOR_ELT(out, 1, xlengthgets(seg, w));
	}
	UNPROTECT(1);
	return out;
}

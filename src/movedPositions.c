#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif

/* at most this many threads, and only from this many bait blocks */
#define MOVED_MAX_THREADS 4
#define MOVED_PAR_MIN 16384

/* R's integer + and - give NA outside [-INT_MAX, INT_MAX] */
#define INT_OK(v) ((v) <= INT_MAX && (v) >= -(int64_t) INT_MAX)

/* as.integer() of one double: NA for NaN and outside the integer range */
static inline int int_of_real(double v)
{
	if (ISNAN(v) || v >= 2147483648.0 || v <= -2147483648.0)
		return NA_INTEGER;
	return (int) v;
}

/* bait_blocks()' inputs and each block's segment, from which each block's
 * start, width and shift are derived as bait_blocks() computes them and its
 * seqlevel code as match(chr_names, seqlevels(y))[bait_seg] gives it: start
 * (si or sd) per block, width (wi or wd) recycled over n_w, target start (ti
 * or td) per block of a replicate, n_ts of them, and seg_code per segment */
typedef struct {
	R_xlen_t n, n_w, n_ts;
	int n_seg;
	const int *seg, *seg_code, *si, *wi, *ti;
	const double *sd, *wd, *td;
} bait_raw;

/* block b's code: NA for a segment index outside [1, n_seg], where R's
 * subset gives NA or drops the element */
static inline int raw_code(const bait_raw *r, R_xlen_t b)
{
	int s = r->seg[b];
	return s >= 1 && s <= r->n_seg ? r->seg_code[s - 1] : NA_INTEGER;
}

/* block b's start, width (index j) and shift (target index t), once
 * raw_block() has accepted them */
static inline int raw_start(const bait_raw *r, R_xlen_t b)
{
	return r->si != NULL ? r->si[b] : (int) r->sd[b];
}

static inline int raw_width(const bait_raw *r, R_xlen_t j)
{
	return r->wi != NULL ? r->wi[j] : (int) r->wd[j];
}

static inline int raw_shift(const bait_raw *r, R_xlen_t t, int sb)
{
	if (r->ti != NULL)
		return r->ti[t] == NA_INTEGER ? NA_INTEGER :
			(int) ((int64_t) r->ti[t] - sb);
	return int_of_real(r->td[t] - (double) sb);
}

/* bait_blocks()' checks of block b: 0 where it declines */
static inline int raw_block(const bait_raw *r, R_xlen_t b, R_xlen_t j,
	R_xlen_t t)
{
	int sb;
	if (r->si != NULL) {
		sb = r->si[b];
		if (sb == NA_INTEGER)
			return 0;
	} else {
		double v = r->sd[b];
		if (ISNAN(v) || v >= 2147483648.0 || v <= -2147483648.0)
			return 0;
		sb = (int) v;
	}
	int wb;
	if (r->wi != NULL) {
		wb = r->wi[j];
		if (wb == NA_INTEGER || wb < 0)
			return 0;
	} else {
		double v = r->wd[j];
		if (ISNAN(v) || v < 0 || v >= 2147483648.0)
			return 0;
		wb = (int) v;
	}
	int64_t e = (int64_t) sb + wb - 1;
	if (!INT_OK(e))
		return 0;
	if (r->ti != NULL && r->ti[t] != NA_INTEGER &&
		!INT_OK((int64_t) r->ti[t] - sb))
		return 0;
	return 1;
}

/* Fixed-width buckets along one sequence's sorted positions: t[k] is the
 * first entry whose bucket is k or later, t[nb] the sequence's end. */
typedef struct {
	int64_t lo;
	int64_t nb;
	int sh;
	R_xlen_t off;
} buckets;

static R_xlen_t plan_buckets(const int *p, const int *seg, int n_lev,
	buckets *b)
{
	R_xlen_t off = 0;
	for (int c = 0; c < n_lev; c++) {
		int cs = seg[c], ce = seg[c + 1];
		b[c].off = off;
		b[c].lo = 0;
		b[c].sh = 0;
		b[c].nb = 0;
		if (ce > cs) {
			int64_t span = (int64_t) p[ce - 1] - p[cs];
			int64_t n = ce - cs;
			int sh = 0;
			while ((span >> sh) + 1 > 2 * n)
				sh++;
			b[c].lo = p[cs];
			b[c].sh = sh;
			b[c].nb = (span >> sh) + 1;
		}
		off += b[c].nb + 1;
	}
	return off;
}

static void fill_buckets(const int *p, const int *seg, int n_lev,
	const buckets *b, int *t_all)
{
	for (int c = 0; c < n_lev; c++) {
		int *t = t_all + b[c].off;
		int64_t k = 0;
		for (int i = seg[c]; i < seg[c + 1]; i++) {
			int64_t bi = ((int64_t) p[i] - b[c].lo) >> b[c].sh;
			while (k <= bi)
				t[k++] = i;
		}
		while (k <= b[c].nb)
			t[k++] = seg[c + 1];
	}
}

/* The end of the entries in [cs, ce) at or before x */
static inline int count_le(const int *p, int cs, int ce, const buckets *b,
	const int *t_all, int64_t x)
{
	if (x < b->lo)
		return cs;
	int64_t k = (x - b->lo) >> b->sh;
	if (k >= b->nb)
		return ce;
	const int *t = t_all + b->off;
	int a = t[k], e = t[k + 1];
	while (e - a > 8) {
		int m = a + (e - a) / 2;
		if (p[m] <= x)
			a = m + 1;
		else
			e = m;
	}
	while (a < e && p[a] <= x)
		a++;
	return a;
}

/* shift_swap_trim_positions() for one hit; 0 where R gives NA */
static inline int move_one(int start, int width, int shift, int code,
	const int *seqlen, const int *bounded, int *new_start, int *new_width)
{
	int64_t ns = (int64_t) start + shift;
	if (!INT_OK(ns))
		return 0;
	int64_t ne = ns + ((int64_t) width - 1);
	if (!INT_OK(ne))
		return 0;
	int64_t w = width;
	if (bounded[code - 1]) {
		int64_t len = seqlen[code - 1];
		if (ns < 1 || ne > len) {
			int64_t os = ns < 1 ? 1 : ns;
			int64_t oe = ne < len ? ne : len;
			int64_t d = oe - os;
			if (!INT_OK(d) || !INT_OK(d + 1))
				return 0;
			w = d + 1 < 0 ? 0 : d + 1;
			ns = os;
		}
	}
	*new_start = (int) ns;
	*new_width = (int) w;
	return 1;
}

/* The first of the intervals [lo, hi) of a sequence, sorted and disjoint,
 * whose end is at or after x; k, the answer for the previous row when it
 * was on the same sequence (same != 0), is kept when it still is. */
static inline int first_end_ge(const int *ie, int lo, int hi, int k,
	int same, int x)
{
	if (same && (k == lo || ie[k - 1] < x) && (k == hi || ie[k] >= x))
		return k;
	int a = lo, b = hi;
	while (a < b) {
		int m = a + (b - a) / 2;
		if (ie[m] < x)
			a = m + 1;
		else
			b = m;
	}
	return a;
}

/* The exclude step for one moved row (code, ns, ne = ns + width - 1) in
 * sequence c's intervals [lo, hi): for mode 1 (drop) 1 when it overlaps
 * none, else 0, as outside_intervals(); for mode 2 (trim) the number of
 * gaps it overlaps, its pieces, as intersect_gaps(). *k is the first
 * interval ending at or after ns. */
static inline int ex_count(int mode, const int *xs, const int *xe, int lo,
	int hi, int *k, int same, int ns, int ne)
{
	int a = first_end_ge(xe, lo, hi, *k, same, ns);
	*k = a;
	if (mode == 1)
		return !(a < hi && xs[a] <= ne);
	int j = a;
	while (j < hi && xs[j] <= ne)
		j++;
	return j - a;
}

/* What both passes read: the bait blocks, y sorted by seqlevel with its
 * bucket tables, the targets and the exclude tables */
typedef struct {
	/* the blocks' codes, starts, widths and shifts, or NULL with 'raw' */
	const int *bc, *bs, *bw, *shift, *tc, *seqlen, *bounded;
	const bait_raw *raw;
	const int *ys, *ye, *yw, *em, *n_short, *seg, *o;
	const int *t_start, *t_max;
	const buckets *bk_start, *bk_max;
	const int *xseg, *xs, *xe;
	int mode, n_lev;
	R_xlen_t n_t, n_t_rep, n_blocks;
	int *first, *last;
	unsigned char *plain;
	/* columns of y gathered at each row's index in y: 8-byte (double) or
	 * 4-byte (integer, logical) elements, copied bit for bit */
	int ncol;
	const unsigned char *col8;
	const void *const *col_in;
	void *const *col_out;
} moved_ctx;

/* Row 'at' of each gathered column, from row 'src' of y */
static inline void put_cols(const moved_ctx *x, R_xlen_t at, int src)
{
	for (int j = 0; j < x->ncol; j++) {
		if (x->col8[j])
			memcpy((double *) x->col_out[j] + at,
				(const double *) x->col_in[j] + src, sizeof(double));
		else
			memcpy((int *) x->col_out[j] + at,
				(const int *) x->col_in[j] + src, sizeof(int));
	}
}

/* Rows at.. of each gathered column for the entries lo..hi - 1 of y
 * sorted by seqlevel */
static inline void put_cols_run(const moved_ctx *x, R_xlen_t at, int lo,
	int hi)
{
	const int *o = x->o;
	for (int j = 0; j < x->ncol; j++) {
		R_xlen_t a = at;
		if (x->col8[j]) {
			double *d = (double *) x->col_out[j];
			const double *s = (const double *) x->col_in[j];
			for (int k = lo; k < hi; k++, a++)
				memcpy(d + a, s + (o == NULL ? k : o[k]), sizeof(double));
		} else {
			int *d = (int *) x->col_out[j];
			const int *s = (const int *) x->col_in[j];
			for (int k = lo; k < hi; k++, a++)
				memcpy(d + a, s + (o == NULL ? k : o[k]), sizeof(int));
		}
	}
}

/* Pass 1 over blocks [b0, b1): each block's candidates, the checks R makes
 * over the hits, and the numbers of hits and rows kept; 0 where R declines */
static int moved_pass1(const moved_ctx *x, R_xlen_t b0, R_xlen_t b1,
	R_xlen_t *hits, R_xlen_t *keeps)
{
	const int mode = x->mode, n_lev = x->n_lev;
	const int *ys = x->ys, *ye = x->ye, *yw = x->yw, *em = x->em;
	const int *xseg = x->xseg, *xs = x->xs, *xe = x->xe;
	R_xlen_t n_hit = 0, n_keep = 0;
	int ex_prev = 0, ex_k = 0;
	const bait_raw *rw = x->raw;
	/* with 'raw', the indices of block b's width and target start */
	R_xlen_t j = rw != NULL ? b0 % rw->n_w : 0;
	R_xlen_t t = rw != NULL ? b0 % rw->n_ts : 0;
	for (R_xlen_t b = b0; b < b1; b++) {
		int bc, bstart, bw;
		R_xlen_t tb = t;
		if (rw != NULL) {
			bc = raw_code(rw, b);
			bstart = raw_start(rw, b);
			bw = raw_width(rw, j);
			if (++j == rw->n_w)
				j = 0;
			if (++t == rw->n_ts)
				t = 0;
		} else {
			bc = x->bc[b];
			bstart = x->bs[b];
			bw = x->bw[b];
		}
		int c = bc - 1, cs = x->seg[c], ce = x->seg[c + 1];
		int lo = count_le(em, cs, ce, x->bk_max + c, x->t_max,
			(int64_t) bstart - 1);
		int hi = count_le(ys, cs, ce, x->bk_start + c, x->t_start,
			(int64_t) bstart + bw - 1);
		x->first[b] = lo;
		x->last[b] = hi;
		x->plain[b] = 0;
		if (lo >= hi)
			continue;
		int sh = rw != NULL ? raw_shift(rw, tb, bstart) : x->shift[b];
		int code = b < x->n_t_rep ? x->tc[b % x->n_t] : NA_INTEGER;
		if (x->n_short[hi] == x->n_short[lo] && sh != NA_INTEGER &&
			code != NA_INTEGER && code >= 1 && code <= n_lev) {
			int64_t s0 = (int64_t) ys[lo] + sh, e1 = (int64_t) em[hi - 1] + sh;
			if (x->bounded[code - 1] ? (s0 >= 1 && e1 <= x->seqlen[code - 1]) :
				(INT_OK(s0) && INT_OK(e1))) {
				int ok = 1;
				if (mode != 0) {
					int xlo = xseg[code - 1], xhi = xseg[code];
					int a = first_end_ge(xe, xlo, xhi, ex_k,
						code == ex_prev, (int) s0);
					ex_prev = code;
					ex_k = a;
					ok = mode == 1 ? (a == xhi || xs[a] > e1) :
						(a < xhi && xs[a] <= s0 && e1 <= xe[a]);
				}
				if (ok) {
					x->plain[b] = 1;
					n_hit += hi - lo;
					n_keep += hi - lo;
					continue;
				}
			}
		}
		for (int k = lo; k < hi; k++) {
			if (ye[k] < bstart)
				continue;
			n_hit++;
			int ns, nw;
			if (sh == NA_INTEGER || code == NA_INTEGER || code < 1 ||
				code > n_lev || !move_one(ys[k], yw[k], sh, code, x->seqlen,
				x->bounded, &ns, &nw))
				return 0;
			if (nw <= 0)
				continue;
			if (mode == 0) {
				n_keep++;
				continue;
			}
			n_keep += ex_count(mode, xs, xe, xseg[code - 1], xseg[code], &ex_k,
				code == ex_prev, ns, ns + (nw - 1));
			ex_prev = code;
		}
	}
	*hits = n_hit;
	*keeps = n_keep;
	return 1;
}

/* Pass 2 over blocks [b0, b1): their rows, written from row 'at' of out
 * (out[0], 's', out[1], the code, and out[4], the replicate, left out when
 * NULL) and of the gathered columns; last[b] then holds the number of rows
 * block b wrote */
static void moved_pass2(const moved_ctx *x, R_xlen_t b0, R_xlen_t b1,
	R_xlen_t at, int *const *out)
{
	const int mode = x->mode;
	const int *ys = x->ys, *ye = x->ye, *yw = x->yw, *o = x->o;
	const int *xseg = x->xseg, *xs = x->xs, *xe = x->xe;
	int *sv = out[0], *cv = out[1], *blk = out[4];
	int ex_prev = 0, ex_k = 0;
	const bait_raw *rw = x->raw;
	R_xlen_t t = rw != NULL ? b0 % rw->n_ts : 0;
	for (R_xlen_t b = b0; b < b1; b++) {
		int lo = x->first[b], hi = x->last[b];
		R_xlen_t tb = t;
		if (rw != NULL && ++t == rw->n_ts)
			t = 0;
		if (lo >= hi) {
			x->last[b] = 0;
			continue;
		}
		R_xlen_t at0 = at;
		int bstart = rw != NULL ? raw_start(rw, b) : x->bs[b];
		int sh = rw != NULL ? raw_shift(rw, tb, bstart) : x->shift[b];
		int code = x->tc[b % x->n_t];
		int rep = (int) (b / x->n_blocks) + 1;
		if (x->plain[b]) {
			if (sv != NULL && cv != NULL) {
				for (int k = lo; k < hi; k++, at++) {
					sv[at] = (o == NULL ? k : o[k]) + 1;
					cv[at] = code;
					out[2][at] = ys[k] + sh;
					out[3][at] = yw[k];
				}
			} else if (cv != NULL) {
				for (int k = lo; k < hi; k++, at++) {
					cv[at] = code;
					out[2][at] = ys[k] + sh;
					out[3][at] = yw[k];
				}
			} else if (sv != NULL) {
				for (int k = lo; k < hi; k++, at++) {
					sv[at] = (o == NULL ? k : o[k]) + 1;
					out[2][at] = ys[k] + sh;
					out[3][at] = yw[k];
				}
			} else {
				for (int k = lo; k < hi; k++, at++) {
					out[2][at] = ys[k] + sh;
					out[3][at] = yw[k];
				}
			}
			if (sv == NULL)
				put_cols_run(x, at0, lo, hi);
			if (blk != NULL)
				for (R_xlen_t i = at0; i < at; i++)
					blk[i] = rep;
			x->last[b] = hi - lo;
			continue;
		}
		int xlo = mode != 0 ? xseg[code - 1] : 0;
		int xhi = mode != 0 ? xseg[code] : 0;
		for (int k = lo; k < hi; k++) {
			if (ye[k] < bstart)
				continue;
			/* pass 1 checked that move_one() succeeds here */
			int ns = 0, nw = 0;
			move_one(ys[k], yw[k], sh, code, x->seqlen, x->bounded, &ns, &nw);
			if (nw <= 0)
				continue;
			int s_k = (o == NULL ? k : o[k]) + 1;
			if (mode == 0) {
				if (sv != NULL)
					sv[at] = s_k;
				put_cols(x, at, s_k - 1);
				if (cv != NULL)
					cv[at] = code;
				out[2][at] = ns;
				out[3][at] = nw;
				at++;
				continue;
			}
			int ne = ns + (nw - 1);
			int a = first_end_ge(xe, xlo, xhi, ex_k, code == ex_prev, ns);
			ex_prev = code;
			ex_k = a;
			if (mode == 1) {
				if (a < xhi && xs[a] <= ne)
					continue;
				if (sv != NULL)
					sv[at] = s_k;
				put_cols(x, at, s_k - 1);
				if (cv != NULL)
					cv[at] = code;
				out[2][at] = ns;
				out[3][at] = nw;
				at++;
				continue;
			}
			for (int j = a; j < xhi && xs[j] <= ne; j++, at++) {
				int ps = ns > xs[j] ? ns : xs[j];
				int pe = ne < xe[j] ? ne : xe[j];
				if (sv != NULL)
					sv[at] = s_k;
				put_cols(x, at, s_k - 1);
				if (cv != NULL)
					cv[at] = code;
				out[2][at] = ps;
				out[3][at] = pe - ps + 1;
			}
		}
		if (blk != NULL)
			for (R_xlen_t i = at0; i < at; i++)
				blk[i] = rep;
		x->last[b] = (int) (at - at0);
	}
}

/* The moved ranges of every bait block in one pass: per block, in order,
 * the ranges of y it overlaps, by increasing index in y, with their index
 * in y ('s'), target seqlevel code, shifted and trimmed start and width,
 * and, with with_block, replicate ('block'), keeping width > 0, and the
 * runs of the code column ('run_code', 'run_length'); without with_block,
 * the number of rows of each replicate ('lens') in place of 'block', as
 * tabulate(block, n_rep) counts them. R_NilValue wherever the R
 * route (bait_overlaps(), shift_swap_trim_positions()) would not take
 * sorted_overlaps() for every block or could decline, or has no hit.
 *
 * With mode 1 or 2 the exclude step is applied to those rows: sequence c
 * holds the sorted, disjoint intervals xseg[c - 1] to xseg[c] - 1 of xs and
 * xe; mode 1 (drop) keeps the rows overlapping none of them, mode 2 (trim)
 * puts in each row's place its pieces clipped to each interval (gap) it
 * overlaps, by gap. FALSE where the exclude tables are malformed or the
 * pieces outnumber an integer index.
 *
 * With 'with_cols', the result has an eighth element, 'cols': when cols_ is
 * a list of double, integer or logical vectors, each of y's length, the
 * rows' values of each, y's at the row's index, in a list named as cols_,
 * and 's' is then NULL; otherwise 'cols' is NULL and 's' is written. The
 * caller passes only vectors without attributes (gather_columns()).
 *
 * Without 'with_code', 'code' is NULL: only its runs are returned.
 *
 * With 'raw' (b_code_, b_start_, b_width_ and shift_ unread), each block's
 * code, start, width and shift are derived from it, and R_NilValue is also
 * returned wherever bait_blocks() declines. */
static SEXP moved_impl(SEXP b_code_, SEXP b_start_, SEXP b_width_,
	SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_,
	int mode, SEXP xseg_, SEXP xs_, SEXP xe_, int with_block, int with_cols,
	SEXP cols_, int with_code, const bait_raw *raw)
{
	if ((raw == NULL && (TYPEOF(b_code_) != INTSXP ||
		TYPEOF(b_start_) != INTSXP || TYPEOF(b_width_) != INTSXP ||
		TYPEOF(shift_) != INTSXP)) || TYPEOF(n_rep_) != INTSXP ||
		TYPEOF(y_code_) != INTSXP || TYPEOF(y_start_) != INTSXP ||
		TYPEOF(y_end_) != INTSXP || TYPEOF(y_width_) != INTSXP ||
		TYPEOF(circular_) != LGLSXP || TYPEOF(seqlen_) != INTSXP ||
		TYPEOF(bounded_) != LGLSXP ||
		TYPEOF(t_code_) != INTSXP || XLENGTH(n_rep_) != 1)
		return R_NilValue;
	R_xlen_t n = raw != NULL ? raw->n : XLENGTH(b_code_);
	R_xlen_t n_y = XLENGTH(y_code_);
	R_xlen_t n_t = XLENGTH(t_code_);
	int n_rep = INTEGER(n_rep_)[0];
	int n_lev = LENGTH(circular_);
	if (n == 0 || n_rep == NA_INTEGER || n_rep < 1 || n % n_rep != 0 ||
		n > INT_MAX || n_y > INT_MAX || (raw == NULL &&
		(XLENGTH(b_start_) != n || XLENGTH(b_width_) != n ||
		XLENGTH(shift_) != n)) ||
		XLENGTH(y_start_) != n_y || XLENGTH(y_end_) != n_y ||
		XLENGTH(y_width_) != n_y || XLENGTH(seqlen_) != n_lev ||
		XLENGTH(bounded_) != n_lev || (int64_t) n_rep * n_lev > INT_MAX)
		return R_NilValue;

	/* the exclude tables: one run of intervals per seqlevel of y, each
	 * sorted and disjoint, as exclude_intervals() and gap_tables() give */
	const int *xseg = NULL, *xs = NULL, *xe = NULL;
	if (mode != 0) {
		if (TYPEOF(xseg_) != INTSXP || TYPEOF(xs_) != INTSXP ||
			TYPEOF(xe_) != INTSXP || XLENGTH(xseg_) != (R_xlen_t) n_lev + 1 ||
			XLENGTH(xe_) != XLENGTH(xs_) || XLENGTH(xs_) > INT_MAX)
			return ScalarLogical(FALSE);
		xseg = INTEGER(xseg_);
		xs = INTEGER(xs_);
		xe = INTEGER(xe_);
		if (xseg[0] != 0 || xseg[n_lev] != XLENGTH(xs_))
			return ScalarLogical(FALSE);
		for (int c = 0; c < n_lev; c++) {
			if (xseg[c + 1] < xseg[c])
				return ScalarLogical(FALSE);
			for (int j = xseg[c]; j < xseg[c + 1]; j++)
				if (xs[j] == NA_INTEGER || xe[j] == NA_INTEGER ||
					xe[j] < xs[j] || (j > xseg[c] && xs[j] <= xe[j - 1]))
					return ScalarLogical(FALSE);
		}
	}
	R_xlen_t n_blocks = n / n_rep;
	const int *yc = INTEGER(y_code_);

	/* y sorted by seqlevel, stably, as order(method = "radix") does */
	int *seg = (int *) R_alloc(n_lev + 1, sizeof(int));
	memset(seg, 0, (n_lev + 1) * sizeof(int));
	int sorted = 1;
	for (R_xlen_t i = 0; i < n_y; i++) {
		int c = yc[i];
		if (c == NA_INTEGER || c < 1 || c > n_lev)
			return R_NilValue;
		seg[c]++;
		if (i > 0 && c < yc[i - 1])
			sorted = 0;
	}
	for (int c = 1; c <= n_lev; c++)
		seg[c] += seg[c - 1];
	const int *ys = INTEGER(y_start_), *ye = INTEGER(y_end_);
	const int *yw = INTEGER(y_width_);
	int *o = NULL;
	if (!sorted) {
		o = (int *) R_alloc(n_y, sizeof(int));
		int *pos = (int *) R_alloc(n_lev, sizeof(int));
		memcpy(pos, seg, n_lev * sizeof(int));
		for (R_xlen_t i = 0; i < n_y; i++)
			o[pos[yc[i] - 1]++] = (int) i;
		int *s2 = (int *) R_alloc(n_y, sizeof(int));
		int *e2 = (int *) R_alloc(n_y, sizeof(int));
		int *w2 = (int *) R_alloc(n_y, sizeof(int));
		for (R_xlen_t k = 0; k < n_y; k++) {
			s2[k] = ys[o[k]];
			e2[k] = ye[o[k]];
			w2[k] = yw[o[k]];
		}
		ys = s2;
		ye = e2;
		yw = w2;
	}
	/* sorted_overlaps(): no empty range, starts in order on each sequence;
	 * end_max, the running maximum of ends along each sequence */
	int *em = (int *) R_alloc(n_y > 0 ? n_y : 1, sizeof(int));
	int *n_short = (int *) R_alloc(n_y + 1, sizeof(int));
	n_short[0] = 0;
	for (int c = 0; c < n_lev; c++) {
		for (int k = seg[c]; k < seg[c + 1]; k++) {
			if (ye[k] < ys[k] || (k > seg[c] && ys[k] < ys[k - 1]))
				return R_NilValue;
			em[k] = (k > seg[c] && em[k - 1] > ye[k]) ? em[k - 1] : ye[k];
			n_short[k + 1] = n_short[k] + (ye[k] < em[k]);
		}
	}

	/* bait_overlaps(): every block on sorted_overlaps()' path, i.e. its
	 * replicate has fewer blocks on its sequence than y has ranges there;
	 * sorted_overlaps(): no empty block and none on a circular sequence;
	 * end(bait), start + width - 1L, not NA. With 'raw', bait_blocks()'
	 * checks first. */
	const int *bc = raw == NULL ? INTEGER(b_code_) : NULL;
	const int *bs = raw == NULL ? INTEGER(b_start_) : NULL;
	const int *bw = raw == NULL ? INTEGER(b_width_) : NULL;
	const int *circ = LOGICAL(circular_);
	int *cnt = (int *) R_alloc(n_lev > 0 ? n_lev : 1, sizeof(int));
	memset(cnt, 0, (n_lev > 0 ? n_lev : 1) * sizeof(int));
	for (R_xlen_t r = 0; r < n_rep; r++) {
		R_xlen_t b0 = r * n_blocks, b1 = b0 + n_blocks;
		R_xlen_t j = raw != NULL ? b0 % raw->n_w : 0;
		for (R_xlen_t b = b0; b < b1; b++) {
			int c, sb, wb;
			if (raw != NULL) {
				/* a replicate's blocks take the targets in order */
				if (!raw_block(raw, b, j, b - b0))
					return R_NilValue;
				c = raw_code(raw, b);
				sb = raw_start(raw, b);
				wb = raw_width(raw, j);
				if (++j == raw->n_w)
					j = 0;
			} else {
				c = bc[b];
				sb = bs[b];
				wb = bw[b];
			}
			int64_t e = (int64_t) sb + wb;
			if (c == NA_INTEGER || c < 1 || c > n_lev || !INT_OK(e) ||
				!INT_OK(e - 1) || e - 1 < sb || circ[c - 1] == TRUE)
				return R_NilValue;
			cnt[c - 1]++;
		}
		for (R_xlen_t b = b0; b < b1; b++) {
			int c = raw != NULL ? raw_code(raw, b) : bc[b];
			if (cnt[c - 1] >= seg[c] - seg[c - 1])
				return R_NilValue;
		}
		for (R_xlen_t b = b0; b < b1; b++)
			cnt[(raw != NULL ? raw_code(raw, b) : bc[b]) - 1] = 0;
	}

	buckets *bk_start = (buckets *) R_alloc(n_lev, sizeof(buckets));
	buckets *bk_max = (buckets *) R_alloc(n_lev, sizeof(buckets));
	int *t_start = (int *) R_alloc(plan_buckets(ys, seg, n_lev, bk_start),
		sizeof(int));
	int *t_max = (int *) R_alloc(plan_buckets(em, seg, n_lev, bk_max),
		sizeof(int));
	fill_buckets(ys, seg, n_lev, bk_start, t_start);
	fill_buckets(em, seg, n_lev, bk_max, t_max);

	/* the blocks in chunks of consecutive blocks: pass 1 counts each
	 * chunk's rows, and pass 2 writes them from the sum of the counts
	 * before it, so the rows land as one pass in block order puts them.
	 * The answer for an exclude interval does not depend on the previous
	 * row's, which each chunk starts without. */
	int n_thr = 1;
#ifdef _OPENMP
	if (n >= MOVED_PAR_MIN) {
		n_thr = omp_get_max_threads();
		if (n_thr > MOVED_MAX_THREADS)
			n_thr = MOVED_MAX_THREADS;
		if (n_thr < 1)
			n_thr = 1;
	}
#endif
	int n_chunk = n_thr == 1 ? 1 : 8 * n_thr;
	R_xlen_t *c_hit = (R_xlen_t *) R_alloc(n_chunk, sizeof(R_xlen_t));
	R_xlen_t *c_keep = (R_xlen_t *) R_alloc(n_chunk, sizeof(R_xlen_t));
	int *c_ok = (int *) R_alloc(n_chunk, sizeof(int));
	moved_ctx x = {
		.bc = bc, .bs = bs, .bw = bw,
		.shift = raw == NULL ? INTEGER(shift_) : NULL, .raw = raw,
		.tc = INTEGER(t_code_), .seqlen = INTEGER(seqlen_),
		.bounded = LOGICAL(bounded_), .ys = ys, .ye = ye, .yw = yw, .em = em,
		.n_short = n_short, .seg = seg, .o = o, .t_start = t_start,
		.t_max = t_max, .bk_start = bk_start, .bk_max = bk_max, .xseg = xseg,
		.xs = xs, .xe = xe, .mode = mode, .n_lev = n_lev, .n_t = n_t,
		.n_t_rep = n_t * (R_xlen_t) n_rep, .n_blocks = n_blocks,
		.first = (int *) R_alloc(n, sizeof(int)),
		.last = (int *) R_alloc(n, sizeof(int)),
		.plain = (unsigned char *) R_alloc(n, 1),
		.ncol = 0, .col8 = NULL, .col_in = NULL, .col_out = NULL
	};

	/* the columns to gather, when every one can be */
	int gather = with_cols && TYPEOF(cols_) == VECSXP;
	R_xlen_t ncol = gather ? XLENGTH(cols_) : 0;
	for (R_xlen_t j = 0; gather && j < ncol; j++) {
		SEXP v = VECTOR_ELT(cols_, j);
		int ty = TYPEOF(v);
		if ((ty != REALSXP && ty != INTSXP && ty != LGLSXP) ||
			XLENGTH(v) != n_y)
			gather = 0;
	}
	if (!gather)
		ncol = 0;

	/* pass 1 */
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_thr) schedule(dynamic, 1) if(n_thr > 1)
#endif
	for (int c = 0; c < n_chunk; c++)
		c_ok[c] = moved_pass1(&x, n * c / n_chunk, n * (c + 1) / n_chunk,
			c_hit + c, c_keep + c);
	R_xlen_t n_hit = 0, n_keep = 0;
	for (int c = 0; c < n_chunk; c++)
		if (!c_ok[c])
			return R_NilValue;
	for (int c = 0; c < n_chunk; c++) {
		R_xlen_t k = c_keep[c];
		n_hit += c_hit[c];
		c_keep[c] = n_keep;
		n_keep += k;
	}
	if (n_hit == 0 || n_hit > INT_MAX)
		return R_NilValue;
	if (n_keep > INT_MAX)
		return ScalarLogical(FALSE);

	/* pass 2: the rows */
	/* element 4 is 'block' (n_keep rows) or 'lens' (n_rep counts) */
	int n_el = with_cols ? 8 : 7;
	SEXP ans = PROTECT(allocVector(VECSXP, n_el));
	SEXP nms = PROTECT(allocVector(STRSXP, n_el));
	const char *names[8] = {"s", "code", "start", "width",
		with_block ? "block" : "lens", "run_code", "run_length", "cols"};
	int *out[5];
	out[0] = NULL;
	out[1] = NULL;
	for (int j = gather ? 1 : 0; j < 5; j++) {
		if (j == 1 && !with_code)
			continue;
		SET_VECTOR_ELT(ans, j, allocVector(INTSXP,
			j < 4 || with_block ? n_keep : (R_xlen_t) n_rep));
		out[j] = INTEGER(VECTOR_ELT(ans, j));
	}
	for (int j = 0; j < n_el; j++)
		SET_STRING_ELT(nms, j, mkChar(names[j]));
	setAttrib(ans, R_NamesSymbol, nms);
	int *lens = NULL;
	if (!with_block) {
		lens = out[4];
		out[4] = NULL;
	}
	if (gather) {
		SEXP cols = allocVector(VECSXP, ncol);
		SET_VECTOR_ELT(ans, 7, cols);
		setAttrib(cols, R_NamesSymbol, getAttrib(cols_, R_NamesSymbol));
		const void **cin = (const void **) R_alloc(ncol > 0 ? ncol : 1,
			sizeof(void *));
		void **cout = (void **) R_alloc(ncol > 0 ? ncol : 1, sizeof(void *));
		unsigned char *c8 = (unsigned char *) R_alloc(ncol > 0 ? ncol : 1, 1);
		for (R_xlen_t j = 0; j < ncol; j++) {
			SEXP v = VECTOR_ELT(cols_, j);
			SEXP w = allocVector(TYPEOF(v), n_keep);
			SET_VECTOR_ELT(cols, j, w);
			c8[j] = TYPEOF(v) == REALSXP;
			if (TYPEOF(v) == REALSXP) {
				cin[j] = REAL_RO(v);
				cout[j] = REAL(w);
			} else if (TYPEOF(v) == INTSXP) {
				cin[j] = INTEGER_RO(v);
				cout[j] = INTEGER(w);
			} else {
				cin[j] = LOGICAL_RO(v);
				cout[j] = LOGICAL(w);
			}
		}
		x.ncol = (int) ncol;
		x.col8 = c8;
		x.col_in = (const void *const *) cin;
		x.col_out = (void *const *) cout;
	}
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_thr) schedule(dynamic, 1) if(n_thr > 1)
#endif
	for (int c = 0; c < n_chunk; c++)
		moved_pass2(&x, n * c / n_chunk, n * (c + 1) / n_chunk, c_keep[c],
			out);

	/* the runs of the code column: each block's rows share its target
	 * sequence, so one run per block that wrote rows, merged with the run
	 * before when the sequence is the same; maximal, as Rle() finds them */
	R_xlen_t n_run = 0;
	int prev = NA_INTEGER;
	for (R_xlen_t b = 0, t = 0; b < n; b++, t = t + 1 == n_t ? 0 : t + 1) {
		if (x.last[b] != 0 && x.tc[t] != prev) {
			n_run++;
			prev = x.tc[t];
		}
	}
	SET_VECTOR_ELT(ans, 5, allocVector(INTSXP, n_run));
	SET_VECTOR_ELT(ans, 6, allocVector(INTSXP, n_run));
	int *rc = INTEGER(VECTOR_ELT(ans, 5)), *rl = INTEGER(VECTOR_ELT(ans, 6));
	R_xlen_t r = -1;
	prev = NA_INTEGER;
	for (R_xlen_t b = 0, t = 0; b < n; b++, t = t + 1 == n_t ? 0 : t + 1) {
		if (x.last[b] == 0)
			continue;
		if (x.tc[t] != prev) {
			prev = x.tc[t];
			rc[++r] = prev;
			rl[r] = 0;
		}
		rl[r] += x.last[b];
	}
	/* each replicate's rows; their sum, n_keep, fits an int */
	if (lens != NULL) {
		for (int rr = 0; rr < n_rep; rr++) {
			R_xlen_t sum = 0;
			const int *nb = x.last + rr * n_blocks;
			for (R_xlen_t b = 0; b < n_blocks; b++)
				sum += nb[b];
			lens[rr] = (int) sum;
		}
	}
	UNPROTECT(2);
	return ans;
}

/* moved_impl() without an exclude step */
SEXP moved_positions(SEXP b_code_, SEXP b_start_, SEXP b_width_, SEXP n_rep_,
	SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_, SEXP circular_,
	SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_)
{
	return moved_impl(b_code_, b_start_, b_width_, n_rep_, y_code_, y_start_,
		y_end_, y_width_, circular_, seqlen_, bounded_, shift_, t_code_, 0,
		R_NilValue, R_NilValue, R_NilValue, 1, 0, R_NilValue, 1, NULL);
}

static SEXP moved_rows_impl(SEXP b_code_, SEXP b_start_, SEXP b_width_,
	SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, int with_block, int with_cols,
	SEXP cols_, int with_code, const bait_raw *raw)
{
	if (TYPEOF(mode_) != INTSXP || XLENGTH(mode_) != 1)
		return ScalarLogical(FALSE);
	int mode = INTEGER(mode_)[0];
	if (mode != 0 && mode != 1 && mode != 2)
		return ScalarLogical(FALSE);
	return moved_impl(b_code_, b_start_, b_width_, n_rep_, y_code_, y_start_,
		y_end_, y_width_, circular_, seqlen_, bounded_, shift_, t_code_, mode,
		xseg_, xs_, xe_, with_block, with_cols, cols_, with_code, raw);
}

/* moved_impl() with the exclude step of mode_ (0 none, 1 drop, 2 trim) */
SEXP moved_rows(SEXP b_code_, SEXP b_start_, SEXP b_width_, SEXP n_rep_,
	SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_, SEXP circular_,
	SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_, SEXP mode_,
	SEXP xseg_, SEXP xs_, SEXP xe_)
{
	return moved_rows_impl(b_code_, b_start_, b_width_, n_rep_, y_code_,
		y_start_, y_end_, y_width_, circular_, seqlen_, bounded_, shift_,
		t_code_, mode_, xseg_, xs_, xe_, 1, 0, R_NilValue, 1, NULL);
}

/* moved_rows() with the rows per replicate ('lens') in place of 'block',
 * gathering the columns 'cols' (see moved_impl()) */
SEXP moved_rows_lens(SEXP b_code_, SEXP b_start_, SEXP b_width_,
	SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_)
{
	return moved_rows_impl(b_code_, b_start_, b_width_, n_rep_, y_code_,
		y_start_, y_end_, y_width_, circular_, seqlen_, bounded_, shift_,
		t_code_, mode_, xseg_, xs_, xe_, 0, 1, cols_, 1, NULL);
}

/* moved_rows_lens() whose per-row code column is written only when
 * with_code_ is not FALSE; otherwise 'code' is NULL and only its runs
 * ('run_code', 'run_length') are returned */
SEXP moved_rows_opt(SEXP b_code_, SEXP b_start_, SEXP b_width_,
	SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_, SEXP with_code_)
{
	return moved_rows_impl(b_code_, b_start_, b_width_, n_rep_, y_code_,
		y_start_, y_end_, y_width_, circular_, seqlen_, bounded_, shift_,
		t_code_, mode_, xseg_, xs_, xe_, 0, 1, cols_,
		asLogical(with_code_) != FALSE, NULL);
}

/* What move_blocks_bulk() builds for moved_positions() from the bait starts
 * of all replicates, the block widths (recycled as rep_len()) and the
 * target starts of one replicate (recycled n_rep times): the integer start
 * and width slots of IRanges(start, width = rep_len(width, n)) and
 * as.integer(rep(target, n_rep) - start), as a list (start, width, shift).
 * R_NilValue wherever IRanges() errors (an NA start or width, a negative
 * width, a start, width or end outside the integer range), integer
 * subtraction overflows (with a warning) or the lengths are not those of
 * one bait block per target per replicate: the caller then builds them in
 * R. */
SEXP bait_blocks(SEXP start_, SEXP width_, SEXP target_, SEXP n_rep_)
{
	int ts = TYPEOF(start_), tw = TYPEOF(width_), tt = TYPEOF(target_);
	if ((ts != INTSXP && ts != REALSXP) || (tw != INTSXP && tw != REALSXP) ||
		(tt != INTSXP && tt != REALSXP) || TYPEOF(n_rep_) != INTSXP ||
		XLENGTH(n_rep_) != 1)
		return R_NilValue;
	R_xlen_t n = XLENGTH(start_), n_w = XLENGTH(width_);
	R_xlen_t n_t = XLENGTH(target_);
	int n_rep = INTEGER(n_rep_)[0];
	if (n == 0 || n_w == 0 || n_t == 0 || n > INT_MAX ||
		n_rep == NA_INTEGER || n_rep < 1 || n % n_rep != 0 ||
		n / n_rep != n_t)
		return R_NilValue;
	SEXP ans = PROTECT(allocVector(VECSXP, 3));
	SEXP nms = PROTECT(allocVector(STRSXP, 3));
	const char *names[3] = {"start", "width", "shift"};
	for (int k = 0; k < 3; k++) {
		SET_VECTOR_ELT(ans, k, allocVector(INTSXP, n));
		SET_STRING_ELT(nms, k, mkChar(names[k]));
	}
	setAttrib(ans, R_NamesSymbol, nms);
	int *s = INTEGER(VECTOR_ELT(ans, 0)), *w = INTEGER(VECTOR_ELT(ans, 1));
	int *sh = INTEGER(VECTOR_ELT(ans, 2));
	const int *si = ts == INTSXP ? INTEGER_RO(start_) : NULL;
	const double *sd = ts == REALSXP ? REAL_RO(start_) : NULL;
	const int *wi = tw == INTSXP ? INTEGER_RO(width_) : NULL;
	const double *wd = tw == REALSXP ? REAL_RO(width_) : NULL;
	const int *ti = tt == INTSXP ? INTEGER_RO(target_) : NULL;
	const double *td = tt == REALSXP ? REAL_RO(target_) : NULL;
	int ok = 1;
	for (R_xlen_t b = 0, j = 0, t = 0; b < n; b++) {
		int sb;
		if (si != NULL) {
			sb = si[b];
			if (sb == NA_INTEGER) {
				ok = 0;
				break;
			}
		} else {
			double v = sd[b];
			if (ISNAN(v) || v >= 2147483648.0 || v <= -2147483648.0) {
				ok = 0;
				break;
			}
			sb = (int) v;
		}
		int wb;
		if (wi != NULL) {
			wb = wi[j];
			if (wb == NA_INTEGER || wb < 0) {
				ok = 0;
				break;
			}
		} else {
			double v = wd[j];
			if (ISNAN(v) || v < 0 || v >= 2147483648.0) {
				ok = 0;
				break;
			}
			wb = (int) v;
		}
		/* IRanges' end, -1 + start + width, as.integer() under warn = 2 */
		int64_t e = (int64_t) sb + wb - 1;
		if (!INT_OK(e)) {
			ok = 0;
			break;
		}
		int shb;
		if (ti != NULL) {
			int tv = ti[t];
			if (tv == NA_INTEGER) {
				shb = NA_INTEGER;
			} else {
				int64_t d = (int64_t) tv - sb;
				if (!INT_OK(d)) {
					ok = 0;
					break;
				}
				shb = (int) d;
			}
		} else {
			shb = int_of_real(td[t] - (double) sb);
		}
		s[b] = sb;
		w[b] = wb;
		sh[b] = shb;
		if (++j == n_w)
			j = 0;
		if (++t == n_t)
			t = 0;
	}
	UNPROTECT(2);
	return ok ? ans : R_NilValue;
}

/* moved_rows_opt() given bait_blocks()' inputs (bait starts, widths and
 * target starts, integer or double) and each block's segment (b_seg_, an
 * index into seg_code_, each segment's seqlevel code) in place of the
 * blocks' codes, starts, widths and shifts, which it derives where it reads
 * them: R_NilValue wherever bait_blocks() declines, a block's code is NA or
 * its segment index is outside [1, length(seg_code_)], or moved_rows_opt()
 * on bait_blocks()' list and match(chr_names, seqlevels(y))[bait_seg]
 * returns R_NilValue; FALSE where it returns FALSE. */
SEXP moved_rows_bait(SEXP b_seg_, SEXP seg_code_, SEXP b_start_,
	SEXP b_width_, SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_,
	SEXP y_width_, SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP t_start_,
	SEXP t_code_, SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_,
	SEXP with_code_)
{
	int ts = TYPEOF(b_start_), tw = TYPEOF(b_width_), tt = TYPEOF(t_start_);
	if (TYPEOF(b_seg_) != INTSXP || TYPEOF(seg_code_) != INTSXP ||
		(ts != INTSXP && ts != REALSXP) || (tw != INTSXP && tw != REALSXP) ||
		(tt != INTSXP && tt != REALSXP) || TYPEOF(n_rep_) != INTSXP ||
		XLENGTH(n_rep_) != 1 || XLENGTH(seg_code_) > INT_MAX)
		return R_NilValue;
	bait_raw r;
	r.n = XLENGTH(b_start_);
	r.n_w = XLENGTH(b_width_);
	r.n_ts = XLENGTH(t_start_);
	int n_rep = INTEGER(n_rep_)[0];
	if (XLENGTH(b_seg_) != r.n || r.n == 0 || r.n_w == 0 || r.n_ts == 0 ||
		r.n > INT_MAX || n_rep == NA_INTEGER || n_rep < 1 ||
		r.n % n_rep != 0 || r.n / n_rep != r.n_ts)
		return R_NilValue;
	r.n_seg = (int) XLENGTH(seg_code_);
	r.seg = INTEGER_RO(b_seg_);
	r.seg_code = INTEGER_RO(seg_code_);
	r.si = ts == INTSXP ? INTEGER_RO(b_start_) : NULL;
	r.sd = ts == REALSXP ? REAL_RO(b_start_) : NULL;
	r.wi = tw == INTSXP ? INTEGER_RO(b_width_) : NULL;
	r.wd = tw == REALSXP ? REAL_RO(b_width_) : NULL;
	r.ti = tt == INTSXP ? INTEGER_RO(t_start_) : NULL;
	r.td = tt == REALSXP ? REAL_RO(t_start_) : NULL;
	return moved_rows_impl(R_NilValue, R_NilValue, R_NilValue, n_rep_,
		y_code_, y_start_, y_end_, y_width_, circular_, seqlen_, bounded_,
		R_NilValue, t_code_, mode_, xseg_, xs_, xe_, 0, 1, cols_,
		asLogical(with_code_) != FALSE, &r);
}

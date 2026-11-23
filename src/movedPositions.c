/* omp.h before the R headers: Rinternals.h defines match, which LLVM
 * libomp's omp.h uses in a pragma */
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include <Rversion.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <string.h>
#include "segDraws.h"

/* at most this many threads, and only from this many bait blocks */
#define MOVED_MAX_THREADS 4
#define MOVED_PAR_MIN 16384

/* seg_draws_moved() runs pass 2 inside the draws' region, into vectors
 * allocated before it with R's resizable vector API (R >= 4.6.0), of
 * ceil(MOVED_P2_SCALE * k0 * n_rep) + MOVED_P2_SLACK rows, k0 the rows
 * replicate 0 keeps; where the rows outgrow them, or without the API,
 * pass 2 runs after the draws */
#if defined(_OPENMP) && R_VERSION >= R_Version(4, 6, 0)
#define MOVED_P2_REGION 1
#else
#define MOVED_P2_REGION 0
#endif
#define MOVED_P2_SCALE 1.125
#define MOVED_P2_SLACK 65536
/* the pieces replicate 0's pass 1 runs in before the region */
#define MOVED_P1_PIECES 128

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

/* one sequence's table, entries cs to ce - 1 of p, planned in b */
static void fill_buckets(const int *p, int cs, int ce, const buckets *b,
	int *t_all)
{
	int *t = t_all + b->off;
	int64_t k = 0;
	for (int i = cs; i < ce; i++) {
		int64_t bi = ((int64_t) p[i] - b->lo) >> b->sh;
		while (k <= bi)
			t[k++] = i;
	}
	while (k <= b->nb)
		t[k++] = ce;
}

/* sorted_overlaps()' checks of one sequence's ranges cs to ce - 1 (no
 * empty range, starts in order) and their end_max, the running maximum of
 * ends; ns[k + 1] counts the ranges from cs to k ending before their
 * end_max. Returns that count over the sequence, or -1 where a check
 * fails. */
static int end_max(const int *ys, const int *ye, int cs, int ce, int *em,
	int *ns)
{
	int m = 0;
	for (int k = cs; k < ce; k++) {
		if (ye[k] < ys[k] || (k > cs && ys[k] < ys[k - 1]))
			return -1;
		em[k] = (k > cs && em[k - 1] > ye[k]) ? em[k - 1] : ye[k];
		m += ye[k] < em[k];
		ns[k + 1] = m;
	}
	return m;
}

/* 1 where ranges cs to ce - 1 (starts in order) are in the order of
 * IRanges' NCList pre-order, start ascending, end descending, index
 * ascending: no equal-start pair whose later range ends later. Then a
 * block's hits by increasing index are findOverlaps()' on either side. */
static unsigned char nclist_order(const int *ys, const int *ye, int cs,
	int ce)
{
	for (int k = cs + 1; k < ce; k++)
		if (ys[k] == ys[k - 1] && ye[k] > ye[k - 1])
			return 0;
	return 1;
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
	/* with 'lens', the row at which each replicate starts, or NULL */
	int *rep_at;
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

/* Pass 1 over blocks [b0, b1): each block's checks (bait_blocks()' with
 * 'raw', and end(bait) not NA and not before its start; the code and the
 * circular check are made before), its candidates, the checks R makes over
 * the hits, and the numbers of hits and rows kept; 0 where R declines */
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
			/* a replicate's blocks take the targets in order */
			if (!raw_block(rw, b, j, t))
				return 0;
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
		int64_t be = (int64_t) bstart + bw;
		if (!INT_OK(be) || !INT_OK(be - 1) || be - 1 < bstart)
			return 0;
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

/* The rows block b wrote, n of them with target code 'code', added to the
 * runs of its chunk, which are kept in first[b0 + k] (code) and last[b0 + k]
 * (rows): never past b, whose entries pass 2 has read */
static inline void add_run(const moved_ctx *x, R_xlen_t b0, R_xlen_t *nr,
	int code, int n)
{
	if (n == 0)
		return;
	if (*nr > 0 && x->first[b0 + *nr - 1] == code) {
		x->last[b0 + *nr - 1] += n;
		return;
	}
	x->first[b0 + *nr] = code;
	x->last[b0 + *nr] = n;
	(*nr)++;
}

/* Pass 2 over blocks [b0, b1): their rows, written from row 'at' of out
 * (out[0], 's', out[1], the code, and out[4], the replicate, left out when
 * NULL) and of the gathered columns; the chunk's maximal runs of the code
 * column then fill first[b0..] and last[b0..], *n_run of them, and rep_at
 * the first row of each replicate starting in it */
static void moved_pass2(const moved_ctx *x, R_xlen_t b0, R_xlen_t b1,
	R_xlen_t at, int *const *out, R_xlen_t *n_run)
{
	const int mode = x->mode;
	const int *ys = x->ys, *ye = x->ye, *yw = x->yw, *o = x->o;
	const int *xseg = x->xseg, *xs = x->xs, *xe = x->xe;
	int *sv = out[0], *cv = out[1], *blk = out[4];
	int ex_prev = 0, ex_k = 0;
	const bait_raw *rw = x->raw;
	R_xlen_t t = rw != NULL ? b0 % rw->n_ts : 0;
	R_xlen_t nr = 0, rb = b0 % x->n_blocks;
	for (R_xlen_t b = b0; b < b1; b++) {
		int lo = x->first[b], hi = x->last[b];
		R_xlen_t tb = t;
		if (rw != NULL && ++t == rw->n_ts)
			t = 0;
		if (rb == 0 && x->rep_at != NULL)
			x->rep_at[b / x->n_blocks] = (int) at;
		if (++rb == x->n_blocks)
			rb = 0;
		if (lo >= hi)
			continue;
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
			add_run(x, b0, &nr, code, hi - lo);
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
		add_run(x, b0, &nr, code, (int) (at - at0));
	}
	*n_run = nr;
}

/* Where moved_setup() allocates: R_alloc() without an arena; with one, raw
 * vectors kept in its list, which outlive an R_tryCatchError() body */
typedef struct {
	SEXP keep;
	int k;
} m_arena;

#define ARENA_SLOTS 40

static void *m_alloc(m_arena *a, size_t n, size_t sz)
{
	if (a == NULL)
		return R_alloc(n, sz);
	SEXP v = allocVector(RAWSXP, (R_xlen_t) (n * sz));
	SET_VECTOR_ELT(a->keep, a->k++, v);
	return RAW(v);
}

/* What moved_impl() builds before its passes */
typedef struct {
	moved_ctx x;
	R_xlen_t n, n_blocks;
	int n_rep, n_lev, n_thr, by_lev;
	/* per seqlevel, 1 where y's ranges there are in NCList order or none */
	const unsigned char *nc_ok;
	int with_block, with_cols, with_code, gather;
	R_xlen_t ncol;
	const int *circ;
	SEXP cols_;
} moved_plan;

/* moved_impl()'s checks of its arguments, y sorted by seqlevel with its
 * end_max and bucket tables, the pass buffers and the thread count (at most
 * n_threads), in pl: 0, or 1 where moved_impl() returns R_NilValue, 2 where
 * it returns FALSE. The checks of each replicate's blocks are
 * moved_rep_check()'s. */
static int moved_setup(moved_plan *pl, SEXP b_code_, SEXP b_start_,
	SEXP b_width_, SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_,
	SEXP y_width_, SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_,
	SEXP t_code_, int mode, SEXP xseg_, SEXP xs_, SEXP xe_, int with_block,
	int with_cols, SEXP cols_, int with_code, const bait_raw *raw,
	m_arena *arena, int n_threads)
{
	if ((raw == NULL && (TYPEOF(b_code_) != INTSXP ||
		TYPEOF(b_start_) != INTSXP || TYPEOF(b_width_) != INTSXP ||
		TYPEOF(shift_) != INTSXP)) || TYPEOF(n_rep_) != INTSXP ||
		TYPEOF(y_code_) != INTSXP || TYPEOF(y_start_) != INTSXP ||
		TYPEOF(y_end_) != INTSXP || TYPEOF(y_width_) != INTSXP ||
		TYPEOF(circular_) != LGLSXP || TYPEOF(seqlen_) != INTSXP ||
		TYPEOF(bounded_) != LGLSXP ||
		TYPEOF(t_code_) != INTSXP || XLENGTH(n_rep_) != 1)
		return 1;
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
		return 1;

	/* the exclude tables: one run of intervals per seqlevel of y, each
	 * sorted and disjoint, as exclude_intervals() and gap_tables() give */
	const int *xseg = NULL, *xs = NULL, *xe = NULL;
	if (mode != 0) {
		if (TYPEOF(xseg_) != INTSXP || TYPEOF(xs_) != INTSXP ||
			TYPEOF(xe_) != INTSXP || XLENGTH(xseg_) != (R_xlen_t) n_lev + 1 ||
			XLENGTH(xe_) != XLENGTH(xs_) || XLENGTH(xs_) > INT_MAX)
			return 2;
		xseg = INTEGER(xseg_);
		xs = INTEGER(xs_);
		xe = INTEGER(xe_);
		if (xseg[0] != 0 || xseg[n_lev] != XLENGTH(xs_))
			return 2;
		for (int c = 0; c < n_lev; c++) {
			if (xseg[c + 1] < xseg[c])
				return 2;
			for (int j = xseg[c]; j < xseg[c + 1]; j++)
				if (xs[j] == NA_INTEGER || xe[j] == NA_INTEGER ||
					xe[j] < xs[j] || (j > xseg[c] && xs[j] <= xe[j - 1]))
					return 2;
		}
	}
	R_xlen_t n_blocks = n / n_rep;
	const int *yc = INTEGER(y_code_);

	/* y sorted by seqlevel, stably, as order(method = "radix") does */
	int *seg = (int *) m_alloc(arena, n_lev + 1, sizeof(int));
	memset(seg, 0, (n_lev + 1) * sizeof(int));
	int sorted = 1;
	for (R_xlen_t i = 0; i < n_y; i++) {
		int c = yc[i];
		if (c == NA_INTEGER || c < 1 || c > n_lev)
			return 1;
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
		o = (int *) m_alloc(arena, n_y, sizeof(int));
		int *pos = (int *) m_alloc(arena, n_lev, sizeof(int));
		memcpy(pos, seg, n_lev * sizeof(int));
		for (R_xlen_t i = 0; i < n_y; i++)
			o[pos[yc[i] - 1]++] = (int) i;
		int *s2 = (int *) m_alloc(arena, n_y, sizeof(int));
		int *e2 = (int *) m_alloc(arena, n_y, sizeof(int));
		int *w2 = (int *) m_alloc(arena, n_y, sizeof(int));
		for (R_xlen_t k = 0; k < n_y; k++) {
			s2[k] = ys[o[k]];
			e2[k] = ye[o[k]];
			w2[k] = yw[o[k]];
		}
		ys = s2;
		ye = e2;
		yw = w2;
	}
	/* the threads of the per-sequence tables and of the two passes: at
	 * most n_threads, the caller's omp_threads() (R/zzz.R) */
	int n_thr = 1;
#ifdef _OPENMP
	if (n >= MOVED_PAR_MIN && n_threads > 1) {
		n_thr = omp_get_max_threads();
		if (n_thr > n_threads)
			n_thr = n_threads;
		if (n_thr > MOVED_MAX_THREADS)
			n_thr = MOVED_MAX_THREADS;
		if (n_thr < 1)
			n_thr = 1;
	}
#else
	(void) n_threads;
#endif

	/* end_max() per sequence; n_short[k + 1], the ranges up to k along
	 * all sequences ending before end_max, is counted per sequence and
	 * then offset by n_base[c], the counts of the sequences before c */
	int *em = (int *) m_alloc(arena, n_y > 0 ? n_y : 1, sizeof(int));
	int *n_short = (int *) m_alloc(arena, n_y + 1, sizeof(int));
	int *n_base = (int *) m_alloc(arena, n_lev > 0 ? n_lev : 1, sizeof(int));
	unsigned char *nc_ok = (unsigned char *) m_alloc(arena,
		n_lev > 0 ? n_lev : 1, 1);
	n_short[0] = 0;
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_thr) schedule(dynamic, 1) if(n_thr > 1)
#endif
	for (int c = 0; c < n_lev; c++) {
		n_base[c] = end_max(ys, ye, seg[c], seg[c + 1], em, n_short);
		nc_ok[c] = nclist_order(ys, ye, seg[c], seg[c + 1]);
	}
	int n_before = 0;
	for (int c = 0; c < n_lev; c++) {
		int m = n_base[c];
		if (m < 0)
			return 1;
		n_base[c] = n_before;
		n_before += m;
	}

	buckets *bk_start = (buckets *) m_alloc(arena, n_lev, sizeof(buckets));
	buckets *bk_max = (buckets *) m_alloc(arena, n_lev, sizeof(buckets));
	int *t_start = (int *) m_alloc(arena,
		plan_buckets(ys, seg, n_lev, bk_start), sizeof(int));
	int *t_max = (int *) m_alloc(arena,
		plan_buckets(em, seg, n_lev, bk_max), sizeof(int));
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_thr) schedule(dynamic, 1) if(n_thr > 1)
#endif
	for (int c = 0; c < n_lev; c++) {
		int cs = seg[c], ce = seg[c + 1], b = n_base[c];
		fill_buckets(ys, cs, ce, bk_start + c, t_start);
		fill_buckets(em, cs, ce, bk_max + c, t_max);
		for (int k = cs; k < ce; k++)
			n_short[k + 1] += b;
	}

	moved_ctx x = {
		.bc = raw == NULL ? INTEGER(b_code_) : NULL,
		.bs = raw == NULL ? INTEGER(b_start_) : NULL,
		.bw = raw == NULL ? INTEGER(b_width_) : NULL,
		.shift = raw == NULL ? INTEGER(shift_) : NULL, .raw = raw,
		.tc = INTEGER(t_code_), .seqlen = INTEGER(seqlen_),
		.bounded = LOGICAL(bounded_), .ys = ys, .ye = ye, .yw = yw, .em = em,
		.n_short = n_short, .seg = seg, .o = o, .t_start = t_start,
		.t_max = t_max, .bk_start = bk_start, .bk_max = bk_max, .xseg = xseg,
		.xs = xs, .xe = xe, .mode = mode, .n_lev = n_lev, .n_t = n_t,
		.n_t_rep = n_t * (R_xlen_t) n_rep, .n_blocks = n_blocks,
		.first = (int *) m_alloc(arena, n, sizeof(int)),
		.last = (int *) m_alloc(arena, n, sizeof(int)),
		.plain = (unsigned char *) m_alloc(arena, n, 1), .rep_at = NULL,
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

	pl->x = x;
	pl->n = n;
	pl->n_blocks = n_blocks;
	pl->n_rep = n_rep;
	pl->n_lev = n_lev;
	pl->n_thr = n_thr;
	pl->nc_ok = nc_ok;
	/* the count test over the seqlevels a replicate's blocks counted when
	 * there are no more of them than blocks */
	pl->by_lev = (R_xlen_t) n_lev <= n_blocks;
	pl->with_block = with_block;
	pl->with_cols = with_cols;
	pl->with_code = with_code;
	pl->gather = gather;
	pl->ncol = ncol;
	pl->circ = LOGICAL(circular_);
	pl->cols_ = cols_;
	return 0;
}

/* bait_overlaps(): every block of replicate r (0-based) on
 * sorted_overlaps()' path, i.e. its replicate has fewer blocks on its
 * sequence than y has ranges there, or on y's NCList where that gives the
 * same hits in the same order (y has no ranges there, or they are in NCList
 * order); sorted_overlaps(): none on a circular sequence. The code, the circular check and the count test; pass 1 makes
 * the checks of each block's start, width and target. 0 where moved_impl()
 * returns R_NilValue. cnt holds n_lev zeros, as it is left on success. */
static int moved_rep_check(const moved_plan *pl, R_xlen_t r, int *cnt)
{
	const bait_raw *raw = pl->x.raw;
	const int *bc = pl->x.bc, *seg = pl->x.seg, *circ = pl->circ;
	int n_lev = pl->n_lev;
	R_xlen_t b0 = r * pl->n_blocks, b1 = b0 + pl->n_blocks;
	for (R_xlen_t b = b0; b < b1; b++) {
		int c = raw != NULL ? raw_code(raw, b) : bc[b];
		if (c == NA_INTEGER || c < 1 || c > n_lev || circ[c - 1] == TRUE)
			return 0;
		cnt[c - 1]++;
	}
	if (pl->by_lev) {
		for (int c = 0; c < n_lev; c++)
			if (cnt[c] != 0 && cnt[c] >= seg[c + 1] - seg[c] &&
				!pl->nc_ok[c])
				return 0;
		memset(cnt, 0, n_lev * sizeof(int));
		return 1;
	}
	for (R_xlen_t b = b0; b < b1; b++) {
		int c = raw != NULL ? raw_code(raw, b) : bc[b];
		if (cnt[c - 1] >= seg[c] - seg[c - 1] && !pl->nc_ok[c - 1])
			return 0;
	}
	for (R_xlen_t b = b0; b < b1; b++)
		cnt[(raw != NULL ? raw_code(raw, b) : bc[b]) - 1] = 0;
	return 1;
}

static void moved_runs_lens(SEXP ans, const moved_ctx *xp, int n_chunk,
	const R_xlen_t *cb, const R_xlen_t *c_run, int *lens, int n_rep,
	R_xlen_t n_keep);

/* After pass 1 over the chunks of blocks [cb[c], cb[c + 1]), which gave
 * c_ok, c_hit and c_keep: pass 2, writing each chunk from the sum of the
 * counts before it, so the rows land as one pass in block order puts
 * them, and the value moved_impl() returns */
static SEXP moved_finish(moved_plan *pl, int n_chunk, const R_xlen_t *cb,
	const R_xlen_t *c_hit, R_xlen_t *c_keep, const int *c_ok)
{
	moved_ctx x = pl->x;
	int n_rep = pl->n_rep, with_block = pl->with_block;
	int gather = pl->gather;
	R_xlen_t ncol = pl->ncol;
	SEXP cols_ = pl->cols_;
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
	int n_el = pl->with_cols ? 8 : 7;
	SEXP ans = PROTECT(allocVector(VECSXP, n_el));
	SEXP nms = PROTECT(allocVector(STRSXP, n_el));
	const char *names[8] = {"s", "code", "start", "width",
		with_block ? "block" : "lens", "run_code", "run_length", "cols"};
	int *out[5];
	out[0] = NULL;
	out[1] = NULL;
	for (int j = gather ? 1 : 0; j < 5; j++) {
		if (j == 1 && !pl->with_code)
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
		x.rep_at = lens;
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
	R_xlen_t *c_run = (R_xlen_t *) R_alloc(n_chunk, sizeof(R_xlen_t));
#ifdef _OPENMP
	int n_thr = pl->n_thr;
#pragma omp parallel for num_threads(n_thr) schedule(dynamic, 1) if(n_thr > 1)
#endif
	for (int c = 0; c < n_chunk; c++)
		moved_pass2(&x, cb[c], cb[c + 1], c_keep[c], out, c_run + c);
	moved_runs_lens(ans, &x, n_chunk, cb, c_run, lens, n_rep, n_keep);
	UNPROTECT(2);
	return ans;
}

/* After pass 2 over every chunk: elements 5 and 6 of ans, the runs of the
 * code column, and lens (or NULL) from the row each replicate starts at
 * to the rows of each */
static void moved_runs_lens(SEXP ans, const moved_ctx *xp, int n_chunk,
	const R_xlen_t *cb, const R_xlen_t *c_run, int *lens, int n_rep,
	R_xlen_t n_keep)
{
	const moved_ctx x = *xp;
	/* the runs of the code column: each block's rows share its target
	 * sequence, so one run per block that wrote rows, merged with the run
	 * before when the sequence is the same; maximal, as Rle() finds them.
	 * Each chunk's runs are maximal, so only its first can join the run
	 * before. */
	R_xlen_t n_run = 0;
	int prev = NA_INTEGER;
	for (int c = 0; c < n_chunk; c++) {
		const int *code = x.first + cb[c];
		for (R_xlen_t k = 0; k < c_run[c]; k++) {
			if (code[k] != prev) {
				n_run++;
				prev = code[k];
			}
		}
	}
	SET_VECTOR_ELT(ans, 5, allocVector(INTSXP, n_run));
	SET_VECTOR_ELT(ans, 6, allocVector(INTSXP, n_run));
	int *rc = INTEGER(VECTOR_ELT(ans, 5)), *rl = INTEGER(VECTOR_ELT(ans, 6));
	R_xlen_t r = -1;
	prev = NA_INTEGER;
	for (int c = 0; c < n_chunk; c++) {
		const int *code = x.first + cb[c];
		const int *rows = x.last + cb[c];
		for (R_xlen_t k = 0; k < c_run[c]; k++) {
			if (code[k] != prev) {
				prev = code[k];
				rc[++r] = prev;
				rl[r] = 0;
			}
			rl[r] += rows[k];
		}
	}
	/* each replicate's rows, from the row it starts at; their sum, n_keep,
	 * fits an int */
	if (lens != NULL) {
		for (int rr = 0; rr < n_rep; rr++)
			lens[rr] = (rr + 1 < n_rep ? lens[rr + 1] : (int) n_keep) -
				lens[rr];
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
	SEXP cols_, int with_code, const bait_raw *raw, int n_threads)
{
	moved_plan pl;
	int st = moved_setup(&pl, b_code_, b_start_, b_width_, n_rep_, y_code_,
		y_start_, y_end_, y_width_, circular_, seqlen_, bounded_, shift_,
		t_code_, mode, xseg_, xs_, xe_, with_block, with_cols, cols_,
		with_code, raw, NULL, n_threads);
	if (st != 0)
		return st == 1 ? R_NilValue : ScalarLogical(FALSE);
	int n_lev = pl.n_lev;
	int *cnt = (int *) R_alloc(n_lev > 0 ? n_lev : 1, sizeof(int));
	memset(cnt, 0, (n_lev > 0 ? n_lev : 1) * sizeof(int));
	for (R_xlen_t r = 0; r < pl.n_rep; r++)
		if (!moved_rep_check(&pl, r, cnt))
			return R_NilValue;

	/* the blocks in chunks of consecutive blocks: pass 1 counts each
	 * chunk's rows, and pass 2 writes them from the prefix sum of the
	 * counts. The answer for an exclude interval does not depend on the
	 * previous row's, which each chunk starts without. */
	int n_thr = pl.n_thr;
	R_xlen_t n = pl.n;
	int n_chunk = n_thr == 1 ? 1 : 8 * n_thr;
	R_xlen_t *cb = (R_xlen_t *) R_alloc(n_chunk + 1, sizeof(R_xlen_t));
	for (int c = 0; c <= n_chunk; c++)
		cb[c] = n * c / n_chunk;
	R_xlen_t *c_hit = (R_xlen_t *) R_alloc(n_chunk, sizeof(R_xlen_t));
	R_xlen_t *c_keep = (R_xlen_t *) R_alloc(n_chunk, sizeof(R_xlen_t));
	int *c_ok = (int *) R_alloc(n_chunk, sizeof(int));

	/* pass 1 */
#ifdef _OPENMP
#pragma omp parallel for num_threads(n_thr) schedule(dynamic, 1) if(n_thr > 1)
#endif
	for (int c = 0; c < n_chunk; c++)
		c_ok[c] = moved_pass1(&pl.x, cb[c], cb[c + 1], c_hit + c, c_keep + c);
	return moved_finish(&pl, n_chunk, cb, c_hit, c_keep, c_ok);
}

/* The thread count a _par entry point was given: 1 (serial) unless one
 * integer, not NA; the caller passes omp_threads() (R/zzz.R), which is 1 in
 * a forked child, where libgomp's threads from the parent do not exist and
 * a parallel region would wait for them */
static int threads_arg(SEXP n_threads_)
{
	if (TYPEOF(n_threads_) == INTSXP && XLENGTH(n_threads_) == 1 &&
		INTEGER(n_threads_)[0] != NA_INTEGER)
		return INTEGER(n_threads_)[0];
	return 1;
}

/* moved_impl() without an exclude step, on one thread */
SEXP moved_positions(SEXP b_code_, SEXP b_start_, SEXP b_width_, SEXP n_rep_,
	SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_, SEXP circular_,
	SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_)
{
	return moved_impl(b_code_, b_start_, b_width_, n_rep_, y_code_, y_start_,
		y_end_, y_width_, circular_, seqlen_, bounded_, shift_, t_code_, 0,
		R_NilValue, R_NilValue, R_NilValue, 1, 0, R_NilValue, 1, NULL, 1);
}

static SEXP moved_rows_impl(SEXP b_code_, SEXP b_start_, SEXP b_width_,
	SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, int with_block, int with_cols,
	SEXP cols_, int with_code, const bait_raw *raw, int n_threads)
{
	if (TYPEOF(mode_) != INTSXP || XLENGTH(mode_) != 1)
		return ScalarLogical(FALSE);
	int mode = INTEGER(mode_)[0];
	if (mode != 0 && mode != 1 && mode != 2)
		return ScalarLogical(FALSE);
	return moved_impl(b_code_, b_start_, b_width_, n_rep_, y_code_, y_start_,
		y_end_, y_width_, circular_, seqlen_, bounded_, shift_, t_code_, mode,
		xseg_, xs_, xe_, with_block, with_cols, cols_, with_code, raw,
		n_threads);
}

/* moved_impl() with the exclude step of mode_ (0 none, 1 drop, 2 trim), on
 * one thread */
SEXP moved_rows(SEXP b_code_, SEXP b_start_, SEXP b_width_, SEXP n_rep_,
	SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_, SEXP circular_,
	SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_, SEXP mode_,
	SEXP xseg_, SEXP xs_, SEXP xe_)
{
	return moved_rows_impl(b_code_, b_start_, b_width_, n_rep_, y_code_,
		y_start_, y_end_, y_width_, circular_, seqlen_, bounded_, shift_,
		t_code_, mode_, xseg_, xs_, xe_, 1, 0, R_NilValue, 1, NULL, 1);
}

/* moved_rows() with the rows per replicate ('lens') in place of 'block',
 * gathering the columns 'cols' (see moved_impl()), on one thread */
SEXP moved_rows_lens(SEXP b_code_, SEXP b_start_, SEXP b_width_,
	SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_)
{
	return moved_rows_impl(b_code_, b_start_, b_width_, n_rep_, y_code_,
		y_start_, y_end_, y_width_, circular_, seqlen_, bounded_, shift_,
		t_code_, mode_, xseg_, xs_, xe_, 0, 1, cols_, 1, NULL, 1);
}

/* moved_rows_lens() whose per-row code column is written only when
 * with_code_ is not FALSE; otherwise 'code' is NULL and only its runs
 * ('run_code', 'run_length') are returned. On one thread; moved_rows_opt_par()
 * on up to n_threads_. */
SEXP moved_rows_opt(SEXP b_code_, SEXP b_start_, SEXP b_width_,
	SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_, SEXP with_code_)
{
	return moved_rows_impl(b_code_, b_start_, b_width_, n_rep_, y_code_,
		y_start_, y_end_, y_width_, circular_, seqlen_, bounded_, shift_,
		t_code_, mode_, xseg_, xs_, xe_, 0, 1, cols_,
		asLogical(with_code_) != FALSE, NULL, 1);
}

SEXP moved_rows_opt_par(SEXP b_code_, SEXP b_start_, SEXP b_width_,
	SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_, SEXP with_code_,
	SEXP n_threads_)
{
	return moved_rows_impl(b_code_, b_start_, b_width_, n_rep_, y_code_,
		y_start_, y_end_, y_width_, circular_, seqlen_, bounded_, shift_,
		t_code_, mode_, xseg_, xs_, xe_, 0, 1, cols_,
		asLogical(with_code_) != FALSE, NULL, threads_arg(n_threads_));
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
 * returns R_NilValue; FALSE where it returns FALSE. On up to n_threads
 * threads. */
static SEXP moved_bait_impl(SEXP b_seg_, SEXP seg_code_, SEXP b_start_,
	SEXP b_width_, SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_,
	SEXP y_width_, SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP t_start_,
	SEXP t_code_, SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_,
	SEXP with_code_, int n_threads)
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
		asLogical(with_code_) != FALSE, &r, n_threads);
}

/* moved_bait_impl() on one thread */
SEXP moved_rows_bait(SEXP b_seg_, SEXP seg_code_, SEXP b_start_,
	SEXP b_width_, SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_,
	SEXP y_width_, SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP t_start_,
	SEXP t_code_, SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_,
	SEXP with_code_)
{
	return moved_bait_impl(b_seg_, seg_code_, b_start_, b_width_, n_rep_,
		y_code_, y_start_, y_end_, y_width_, circular_, seqlen_, bounded_,
		t_start_, t_code_, mode_, xseg_, xs_, xe_, cols_, with_code_, 1);
}

/* moved_bait_impl() on up to n_threads_ threads */
SEXP moved_rows_bait_par(SEXP b_seg_, SEXP seg_code_, SEXP b_start_,
	SEXP b_width_, SEXP n_rep_, SEXP y_code_, SEXP y_start_, SEXP y_end_,
	SEXP y_width_, SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP t_start_,
	SEXP t_code_, SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_,
	SEXP with_code_, SEXP n_threads_)
{
	return moved_bait_impl(b_seg_, seg_code_, b_start_, b_width_, n_rep_,
		y_code_, y_start_, y_end_, y_width_, circular_, seqlen_, bounded_,
		t_start_, t_code_, mode_, xseg_, xs_, xe_, cols_, with_code_,
		threads_arg(n_threads_));
}

/* seg_draws_moved()' work before GetRNGstate(), run under
 * R_tryCatchError(): moved_setup() given the draws' vectors, which it
 * does not read, and the buffers of the units, all in the arena */
typedef struct {
	moved_plan *pl;
	m_arena *arena;
	const bait_raw *raw;
	SEXP n_rep_, y_code_, y_start_, y_end_, y_width_, circular_, seqlen_;
	SEXP bounded_, t_code_, xseg_, xs_, xe_, cols_;
	int mode, with_code, n_threads, q, n_unit, status;
	R_xlen_t *cb, *c_hit, *c_keep;
	int *c_ok, *rep_state, *cnt;
} fuse_setup;

static SEXP fuse_setup_body(void *d)
{
	fuse_setup *f = (fuse_setup *) d;
	f->status = moved_setup(f->pl, R_NilValue, R_NilValue, R_NilValue,
		f->n_rep_, f->y_code_, f->y_start_, f->y_end_, f->y_width_,
		f->circular_, f->seqlen_, f->bounded_, R_NilValue, f->t_code_, f->mode,
		f->xseg_, f->xs_, f->xe_, 0, 1, f->cols_, f->with_code, f->raw,
		f->arena, f->n_threads);
	if (f->status != 0)
		return R_NilValue;
	moved_plan *pl = f->pl;
	R_xlen_t nb = pl->n_blocks;
	int n_rep = pl->n_rep, n_lev1 = pl->n_lev > 0 ? pl->n_lev : 1;
	/* units of consecutive blocks within one replicate: q per replicate,
	 * at least 8 per thread in all where there are blocks enough */
	int q = 1;
	if (pl->n_thr > 1) {
		int want = 8 * pl->n_thr;
		q = (want + n_rep - 1) / n_rep;
		if (q > nb)
			q = (int) nb;
		if (q < 1)
			q = 1;
	}
	int n_unit = pl->n_thr > 1 ? n_rep * q : 1;
	f->q = q;
	f->n_unit = n_unit;
	f->cb = (R_xlen_t *) m_alloc(f->arena, n_unit + 1, sizeof(R_xlen_t));
	f->c_hit = (R_xlen_t *) m_alloc(f->arena, n_unit, sizeof(R_xlen_t));
	f->c_keep = (R_xlen_t *) m_alloc(f->arena, n_unit, sizeof(R_xlen_t));
	f->c_ok = (int *) m_alloc(f->arena, n_unit, sizeof(int));
	f->rep_state = (int *) m_alloc(f->arena, n_rep, sizeof(int));
	f->cnt = (int *) m_alloc(f->arena, (size_t) pl->n_thr * n_lev1,
		sizeof(int));
	if (pl->n_thr > 1) {
		for (int r = 0; r < n_rep; r++)
			for (int i = 0; i < q; i++)
				f->cb[(R_xlen_t) r * q + i] = r * nb + nb * i / q;
	} else
		f->cb[0] = 0;
	f->cb[n_unit] = pl->n;
	memset(f->c_ok, 0, n_unit * sizeof(int));
	memset(f->rep_state, 0, n_rep * sizeof(int));
	memset(f->cnt, 0, (size_t) pl->n_thr * n_lev1 * sizeof(int));
	return R_NilValue;
}

typedef struct {
	moved_plan *pl;
	int n_chunk;
	const R_xlen_t *cb, *c_hit;
	R_xlen_t *c_keep;
	const int *c_ok;
} fuse_finish;

static SEXP fuse_finish_body(void *d)
{
	fuse_finish *f = (fuse_finish *) d;
	return moved_finish(f->pl, f->n_chunk, f->cb, f->c_hit, f->c_keep,
		f->c_ok);
}

static SEXP fuse_error(SEXP cond, void *d)
{
	*(int *) d = 1;
	return R_NilValue;
}

#if MOVED_P2_REGION
/* Pass 2 in the draws' region: units below the frontier 'front' (the
 * first unit whose pass 1 is not done, done1) have their row offsets in
 * off, and are written in order from 'next' while their rows end at most
 * maxlen rows in; 'stop' once one would not. done2 marks the units
 * written. front, next and stop change only under the critical section. */
typedef struct {
	const moved_ctx *x;
	const R_xlen_t *cb, *c_keep;
	R_xlen_t *off, *c_run, maxlen;
	int *done1, *done2;
	int *const *out;
	int n_unit, front, next, stop;
} p2_sched;

/* One unit of pass 2, where one is ready: 1 when it wrote one */
static int p2_take(p2_sched *s)
{
	int nx, fr, st, d = 1;
#pragma omp atomic read
	st = s->stop;
	if (st)
		return 0;
#pragma omp atomic read
	nx = s->next;
#pragma omp atomic read
	fr = s->front;
	if (nx >= fr) {
		if (fr >= s->n_unit)
			return 0;
#pragma omp atomic read
		d = s->done1[fr];
		if (!d)
			return 0;
	}
	int p = -1;
#pragma omp critical(moved_p2)
	{
		int f = s->front;
		while (f < s->n_unit) {
#pragma omp atomic read
			d = s->done1[f];
			if (!d)
				break;
#pragma omp flush
			s->off[f + 1] = s->off[f] + s->c_keep[f];
			f++;
		}
#pragma omp atomic write
		s->front = f;
		if (!s->stop && s->next < f) {
			if (s->off[s->next + 1] <= s->maxlen) {
				p = s->next;
#pragma omp atomic write
				s->next = p + 1;
			} else {
#pragma omp atomic write
				s->stop = 1;
			}
		}
	}
	if (p < 0)
		return 0;
	moved_pass2(s->x, s->cb[p], s->cb[p + 1], s->off[p], s->out,
		s->c_run + p);
	s->done2[p] = 1;
	return 1;
}

/* 0 once no unit is left for pass 2 in the region */
static int p2_more(p2_sched *s, int *fail)
{
	int f, st, nx;
#pragma omp atomic read
	f = *fail;
#pragma omp atomic read
	st = s->stop;
#pragma omp atomic read
	nx = s->next;
	return !f && !st && nx < s->n_unit;
}

/* Before the region, under R_tryCatchError(): moved_finish()'s value with
 * every vector pass 2 writes allocated resizable at maxlen rows (lens at
 * n_rep), pass 2's context x2 writing into them, and the units' buffers */
typedef struct {
	moved_plan *pl;
	m_arena *arena;
	R_xlen_t maxlen;
	int n_unit;
	moved_ctx x2;
	int *out[5];
	int *lens;
	R_xlen_t *off, *c_run;
	int *done1, *done2;
} fuse_alloc;

static SEXP fuse_alloc_body(void *d)
{
	fuse_alloc *f = (fuse_alloc *) d;
	moved_plan *pl = f->pl;
	int with_block = pl->with_block, gather = pl->gather;
	R_xlen_t ncol = pl->ncol, maxlen = f->maxlen;
	SEXP cols_ = pl->cols_;
	int n_unit = f->n_unit;
	f->off = (R_xlen_t *) m_alloc(f->arena, n_unit + 1, sizeof(R_xlen_t));
	f->c_run = (R_xlen_t *) m_alloc(f->arena, n_unit, sizeof(R_xlen_t));
	f->done1 = (int *) m_alloc(f->arena, n_unit, sizeof(int));
	f->done2 = (int *) m_alloc(f->arena, n_unit, sizeof(int));
	memset(f->done1, 0, n_unit * sizeof(int));
	memset(f->done2, 0, n_unit * sizeof(int));
	f->x2 = pl->x;
	int n_el = pl->with_cols ? 8 : 7;
	SEXP ans = PROTECT(allocVector(VECSXP, n_el));
	SEXP nms = PROTECT(allocVector(STRSXP, n_el));
	const char *names[8] = {"s", "code", "start", "width",
		with_block ? "block" : "lens", "run_code", "run_length", "cols"};
	f->out[0] = NULL;
	f->out[1] = NULL;
	for (int j = gather ? 1 : 0; j < 5; j++) {
		if (j == 1 && !pl->with_code)
			continue;
		SEXP v = j < 4 || with_block ?
			R_allocResizableVector(INTSXP, maxlen) :
			allocVector(INTSXP, (R_xlen_t) pl->n_rep);
		SET_VECTOR_ELT(ans, j, v);
		f->out[j] = INTEGER(v);
	}
	for (int j = 0; j < n_el; j++)
		SET_STRING_ELT(nms, j, mkChar(names[j]));
	setAttrib(ans, R_NamesSymbol, nms);
	f->lens = NULL;
	if (!with_block) {
		f->lens = f->out[4];
		f->out[4] = NULL;
		f->x2.rep_at = f->lens;
	}
	if (gather) {
		SEXP cols = allocVector(VECSXP, ncol);
		SET_VECTOR_ELT(ans, 7, cols);
		setAttrib(cols, R_NamesSymbol, getAttrib(cols_, R_NamesSymbol));
		const void **cin = (const void **) m_alloc(f->arena,
			ncol > 0 ? ncol : 1, sizeof(void *));
		void **cout = (void **) m_alloc(f->arena, ncol > 0 ? ncol : 1,
			sizeof(void *));
		unsigned char *c8 = (unsigned char *) m_alloc(f->arena,
			ncol > 0 ? ncol : 1, 1);
		for (R_xlen_t j = 0; j < ncol; j++) {
			SEXP v = VECTOR_ELT(cols_, j);
			SEXP w = R_allocResizableVector(TYPEOF(v), maxlen);
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
		f->x2.ncol = (int) ncol;
		f->x2.col8 = c8;
		f->x2.col_in = (const void *const *) cin;
		f->x2.col_out = (void *const *) cout;
	}
	UNPROTECT(2);
	return ans;
}

/* After the region, under R_tryCatchError(): moved_finish()'s value from
 * the rows pass 2 wrote in the region, the units it did not write written
 * now, each vector cut to the rows; moved_finish() itself where the rows
 * outgrew the vectors */
typedef struct {
	moved_plan *pl;
	fuse_alloc *fa;
	SEXP ans;
	int n_chunk, stop;
	const R_xlen_t *cb, *c_hit;
	R_xlen_t *c_keep;
	const int *c_ok;
} fuse_finish2;

static SEXP fuse_finish2_body(void *d)
{
	fuse_finish2 *f = (fuse_finish2 *) d;
	moved_plan *pl = f->pl;
	fuse_alloc *fa = f->fa;
	int n_chunk = f->n_chunk, n_thr = pl->n_thr;
	const R_xlen_t *cb = f->cb;
	R_xlen_t n_hit = 0, n_keep = 0;
	for (int c = 0; c < n_chunk; c++) {
		if (!f->c_ok[c])
			return R_NilValue;
		n_hit += f->c_hit[c];
		n_keep += f->c_keep[c];
	}
	if (n_hit == 0 || n_hit > INT_MAX)
		return R_NilValue;
	if (f->stop || n_keep > fa->maxlen) {
		/* pass 2 overwrote pass 1's entries of the blocks it wrote */
#pragma omp parallel for num_threads(n_thr) schedule(dynamic, 1) if(n_thr > 1)
		for (int c = 0; c < n_chunk; c++) {
			R_xlen_t h, k;
			if (fa->done2[c])
				moved_pass1(&pl->x, cb[c], cb[c + 1], &h, &k);
		}
		return moved_finish(pl, n_chunk, cb, f->c_hit, f->c_keep, f->c_ok);
	}
	R_xlen_t *off = fa->off;
	off[0] = 0;
	for (int c = 0; c < n_chunk; c++)
		off[c + 1] = off[c] + f->c_keep[c];
#pragma omp parallel for num_threads(n_thr) schedule(dynamic, 1) if(n_thr > 1)
	for (int c = 0; c < n_chunk; c++)
		if (!fa->done2[c])
			moved_pass2(&fa->x2, cb[c], cb[c + 1], off[c], fa->out,
				fa->c_run + c);
	SEXP ans = f->ans;
	for (int j = 0; j < 5; j++) {
		SEXP v = VECTOR_ELT(ans, j);
		if (v != R_NilValue && (j < 4 || pl->with_block))
			R_resizeVector(v, n_keep);
	}
	if (pl->gather) {
		SEXP cols = VECTOR_ELT(ans, 7);
		for (R_xlen_t j = 0; j < XLENGTH(cols); j++)
			R_resizeVector(VECTOR_ELT(cols, j), n_keep);
	}
	moved_runs_lens(ans, &fa->x2, n_chunk, cb, fa->c_run, fa->lens,
		pl->n_rep, n_keep);
	return ans;
}
#endif

/* seg_draws() and moved_rows_bait() on its draws in one call: pass 1 of
 * moved_impl() (each replicate's checks, then each block's) runs on the
 * other threads on each replicate as soon as it is drawn, while this
 * thread, R's, draws the next, exactly as seg_draws() does. Pass 2 of a
 * unit runs in the same region once pass 1 of it and of every unit before
 * it is done, into vectors fuse_alloc_body() allocates before the region
 * from replicate 0's rows (drawn and passed before it), and otherwise
 * after the draws. The arguments are seg_draws()' and then those of
 * moved_rows_bait() but b_seg_, b_start_ and n_rep_, which the draws
 * give. Returns NULL where seg_draws() does, else list(start, seg, len,
 * moved): the draws as seg_draws() returns them, and moved_rows_bait()'s
 * list on them, or NULL where it gives no list, where a replicate holds
 * other than length(t_start_) blocks, or where an allocation for the
 * moved rows fails. The region has at most n_threads threads. */
static SEXP seg_draws_moved_impl(SEXP s_n_rep, SEXP prob, SEXP rmin,
	SEXP rmax, SEXP gidx, SEXP off, SEXP n_draw, SEXP st, SEXP rearr,
	SEXP seg_code_, SEXP b_width_, SEXP y_code_, SEXP y_start_, SEXP y_end_,
	SEXP y_width_, SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP t_start_,
	SEXP t_code_, SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_,
	SEXP with_code_, int n_threads)
{
	seg_plan sp;
	if (!seg_plan_init(&sp, s_n_rep, prob, rmin, rmax, gidx, off, n_draw,
		st, rearr))
		return R_NilValue;
	SEXP out = PROTECT(seg_plan_out(&sp, 4));
	SEXP nms = PROTECT(allocVector(STRSXP, 4));
	const char *names[4] = {"start", "seg", "len", "moved"};
	for (int k = 0; k < 4; k++)
		SET_STRING_ELT(nms, k, mkChar(names[k]));
	int n_rep = sp.n_rep;

	/* moved_rows_bait()'s and moved_rows_impl()'s checks, with the draws'
	 * starts and segments as b_start_ and b_seg_, all replicates full */
	int fuse = 1;
	int tw = TYPEOF(b_width_), tt = TYPEOF(t_start_);
	if (TYPEOF(seg_code_) != INTSXP || (tw != INTSXP && tw != REALSXP) ||
		(tt != INTSXP && tt != REALSXP) || XLENGTH(seg_code_) > INT_MAX ||
		TYPEOF(mode_) != INTSXP || XLENGTH(mode_) != 1)
		fuse = 0;
	int mode = fuse ? INTEGER(mode_)[0] : 0;
	if (mode != 0 && mode != 1 && mode != 2)
		fuse = 0;
	bait_raw r;
	r.n = sp.total;
	r.n_w = XLENGTH(b_width_);
	r.n_ts = XLENGTH(t_start_);
	if (r.n == 0 || r.n_w == 0 || r.n_ts == 0 || r.n > INT_MAX ||
		sp.per_rep != r.n_ts)
		fuse = 0;
	moved_plan pl;
	fuse_setup fs;
#if MOVED_P2_REGION
	fuse_alloc fa;
	p2_sched ps;
	SEXP p2_ans = R_NilValue;
#endif
	SEXP keep = R_NilValue;
	if (fuse) {
		r.n_seg = (int) XLENGTH(seg_code_);
		r.seg = sp.og;
		r.seg_code = INTEGER_RO(seg_code_);
		r.si = sp.os;
		r.sd = NULL;
		r.wi = tw == INTSXP ? INTEGER_RO(b_width_) : NULL;
		r.wd = tw == REALSXP ? REAL_RO(b_width_) : NULL;
		r.ti = tt == INTSXP ? INTEGER_RO(t_start_) : NULL;
		r.td = tt == REALSXP ? REAL_RO(t_start_) : NULL;
		keep = allocVector(VECSXP, ARENA_SLOTS);
	}
	PROTECT(keep);
	m_arena arena = {keep, 0};
	if (fuse) {
		fs.pl = &pl;
		fs.arena = &arena;
		fs.raw = &r;
		fs.n_rep_ = s_n_rep;
		fs.y_code_ = y_code_;
		fs.y_start_ = y_start_;
		fs.y_end_ = y_end_;
		fs.y_width_ = y_width_;
		fs.circular_ = circular_;
		fs.seqlen_ = seqlen_;
		fs.bounded_ = bounded_;
		fs.t_code_ = t_code_;
		fs.xseg_ = xseg_;
		fs.xs_ = xs_;
		fs.xe_ = xe_;
		fs.cols_ = cols_;
		fs.mode = mode;
		fs.with_code = asLogical(with_code_) != FALSE;
		fs.n_threads = n_threads;
		fs.status = 1;
		int err = 0;
		R_tryCatchError(fuse_setup_body, &fs, fuse_error, &err);
		if (err || fs.status != 0)
			fuse = 0;
	}

	GetRNGstate();
	sp.kind = R_sample_kind();
	int n_thr = fuse ? pl.n_thr : 1;
	if (!fuse || n_thr == 1) {
		for (int rr = 0; rr < n_rep; rr++)
			seg_plan_rep(&sp, rr);
	}
#ifdef _OPENMP
	else {
		R_xlen_t nb = pl.n_blocks;
		int q = fs.q, n_unit = fs.n_unit, n_lev1 = pl.n_lev > 0 ? pl.n_lev : 1;
		int next_u = 0, fail = 0, r0 = 0;
#if MOVED_P2_REGION
		int p2 = 0;
#endif
		int *rep_state = fs.rep_state;
		const R_xlen_t *cb = fs.cb;
#if MOVED_P2_REGION
		/* replicate 0 before the region: its draws, its checks and pass 1
		 * of its units, whose rows size the vectors of pass 2 */
		seg_plan_rep(&sp, 0);
		r0 = 1;
		int ok0 = sp.ol[0] == nb && moved_rep_check(&pl, 0, fs.cnt);
		if (ok0) {
			/* each unit in n_s pieces on the threads; pass 1 gives the
			 * same entries and counts however the blocks are cut */
			int n_s = n_thr, n_p = q * n_s;
			int p_ok[MOVED_P1_PIECES];
			R_xlen_t p_hit[MOVED_P1_PIECES], p_keep[MOVED_P1_PIECES];
			if (n_p > MOVED_P1_PIECES) {
				n_s = 1;
				n_p = q;
			}
#pragma omp parallel for num_threads(n_thr) schedule(static)
			for (int i = 0; i < n_p; i++) {
				int u = i / n_s, k = i % n_s;
				R_xlen_t w = cb[u + 1] - cb[u];
				R_xlen_t a = cb[u] + w * k / n_s, b = cb[u] + w * (k + 1) / n_s;
				p_hit[i] = 0;
				p_keep[i] = 0;
				p_ok[i] = a < b ? moved_pass1(&pl.x, a, b, p_hit + i,
					p_keep + i) : 1;
			}
			for (int u = 0; u < q; u++) {
				int ok = 1;
				R_xlen_t h = 0, k = 0;
				for (int i = u * n_s; i < (u + 1) * n_s; i++) {
					ok = ok && p_ok[i];
					h += p_hit[i];
					k += p_keep[i];
				}
				fs.c_ok[u] = ok;
				fs.c_hit[u] = h;
				fs.c_keep[u] = k;
				ok0 = ok0 && ok;
			}
		}
		if (!ok0) {
			fuse = 0;
			for (int rr = 1; rr < n_rep; rr++)
				seg_plan_rep(&sp, rr);
		} else {
			rep_state[0] = 2;
			next_u = q;
			R_xlen_t k0 = 0;
			for (int u = 0; u < q; u++)
				k0 += fs.c_keep[u];
			double want = ceil(MOVED_P2_SCALE * (double) k0 * n_rep) +
				MOVED_P2_SLACK;
			if (want <= INT_MAX) {
				fa.pl = &pl;
				fa.arena = &arena;
				fa.maxlen = (R_xlen_t) want;
				fa.n_unit = n_unit;
				int err = 0;
				SEXP a = R_tryCatchError(fuse_alloc_body, &fa, fuse_error,
					&err);
				if (!err && TYPEOF(a) == VECSXP &&
					arena.k < ARENA_SLOTS) {
					SET_VECTOR_ELT(arena.keep, arena.k++, a);
					p2_ans = a;
					p2 = 1;
				}
			}
		}
		if (p2) {
			ps.x = &fa.x2;
			ps.cb = cb;
			ps.c_keep = fs.c_keep;
			ps.off = fa.off;
			ps.c_run = fa.c_run;
			ps.maxlen = fa.maxlen;
			ps.done1 = fa.done1;
			ps.done2 = fa.done2;
			ps.out = fa.out;
			ps.n_unit = n_unit;
			ps.front = q;
			ps.next = 0;
			ps.stop = 0;
			fa.off[0] = 0;
			for (int u = 0; u < q; u++) {
				fa.off[u + 1] = fa.off[u] + fs.c_keep[u];
				fa.done1[u] = 1;
			}
		}
#endif
		if (fuse) {
#pragma omp parallel num_threads(n_thr)
			{
				int tid = omp_get_thread_num();
				int *cnt = fs.cnt + (R_xlen_t) tid * n_lev1;
				if (tid == 0) {
					/* R's thread: the draws, published replicate by
					 * replicate; one of another length moves the blocks of
					 * every later one */
					int bad = 0;
					for (int rr = r0; rr < n_rep; rr++) {
						seg_plan_rep(&sp, rr);
						if (sp.ol[rr] != nb)
							bad = 1;
#pragma omp flush
#pragma omp atomic write
						rep_state[rr] = bad ? -1 : 1;
#pragma omp flush
					}
				}
				/* the units in order: a replicate's first checks it (state
				 * 1 -> 2), the others wait for that; between them, and while
				 * waiting, pass 2 of the units below the frontier */
				for (;;) {
					int u, f, s;
#if MOVED_P2_REGION
					if (p2 && p2_take(&ps))
						continue;
#endif
#pragma omp atomic capture
					u = next_u++;
					if (u >= n_unit)
						break;
#pragma omp atomic read
					f = fail;
					if (f)
						continue;
					int rr = u / q, need = u % q == 0 ? 1 : 2;
					for (;;) {
#pragma omp flush
#pragma omp atomic read
						s = rep_state[rr];
#pragma omp atomic read
						f = fail;
						if (s >= need || s < 0 || f)
							break;
#if MOVED_P2_REGION
						if (p2)
							p2_take(&ps);
#endif
					}
#pragma omp flush
					if (s < need) {
#pragma omp atomic write
						fail = 1;
						continue;
					}
					if (need == 1) {
						int ok = moved_rep_check(&pl, rr, cnt);
#pragma omp flush
#pragma omp atomic write
						rep_state[rr] = ok ? 2 : -1;
#pragma omp flush
						if (!ok) {
#pragma omp atomic write
							fail = 1;
							continue;
						}
					}
					fs.c_ok[u] = moved_pass1(&pl.x, cb[u], cb[u + 1],
						fs.c_hit + u, fs.c_keep + u);
					if (!fs.c_ok[u]) {
#pragma omp atomic write
						fail = 1;
					}
#if MOVED_P2_REGION
					else if (p2) {
#pragma omp flush
#pragma omp atomic write
						ps.done1[u] = 1;
					}
#endif
				}
#if MOVED_P2_REGION
				/* every pass-1 unit taken: pass 2 of the rest as their
				 * pass 1 ends */
				while (p2 && p2_more(&ps, &fail))
					p2_take(&ps);
#endif
			}
			if (fail)
				fuse = 0;
		}
	}
#endif
	PutRNGstate();
	seg_plan_trim(&sp, out);

	if (fuse && n_thr == 1) {
		/* moved_impl() on one thread */
		if (sp.w != sp.total)
			fuse = 0;
		for (R_xlen_t rr = 0; fuse && rr < n_rep; rr++)
			if (!moved_rep_check(&pl, rr, fs.cnt))
				fuse = 0;
		if (fuse)
			fs.c_ok[0] = moved_pass1(&pl.x, 0, pl.n, fs.c_hit, fs.c_keep);
	}
	SEXP moved = R_NilValue;
	if (fuse) {
		int err = 0;
#if MOVED_P2_REGION
		if (p2_ans != R_NilValue) {
			fuse_finish2 f2 = {&pl, &fa, p2_ans, fs.n_unit, ps.stop, fs.cb,
				fs.c_hit, fs.c_keep, fs.c_ok};
			moved = R_tryCatchError(fuse_finish2_body, &f2, fuse_error, &err);
		} else
#endif
		{
			fuse_finish ff = {&pl, fs.n_unit, fs.cb, fs.c_hit, fs.c_keep,
				fs.c_ok};
			moved = R_tryCatchError(fuse_finish_body, &ff, fuse_error, &err);
		}
		if (err || TYPEOF(moved) != VECSXP)
			moved = R_NilValue;
	}
	SET_VECTOR_ELT(out, 3, moved);
	setAttrib(out, R_NamesSymbol, nms);
	UNPROTECT(3);
	return out;
}

/* seg_draws_moved_impl() on one thread */
SEXP seg_draws_moved(SEXP s_n_rep, SEXP prob, SEXP rmin, SEXP rmax,
	SEXP gidx, SEXP off, SEXP n_draw, SEXP st, SEXP rearr, SEXP seg_code_,
	SEXP b_width_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP t_start_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_, SEXP with_code_)
{
	return seg_draws_moved_impl(s_n_rep, prob, rmin, rmax, gidx, off, n_draw,
		st, rearr, seg_code_, b_width_, y_code_, y_start_, y_end_, y_width_,
		circular_, seqlen_, bounded_, t_start_, t_code_, mode_, xseg_, xs_,
		xe_, cols_, with_code_, 1);
}

/* seg_draws_moved_impl() on up to n_threads_ threads */
SEXP seg_draws_moved_par(SEXP s_n_rep, SEXP prob, SEXP rmin, SEXP rmax,
	SEXP gidx, SEXP off, SEXP n_draw, SEXP st, SEXP rearr, SEXP seg_code_,
	SEXP b_width_, SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_,
	SEXP circular_, SEXP seqlen_, SEXP bounded_, SEXP t_start_, SEXP t_code_,
	SEXP mode_, SEXP xseg_, SEXP xs_, SEXP xe_, SEXP cols_, SEXP with_code_,
	SEXP n_threads_)
{
	return seg_draws_moved_impl(s_n_rep, prob, rmin, rmax, gidx, off, n_draw,
		st, rearr, seg_code_, b_width_, y_code_, y_start_, y_end_, y_width_,
		circular_, seqlen_, bounded_, t_start_, t_code_, mode_, xseg_, xs_,
		xe_, cols_, with_code_, threads_arg(n_threads_));
}

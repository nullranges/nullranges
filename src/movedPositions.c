#include <R.h>
#include <Rinternals.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>

/* R's integer + and - give NA outside [-INT_MAX, INT_MAX] */
#define INT_OK(v) ((v) <= INT_MAX && (v) >= -(int64_t) INT_MAX)

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

/* The moved ranges of every bait block in one pass: per block, in order,
 * the ranges of y it overlaps, by increasing index in y, with their index
 * in y ('s'), target seqlevel code, shifted and trimmed start and width,
 * and replicate ('block'), keeping width > 0. R_NilValue wherever the R
 * route (bait_overlaps(), shift_swap_trim_positions()) would not take
 * sorted_overlaps() for every block or could decline, or has no hit. */
SEXP moved_positions(SEXP b_code_, SEXP b_start_, SEXP b_width_, SEXP n_rep_,
	SEXP y_code_, SEXP y_start_, SEXP y_end_, SEXP y_width_, SEXP circular_,
	SEXP seqlen_, SEXP bounded_, SEXP shift_, SEXP t_code_)
{
	if (TYPEOF(b_code_) != INTSXP || TYPEOF(b_start_) != INTSXP ||
		TYPEOF(b_width_) != INTSXP || TYPEOF(n_rep_) != INTSXP ||
		TYPEOF(y_code_) != INTSXP || TYPEOF(y_start_) != INTSXP ||
		TYPEOF(y_end_) != INTSXP || TYPEOF(y_width_) != INTSXP ||
		TYPEOF(circular_) != LGLSXP || TYPEOF(seqlen_) != INTSXP ||
		TYPEOF(bounded_) != LGLSXP || TYPEOF(shift_) != INTSXP ||
		TYPEOF(t_code_) != INTSXP || XLENGTH(n_rep_) != 1)
		return R_NilValue;
	R_xlen_t n = XLENGTH(b_code_);
	R_xlen_t n_y = XLENGTH(y_code_);
	R_xlen_t n_t = XLENGTH(t_code_);
	int n_rep = INTEGER(n_rep_)[0];
	int n_lev = LENGTH(circular_);
	if (n == 0 || n_rep == NA_INTEGER || n_rep < 1 || n % n_rep != 0 ||
		n > INT_MAX || n_y > INT_MAX || XLENGTH(b_start_) != n ||
		XLENGTH(b_width_) != n || XLENGTH(shift_) != n ||
		XLENGTH(y_start_) != n_y || XLENGTH(y_end_) != n_y ||
		XLENGTH(y_width_) != n_y || XLENGTH(seqlen_) != n_lev ||
		XLENGTH(bounded_) != n_lev || (int64_t) n_rep * n_lev > INT_MAX)
		return R_NilValue;
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
	 * end(bait), start + width - 1L, not NA */
	const int *bc = INTEGER(b_code_), *bs = INTEGER(b_start_);
	const int *bw = INTEGER(b_width_), *circ = LOGICAL(circular_);
	int *cnt = (int *) R_alloc(n_lev > 0 ? n_lev : 1, sizeof(int));
	memset(cnt, 0, (n_lev > 0 ? n_lev : 1) * sizeof(int));
	for (R_xlen_t r = 0; r < n_rep; r++) {
		R_xlen_t b0 = r * n_blocks, b1 = b0 + n_blocks;
		for (R_xlen_t b = b0; b < b1; b++) {
			int c = bc[b];
			int64_t e = (int64_t) bs[b] + bw[b];
			if (c == NA_INTEGER || c < 1 || c > n_lev || !INT_OK(e) ||
				!INT_OK(e - 1) || e - 1 < bs[b] || circ[c - 1] == TRUE)
				return R_NilValue;
			cnt[c - 1]++;
		}
		for (R_xlen_t b = b0; b < b1; b++)
			if (cnt[bc[b] - 1] >= seg[bc[b]] - seg[bc[b] - 1])
				return R_NilValue;
		for (R_xlen_t b = b0; b < b1; b++)
			cnt[bc[b] - 1] = 0;
	}

	buckets *bk_start = (buckets *) R_alloc(n_lev, sizeof(buckets));
	buckets *bk_max = (buckets *) R_alloc(n_lev, sizeof(buckets));
	int *t_start = (int *) R_alloc(plan_buckets(ys, seg, n_lev, bk_start),
		sizeof(int));
	int *t_max = (int *) R_alloc(plan_buckets(em, seg, n_lev, bk_max),
		sizeof(int));
	fill_buckets(ys, seg, n_lev, bk_start, t_start);
	fill_buckets(em, seg, n_lev, bk_max, t_max);

	/* pass 1: each block's candidates, the checks R makes over the hits,
	 * and the number of rows kept */
	const int *shift = INTEGER(shift_), *tc = INTEGER(t_code_);
	const int *seqlen = INTEGER(seqlen_), *bounded = LOGICAL(bounded_);
	int *first = (int *) R_alloc(n, sizeof(int));
	int *last = (int *) R_alloc(n, sizeof(int));
	/* a block whose candidates are all hits, moved without overflow or
	 * trimming: every row kept as it is */
	unsigned char *plain = (unsigned char *) R_alloc(n, 1);
	R_xlen_t n_t_rep = n_t * (R_xlen_t) n_rep;
	R_xlen_t n_hit = 0, n_keep = 0;
	for (R_xlen_t b = 0; b < n; b++) {
		int c = bc[b] - 1, cs = seg[c], ce = seg[c + 1];
		int lo = count_le(em, cs, ce, bk_max + c, t_max, (int64_t) bs[b] - 1);
		int hi = count_le(ys, cs, ce, bk_start + c, t_start,
			(int64_t) bs[b] + bw[b] - 1);
		first[b] = lo;
		last[b] = hi;
		plain[b] = 0;
		if (lo >= hi)
			continue;
		int sh = shift[b];
		int code = b < n_t_rep ? tc[b % n_t] : NA_INTEGER;
		int bstart = bs[b];
		if (n_short[hi] == n_short[lo] && sh != NA_INTEGER &&
			code != NA_INTEGER && code >= 1 && code <= n_lev) {
			int64_t s0 = (int64_t) ys[lo] + sh, e1 = (int64_t) em[hi - 1] + sh;
			if (bounded[code - 1] ? (s0 >= 1 && e1 <= seqlen[code - 1]) :
				(INT_OK(s0) && INT_OK(e1))) {
				plain[b] = 1;
				n_hit += hi - lo;
				n_keep += hi - lo;
				continue;
			}
		}
		for (int k = lo; k < hi; k++) {
			if (ye[k] < bstart)
				continue;
			n_hit++;
			int ns, nw;
			if (sh == NA_INTEGER || code == NA_INTEGER || code < 1 ||
				code > n_lev ||
				!move_one(ys[k], yw[k], sh, code, seqlen, bounded, &ns, &nw))
				return R_NilValue;
			n_keep += nw > 0;
		}
	}
	if (n_hit == 0 || n_hit > INT_MAX)
		return R_NilValue;

	/* pass 2: the rows */
	SEXP ans = PROTECT(allocVector(VECSXP, 5));
	SEXP nms = PROTECT(allocVector(STRSXP, 5));
	const char *names[5] = {"s", "code", "start", "width", "block"};
	int *out[5];
	for (int j = 0; j < 5; j++) {
		SET_VECTOR_ELT(ans, j, allocVector(INTSXP, n_keep));
		SET_STRING_ELT(nms, j, mkChar(names[j]));
		out[j] = INTEGER(VECTOR_ELT(ans, j));
	}
	setAttrib(ans, R_NamesSymbol, nms);
	R_xlen_t at = 0;
	for (R_xlen_t b = 0; b < n; b++) {
		int lo = first[b], hi = last[b];
		if (lo >= hi)
			continue;
		int sh = shift[b];
		int code = tc[b % n_t];
		int bstart = bs[b];
		int rep = (int) (b / n_blocks) + 1;
		if (plain[b]) {
			for (int k = lo; k < hi; k++, at++) {
				out[0][at] = (o == NULL ? k : o[k]) + 1;
				out[1][at] = code;
				out[2][at] = ys[k] + sh;
				out[3][at] = yw[k];
				out[4][at] = rep;
			}
			continue;
		}
		for (int k = lo; k < hi; k++) {
			if (ye[k] < bstart)
				continue;
			int ns, nw;
			move_one(ys[k], yw[k], sh, code, seqlen, bounded, &ns, &nw);
			if (nw <= 0)
				continue;
			out[0][at] = (o == NULL ? k : o[k]) + 1;
			out[1][at] = code;
			out[2][at] = ns;
			out[3][at] = nw;
			out[4][at] = rep;
			at++;
		}
	}
	UNPROTECT(2);
	return ans;
}

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Random.h>
#include <float.h>
#include <math.h>
#include <limits.h>
#include <string.h>

/* buckets per distinct focal score in the index over u */
#define NN_BUCKETS 16

/* finite: not NaN, NA or infinite, without a call into R */
#define NN_FINITE(x) (fabs(x) <= DBL_MAX)

/* An index over the m sorted distinct focal scores u: pool score x lies in
 * interval k, u[k - 1] < x <= u[k] (k = 0 at or below u[0], k = m above
 * u[m - 1]). start[b] is the number of u whose bucket is below b, so the
 * interval of an x in bucket b is at least start[b] (the bucket is monotone
 * in x). */
typedef struct {
	int m, B, use_b;
	double lo, hi, scale;
	int *start;
} nn_index;

static void nn_index_build(nn_index *ix, const double *u, int m)
{
	double lo = u[0], hi = u[m - 1], span = hi - lo;
	int B = m > 1 ? NN_BUCKETS * m : 1;
	double scale = m > 1 ? (double) B / span : 0;
	int use_b = m > 1 && R_FINITE(span) && span > 0 && R_FINITE(scale);
	int *start = NULL;
	if (use_b) {
		start = (int *) R_alloc((size_t) B + 1, sizeof(int));
		memset(start, 0, ((size_t) B + 1) * sizeof(int));
		for (int j = 0; j < m; j++) {
			int b = (int) ((u[j] - lo) * scale);
			if (b > B - 1)
				b = B - 1;
			start[b + 1]++;
		}
		for (int b = 1; b <= B; b++)
			start[b] += start[b - 1];
	}
	ix->m = m;
	ix->B = B;
	ix->use_b = use_b;
	ix->lo = lo;
	ix->hi = hi;
	ix->scale = scale;
	ix->start = start;
}

/* the interval of a finite x */
static inline int nn_interval(double x, const double *u, int m, double lo,
	double hi, int use_b, const int *start, int B, double scale)
{
	if (x <= lo)
		return 0;
	if (x > hi)
		return m;
	if (use_b) {
		int b = (int) ((x - lo) * scale);
		if (b > B - 1)
			b = B - 1;
		/* two steps without a branch, as most scans are short; k stays
		 * below m, since u[m - 1] = hi >= x */
		int k = start[b];
		k += u[k] < x;
		k += u[k] < x;
		while (u[k] < x)
			k++;
		return k;
	}
	/* the first k with u[k] >= x; u[m - 1] = hi >= x */
	int a = 0, z = m - 1;
	while (a < z) {
		int mid = a + (z - a) / 2;
		if (u[mid] < x)
			a = mid + 1;
		else
			z = mid;
	}
	return a;
}

/* nnMatch's candidate pool rows. The rolling nearest join gives each
 * distinct focal score u the rows of the largest pool score <= u or of the
 * smallest >= u, and each of those is its interval's max or min. Returns
 * the 1-based pool rows whose score equals its interval's min or max, in
 * increasing order. NULL when fps or pps is not a double vector, either is
 * empty or longer than INT_MAX, or either holds a value that is not
 * finite. */
SEXP nn_keep(SEXP fps, SEXP pps)
{
	if (TYPEOF(fps) != REALSXP || TYPEOF(pps) != REALSXP)
		return R_NilValue;
	R_xlen_t nf_ = XLENGTH(fps), np_ = XLENGTH(pps);
	if (nf_ < 1 || np_ < 1 || nf_ > INT_MAX / (2 * NN_BUCKETS) ||
	    np_ > INT_MAX)
		return R_NilValue;
	int nf = (int) nf_, np = (int) np_;
	const double *f = REAL_RO(fps), *p = REAL_RO(pps);

	double *u = (double *) R_alloc(nf, sizeof(double));
	for (int i = 0; i < nf; i++) {
		if (!NN_FINITE(f[i]))
			return R_NilValue;
		u[i] = f[i];
	}
	R_qsort(u, 1, (size_t) nf);
	int m = 1;
	for (int i = 1; i < nf; i++)
		if (u[i] != u[m - 1])
			u[m++] = u[i];

	nn_index ix;
	nn_index_build(&ix, u, m);
	const double lo = ix.lo, hi = ix.hi, scale = ix.scale;
	const int B = ix.B, use_b = ix.use_b;
	const int *start = ix.start;

	int *kk = (int *) R_alloc(np, sizeof(int));
	double *mn = (double *) R_alloc((size_t) m + 1, sizeof(double));
	double *mx = (double *) R_alloc((size_t) m + 1, sizeof(double));
	for (int k = 0; k <= m; k++) {
		mn[k] = R_PosInf;
		mx[k] = R_NegInf;
	}
	for (int i = 0; i < np; i++) {
		double x = p[i];
		if (!NN_FINITE(x))
			return R_NilValue;
		int k = nn_interval(x, u, m, lo, hi, use_b, start, B, scale);
		kk[i] = k;
		double a = mn[k], z = mx[k];
		mn[k] = x < a ? x : a;
		mx[k] = x > z ? x : z;
	}

	/* kept rows written over kk, which is read ahead of the writes */
	int w = 0;
	for (int i = 0; i < np; i++) {
		int k = kk[i];
		if (p[i] == mn[k] || p[i] == mx[k])
			kk[w++] = i + 1;
	}
	SEXP out = PROTECT(allocVector(INTSXP, w));
	if (w)
		memcpy(INTEGER(out), kk, (size_t) w * sizeof(int));
	UNPROTECT(1);
	return out;
}

/* nnMatch with replace = TRUE: the matched pool row of every focal row, as
 * dt[.(uniq.fps), roll = 'nearest', mult = 'all'] on the pool keyed by
 * score (stable, so tied rows in increasing order) and sample_groups give
 * it. For each distinct focal score u, increasing: when pool scores equal
 * u, all of those rows; otherwise one row, the last of the largest pool
 * score lo below u when there is no score above u or when u - lo <= hi - u
 * (in doubles), else the first of the smallest score hi above u. Then one
 * R_unif_index(number of rows) draw per focal row with score u, in
 * increasing focal row order. Returns the 1-based pool rows, or NULL where
 * nn_keep does and when a focal score is 0 (nnMatch's grouped fps column
 * then holds the group's first zero, whose sign may differ). */
SEXP nn_draws(SEXP fps, SEXP pps)
{
	if (TYPEOF(fps) != REALSXP || TYPEOF(pps) != REALSXP)
		return R_NilValue;
	R_xlen_t nf_ = XLENGTH(fps), np_ = XLENGTH(pps);
	if (nf_ < 1 || np_ < 1 || nf_ > INT_MAX / (2 * NN_BUCKETS) ||
	    np_ > INT_MAX)
		return R_NilValue;
	int nf = (int) nf_, np = (int) np_;
	const double *f = REAL_RO(fps), *p = REAL_RO(pps);

	/* the focal scores sorted with their rows, tied rows in increasing
	 * order; then u[j] (over v) the distinct scores, and the rows of u[j]
	 * at ord[gs[j]] to ord[gs[j + 1] - 1] */
	double *v = (double *) R_alloc((size_t) nf + 1, sizeof(double));
	int *ord = (int *) R_alloc(nf, sizeof(int));
	int *gs = (int *) R_alloc((size_t) nf + 1, sizeof(int));
	for (int i = 0; i < nf; i++) {
		if (!NN_FINITE(f[i]) || f[i] == 0)
			return R_NilValue;
		v[i] = f[i];
		ord[i] = i;
	}
	R_qsort_I(v, ord, 1, nf);
	int m = 0;
	for (int i = 0; i < nf;) {
		int e = i + 1;
		while (e < nf && v[e] == v[i])
			e++;
		if (e - i > 1)
			R_isort(ord + i, e - i);
		v[m] = v[i];
		gs[m++] = i;
		i = e;
	}
	gs[m] = nf;
	/* no pool score equals the sentinel above u[m - 1] */
	v[m] = R_PosInf;
	const double *u = v;

	nn_index ix;
	nn_index_build(&ix, u, m);
	const double lo = ix.lo, hi = ix.hi, scale = ix.scale;
	const int B = ix.B, use_b = ix.use_b;
	const int *start = ix.start;

	/* per interval: min and its first row, max and its last row (each
	 * stored through an index that points at the spare slot m + 1 when the
	 * row does not set it, which is cheaper here than a load and a select);
	 * and the rows equal to u[k], listed with k */
	double *mn = (double *) R_alloc((size_t) m + 1, sizeof(double));
	double *mx = (double *) R_alloc((size_t) m + 1, sizeof(double));
	int *first = (int *) R_alloc((size_t) m + 2, sizeof(int));
	int *last = (int *) R_alloc((size_t) m + 2, sizeof(int));
	for (int k = 0; k <= m; k++) {
		mn[k] = R_PosInf;
		mx[k] = R_NegInf;
		first[k] = -1;
		last[k] = -1;
	}
	int *exi = (int *) R_alloc(np, sizeof(int));
	int *exk = (int *) R_alloc(np, sizeof(int));
	int w = 0;
	for (int i = 0; i < np; i++) {
		double x = p[i];
		if (!NN_FINITE(x))
			return R_NilValue;
		int k = nn_interval(x, u, m, lo, hi, use_b, start, B, scale);
		double a = mn[k], z = mx[k];
		int lt = x < a, ge = x >= z;
		mn[k] = lt ? x : a;
		mx[k] = ge ? x : z;
		first[lt ? k : m + 1] = i;
		last[ge ? k : m + 1] = i;
		/* written always, kept when x equals u[k] */
		exi[w] = i;
		exk[w] = k;
		w += x == u[k];
	}

	/* the rows equal to u[k] grouped by k, in increasing row order */
	int *cnt = (int *) R_alloc((size_t) m + 1, sizeof(int));
	memset(cnt, 0, ((size_t) m + 1) * sizeof(int));
	for (int t = 0; t < w; t++)
		cnt[exk[t]]++;
	int *off = (int *) R_alloc((size_t) m + 1, sizeof(int));
	int *pos = (int *) R_alloc((size_t) m + 1, sizeof(int));
	off[0] = 0;
	for (int k = 0; k < m; k++)
		off[k + 1] = off[k] + cnt[k];
	memcpy(pos, off, ((size_t) m + 1) * sizeof(int));
	int *rows = (int *) R_alloc(w > 0 ? (size_t) w : 1, sizeof(int));
	for (int t = 0; t < w; t++)
		rows[pos[exk[t]]++] = exi[t] + 1;

	/* nxt[k]: the first non-empty interval at or above k, -1 if none */
	int *nxt = (int *) R_alloc((size_t) m + 2, sizeof(int));
	nxt[m + 1] = -1;
	for (int k = m; k >= 0; k--)
		nxt[k] = last[k] >= 0 ? k : nxt[k + 1];

	SEXP out = PROTECT(allocVector(INTSXP, nf));
	int *o = INTEGER(out);
	int lo_k = -1;
	GetRNGstate();
	for (int j = 0; j < m; j++) {
		if (last[j] >= 0)
			lo_k = j;
		const int *cand;
		int nc, one;
		if (cnt[j] > 0) {
			cand = rows + off[j];
			nc = cnt[j];
		} else {
			/* lo: interval lo_k's max, below u[j]; hi: interval hi_k's
			 * min, above u[j]; at least one exists, as np >= 1 */
			int hi_k = nxt[j + 1];
			if (lo_k >= 0 && (hi_k < 0 ||
			    u[j] - mx[lo_k] <= mn[hi_k] - u[j]))
				one = last[lo_k] + 1;
			else
				one = first[hi_k] + 1;
			cand = &one;
			nc = 1;
		}
		double dn = (double) nc;
		for (int t = gs[j]; t < gs[j + 1]; t++)
			o[ord[t]] = cand[(R_xlen_t) R_unif_index(dn)];
	}
	PutRNGstate();
	UNPROTECT(1);
	return out;
}

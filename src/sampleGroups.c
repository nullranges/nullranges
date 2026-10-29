#include <R.h>
#include <Rinternals.h>
#include <R_ext/Random.h>

/* Groups are consecutive blocks of `index`, block g holding n_rows[g] rows.
 * From each block in turn, draw size[g] rows with replacement, exactly as
 * x[sample.int(n_rows[g], size[g], replace = TRUE)] does: base R's
 * do_sample draws one R_unif_index(n) per element. */
SEXP sample_groups(SEXP index, SEXP n_rows, SEXP size)
{
	if (!isInteger(index) || !isInteger(n_rows) || !isInteger(size) ||
	    XLENGTH(n_rows) != XLENGTH(size))
		error("sample_groups: invalid arguments");
	R_xlen_t ng = XLENGTH(n_rows), rows = 0, total = 0;
	const int *ix = INTEGER(index), *k = INTEGER(n_rows), *n = INTEGER(size);
	for (R_xlen_t g = 0; g < ng; g++) {
		if (k[g] < 1 || n[g] < 0)
			error("sample_groups: invalid group size");
		rows += k[g];
		total += n[g];
	}
	if (rows != XLENGTH(index))
		error("sample_groups: group sizes do not sum to the index length");

	SEXP out = PROTECT(allocVector(INTSXP, total));
	int *o = INTEGER(out);
	R_xlen_t start = 0, pos = 0;
	GetRNGstate();
	for (R_xlen_t g = 0; g < ng; g++) {
		double dn = (double) k[g];
		for (int j = 0; j < n[g]; j++)
			o[pos++] = ix[start + (R_xlen_t) R_unif_index(dn)];
		start += k[g];
	}
	PutRNGstate();
	UNPROTECT(1);
	return out;
}

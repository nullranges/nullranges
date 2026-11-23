#' @useDynLib nullranges, .registration = TRUE, .fixes = "C_"
NULL

# bootRanges() and matchRanges() have compiled implementations, used when
# their argument compiledCode is TRUE (the default). Each returns NULL where
# it declines a call, and the function then runs its R code, so that the
# draws, warnings and errors are those of that code.

# Runs fun(...) and returns its value. Returns NULL, with the random number
# generator as it was before the call, where the generator is user-supplied
# or where fun signals a warning or an error, so that the R code reproduces
# the condition with its own call, at its own point of the random stream.
#
# @param fun a compiled implementation
# @param ... the arguments of the call
run_or_decline <- function(fun, ...) {
  if (user_rng()) {
    return(NULL)
  }
  mark <- rng_mark()
  out <- tryCatch(
    fun(...),
    warning = function(w) NULL, error = function(e) NULL
  )
  if (is.null(out)) {
    rng_restore(mark)
  }
  out
}

# TRUE when the generator is user-supplied (RNGkind()), whose state need
# not be in .Random.seed, so that rng_restore() could not rewind it
user_rng <- function() {
  seed <- get0(".Random.seed", envir = globalenv(), inherits = FALSE)
  is.integer(seed) && length(seed) >= 1L && !is.na(seed[1L]) &&
    seed[1L] %% 100L == 5L
}

# .Random.seed as it stands, NULL when there is none, for rng_restore()
rng_mark <- function() {
  list(seed = get0(".Random.seed", envir = globalenv(), inherits = FALSE))
}

# Puts back the .Random.seed of rng_mark(): the next draw continues from
# it, or, when there was none, seeds the generator afresh, as the first
# draw of a session does
#
# @param mark the value of rng_mark()
rng_restore <- function(mark) {
  if (!is.null(mark$seed)) {
    assign(".Random.seed", mark$seed, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
}

# Helpers for the compiled bootRanges() and matchRanges(): the metadata
# columns of a GRanges that the compiled passes gather.

# The ranges and metadata columns part of plain_rows()
#
# @param y the input GRanges
plain_ranges_mcols <- function(y) {
  rng <- y@ranges
  mc <- y@elementMetadata
  identical(as.vector(class(rng)), "IRanges") &&
    is.null(rng@NAMES) && is.null(rng@elementMetadata) &&
    identical(as.vector(class(mc)), "DFrame") &&
    is.null(mc@rownames) && is.null(mc@elementMetadata) &&
    all(vapply(mc@listData, function(v) {
      is.atomic(v) && is.null(attributes(v))
    }, NA))
}

# The metadata columns of 'y' (its DFrame's listData) when moved_ranges()
# takes its plain branch for them and each is a double, integer or logical
# vector of y's length without attributes, which the compiled pass gathers
# bit for bit as `[` does; NULL otherwise.
#
# @param y the input GRanges
gather_columns <- function(y) {
  if (!identical(as.vector(class(y)), "GRanges") ||
      nrun(y@strand) != 1L || !plain_ranges_mcols(y)) {
    return(NULL)
  }
  cols <- y@elementMetadata@listData
  n <- length(y)
  if (!all(vapply(cols, function(v) {
    typeof(v) %in% c("double", "integer", "logical") && length(v) == n
  }, NA))) {
    return(NULL)
  }
  cols
}

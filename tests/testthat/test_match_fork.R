library(nullranges)

# matchRanges in a forked child (parallel::mcparallel, as mclapply forks it)
# after the parent has run an OpenMP region. The parent's OpenMP threads do
# not exist in the child, so a parallel region there can wait for them
# forever: each child is collected with a timeout and killed if it has not
# finished, so a hang fails the test instead of stalling the suite.

# 22,000 rows, enough for the propensity fit's threaded loops
fork_match_input <- function() {
  set.seed(40)
  list(focal = data.frame(x = rnorm(2000, 1), y = runif(2000)),
       pool = data.frame(x = rnorm(20000), y = runif(20000)))
}

fork_match <- function(d) {
  force(d)
  set.seed(41)
  matchRanges(d$focal, d$pool, covar = ~x + y,
              method = "nearest", replace = TRUE)
}

# a value as it comes back from a child: serialize() drops the data.table
# self-reference pointer in matchedData, so only such copies are identical()
via_child <- function(x) unserialize(serialize(x, NULL))

# the child's value, or NULL if it did not finish within 'timeout' seconds
in_match_child <- function(f, timeout = 30) {
  job <- parallel::mcparallel(f(), silent = TRUE)
  res <- parallel::mccollect(job, wait = FALSE, timeout = timeout)
  if (is.null(res)) {
    tools::pskill(job$pid, tools::SIGKILL)
    suppressWarnings(parallel::mccollect(job))
    return(NULL)
  }
  res[[1]]
}

test_that("matchRanges finishes in a forked child after a threaded data.table sort", {
  skip_on_os("windows")
  d <- fork_match_input()
  dt <- data.table::data.table(a = runif(1e6))
  data.table::setorder(dt, a)
  got <- in_match_child(function() fork_match(d))
  expect_false(is.null(got))
  expect_identical(got, via_child(fork_match(d)))
})

test_that("matchRanges finishes in a forked child after its own threaded call", {
  skip_on_os("windows")
  d <- fork_match_input()
  ref <- fork_match(d)
  got <- in_match_child(function() fork_match(d))
  expect_false(is.null(got))
  expect_identical(got, via_child(ref))
})

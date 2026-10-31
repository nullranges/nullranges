library(nullranges)
library(GenomicRanges)

# bootRanges in a forked child (parallel::mcparallel, as mclapply forks it)
# after the parent has run an OpenMP region. The parent's OpenMP threads do
# not exist in the child, so a parallel region there can wait for them
# forever: each child is collected with a timeout and killed if it has not
# finished, so a hang fails the test instead of stalling the suite.

fork_input <- function() {
  set.seed(30)
  lens <- c(chr1 = 1e7, chr2 = 1e7, chr3 = 1e7)
  n <- 5000
  starts <- round(runif(n * length(lens), 1, 1e7 - 200))
  y <- GRanges(rep(names(lens), each = n),
               IRanges(starts, width = round(runif(n * length(lens), 50, 200))),
               seqlengths = lens)
  sort(y)
}

# about 3,000 blocks per replicate, enough for bootRanges' threaded path
fork_boot <- function(y) {
  force(y)
  set.seed(31)
  bootRanges(y, blockLength = 1e4, R = 8)
}

# the child's value, or NULL if it did not finish within 'timeout' seconds
in_child <- function(f, timeout = 30) {
  job <- parallel::mcparallel(f(), silent = TRUE)
  res <- parallel::mccollect(job, wait = FALSE, timeout = timeout)
  if (is.null(res)) {
    tools::pskill(job$pid, tools::SIGKILL)
    suppressWarnings(parallel::mccollect(job))
    return(NULL)
  }
  res[[1]]
}

test_that("bootRanges finishes in a forked child after a threaded data.table sort", {
  skip_on_os("windows")
  y <- fork_input()
  dt <- data.table::data.table(a = runif(1e6))
  data.table::setorder(dt, a)
  got <- in_child(function() fork_boot(y))
  expect_false(is.null(got))
  expect_identical(got, fork_boot(y))
})

test_that("bootRanges finishes in a forked child after its own threaded call", {
  skip_on_os("windows")
  y <- fork_input()
  ref <- fork_boot(y)
  got <- in_child(function() fork_boot(y))
  expect_false(is.null(got))
  expect_identical(got, ref)
})

test_that("bootRanges finishes in a forked child that raised data.table's thread count", {
  skip_on_os("windows")
  y <- fork_input()
  ref <- fork_boot(y)
  got <- in_child(function() {
    data.table::setDTthreads(4)
    fork_boot(y)
  })
  expect_false(is.null(got))
  expect_identical(got, ref)
})

test_that("bootRanges finishes in a forked child that loaded nullranges after the fork", {
  skip_on_os("windows")
  y <- fork_input()
  ref <- fork_boot(y)
  # the state loading data.table and nullranges in the child leaves:
  # data.table's thread count from the environment, and the child's id as
  # the process nullranges was loaded in
  got <- in_child(function() {
    data.table::setDTthreads(4)
    ns <- asNamespace("nullranges")
    if (exists("load_state", envir = ns, inherits = FALSE)) {
      assign("pid", Sys.getpid(), envir = ns$load_state)
    }
    fork_boot(y)
  })
  expect_false(is.null(got))
  expect_identical(got, ref)
})

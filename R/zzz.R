# suppress R CMD CHECK note on globals
if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))

# the process nullranges was loaded in
load_state <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  load_state$pid <- Sys.getpid()
}

# Threads for the OpenMP regions in src/: 1 in a forked child (mclapply,
# mcparallel), where a parallel region waits forever for the parent's
# threads. getDTthreads() is 1 there only if data.table was loaded before
# the fork and the child left it so; parallel's flag also covers a child
# that loaded nullranges after the fork.
omp_threads <- function() {
  if (!identical(Sys.getpid(), load_state$pid) || forked_child()) {
    return(1L)
  }
  getDTthreads()
}

forked_child <- function() {
  isNamespaceLoaded("parallel") &&
    isTRUE(tryCatch(parallel:::isChild(), error = function(e) FALSE))
}

## Internal functions for matchRanges ------------------------------------------------------

#' Validate covariate formula
#' @param f A formula object.
#' @return A validated character vector of covariates
#' @noRd
parseFormula <- function(f) {

  ## Check for valid formula
  if (length(f) != 0) {

    lhs <- rlang::f_lhs(f)
    rhs <- rlang::f_rhs(f)

    if (!is.null(lhs)) {
      msg <- paste0('covar only accepts rhs of formula.')
      stop(msg)
    }

  } else {
    msg <- paste0('covar cannot be an empty formula.')
    stop(msg)
  }

  ## Return covariates
  return(labels(terms(f)))

}

#' Define sample.vec to handle vectors of varying length
#' @noRd
sample.vec <- function(x, ...) x[sample(length(x), ...)]

#' Nearest neighbor matching with replacement (data.table)
#' @param fps A numeric vector of `focal` propensity scores.
#' @param pps A numeric vector of `pool` propensity scores.
#' @inheritParams matchRanges
#' @return A matched `data.table`.
#' @noRd
nnMatch <- function(fps, pps, replace) {

  ## Suppress R CMD CHECK NOTE
  ppsIndex <- NULL
  fpsIndex <- NULL
  N <- NULL

  ## Create data table with original pps index and setkey to sort
  dt <- data.table(pps, val = pps, ppsIndex = seq_along(pps))
  setkey(dt, pps)

  ## Map ids to unique propensity scores
  fpsMap <- data.table(fps = fps, fpsIndex = seq_along(fps))

  ## Order fpsMap by fps (for to adding index later)
  fpsMap <- fpsMap[order(fps)]

  ## Create data table with each unique fps and the number of occurrences
  uniq.fps <- fpsMap[, .(fps, .N), by = fps]

  ## Find all nearest matches for each unique fps
  amdt <- dt[.(uniq.fps), roll = 'nearest', mult = 'all']

  ## Each row of uniq.fps gives one block of rows of amdt; when the blocks
  ## are the fps groups, draw from every block in one compiled loop, which
  ## makes the same draws as sample.vec on each group in turn (an empty
  ## pool gives a malformed amdt, left to the grouped form's error)
  grp <- if (isTRUE(replace) && nrow(amdt) > 0L &&
             is.integer(amdt$ppsIndex) &&
             length(amdt$ppsIndex) == nrow(amdt))
    amdt[, .N, by = fps]
  if (!is.null(grp) && nrow(grp) == nrow(uniq.fps)) {
    size <- amdt$N[cumsum(c(1L, grp$N[-nrow(grp)]))]
    mdt <- data.table(fps = rep(grp$fps, size),
                      ppsIndex = .Call(C_sample_groups, amdt$ppsIndex,
                                       grp$N, size))
  } else {
    ## Randomly subsample with replacement among duplicates;
    ## Check for categorical variables during sampling
    mdt <-
      amdt[, {
        if (length(N) > 1) {
          stopifnot(length(unique(N)) == 1)
          .(ppsIndex = sample.vec(ppsIndex, N[1], replace = replace))
        } else {
          .(ppsIndex = sample.vec(ppsIndex, N, replace = replace))
        }
      }, by = fps]
  }

  ## Add fpsIndex to mdt
  stopifnot(mdt$fps == fpsMap$fps)
  mdt$fpsIndex <- fpsMap$fpsIndex

  ## Reorder by fpsIndex (to match fps)
  mdt <- mdt[order(fpsIndex)]

  ## Return matched data table
  return(mdt)

}

#' Perform rejection sampling
#' @inheritParams nnMatch
#' @return A binary numeric vector equal to the length
#'   of pps representing whether or not to accept each
#'   pps based on fps distribution.
#' @noRd
rejectSample <- function(fps, pps) {

  if (!requireNamespace("ks", quietly = TRUE)) {
    stop("rejection sampling requires installing the CRAN package 'ks'")
  }
  
  ## Ensure fps <= pps
  if (length(fps) > length(pps)) {
    stop("focal must be <= pool for method = 'rejection'.")
  }

  ## Kernal density estimates for focal and pool
  ## (compute.cont = FALSE: the contour levels are never read)
  df <- ks::kde(fps, compute.cont = FALSE)
  dp <- ks::kde(pps, compute.cont = FALSE)

  ## Set scale by finding the highest point of density ratios (focal/pool)
  ## This ensures that pool covers focal at all points
  grid <- seq(from=quantile(pps, .001), to=quantile(pps, .999), length=1000)
  fgrid <- predict(df, x=grid)
  pgrid <- predict(dp, x=grid)
  if (any(fgrid < 0 | pgrid < 0)) {
    stop("kernel density estimates by ks::kde are negative, cannot perform rejection sampling")
  }
  scale <- max(fgrid/pgrid)
  if (scale > 1e3) {
    stop("scaling factor for density of the PS for pool is > 1e3, could lead to instability")
  }

  ## Calculate the probability of accepting each pool
  thresh <- function(x) ifelse(x > 1e-3, x, 0)
  accept_prob <- pmin(thresh(predict(df, x=pps))/(scale * predict(dp, x=pps)), 1)

  ## Randomly select pps according to accept probability
  accept <- rbinom(length(pps), size=1, prob=accept_prob)

  ## Return random selection equal to the length of pps
  return(accept)

}

#' Rejection sampling and matching for discrete and
#' continuous distributions.
#' @inherit nnMatch params return
#' @noRd
rsMatch <- function(fps, pps, replace) {

  ## Detect descrete distribution
  accept <- tryCatch(expr = {

    ## Rejection sampling
    rejectSample(fps, pps)

  }, error = function(e) {

    ## Introduce random noise to convert discrete into continuous
    stdev <- seq(1e-06, 0.1, length.out = 100)

    ## Iterate over possible noise levels
    for(i in seq_along(stdev)) {

      ## Add noise to fps and pps
      nfps <- fps + (rnorm(length(fps), mean = 0, sd = stdev[i]))
      npps <- pps + (rnorm(length(pps), mean = 0, sd = stdev[i]))

      ## Bound by 0 and 1
      nfps <- pmax(0, pmin(nfps, 1))
      npps <- pmax(0, pmin(npps, 1))

      ## Try rejection sampling
      # Mike: I'm taking off the suppressWarnings,
      #   we can instead catch these upstream?
      accept <- rejectSample(nfps, npps)

      ## Keep trying until there are enough options
      if (any(is.na(accept))) next
      if (sum(accept) >= length(fps)) break

    }

    return(accept)

  })

  ## Ensure there are adequate options
  stopifnot(sum(accept) >= length(fps))

  ## Select matching indices irrespective of data order
  ppsIndex <- sample(which(accept == 1), length(fps), replace = replace)

  ## Assemble matched data table
  mdt <- data.table(fps, ppsIndex, fpsIndex = seq_along(fps))

  ## Return matched data table
  return(mdt)

}

#' `.bincode(x, br, include.lowest = TRUE)` for breaks spaced `by`
#' apart from `mn`
#' @param x A numeric vector.
#' @param br The breaks, `seq(from = mn, to = mx, by = by)`.
#' @param mn,by The first break and the spacing of `br`.
#' @return The integer bin codes `.bincode` returns.
#' @noRd
binCode <- function(x, br, mn, by) {
  nb <- length(br)
  if (nb < 2L || anyNA(br) || !is.finite(mn) || !is.finite(by) || by <= 0 ||
      is.unsorted(br))
    return(.bincode(x, br, include.lowest = TRUE))

  ## For sorted breaks the one g with br[g] < x <= br[g + 1] is the
  ## code .bincode returns; a guess from the spacing finds it without a
  ## binary search, and the misses (x == br[1], outside the breaks,
  ## rounding) go to .bincode. NA stays NA.
  g <- as.integer(pmin(pmax(ceiling((x - mn) / by), 1), nb - 1L))
  bad <- which(!(br[g] < x & x <= br[g + 1L]))
  if (length(bad))
    g[bad] <- .bincode(x[bad], br, include.lowest = TRUE)
  g
}

#' Stratifies propensity scores by bins
#' @param f A numeric vector of the remaining `focal` propensity scores.
#' @param p A numeric vector of the remaining `pool` propensity scores.
#' @param n The number of bins used for stratifying
#'   propensity scores.
#' @return `NULL` when no bin can be sampled, otherwise a list
#'   with `fpsIndices` and `ppsIndices`, the positions in `f` and
#'   `p` of each bin holding at least one focal and no fewer
#'   pool than focal, in the order the bins first appear in `p`.
#' @noRd
stratify <- function(f, p, n) {

  ## Define breaks using fps and pps (min(c(f, p)) warns once
  ## when only one of them is empty)
  if (length(f) && length(p)) {
    mn <- min(min(f), min(p))
    mx <- max(max(f), max(p))
  } else {
    mn <- min(c(f, p))
    mx <- max(c(f, p))
  }
  by <- (mx-mn)/n
  br <- seq(from=mn, to=mx, by=by)

  ## Assign fps and pps to bin; values outside the breaks (NA)
  ## form one more bin, as NA matches NA in a data.table join
  fbin <- binCode(f, br, mn, by)
  pbin <- binCode(p, br, mn, by)
  nb <- length(br)
  fbin[is.na(fbin)] <- nb
  pbin[is.na(pbin)] <- nb

  ## Bins that can be sampled
  fpsN <- tabulate(fbin, nb)
  ppsN <- tabulate(pbin, nb)
  ok <- fpsN > 0L & fpsN <= ppsN
  if (!any(ok)) return(NULL)

  ## Gather positions of those bins only, keeping row order
  pq <- which(ok[pbin])
  fq <- which(ok[fbin])
  bins <- unique(pbin[pq])
  list(fpsIndices = split(fq, factor(fbin[fq], levels = bins)),
       ppsIndices = split(pq, factor(pbin[pq], levels = bins)))

}

#' One stratify-and-match pass of ssMatch on its compiled state
#' @param st The state from `C_ss_init`, updated in place.
#' @param n The number of bins used for stratifying.
#' @param replace Draw with replacement within a bin.
#' @param count Count the matched rows out of nLeft.
#' @return `c(matched, nLeft)`: matched is 0 when no bin qualifies (as
#'   `stratify` returns `NULL`) and `NA` when the breaks fail `binCode`'s
#'   guard, in which case the state is unchanged and the pass is to be
#'   made in R.
#' @noRd
ssPass <- function(st, n, replace, count) {
  rng <- .Call(C_ss_range, st)
  mn <- rng[1L]
  mx <- rng[2L]
  by <- (mx-mn)/n
  br <- seq(from=mn, to=mx, by=by)
  .Call(C_ss_pass, st, br, mn, by, replace, count)
}

#' Perform iterative stratified sampling
#' @inherit nnMatch params return
#' @noRd
ssMatch <- function(fps, pps, replace) {
  ## The remaining focal and pool: scores and indices
  nf <- length(fps)
  f <- fps
  fIdx <- seq_along(fps)
  p <- pps
  pIdx <- seq_along(pps)

  ## The number of distinct scores remaining, kept as rows leave:
  ## code[j] is the first position of row j's score in c(fps, pps)
  ## (match() and unique() share one notion of equality)
  x <- c(fps, pps)
  code <- match(x, x)
  cnt <- tabulate(code, length(x))
  nLeft <- sum(cnt > 0L)
  rm(x)

  ## Double scores make each pass in C (st), drawing exactly as
  ## sample() does below the 1e7 rows where sample.int would hash
  st <- NULL
  if (is.double(fps) && is.double(pps) && length(fps) <= 1e7 &&
      length(pps) <= 1e7 && (isTRUE(replace) || isFALSE(replace)))
    st <- .Call(C_ss_init, fps, pps, cnt, code, nLeft)

  ## Matched indices, one element per pass
  resF <- resP <- list(integer())
  nDone <- 0L

  ## Set flags and iteration count
  skip <- FALSE
  i <- 1

  ## Start progress bar
  pb <- progress::progress_bar$new(
    format = "  :step [:bar] :percent elapsed: :elapsedfull",
    clear = FALSE, total = length(fps) + 1)
  pb$tick(0)

  while (nDone != nf) {
    ## Update n
    if (skip)  n <- floor(n/2)
    if (!skip) n <- nLeft

    if (!is.null(st)) {
      ## Stratify, match and remove the matched rows in C
      r <- ssPass(st, n, replace, !skip)
      while (identical(r[1L], 0L)) {
        skip <- TRUE
        n <- floor(n/2)
        r <- ssPass(st, n, replace, FALSE)
      }
      if (!is.na(r[1L])) {
        pb$update(tokens=list(step=sprintf('Iteration %s, %s bin(s)', i, n)),
                  ratio = nDone/nf)
        i <- i + 1
        nDone <- nDone + r[1L]
        nLeft <- r[2L]
        next
      }

      ## The breaks fail binCode's guard: go on in R from this state
      if (nDone > 0L) {
        lens <- st[[7L]]
        fIdx <- st[[2L]][seq_len(lens[1L])]
        f <- fps[fIdx]
        pIdx <- st[[4L]][seq_len(lens[2L])]
        p <- pps[pIdx]
        cnt <- st[[5L]]
        resF <- list(st[[8L]][seq_len(nDone)])
        resP <- list(st[[9L]][seq_len(nDone)])
      }
      st <- NULL
    }

    ## Stratify ps by bins and match focal and pool
    strata <- stratify(f, p, n)
    while (is.null(strata)) {
      ## Enter faster n searching
      skip <- TRUE
      n <- floor(n/2)

      ## Stratify ps by bins and match focal and pool
      strata <- stratify(f, p, n)
    }
    ## Update progress
    pb$update(tokens=list(step=sprintf('Iteration %s, %s bin(s)', i, n)),
              ratio = nDone/nf)
    i <- i + 1

    ## Assign positions that can be sampled, bin by bin
    ppsSampled <- lapply(seq_along(strata$ppsIndices), function(k)
      sample.vec(strata$ppsIndices[[k]], length(strata$fpsIndices[[k]]),
                 replace = replace))
    fq <- unlist(strata$fpsIndices, use.names = FALSE)
    pq <- unique(unlist(ppsSampled, use.names = FALSE))
    resF[[i]] <- fIdx[fq]
    resP[[i]] <- pIdx[unlist(ppsSampled, use.names = FALSE)]
    nDone <- nDone + length(fq)

    ## Count out the scores of the matched rows (n is only read
    ## again while skip is FALSE)
    if (!skip) {
      rc <- code[c(resF[[i]], nf + pIdx[pq])]
      ru <- unique(rc)
      cnt[ru] <- cnt[ru] - tabulate(match(rc, ru), length(ru))
      nLeft <- nLeft - sum(cnt[ru] == 0L)
    }

    ## Remove matched rows from the remaining focal and pool
    f <- f[-fq]
    fIdx <- fIdx[-fq]
    p <- p[-pq]
    pIdx <- pIdx[-pq]
  }

  ## Close progress bar
  pb$update(tokens=list(step=sprintf('Iteration %s, %s bin(s), done!', i, n)),
            ratio = nDone/nf)
  if(pb$finished) pb$terminate()

  ## Reorder by fpsIndex
  if (!is.null(st)) {
    fpsIndex <- st[[8L]]
    ppsIndex <- st[[9L]]
  } else {
    fpsIndex <- unlist(resF, use.names = FALSE)
    ppsIndex <- unlist(resP, use.names = FALSE)
  }
  o <- order(fpsIndex)

  ## Assemble matched data table
  mdt <- data.table(fps = fps[fpsIndex[o]],
                    ppsIndex = ppsIndex[o],
                    fpsIndex = fpsIndex[o])

  ## Return matched data table
  return(mdt)
}

#' glm's fitting method for propensityMatch
#'
#' glm.fit for binomial("logit") with the IRLS loop in compiled code
#' (src/logitIRLS.c), which computes glm.fit's fitted values bitwise. It
#' covers the designs propensityMatch makes: a factor response, prior weights
#' of 1, no offset or starting values and glm's default control. glm.fit
#' itself runs for any other design and wherever the loop would stop, warn
#' inside it, halve a step or leave the BLAS product.
#' @param x,y,weights,start,etastart,mustart,offset,family,control,intercept,singular.ok
#'   As glm passes them to its `method`.
#' @return A list holding glm.fit's `fitted.values`, named by `y`'s names,
#'   which is all that `predict(type = "response")` reads; glm.fit's own
#'   value when glm.fit runs.
#' @noRd
logitGlmFit <- function(x, y, weights = NULL, start = NULL, etastart = NULL,
                        mustart = NULL, offset = NULL, family = binomial(),
                        control = list(), intercept = TRUE,
                        singular.ok = TRUE) {
  ## glm applies glm.control only when method is "glm.fit"
  control <- do.call("glm.control", control)
  fit <- NULL
  if (is.null(weights) && is.null(start) && is.null(etastart) &&
      is.null(mustart) && is.null(offset) && isTRUE(singular.ok) &&
      identical(control, list(epsilon = 1e-08, maxit = 25, trace = FALSE)) &&
      identical(family$family, "binomial") && identical(family$link, "logit") &&
      is.matrix(x) && is.double(x) && is.factor(y) && is.null(dim(y)) &&
      nrow(x) == length(y) && !anyNA(levels(y)) &&
      identical(getOption("matprod"), "default") &&
      isTRUE(capabilities("long.double")))
    fit <- .Call(C_logit_irls, x, y, control$epsilon, control$maxit,
                 min(1e-07, control$epsilon/1000))
  if (is.null(fit))
    return(eval(call("glm.fit", x = x, y = y, weights = weights,
                     start = start, etastart = etastart, mustart = mustart,
                     offset = offset, family = family, control = control,
                     intercept = intercept, singular.ok = singular.ok)))

  ## glm.fit's warnings after its loop, in its order and translation
  if (!fit[[2L]])
    warning(gettext("glm.fit: algorithm did not converge", domain = "R-stats"),
            call. = FALSE, domain = NA)
  if (fit[[3L]])
    warning(gettext("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                    domain = "R-stats"), call. = FALSE, domain = NA)
  mu <- fit[[1L]]
  names(mu) <- if (is.matrix(y)) rownames(y) else names(y)
  list(fitted.values = mu)
}

#' Calculate propensity scores, implement matching
#' methods, and construct Matched object
#' @param covarData A `data.table` with covariate data.
#'   Contains id column for designating whether covariates
#'   belong in focal or pool (1 or 0, respectively), and
#'   columns for the values of each covariate.
#' @param nf The number of focal rows: rows 1..nf of `covarData` are
#'   focal and the rest are pool.
#' @inheritParams matchRanges
#' @return A `Matched` object.
#' @noRd
propensityMatch <- function(covarData, covars, method, replace, nf) {

  ## Assemble covariate formula
  f <- as.formula(paste("id ~", paste(covars, collapse = "+")))

  ## Run glm model

  # April 2023 - speedglm is removed from CRAN so we cannot use it here...
  #model <- speedglm(formula = f, data = covarData,
  #                  family = binomial("logit"), fitted = TRUE, model = TRUE)
  ## binomial("logit") without the work nothing here reads: glm.fit calls
  ## aic once, after the IRLS loop, only for model$aic; validmu gives the
  ## same TRUE/FALSE as binomial's without its two logical vectors
  family <- binomial("logit")
  family$aic <- function(y, n, mu, wt, dev) NA_real_
  family$validmu <- function(mu)
    !anyNA(mu) && (length(mu) == 0L || (min(mu) > 0 && max(mu) < 1))
  model <- glm(formula = f, data = covarData, family = family,
               method = logitGlmFit)
  
  ## Get propensity scores of focal and pool sets as vectors
  id <- NULL

  ## Unnamed first: data.table's copy would expand glm's deferred row-name
  ## strings, one per row, only to drop them
  ps <- predict(model, type = "response")
  names(ps) <- NULL
  psData <- data.table(ps = ps, id = covarData$id)
  focalRows <- seq_len(nf)
  poolRows <- nf + seq_len(nrow(covarData) - nf)
  fps <- psData$ps[focalRows]
  pps <- psData$ps[poolRows]

  ## Add propensity scores to covarData (by reference: covarData is the
  ## caller's temporary)
  set(covarData, j = "ps", value = psData$ps)

  ## Implement matching method w or w/o replacement
  if (method == 'nearest' & isTRUE(replace))
    mdt <- nnMatch(fps, pps, replace = TRUE)

  if (method == 'nearest' & isFALSE(replace))
    stop("nearest neighbor matching without replacement not available.")

  if (method == 'rejection' & isTRUE(replace))
    mdt <- rsMatch(fps, pps, replace = TRUE)

  if (method == 'rejection' & isFALSE(replace))
    mdt <- rsMatch(fps, pps, replace = FALSE)

  if (method == 'stratified' & isTRUE(replace))
    mdt <- ssMatch(fps, pps, replace = TRUE)

  if (method == 'stratified' & isFALSE(replace))
    mdt <- ssMatch(fps, pps, replace = FALSE)


  ## Assemble information by set: focal, matched (draw order, repeats
  ## kept), pool, and unmatched (pool order, rows never drawn), as one
  ## positional subset. The rbind stays for a covariate named "set" (it
  ## appends a second "set" column) and a factor with an NA level (rbind
  ## drops the level and turns its codes NA).
  ppsIndex <- mdt$ppsIndex
  np <- length(poolRows)
  if (!("set" %in% names(covarData)) && is.integer(ppsIndex) &&
      !anyNA(ppsIndex) && (length(ppsIndex) == 0L ||
      (min(ppsIndex) >= 1L && max(ppsIndex) <= np)) &&
      !any(vapply(covarData, function(x) is.factor(x) && anyNA(levels(x)),
                  logical(1)))) {
    drawn <- logical(np)
    drawn[ppsIndex] <- TRUE
    unmatchedRows <- which(!drawn)
    matchedData <- covarData[c(focalRows, nf + ppsIndex, poolRows,
                               nf + unmatchedRows)]
    set(matchedData, j = "set",
        value = rep(c("focal", "matched", "pool", "unmatched"),
                    c(nf, length(ppsIndex), np, length(unmatchedRows))))
  } else {
    poolData <- covarData[poolRows]
    matchedData <- rbind(
      covarData[focalRows, c(.SD, set = 'focal')],
      poolData[mdt$ppsIndex, c(.SD, set = 'matched')],
      poolData[, c(.SD, set = 'pool')],
      poolData[!mdt$ppsIndex, c(.SD, set = 'unmatched')]
    )
  }

  ## Matched indicies
  matchedIndex <- mdt$ppsIndex

  ## Combine information in Matched object
  obj <- Matched(matchedData = matchedData,
                 matchedIndex = matchedIndex,
                 covar = covars,
                 method = method,
                 replace = replace)

  return(obj)
}

#' Covariate data.table: the id column, then each covariate's focal and pool
#' values, equal to as.data.table(cbind(id = ids, rbind(focal[covars],
#' pool[covars]))). Columns are concatenated directly when every covariate is
#' an attribute-free logical, integer, double or character vector, or a plain
#' factor, of one type (and levels) in both sets; otherwise the DataFrame route
#' does the coercions (its rbind drops an NA factor level, for one).
#' @param focal,pool `DataFrame` objects holding the covariates.
#' @param covars Character vector of covariate column names.
#' @param ids Factor of length nrow(focal) + nrow(pool).
#' @return A `data.table`.
#' @noRd
covarTable <- function(focal, pool, covars, ids) {
  direct <- identical(make.names(c("id", covars), unique = TRUE),
                      c("id", covars))
  cols <- vector("list", length(covars))
  for (i in seq_along(covars)) {
    if (!direct) break
    x <- focal[[covars[i]]]
    y <- pool[[covars[i]]]
    ax <- attributes(x)
    ay <- attributes(y)
    if (is.null(ax) && is.null(ay)) {
      direct <- typeof(x) == typeof(y) &&
        typeof(x) %in% c("logical", "integer", "double", "character")
    } else {
      direct <- identical(class(x), "factor") &&
        identical(class(y), "factor") &&
        setequal(names(ax), c("levels", "class")) && length(ax) == 2L &&
        setequal(names(ay), c("levels", "class")) && length(ay) == 2L &&
        identical(ax$levels, ay$levels) && !anyNA(ax$levels)
    }
    if (direct) cols[[i]] <- c(x, y)
  }
  if (!direct)
    return(as.data.table(cbind(id = ids, rbind(focal[covars], pool[covars]))))
  names(cols) <- covars
  setDT(c(list(id = ids), cols))
}

#' Core matching for DataFrames/GRanges/GInteractions
#' @param focal A `DataFrame` object containing
#'   the focal data to match.
#' @param pool A `DataFrame` object containing
#'   the pool from which to select matches.
#' @inheritParams matchRanges
#' @return A `Matched` object.
#' @noRd
matchRanges_Core <- function(focal, pool, covar, method, replace) {

  ## Extract covariates from formula as character vector
  covars <- parseFormula(covar)

  ## Check that all covariates are in both focal and pool
  if (!(all(covars %in% colnames(focal)) &
        all(covars %in% colnames(pool)))) {
    stop("All variables in covar must be columns in both focal and pool.")
  }

  ## Check method and replace arguments
  method <- match.arg(method, choices = c('nearest', 'rejection', 'stratified'))

  if (isFALSE(replace) & nrow(focal) >= nrow(pool))
    stop("focal must be <= pool when replace = FALSE.")

  if (method == 'nearest' & isFALSE(replace))
    stop("nearest neighbor matching without replacement not available.")

  if (method == 'stratified' & nrow(focal) >= nrow(pool))
    stop("focal must be <= pool for stratified sampling.")

  ## id: focal "1", pool "0", equal to factor(c(rep(1, nf), rep(0, np))),
  ## built directly because factor() converts a double vector to character
  nf <- nrow(focal)
  np <- nrow(pool)
  if (nf >= 1L && np >= 1L) {
    ids <- structure(c(rep(2L, nf), rep(1L, np)),
                     levels = c("0", "1"), class = "factor")
  } else {
    ids <- factor(c(rep(1, nf), rep(0, np)))
  }

  ## Create data table with covariate data
  covarData <- covarTable(focal, pool, covars, ids)

  ## Calculate propensity scores and match data
  md <- propensityMatch(covarData, covars, method, replace, nf = nrow(focal))

  return(md)
}


## Matched subclass methods for matchRanges ----------------------------------------------

#' MatchedDataFrame method for matchRanges
#' @inheritParams matchRanges_Core
#' @return A `MatchedDataFrame` object.
#' @noRd
matchRanges_MatchedDataFrame <- function(focal, pool, covar, method, replace) {

  ## Convert focal and pool to DataFrames
  f <- DataFrame(focal)
  p <- DataFrame(pool)

  ## Execute shared GRanges/GInteractions matching code
  md <- matchRanges_Core(f, p, covar, method, replace)
  idx <- indices(md)

  ## Combine matched data into MatchedDataFrame class
  obj <- MatchedDataFrame(focal = f,
                          pool  = p,
                          matchedData = matchedData(md),
                          matchedIndex = idx,
                          covar = covariates(md),
                          method = method(md),
                          replace = withReplacement(md),
                          delegate = pool[idx,])

  ## Return MatchedDataFrame object
  return(obj)
}

#' Generate a covariate-matched control set of ranges
#' 
#' `matchRanges()` uses a propensity score-based method to
#' generate a covariate-matched control set of DataFrame,
#' GRanges, or GInteractions objects.
#' 
#' Available inputs for `focal` and `pool` include `data.frame`,
#' `data.table`, `DataFrame`, `GRanges`, or `GInteractions`.
#' `data.frame` and `data.table` inputs are coerced to `DataFrame`
#' objects and returned as `MatchedDataFrame` while `GRanges` and
#' `GInteractions` objects are returned as `MatchedGRanges` or
#' `MatchedGInteractions`, respectively.
#' 
#' @section Methodology:
#' `matchRanges` uses 
#' [propensity scores](https://en.wikipedia.org/wiki/Propensity_score_matching)
#' to perform subset selection on the `pool` set such that the resulting `matched`
#' set contains similar distributions of covariates to that of the `focal` set.
#' A propensity score is the conditional probability of assigning an element
#' (in our case, a genomic range) to a particular outcome (`Y`) given a set of
#' covariates. Propensity scores are estimated using a logistic regression model
#' where the outcome `Y=1` for `focal` and `Y=0` for `pool`, over the provided 
#' covariates `covar`.
#' 
#' @section Matching methods:
#' \itemize{
#'   \item `method = 'nearest'`: Nearest neighbor matching 
#'   with replacement. Finds the nearest neighbor by using a 
#'   rolling join with `data.table`. Matching without replacement
#'   is not currently supported.
#'   \item `method = 'rejection'`: (Default) Rejection sampling
#'   with or without replacement. Uses a probability-based approach
#'   to select options in the `pool` that match the `focal` distribition.
#'   \item `method = 'stratified'`: Iterative stratified sampling
#'   with or without replacement. Bins `focal` and `pool` propensity
#'   scores by value and selects matches within bins until all `focal`
#'   items have a corresponding match in `pool`.
#' }
#'
#' @param focal A DataFrame, GRanges, or GInteractions object containing
#'   the focal data to match.
#' @param pool A DataFrame, GRanges, or GInteractions object containing
#'   the pool from which to select matches.
#' @param covar A rhs formula with covariates on which to match.
#' @param method A character describing which matching method to use.
#'   supported options are either 'nearest', 'rejection', or 'stratified'.
#' @param replace TRUE/FALSE describing whether to select matches with or without
#'   replacement.
#' @param ... Additional arguments.
#'
#' @return A covariate-matched control set of data.
#'
#' @references
#'
#' matchRanges manuscript:
#'
#' Eric S. Davis, Wancen Mu, Stuart Lee, Mikhail G. Dozmorov,
#' Michael I. Love, Douglas H. Phanstiel. 2023.
#' "matchRanges: Generating null hypothesis genomic ranges
#' via covariate-matched sampling."
#' Bioinformatics. doi: 10.1093/bioinformatics/btad197
#' 
#' @examples 
#' ## Match with DataFrame
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = 'DataFrame')
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3)
#' 
#' ## Match with GRanges
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = "GRanges")
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3)
#' 
#' ## Match with GInteractions
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = "GInteractions")
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3)
#' 
#' ## Nearest neighbor matching with replacement
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = 'DataFrame')
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3,
#'             method = 'nearest',
#'             replace = TRUE)
#' 
#' ## Rejection sampling without replacement
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = 'DataFrame')
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3,
#'             method = 'rejection',
#'             replace = FALSE)
#' 
#' ## Stratified sampling without replacement
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(type = 'DataFrame')
#' matchRanges(focal = x[x$feature1,],
#'             pool = x[!x$feature1,],
#'             covar = ~feature2 + feature3,
#'             method = 'stratified',
#'             replace = FALSE)
#'
#' @rdname matchRanges
#'
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @importFrom rlang f_lhs f_rhs
#' @importFrom stats glm glm.control glm.fit
#' @import S4Vectors
#' 
#' @export
setMethod("matchRanges",
          signature = signature(focal   = "DF_OR_df_OR_dt",
                                pool    = "DF_OR_df_OR_dt",
                                covar   = "formula",
                                method  = "character_OR_missing",
                                replace = "logical_OR_missing"),
          definition = matchRanges_MatchedDataFrame)

#' MatchedGRanges method for matchRanges
#' @param focal A `GRanges` object containing
#'   the focal data to match.
#' @param pool A `GRanges` object containing
#'   the pool from which to select matches.
#' @inheritParams matchRanges
#' @return A `MatchedGRanges` object
#' @noRd
matchRanges_MatchedGRanges <- function(focal, pool, covar, method, replace) {

  ## Extract DataFrame from GRanges/GInteractions objects
  f <- mcols(focal)
  p <- mcols(pool)

  ## Execute matching code to get a Matched object
  md <- matchRanges_Core(f, p, covar, method, replace)
  idx <- indices(md)

  ## Combine matched data into MatchedGRanges class
  obj <- MatchedGRanges(focal = focal,
                        pool = pool,
                        matchedData = matchedData(md),
                        matchedIndex = idx,
                        covar = covariates(md),
                        method = method(md),
                        replace = withReplacement(md),
                        delegate = pool[idx])

  ## Return MatchedGRanges object
  return(obj)
}

#' @rdname matchRanges
#' @export
setMethod("matchRanges",
          signature = signature(focal   = "GRanges",
                                pool    = "GRanges",
                                covar   = "formula",
                                method  = 'character_OR_missing',
                                replace = 'logical_OR_missing'),
          definition = matchRanges_MatchedGRanges)

#' MatchedGInteractions method for matchRanges
#' @param focal A `GInteractions` object containing
#'   the focal data to match.
#' @param pool A `GInteractions` object containing
#'   the pool from which to select matches.
#' @inheritParams matchRanges
#' @return A `MatchedGInteractions` object
#' @noRd
matchRanges_MatchedGInteractions <- function(focal, pool, covar, method, replace) {

  ## Extract DataFrame from GRanges/GInteractions objects
  f <- mcols(focal)
  p <- mcols(pool)

  ## Execute matching code to get a Matched object
  md <- matchRanges_Core(f, p, covar, method, replace)
  idx <- indices(md)

  ## Combine matched data into MatchedGInteractions class
  obj <- MatchedGInteractions(focal = focal,
                              pool = pool,
                              matchedData = matchedData(md),
                              matchedIndex = idx,
                              covar = covariates(md),
                              method = method(md),
                              replace = withReplacement(md),
                              delegate = pool[idx])

  ## Return MatchedGInteractions object
  return(obj)
}

#' @rdname matchRanges
#' @export
setMethod("matchRanges",
          signature = signature(focal   = "GInteractions",
                                pool    = "GInteractions",
                                covar   = "formula",
                                method  = 'character_OR_missing',
                                replace = 'logical_OR_missing'),
          definition = matchRanges_MatchedGInteractions)

## Define accessor functions for matchGRanges class --------------------------------------

#' Get focal set from a Matched object
#' 
#' @param x A `MatchedDataFrame`, `MatchedGRanges`,
#'   or `MatchedGInteractions` object.
#' @param ... Additional options.
#' 
#' @return An object of the same class as `x` representing
#'   the focal set.
#' 
#' @examples 
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' focal(x)
#' 
#' @rdname focal
#' @export
setMethod("focal", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@focal
})

#' Get pool set from a Matched object
#' 
#' @inheritParams focal
#' 
#' @return An object of the same class as `x` representing
#'   the pool set.
#' 
#' @examples 
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' pool(x)
#' 
#' @rdname pool
#' @export
setMethod("pool", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@pool
})

#' Get matched set from a Matched object
#' 
#' @inheritParams focal
#' 
#' @return An object of the same class as `x` representing
#'   the matched set.
#' 
#' @examples 
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' matched(x)
#' 
#' @rdname matched
#' @export
setMethod("matched", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@delegate
})

#' Get unmatched set from a Matched object
#' 
#' @inheritParams focal
#' 
#' @return An object of the same class as `x` representing
#'   the unmatched set.
#' 
#' @examples 
#' set.seed(123)
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' unmatched(x)
#' 
#' @rdname unmatched
#' @export
setMethod("unmatched", "MDF_OR_MGR_OR_MGI", function(x, ...) {
  x@pool[indices(x, set = "unmatched"),]
})

## Coercion from "matchit" object to Matched* class ----------------------------

#' Coerce `matchit` to `MatchedDataFrame`
#' @inheritParams matchitToMatched
#' @importFrom data.table data.table
#' @returns A `Matched` object.
#' @noRd
as_Matched <- function(x) {
  
  ## Suppress R CMD CHECK NOTE
  id <- NULL
  
  if (!requireNamespace("MatchIt", quietly = TRUE)) {
    stop("The 'MatchIt' package is required for this function.")
  }
  
  ## Get matched indices from matchit object
  mi <- MatchIt::match.data(object=x, group='control') |>
    rownames() |>
    as.integer()
  
  ## Put these in terms of the pool set
  pmi <- which(which(x$treat == 0) %in% mi)
  
  ## Assemble covarData and matchedData by set
  covarData <- x$X
  covarData <- data.table(cbind(id=x$treat, covarData, ps=x$distance))
  matchedData <- rbind(
    covarData[id==1, c(.SD, set='focal')],
    covarData[mi, c(.SD, set='matched')],
    covarData[id==0, c(.SD, set='pool')],
    covarData[!mi, c(.SD, set='unmatched')]
  )
  
  ## Construct Matched object
  obj <- Matched(matchedData=matchedData,
                 matchedIndex=pmi,
                 covar=labels(terms(x$formula)),
                 method=paste0("MatchIt_", x$info$method),
                 replace=x$info$replace)
  
  return(obj)
}

#' Coerce `matchit` to `MatchedDataFrame`
#' @inheritParams matchitToMatched
#' @importFrom S4Vectors DataFrame
#' @noRd
as_MatchedDataFrame <- function(x) {
  
  ## Coerce to Matched object
  md <- as_Matched(x)
  d <- x$model$call$data
  
  ## Define focal & pool as DataFrames
  focal <- DataFrame(d[x$treat == 1,])
  pool <- DataFrame(d[x$treat == 0,])
  
  ## Combine matched data into MatchedDataFrame class
  obj <- MatchedDataFrame(focal = focal,
                          pool  = pool,
                          matchedData = matchedData(md),
                          matchedIndex = indices(md),
                          covar = covariates(md),
                          method = method(md),
                          replace = withReplacement(md),
                          delegate = DataFrame(pool[indices(md),]))
  
  ## Return MatchedDataFrame object
  return(obj)
}

#' Coerce to ranges/interactions
#' @param d data from `MatchIt`
#' @param ranges `GRanges` or `GInteractions` object
#' @param keep_mcols boolean, whether to keep metadata columns
#' @param type either "GRanges" or "GInteractions"
#' @return A `GRanges` or `GInteractions` object depending on `type`
#' @noRd
coerceToRanges <- function(d, ranges, keep_mcols, type="GRanges") {
  
  if (!requireNamespace("mariner", quietly = TRUE)) {
    stop("The 'mariner' package is required for this function.")
  }
  
  ## Set coercion function
  if (type=="GRanges") {
    FUN <- \(d, mc) as_granges(d, keep_mcols=mc)
  }
  if (type=="GInteractions") {
    FUN <- \(d, mc) mariner::as_ginteractions(d, keep.extra.columns=mc)
  }
  
  ## Coerce to GRanges/GInteractions or add supplied
  if (is.null(ranges)) {
    g <- FUN(d, keep_mcols)
  } else {
    stopifnot(length(ranges) == nrow(d)) # must be same length
    if (keep_mcols) {
      mcols(ranges) <- cbind(mcols(ranges), d)
    } else {
      mcols(ranges) <- d
    }
    g <- ranges
  }
  
  return(g)
}


#' Coerce `matchit` to `MatchedGRanges`
#' @inheritParams matchitToMatched
#' @importFrom plyranges as_granges
#' @noRd
as_MatchedGRanges <- function(x, ranges, keep_mcols) {

  ## Coerce to Matched object
  md <- as_Matched(x)
  d <- x$model$call$data
  
  ## Coerce to GRanges or add supplied ranges
  gr <- coerceToRanges(d, ranges, keep_mcols, "GRanges")
  
  ## Define focal & pool
  focal <- gr[x$treat == 1,]
  pool <- gr[x$treat == 0,]
  
  ## Combine matched data into MatchedGRanges class
  obj <- MatchedGRanges(focal = focal,
                        pool = pool,
                        matchedData = matchedData(md),
                        matchedIndex = indices(md),
                        covar = covariates(md),
                        method = method(md),
                        replace = withReplacement(md),
                        delegate = pool[indices(md)])
  
  ## Return MatchedGRanges object
  return(obj)
  
}

#' Coerce `matchit` to `MatchedGInteractions`
#' @inheritParams matchitToMatched
#' @param interactions A `GInteractions` object corresponding
#'  to the rows of data input to `MatchIt`.
#' @noRd
as_MatchedGInteractions <- function(x, interactions, keep_mcols) {
  
  ## Coerce to Matched object
  md <- as_Matched(x)
  d <- x$model$call$data
  
  ## Coerce to GInteractions or add supplied interactions
  gi <- coerceToRanges(d, interactions, keep_mcols, "GInteractions")
  
  ## Define focal & pool
  focal <- gi[x$treat == 1,]
  pool <- gi[x$treat == 0,]
  
  ## Combine matched data into MatchedGInteractions class
  obj <- MatchedGInteractions(focal = focal,
                              pool = pool,
                              matchedData = matchedData(md),
                              matchedIndex = indices(md),
                              covar = covariates(md),
                              method = method(md),
                              replace = withReplacement(md),
                              delegate = pool[indices(md),])
  
  ## Return MatchedGInteractions object
  return(obj)
}

#' Execute null method
#' @inheritParams matchitToMatched
#' @noRd
.matchitToMatched <- function(x, ranges, keep_mcols) {
  
  ## Dispatch on ranges
  if (is(ranges, "GRanges")) {
    as_MatchedGRanges(x, ranges, keep_mcols)
  }
  
  if (is(ranges, "GInteractions")) {
    as_MatchedGInteractions(x, ranges, keep_mcols)
  }
 
  ## Try each method if null
  ans <- tryCatch(as_MatchedGRanges(x, ranges, keep_mcols),
                  error=function(e) NULL)
  if (!is.null(ans)) return(ans)
  
  ans <- tryCatch(as_MatchedGInteractions(x, ranges, keep_mcols),
                  error=function(e) NULL)
  if (!is.null(ans)) return(ans)
  
  return(as_MatchedDataFrame(x))
}

#' Coerce `matchit` to `Matched` object
#' @param x A `matchit` object from `MatchIt` package.
#' @param ranges A `GRanges` or `GInteractions` object
#'  that is the same length as the data used to create
#'  the `matchit` object, `x`.
#' @param keep_mcols Logical whether to keep or existing
#'  mcols of the object supplied in `ranges` or not.
#'  Default is TRUE.
#' @rdname matchitToMatched
#' @return A `Matched` subclass object, depending on the
#'  value of `ranges`. If `is(ranges, "GRanges")` then a
#'  `MatchedGRanges` object is returned. If 
#'  `is(ranges, "GInteractions")` then a `MatchedGInteractions`
#'  object is returned. If `is.null(ranges)` then a
#'  the function first attempts to coerce into the other
#'  classes before coercing to `MatchedDataFrame`.
#'  
#' @examples
#' \donttest{
## Create example data.frame dataset
#' if (!requireNamespace("MatchIt", quietly=TRUE)) {
#'  set.seed(123)
#'  x <- makeExampleMatchedDataSet(type="GRanges")
#'  
#'  ## Convert GRanges to data.frame, pass to matchit,
#'  ## and convert to MatchedGRanges object
#'  set.seed(123)
#'  mgr <- 
#'    as.data.frame(x) |>
#'    matchit(formula=feature1 ~ feature2 + feature3,
#'            data=_,
#'            method='nearest',
#'            replace=FALSE) |>
#'    matchitToMatched()
#'  
#'  ## Compatible with GRanges & Matched functions
#'  mgr
#'  plotCovariate(mgr)
#' } else {
#'  message("The 'MatchIt' package is required to run this example.")
#' }
#' } 
#' 
#' @export
setMethod("matchitToMatched",
          signature(x="matchit"),
          definition=.matchitToMatched)
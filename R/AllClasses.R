#' @import methods
#' @export
bootRanges <- setClass(
  Class = "bootRanges",
  contains = "CompressedGRangesList"
)

setValidity(
  Class = "bootRanges",
  method = function(object) {
    stopifnot(all(c("blockLength", "iter") %in% names(mcols(object))))
  }
)

#' Matched objects
#'
#' The Matched class is a container for attributes of covariate-matched
#' data resulting from `matchRanges()`.
#'
#' @examples
#' ## Make Matched example
#' x <- makeExampleMatchedDataSet(matched = TRUE)
#' @template matched-class-slots
#' @template matched-class-details
#'
#' @rdname matchedClass
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @export
Matched <- setClass(
  Class = "Matched",
  slots = list(
    matchedData = "data.table",
    matchedIndex = "integer",
    covar = "character",
    method = "character",
    replace = "logical"
  )
)

setValidity(
  Class = "Matched",
  method = function(object) {
    stopifnot(all(c("id", "ps", "set") %in%
      colnames(object@matchedData)))
  }
)

#' MatchedOverview object
#'
#' Internal container for holding overview summary
#' data for Matched objects.
#'
#' Necessary to enable separate show and extracting
#' methods after using the `overview` method on
#' `Matched` classes and subclasses.
#'
#' @slot aggData A `data.table` with aggregated means
#'   and s.d. for covariates and propensity scores.
#' @slot diffData A `data.table` with the aggregated
#'   means and s.d. between focal and matched sets.
#' @slot quality A numeric representing the absolute
#'   value of the mean difference in propensity scores
#'   between focal and matched sets. Lower values
#'   indicate higher matching quality.
#'
#' @noRd
MatchedOverview <- setClass(
  Class = "MatchedOverview",
  slots = list(
    MatchedClass = "character",
    aggData = "data.table",
    diffData = "data.table",
    quality = "numeric"
  )
)

#' Class union for "DataFrame-like" objects
#' @noRd
setClassUnion("DF_OR_df_OR_dt", c("DFrame", "data.frame", "data.table"))

#' MatchedDataFrame objects
#'
#' The `MatchedDataFrame` class is a subclass of both
#' `Matched` and `DFrame`. Therefore, it contains slots
#' and methods for both of these classes.
#'
#' The `MatchedDataFrame` class uses a delegate object
#' during initialization to assign its `DFrame` slots.
#' `MatchedDataFrame` behaves as a `DataFrame` but also
#' includes additional `Matched` object functionality
#' (see `?Matched`). For more information about
#' `DataFrame` see `?S4Vectors::DataFrame`.
#'
#' @slot focal A `DataFrame` object containing the focal
#'  data to match.
#' @slot pool A `DataFrame` object containing the pool
#'  from which to select matches.
#' @slot delegate A `DataFrame` object used to initialize
#'  `DataFrame`-specific slots. `matchRanges()` assigns
#'  the matched set to the slot.
#' @template matched-class-slots
#' @slot rownames `rownames(delegate)`
#' @slot nrows `nrows(delegate)`
#' @slot listData `as.list(delegate)`
#' @slot elementType `elementType(delegate)`
#' @slot elementMetadata `elementMetadata(delegate)`
#' @slot metadata `metadata(delegate)`
#'
#' @seealso [S4Vectors::DataFrame]
#'
#' @examples
#' ## Constructing MatchedDataFrame with matchRanges
#' ## data.frame
#' x <- makeExampleMatchedDataSet(type = "data.frame")
#' mx <- matchRanges(
#'   focal = x[x$feature1, ],
#'   pool = x[!x$feature1, ],
#'   covar = ~ feature2 + feature3,
#'   method = "rejection",
#'   replace = FALSE
#' )
#' class(mx)
#'
#' ## data.table
#' x <- makeExampleMatchedDataSet(type = "data.table")
#' mx <- matchRanges(
#'   focal = x[x$feature1],
#'   pool = x[!x$feature1],
#'   covar = ~ feature2 + feature3,
#'   method = "rejection",
#'   replace = FALSE
#' )
#' class(mx)
#'
#' ## DataFrame
#' x <- makeExampleMatchedDataSet(type = "DataFrame")
#' mx <- matchRanges(
#'   focal = x[x$feature1, ],
#'   pool = x[!x$feature1, ],
#'   covar = ~ feature2 + feature3,
#'   method = "rejection",
#'   replace = FALSE
#' )
#' class(mx)
#'
#' ## Make MatchedDataFrame example
#' x <- makeExampleMatchedDataSet(type = "DataFrame", matched = TRUE)
#' @template matched-class-details
#' @template matched-subclass-details
#'
#' @rdname MatchedDataFrame
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @export
MatchedDataFrame <- setClass(
  Class = "MatchedDataFrame",
  contains = c("Matched", "DFrame"),
  slots = list(
    focal = "DF_OR_df_OR_dt",
    pool = "DF_OR_df_OR_dt",
    delegate = "DF_OR_df_OR_dt"
  )
)

setMethod(
  "initialize", "MatchedDataFrame",
  function(.Object, ..., delegate = DataFrame()) {

    ## Convert delegate to DataFrame
    delegate <- DataFrame(delegate)

    ## Use delegate to set other attributes in MatchedDataFrame
    .Object@delegate <- delegate
    .Object@rownames <- rownames(delegate)
    .Object@nrows <- nrow(delegate)
    .Object@listData <- as.list(delegate)
    .Object@elementType <- elementType(delegate)
    .Object@elementMetadata <- elementMetadata(delegate)
    .Object@metadata <- metadata(delegate)

    ## Use default class generator function for args in MatchedDataFrame and Matched
    .Object <- callNextMethod(.Object, ...)
    .Object
  }
)

#' MatchedGRanges objects
#'
#' The `MatchedGRanges` class is a subclass of both
#' `Matched` and `GRanges`. Therefore, it contains slots
#' and methods for both of these classes.
#'
#' The `MatchedGRanges` class uses a delegate object
#' during initialization to assign its `GRanges` slots.
#' `MatchedGRanges` behaves as a `GRanges` but also
#' includes additional `Matched` object functionality
#' (see `?Matched`). For more information about
#' `GRanges` see `?GenomicRanges::GRanges`.
#'
#' @slot focal A `GRanges` object containing the focal
#'  data to match.
#' @slot pool A `GRanges` object containing the pool
#'  from which to select matches.
#' @slot delegate A `GRanges` object used to initialize
#'  `GRanges`-specific slots. `matchRanges()` assigns
#'  the matched set to the slot.
#' @template matched-class-slots
#' @slot seqnames `seqnames(delegate)`
#' @slot ranges `ranges(delegate)`
#' @slot strand `strand(delegate)`
#' @slot seqinfo `seqinfo(delegate)`
#' @slot elementMetadata `elementMetadata(delegate)`
#' @slot elementType `elementType(delegate)`
#' @slot metadata `metadata(delegate)`
#'
#' @seealso [GenomicRanges::GRanges]
#'
#' @examples
#' ## Contructing MatchedGRanges with matchRanges
#' gr <- makeExampleMatchedDataSet(type = "GRanges")
#' mgr <- matchRanges(
#'   focal = gr[gr$feature1, ],
#'   pool = gr[!gr$feature1, ],
#'   covar = ~ feature2 + feature3,
#'   method = "rejection",
#'   replace = FALSE
#' )
#' class(mgr)
#'
#' ## Make MatchedGRanges example
#' x <- makeExampleMatchedDataSet(type = "GRanges", matched = TRUE)
#' @template matched-class-details
#' @template matched-subclass-details
#'
#' @rdname MatchedGRanges
#' @import GenomicRanges
#' @export
MatchedGRanges <- setClass(
  Class = "MatchedGRanges",
  contains = c("Matched", "GRanges"),
  slots = list(
    focal = "GRanges",
    pool = "GRanges",
    delegate = "GRanges"
  )
)

setMethod(
  "initialize", "MatchedGRanges",
  function(.Object, ..., delegate = GRanges()) {

    ## Use delegate to set other attributes in GenomicRanges
    .Object@delegate <- delegate

    for (x in slotNames(delegate)) {
      `slot<-`(
        object = .Object,
        name = x,
        value = do.call(x, list(delegate))
      )
    }

    ## Use default class generator function for args in MatchedGRanges and Matched
    .Object <- callNextMethod(.Object, ...)
    .Object
  }
)

#' MatchedGInteractions objects
#'
#' The `MatchedGInteractions` class is a subclass of both
#' `Matched` and `GInteractions`. Therefore, it contains slots
#' and methods for both of these classes.
#'
#' The `MatchedGInteractions` class uses a delegate object
#' during initialization to assign its `GInteractions` slots.
#' `MatchedGInteractions` behaves as a `GInteractions` but also
#' includes additional `Matched` object functionality
#' (see `?Matched`). For more information about
#' `GInteractions` see `?InteractionSet::GInteractions`.
#'
#' @slot focal A `GInteractions` object containing the focal
#'  data to match.
#' @slot pool A `GInteractions` object containing the pool
#'  from which to select matches.
#' @slot delegate A `GInteractions` object used to initialize
#'  `GInteractions`-specific slots. `matchRanges()` assigns
#'  the matched set to the slot.
#' @template matched-class-slots
#' @slot anchor1 `anchorIds(delegate)$first`
#' @slot anchor2 `anchorIds(delegate)$second`
#' @slot regions `regions(delegate)`
#' @slot NAMES `names(delegate)`
#' @slot elementMetadata `elementMetadata(delegate)`
#' @slot metadata `metadata(delegate)`
#'
#' @seealso [InteractionSet::GInteractions]
#'
#' @examples
#' ## Constructing MatchedGInteractions with matchRanges
#' gi <- makeExampleMatchedDataSet(type = "GInteractions")
#' mgi <- matchRanges(
#'   focal = gi[gi$feature1, ],
#'   pool = gi[!gi$feature1, ],
#'   covar = ~ feature2 + feature3,
#'   method = "rejection",
#'   replace = FALSE
#' )
#' class(mgi)
#'
#' ## Make MatchedGInteractions example
#' x <- makeExampleMatchedDataSet(type = "GInteractions", matched = TRUE)
#' @template matched-class-details
#' @template matched-subclass-details
#'
#' @rdname MatchedGInteractions
#' @import InteractionSet
#' @export
MatchedGInteractions <- setClass(
  Class = "MatchedGInteractions",
  contains = c("Matched", "GInteractions"),
  slots = list(
    focal = "GInteractions",
    pool = "GInteractions",
    delegate = "GInteractions"
  )
)

setMethod(
  "initialize", "MatchedGInteractions",
  function(.Object, ..., delegate = GInteractions()) {

    ## Use GInteractions object to set other attributes
    .Object@delegate <- delegate
    .Object@anchor1 <- anchorIds(delegate)$first
    .Object@anchor2 <- anchorIds(delegate)$second
    .Object@regions <- regions(delegate)
    .Object@NAMES <- names(delegate)
    .Object@elementMetadata <- elementMetadata(delegate)
    .Object@metadata <- metadata(delegate)

    ## Use default class generator function for args in MatchedGInteractions and Matched
    .Object <- callNextMethod(.Object, ...)
    .Object
  }
)

#' Class union for Matched subtypes
#' @noRd
setClassUnion(
  "MDF_OR_MGR_OR_MGI",
  c("MatchedDataFrame", "MatchedGRanges", "MatchedGInteractions")
)

#' Class unions for general types
#' @noRd
setClassUnion("character_OR_missing", c("character", "missing"))
#' @noRd
setClassUnion("logical_OR_missing", c("logical", "missing"))
#' @noRd
setClassUnion("numeric_OR_missing", c("numeric", "missing"))

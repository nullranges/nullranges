#' @rdname Matched
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @import methods
#' @export
Matched <- setClass(Class = "Matched",
                    slots = list(matchedData = "data.table",
                                 matchedIndex = "integer",
                                 covar = "character"))

setValidity(Class = "Matched",
            method = function(object){
              
              stopifnot(all(c('id', 'ps', 'group') %in%
                              colnames(object@matchedData)))
              
            })

#' @rdname MatchedDataFrame
#' @import S4Vectors
#' @rawNamespace import(data.table, except = c(between, shift, first, second, indices))
#' @import methods
#' @export
setClassUnion("DF_OR_df_OR_dt", c("DFrame", "data.frame", "data.table"))

MatchedDataFrame <- setClass(Class = "MatchedDataFrame",
                             contains = c("Matched", "DFrame"),
                             slots = list(focal = "DF_OR_df_OR_dt",
                                          pool  = "DF_OR_df_OR_dt",
                                          delegate = "DF_OR_df_OR_dt"))

setMethod("initialize", "MatchedDataFrame",
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
            
          })


#' @rdname MatchedGRanges
#' @import GenomicRanges
#' @import methods
#' @export
MatchedGRanges <- setClass(Class = "MatchedGRanges",
                           contains = c("Matched", "DelegatingGenomicRanges"),
                           slots = list(focal = "GRanges",
                                        pool  = "GRanges"))

setMethod("initialize", "MatchedGRanges",
          function(.Object, ..., delegate = GRanges()) {
            
            ## Use delegate to set other attributes in DelegatingGenomicRanges
            .Object@delegate <- delegate
            .Object@elementMetadata <- elementMetadata(delegate)
            .Object@elementType <- elementType(delegate)
            .Object@metadata <- metadata(delegate)
            
            ## Use default class generator function for args in MatchedGRanges and Matched
            .Object <- callNextMethod(.Object, ...)
            .Object
            
          })


#' @rdname MatchedGInteractions
#' @import InteractionSet
#' @import methods
#' @export
MatchedGInteractions <- setClass(Class = "MatchedGInteractions",
                                 contains = c("Matched", "GInteractions"),
                                 slots = list(focal = "GInteractions",
                                              pool  = "GInteractions",
                                              delegate = "GInteractions"))

setMethod("initialize", "MatchedGInteractions",
          function(.Object, ..., delegate = GInteractions()) {
            
            ## Use GInteractions object to set other attributes
            .Object@delegate <- delegate
            .Object@anchor1 <- anchorIds(delegate)$first
            .Object@anchor2 <- anchorIds(delegate)$second
            .Object@regions <- regions(delegate)
            .Object@elementMetadata <- elementMetadata(delegate)
            .Object@metadata <- metadata(delegate)
            
            ## Use default class generator function for args in MatchedGInteractions and Matched
            .Object <- callNextMethod(.Object, ...)
            .Object
            
          })

## Class Union for shared MatchedDataFrame/MatchedGRanges/MatchedGInteractions methods
setClassUnion("MDF_OR_MGR_OR_MGI",
              c("MatchedDataFrame", "MatchedGRanges", "MatchedGInteractions"))

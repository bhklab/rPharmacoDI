library(data.table)
library(BiocParallel)
library(PharmacoGx)

#' Read in all the gene signatures for a PharmacoSet and assemble into a single `data.table`
#' 
#' @param filePath [`character`] Where to search for the `pattern` argument to get the file paths for the
#'     PSet gene signatures.
#' @param pattern [`character`] Passed to `list.files` to match all .rds files containing `pattern`. This
#'     should usually be the PSet name as it will also be used as the file name if `save` is `TRUE`.
#' @param save [`logical`] Should the `data.table` be saved to disk as a .csv file? Default is `TRUE`.
#' 
#' @return [`data.table`] Optionally returns the gene signature data.table for the specified PSet
#' 
#' @importFrom BiocParallel bplapply
#' @importFrom data.table fread data.table
#' @importFrom PharmacoGx SignatureClass
#' @export
getPSetGeneSigs <- function(filePath, pattern, save=TRUE) {
    
    # read in gene sigs and convert to DTs
    sigFiles <- list.files(filePath, paste0(pattern, ".*rds$"), full.names=TRUE)
    geneSigs <- bplapply(sigFiles, readRDS)
    geneSigDTs <- lapply(geneSigs, function(sig) data.table(sig[,,], keep.rownames='gene_id'))
    
    # extract and format metadata
    .getArgument <- function(geneSig, argument) unlist(geneSig@Arguments[argument])
    drugs <- unlist(lapply(geneSigs, .getArgument, argument='drugs'))
    mDataTypes <- unlist(lapply(geneSigs, .getArgument, argument='mDataType'))
    tissues <- unlist(lapply(sigFiles, gsub, pattern='.*\\_|.rds$', replacement=''))

    # add metadata to DTs then rbind into a single DT
    .refAssignMetadata <- function(DT, drug, mDataType, tissue) {
        DT[, `:=`('drug'=drug, 'mDataType'=mDataType, 'tissue'=tissue)]
    }

    geneSigMetaDTs <- mapply(FUN=.refAssignMetadata,
                             DT=geneSigDTs, drug=drugs, mDataType=mDataTypes, tissue=tissues, 
                             SIMPLIFY=FALSE)
    geneSigPSetDT <- rbindlist(geneSigMetaDTs)

    # save PSet DT
    if (save) fwrite(geneSigPSetDT, file=file.path(filePath, paste0(pattern, '.csv')))
    return(geneSigPSetDT)
}

#' Take a SignatureClass object, convert it to a data.table and optionally save it to a .csv
#'
#' @param geneSig [`PharmacoSig`] object to convert to a `data.table`
#'
#' @return [data.table] containing the results in `geneSig` as a `data.table OR [`NULL`] if you specify `saveDir` and
#'     `fileName`, in which case the file is saved to disk instead of returned.
#'
#' @import data.table
#' @import PharmacoGx
#' @export
convertGeneSignatureToDT <- function(geneSig) {
    geneSigDT <- data.table(geneSig@.Data[,,], keep.rownames='gene')
    geneSigDT$pSet <- rep(geneSig@PSetName, nrow(geneSigDT))
    geneSigDT$drug <- geneSig@Arguments$drugs
    geneSigDT$mDataType <- geneSig@Arguments$mDataType
    return(geneSigDT)
}


#' Takes in the path to a set of gene signature .rds files, selects those pathing `pSetPattern` and converts them
#'
#' @param dataDir [`character`] vector specifying the path to the gene signature .rds files
#' @param pSetPattern [`character`] an identifier for the Pset you wish to read. In general this will be the
#'    pSet indentifier (e.g., GDSC_v1). '.*' is automatically prepended and appended to the pattern, so adding
#'    this is not necssary.
#' @param BPPARAM [`BiocParam`] object specifying what backend to use for parallelization. This defaults to your
#'    registered back-end, retrieved via the `bpparam()` function.
#'
#' @return [`list`] of `PharmacoSig` objects for the selected pSet
#'
#' @import data.table
#' @importFrom BiocParallel bplapply bpparam
#' @export
readGeneSigsForPSet <- function(dataDir, pSetPattern, BPPARAM=BiocParallel::bpparam()) {
    pSetFiles <- list.files(dataDir, pattern=paste0(".*", pSetPattern, ".*"), full.names=TRUE)
    pSetGeneSigs <- bplapply(pSetFiles, FUN=readRDS, BPPARAM=BPPARAM)
    names(pSetGeneSigs) <- gsub("^.*/|.rds$", "", pSetFiles)
    return(pSetGeneSigs)
}


#' Merge the results data for a list of `PharmacoSig` objects containing different gene signatures into a single,
#'    fully annotated long `data.table` and optionally save to disk as a .csv file
#'
#' @param geneSigL [`list`] of `PharmacoSig` objects, as returned by `readGeneSigsForPSet`
#' @param saveDir Optional [`character`] vector specifying the path at which to save a .csv of the `data.table`
#' @param fileName Optional [`character`] vector specifying the name for the .csv file. If you specify one of `saveDir`
#'     or `fileName`, you must specify both.
#'
#' @return A long [`data.table`] containing the gene signature statistics for each drug x gene x tissue combination
#'
#' @importFrom BiocParallel bplapply bpparam
#' @import data.table
#' @export
mergePSetGeneSigsToDT <- function(geneSigL, saveDir, fileName, BPPARAM=bpparam()) {
    geneSigDTs <- bplapply(geneSigL, convertGeneSignatureToDT, BPPARAM=BPPARAM)
    tissues <- gsub("^.*_", "", names(geneSigL))

    .annotateTissue <- function(DT, tissueName) {DT[["tissue"]] <- rep(tissueName, nrow(DT)); return(DT)}
    # TODO:: Determine if I change geneSigL to geneSigE (environment), can I modify each DT by reference without needing to copy the environment?
    geneSigDTs <- mapply(FUN=.annotateTissue,
                         DT=geneSigDTs, tissueName=tissues,
                         SIMPLIFY=FALSE)

    longGeneSigDT <- rbindlist(geneSigDTs, fill=TRUE)

    # Save to disk or return
    ## TODO:: Add try catch to deal with not existent directories
    if (!missing(saveDir) && !missing(fileName)) {
        fwrite(longGeneSigDT, file=file.path(saveDir, fileName))
    } else if (xor(!missing(saveDir), !missing(fileName))) {
        stop("If you specify one of `saveDir` or `fileName`, you must specify both")
    } else {
        return(longGeneSigDT)
    }

    return(longGeneSigDT)
}
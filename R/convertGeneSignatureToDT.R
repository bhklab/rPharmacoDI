#' Take a SignatureClass object, convert it to a data.table and optionally save 
#'   it to a .csv
#'
#' @param geneSig [`PharmacoSig`] object to convert to a `data.table`
#'
#' @return [data.table] containing the results in `geneSig` as a `data.table 
#'   OR [`NULL`] if you specify `saveDir` and `fileName`, in which case the 
#'   file is saved to disk instead of returned.
#'
#' @import data.table
#' @import PharmacoGx
#' @export
convertGeneSignatureToDT <- function(geneSig) {
    geneSigDT <- data.table(geneSig@.Data[, , ], keep.rownames='gene')
    geneSigDT$dataset <- rep(geneSig@PSetName, nrow(geneSigDT))
    geneSigDT$compound <- geneSig@Arguments$drugs
    geneSigDT$mDataType <- geneSig@Arguments$mDataType
    tissues <- geneSig@Arguments$tissues
    geneSigDT$tissue <- if (length(tissues) == 1) tissues else NA_character_
    return(geneSigDT)
}



#' Takes in the path to a set of gene signature .rds files, selects those 
#'   path matching `pSetPattern` and converts them into a list of PharmacoSig
#'   objects.
#'
#' @param dataDir [`character`] vector specifying the path to the gene signature .rds files
#' @param pSetPattern [`character`] an identifier for the Pset you wish to read. In general this will be the
#'    pSet indentifier (e.g., GDSC_v1). '.*' is automatically prepended and appended to the pattern, so adding
#'    this is not necssary.
#' @param mDataTypes [`character`] directory names for each molecular datatype 
#'    to read gene signatures for.
#' @param BPPARAM [`BiocParam`] object specifying what backend to use for 
#'   parallelization. This defaults to your registered back-end, retrieved via 
#'   the `BiocParallel::bpparam()` function.
#'
#' @return [`list`] of `PharmacoSig` objects for the selected pSet
#'
#' @import data.table
#' @importFrom BiocParallel bplapply bpparam
#' @export
readGeneSigsForPSet <- function(dataDir, pSetPattern, mDataTypes, 
    BPPARAM=BiocParallel::bpparam()) 
{
    dataDirL <- vapply(mDataTypes, function(x, y) paste0(y, '/', x), 
        y=dataDir, FUN.VALUE=character(1))
    pSetFilesL <- lapply(dataDirL, FUN=list.files, 
        pattern=paste0(".*", pSetPattern, ".*"), full.names=TRUE)
    mDTlengths <- vapply(pSetFilesL, length, numeric(1))
    fileMDataTypes <- unlist(mapply(rep, mDataTypes, each=mDTlengths))
    pSetFiles <- unlist(pSetFilesL)
    pSetGeneSigs <- bplapply(pSetFiles, FUN=readGeneSig, BPPARAM=BPPARAM)

    names(pSetGeneSigs) <- paste(fileMDataTypes, 
        gsub("^.*/[^/]*/|.rds$|.qs$", "", pSetFiles), sep='_')
    return(pSetGeneSigs)
}

#' Read a gene signature from a binary file format
#'
#' @details
#' Currently supports .rds (readRDS) or .qs (qread) files.
#'
#' @param path Path to read the file from. Uses regex to extract the file type
#'   and use the correct reader.
#' @param reader An optional function to read files. Use this for unsupported
#'   file formats. Please pass the function directly, not a as a character
#'   string.
#'
#'
#'
#' @export
readGeneSig <- function(path, reader=NA) {
    if (!is.na(reader)) return(reader(path))
    if (grepl('.rds$', path)) readRDS(path, reader)
    else if (grel('.qs$', path)) qread(path)
    else .error("[rPharmacoDI::readGeneSig] Unsupported file format in path!")
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

    .annotateTissue <- function(DT, tissueName) {
        DT[["tissue"]] <- rep(tissueName, nrow(DT)); return(DT) 
    }
    # TODO:: Determine if I change geneSigL to geneSigE (environment), 
    #>can I modify each DT by reference without needing to copy the environment?
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
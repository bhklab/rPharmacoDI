#' @export
setGeneric('stripEnsemblVersion', function(object, ...) 
    standardGeneric("stripEnsemblVersion"))

#' Remove version numbers from 
#' 
#' @param object `PharmacoSet` 
#' 
#' @return The `PharmacoSet` object with all Ensembl ids regexed to remove
#'   the version number.
#' 
#' @details 
#' More specifically, this function will regex the `gene_id` column and the
#' `rownames()` of `rowData()` for all `SummarizedExperiment`s in the 
#' `@molecularProfiles` slot.
#' 
#' @md
#' @importFrom PharmacoGx molecularProfilesSlot
#' @importFrom SummarizedExperiment rowData
#' @import S4Vectors
#' @import BiocGenerics
#' @export
setMethod('stripEnsemblVersion', signature(object='PharmacoSet'), 
    function(object) 
{
    # -- molecularProfiles
    # extract the required data
    SE_list <- molecularProfilesSlot(object)
    rowDataL <- lapply(SE_list, FUN=rowData)

    # check for empty SummarizedExperiments and early return
    if (!all(vapply(rowDataL, nrow, numeric(1)) > 1)) {
        return(object)
    }

    # add gene_id column to rowData if it is missing
    rowDataL <- lapply(rowDataL, within, 
        { if (!exists('gene_id')) { 
            if (exists('EnsemblGeneId')) {
                gene_id <- EnsemblGeneId 
            } else {
                gene_id <- rownames
            }
        }})
    geneIdL <- lapply(rowDataL, `[[`, 'gene_id')
    SE_list <- mapply(`rowData<-`, SE_list, value=rowDataL)

    # make rownames ensembl id then regex off version
    rownamesL <- lapply(SE_list, rownames)
    .all_grepl <- function(x) all(grepl('^ENSG.*', x)) # ENS leaves transcripts
    ensIsRownames <- vapply(rownamesL, FUN=.all_grepl, logical(1))
    rownamesL[!ensIsRownames] <- geneIdL[!ensIsRownames]
    rownamesL <- lapply(rownamesL, FUN=.removeEnsemblVersion)
    SE_list <- mapply(`rownames<-`, x=SE_list, value=rownamesL)

    # drop missing and duplicated ensembl gene id rows
    rownamesL <- lapply(SE_list, rownames)
    keepRowsL <- lapply(rownamesL, function(x) !(is.na(x) | duplicated(x)))
    SE_list <- mapply(`[`, SE_list, i=keepRowsL)

    # update the object
    molecularProfilesSlot(object) <- SE_list
    return(object)
    
})

#' Regex off Ensembl identifier version numbers
#' 
#' @description 
#' Converts Ensembl identifiers from 'ENS*.XX' to 'ENS*', where 'XX' is any two
#' digit version number and '*' is the rest of the identifier text.
#' 
#' @param x `character()` vector of Ensembl gene identifiers
#'   to remove the version numbers from.
#' 
#' @return `character()` vector of Ensembl gene identifiers with the version
#'   numbers removed.
#' 
#' @details
#' Specific regex pattern removed is '\\.[0-9]*$', so it will technically
#' remove more than two numbers. Only does regex for string matcing 'ENS.*'.
#' 
#' @md
.removeEnsemblVersion <- function(x) {
    isMatch <- grepl('ENS.*', x)
    x[isMatch] <- gsub('\\.[0-9]*$', '', x[isMatch])
    x
}
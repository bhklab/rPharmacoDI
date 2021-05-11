# -------------------------------------------------------------------------------------------------
# ---- A collection of misfit functions that didn't make it past QC but which may be  -------------
# ---- recycled, in part or whole, for use in other functions/projects ----------------------------
# -------------------------------------------------------------------------------------------------

# ---- Make sensitivity slot into a single long data.table; used to much disk space to be ---------
# ---- practical ----------------------------------------------------------------------------------

#' Merge all the data in the `sensitivity` slot into a single long format `data.table`
#'
#' @param sensSlot [`list`] PharmacoSet sensitivity slot, as returned by `sensitivitySlot`.
#'
#' @return [`data.table`] Long format of data in sensitivity slot
#'
#' @noRd
.sensSlotToLong <- function(sensSlot) {

    # -- sensitivityRaw
    rawDTdose <- data.table(sensSlot$raw[,,1], keep.rownames='.exp_id')
    rawDTviab <- data.table(sensSlot$raw[,,2], keep.rownames='.exp_id')

    moltenRawDose <- melt.data.table(rawDTdose, id.vars='.exp_id', value.vars=colnames(sensSlot$raw[,,1]),
                                     variable.name='.sample', value.name='concentration')
    moltenRawViab <- melt.data.table(rawDTviab, id.vars='.exp_id', value.vars=colnames(sensSlot[,,2]),
                                     variable.name='.sample', value.name='viability')

    # Memory manage - not triggering gc() because it is slow, could if we run out of system memory
    rm(rawDTdose, rawDTviab)

    # Set keys then join on them
    setkey(moltenRawDose, .exp_id, .sample)
    setkey(moltenRawViab, .exp_id, .sample)
    rawDT <- merge.data.table(moltenRawDose, moltenRawViab)

    # Memory manage
    rm(moltenRawDose, moltenRawViab)

    # -- sensitivityInfo
    infoDT <- data.table(sensSlot$info, keep.rownames=".exp_id")
    colnames(infoDT)[2:ncol(infoDT)] <- paste0('info_', colnames(infoDT)[2:ncol(infoDT)])

    # -- sensitivityProfiles
    profDT <- data.table(sensSlot$profiles, keep.rownames=".exp_id")
    colnames(profDT)[2:ncol(profDT)] <- paste0('prof_', colnames(profDT)[2:ncol(profDT)])


    # -- sensNumber
    numDT <- data.table(sensSlot$n, keep.rownames='.cancer_type')
    moltenNumDT <- melt(numDT, id.vars='.cancer_type', measure.vars=colnames(sensSlot$n),
                        variable.name=".drug_id", value.name="n")
    rm(numDT)

    # -- longSensDT

    # Join sensitivityInfo with sensitivityProfiles
    setkey(infoDT, .exp_id)
    setkey(profDT, .exp_id)
    annotDT <- merge.data.table(infoDT, profDT)

    # Memory manage
    rm(infoDT, profDT)

    # Join annotDT with numDT
    setkey(annotDT, info_cellid, info_drugid)
    setkey(moltenNumDT, .cancer_type, .drug_id)
    metaDT <- merge.data.table(annotDT, moltenNumDT,
                                 by.x=c('info_cellid', 'info_drugid'),
                                 by.y=c('.cancer_type', '.drug_id'))

    # Memory manage
    rm(annotDT, moltenNumDT)

    # Join annotsDT with rawDT
    setkey(metaDT, .exp_id)
    setkey(rawDT, .exp_id)
    longSensDT <- merge.data.table(rawDT, metaDT)

    return(longSensDT)
}


# ---- Implementation of `as` function for converting `SummarizedExperiment` objects to -----------
# ---- a single long `data.table`; to0 much disk usage to be practical -----------------------------

#' Coerce a SummarizedExperiment object to a long data.table, retaining the data in all assays, rowData, colData
#'
#' @param from [`SummarizedExperiment`]
#' @param to [`character`]
#'
#' @return ['data.table]
#'
#' @importFrom SummarizedExperiment colData rowData
#' @import data.table
#'
#' @noRd
setMethod("coerce",
    signature(from="SummarizedExperiment", to="data.table"),
    function(from, to) {

    # Get the data to join on
    assaysDT <- .assaysToLongDT(assays(from), names(assays(from)))

    rowDataDT <- data.table(as(rowData(from), 'data.frame'), keep.rownames='features')
    colnames(rowDataDT)[-1] <- paste0('row_', colnames(rowDataDT)[-1])

    colDataDT <- data.table(as(colData(from), 'data.frame'), keep.rownames='samples')


    ## TODO:: Determine if there are any other items we need from `metadata`
    if('protocalData' %in% names(metadata(from))) {
                protocolDataDT <- data.table(as(metadata(from)$protocolData, 'data.frame'),
                                             keep.rownames='samples')
        colDataDT <- merge.data.table(colDataDT, protocolDataDT, by='samples')
    }

    colnames(colDataDT)[-1] <- paste0('col_', colnames(colDataDT)[-1])

    # Join colData and rowData to the long assaysDT
    setkey(assaysDT, features)
    setkey(rowDataDT, features)
    longSummarizedExperiment <- merge.data.table(assaysDT, rowDataDT,
                                                 allow.cartesion=TRUE, fill.missing=TRUE)
    setkey(longSummarizedExperiment, samples)
    setkey(colDataDT, samples)
    longSummarizedExperiment <- merge.data.table(longSummarizedExperiment, colDataDT,
                                                 allow.cartesian=TRUE, fill.missing=TRUE)

    return(longSummarizedExperiment)
})

# Helpers - coerce method

#' Converts each assay in a list of assays to a data.table, then iteratively merged the data.tables by the shared
#'    feature (rownames) and samples (colnames) columns.
#'
#' @param assays [`list`] A list of `matrix` objects, as returned by the `SummarizedExperiment::assays` function.
#' @param assayNames [`character`] Names for each assay, e.g. names of the list returned by
#'     `SummarizedExperiment::assays`
#'
#' @import data.table
#' @importFrom SummarizedExperiment assays
#'
#' @noRd
.assaysToLongDT <- function(assays, assayNames) {
        # Convert each assay to a data.table, return a list
        assaysDtL <- mapply(FUN=.assayToLongDT,
                            assays, assayNames,
                            SIMPLIFY=FALSE)

        # Metaprogram some R code to join an unknown number of data.tables by key
        codeDataTables <- paste0('assaysDtL[[', seq_along(assaysDtL), ']]')
        if (length(codeDataTables) - 2 > 0) {
            codeOperators <- c('[', rep('][', length(codeDataTables) - 2), ']')
        } else {
            codeOperators <- c('[', ']')
        }

        zippedStrings <- unlist(mapply(c, codeDataTables, codeOperators, SIMPLIFY=FALSE))
        chainedJoinExpression <- parse(text=paste0(zippedStrings, collapse=''))

        # Evaluate the R code, this allows merging n data.tables without using inefficient Reduce function
        assayDT <- eval(chainedJoinExpression)

        return(assayDT)
}

#' Coerce an assay matrix from a SummarizedExperiment to a long data.table
#'
#' @param assay [`matrix`] A numeric matrix with features as rows and samples as columns, as returend by
#'     `SummarizedExperiment::assay`
#' @param assayName [`character`] The name of the assay, this becomes the column name for the assay values
#'
#' @return [`data.table`] where the columns of the assay values are in column named `assayName` and the
#'     column names are in the `samples` column
#'
#' @noRd
.assayToLongDT <- function(assay, assayName) {
    DT <- data.table(assay, keep.rownames='features')
    DTmolten = melt.data.table(DT,  id.vars='features', value.vars=rownames(assay),
                                variable.name='samples', value.name=paste0('assay_', assayName)) # So we can identify assays vs annotations
    data.table::setkey(DTmolten, features, samples) # Set keys for joins
    return(DTmolten)
}

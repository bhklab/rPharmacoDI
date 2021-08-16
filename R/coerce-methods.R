#' @title Coerce a SummarizedExperiment to a data.table
#' 
#' @examples 
#' SE <- molecularProfilesSlot(GDSCsmall)[[1]]
#' as(SE, 'data.table')
#' 
#' @param from `SummarizedExperiment` object.
#' 
#' @return `data.table` with long format of data in `from`
#' 
#' @importFrom data.table as.data.table melt.data.table merge.data.table
#' @export
setAs(from='SummarizedExperiment', to='data.table', function(from) {

    # -- extract sample metadata
    colDT <- as.data.table(colData(from), keep.rownames='.sample')

    # -- extract feature metadata
    rowDT <- as.data.table(rowData(from), keep.rownames='.feature')

    # -- extract and process assays
    assayL <- assays(from)
    assayDtL <- lapply(assayL, as.data.table, keep.rownames='.feature')
    meltDtL <- lapply(assayDtL, melt, id.vars='.feature', 
        variable.name='.sample', variable.factor=FALSE)
    assayDT <- meltDtL[[1]][, .(.sample, .feature)]
    for (i in seq_along(meltDtL)) 
        assayDT[[names(assayL)[[i]]]] <- meltDtL[[i]][['value']]

    # -- merge into a single long format table
    DT <- merge.data.table(assayDT, colDT, by='.sample')
    DT <- merge.data.table(DT, rowDT, by='.feature')

    # -- add metadata
    metadata <- metadata(from)
    notS4 <- !vapply(metadata, isS4, logical(1))
    if (!all(notS4)) .warning('Dropped S4 metadata during coercion to data.table!')
    for (name in names(metadata)[notS4]) DT[[name]] <- metadata[[name]]

    return(DT)
})

#' @title Coerce a SummarizedExperiment to a data.frame
#' 
#' @examples 
#' SE <- molecularProfileSlot(GDSCsmall)[[1]]
#' as(SE, 'data.frame')
#' 
#' @param from `SummarizedExperiment` object.
#' 
#' @return `data.frame` with long format of data in `from`.
#' 
#' @importFrom data.table as.data.table melt.data.table merge.data.table
#' @export
setAs(from='SummarizedExperiment', to='data.frame', function(from) {
    setDF(as(from, 'data.table'))
})
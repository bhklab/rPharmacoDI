#' @include writeToCsv-methods.R
NULL

#' Convert all data in an `PharmacoSet` object to disk as .parquet files,
#'   named according to the nested location in the PharmacoSet
#'
#' @section Note:
#'   This function will create a directory with the same name as the `object`
#'   being written to .parquet, unless one already exists. The 
#    `filePath` is automatically appended with the object name, therefore it is 
#'   not necssary to specify the object name in `filePath`
#'
#' @param object [`S4`] An object inheriting from the `CoreGx::CoreSet` class, such as a
#'  `PharmacoSet`, `RadioSet`, `ToxicoSet` or `XevaSet`.
#' @param filePath [`character`] The path to which the .parquet files will be written
#' @param ... `pairlist` Fall through arguments to all instances of 
#'  `arrow::write_parquet`
#'
#' @import arrow
#' @import PharmacoGx
#' @export
setMethod("writeToParquet", signature(object="PharmacoSet"), function(object, filePath, ...) {

    objectName <- name(object)
    pSetDir <- paste0(objectName, '_PSet')

    tryCatch({ dir.create(file.path(filePath, pSetDir)) },
        warning=function(w) message(paste0('\n', pSetDir,
            ' directory already exists, writing .parquet files there\n'))
    )

    if (!grepl(paste0('.*', objectName, '$'), filePath)) {
        filePath <- file.path(filePath, pSetDir)
    } else {
        message(paste0('\nFYI: It is not necessay to specify the ',
            pSetDir, 'directroy in `filePath`\n. We have
            already do that for you!\n'))
    }

    message(paste0('Writing ', objectName, ' to parquet in: \n\t', filePath, '\n'))

    message("       --> Writing molecularProfiles slot to disk\n")
    .writeMolecularProfilesToParquet(molecularProfilesSlot(object), filePath, objectName, ...)

    # annotation and datasetType
    message("       --> Writing annotations slot to disk\n")
    .writeAnnotationToTxt(annotation(object), datasetType(object), filePath, objectName)

    if (datasetType(object) %in% c('sensitivity', 'both'))
        message("       --> Writing sensitivity slot to disk\n")
        .writeSensitivityToParquet(sensitivitySlot(object), filePath, objectName, ...)

    if (datasetType(object) %in% c('perturbation', 'both'))
        ## TODO:: implement this for perturbation dataset!
        message("       SORRY: The method for writing the perturbation slot to .csv has not been implemented yet!\n")
        #message("   Writing perturbation slot to disk")
        #.writePerturbationToCsv(slot(object, 'perturbation'), filePath, objectName)

    .writeDFslotToParquet <- function(slotDF, filePath, objectName, fileSuffix)
        arrow::write_parquet(data.table(slotDF, keep.rownames='rownames') ,
            sink=file.path(filePath, paste0(objectName, fileSuffix)), ...)

    # -- drug
    message("       --> Writing drug slot to disk\n")
    .writeDFslotToParquet(slot(object, 'drug'), filePath, objectName, '@drug.parquet')

    # -- curation
    .writeCurationSlotToParquet <- function(curationSlot, filePath, objectName) {
        curationNames <- names(curationSlot)
        for (i in seq_along(curationSlot)) {
            .writeDFslotToParquet(curationSlot[[i]], filePath, objectName, 
                paste0('@curation$', curationNames[[i]], '.parquet'))
        }
    }

    message("       --> Writing curation slot to disk\n")
    .writeCurationSlotToParquet(curation(object), filePath, objectName)

    message("       --> Writing cell slot to disk\n")
    .writeDFslotToParquet(cellInfo(object), filePath, objectName, '@cell.parquet')

    message('SUCCESS: Your files are in ', filePath, '\n\n')
})

#' Convert each SummarizedExperiment in a CSet into a list of data.tables, then save to disk as .csv files
#'
#' @param SElist A [`list`] of `SummarizedExperiment` objects, as returned by `molecularProfilesSlot`
#' @param filePath [`character`] Path to save the files .csv files
#' @param objectName [`character`] The name of the `CSet` object being written to disk
#' @param ... `pairlist` Fall through arguments to `arrow::write_parquet`. 
#'
#' @import arrow
#' @export
.writeMolecularProfilesToParquet <- function(SElist, filePath, objectName, ...) {
    sumExperDtList <- lapply(SElist, .convertSEToDataTableList)
    sumExperNames <- names(sumExperDtList)

    # -- loop over summarized experiments
    for (i in seq_along(sumExperDtList)) {
        sumExperName <- sumExperNames[i]
        dataTableList <- sumExperDtList[[i]]
        dataTableNames <- names(dataTableList)

        # loop over each data.table, save to disk with appropraite name
        for (j in seq_along(dataTableNames)) {
            arrow::write_parquet(dataTableList[[j]],
                sink=file.path(filePath,
                    paste0(objectName, '@molecularProfiles$', # pset + slot
                        sumExperName, '$', # summarized experiment name
                        dataTableNames[j], '.parquet') # table name
                    ),
                ...
                )
        }
    }
}

#' Convert a PSets' sensitivty slot into a long data.table and save to disk as a .csv
#'
#' @param sensSlot [`list`] PharmacoSet sensitivity slot, as returned by `sensitivitySlot`.
#'
#' @keywords internal
#' @export
.writeSensitivityToParquet <- function(sensSlot, filePath, objectName, ...) {
    sensSlotDTs <- .sensSlotToDataTables(sensSlot)
    for (i in seq_along(sensSlotDTs))
        arrow::write_parquet(sensSlotDTs[[i]],
            sink=file.path(filePath,
                paste0(objectName, '@sensitivity$', names(sensSlotDTs)[i], 
                    '.parquet')),
            ...)
}
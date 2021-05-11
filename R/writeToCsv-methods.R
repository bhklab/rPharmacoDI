#' Convert all data in an `PharmacoSet` object to disk as .csv files, named per
#'
#' @section Note:
#'   This function will create a directory with the same name as the `object`
#'   being written to .csv, unless one already exists. The `filePath` is automatically
#'   appended with the object name, therefore it is not necssary to specify the object
#'   name in `filePath`
#'
#' @param object [`S4`] An object inheriting from the `CoreGx::CoreSet` class, such as a
#'  `PharmacoSet`, `RadioSet`, `ToxicoSet` or `XevaSet`.
#' @param filePath [`character`] The path to which the .csv files will be written
#'
#' @import data.table
#' @import PharmacoGx
#' @export
setMethod("writeToCsv", signature(object="PharmacoSet"), function(object, filePath) {

    objectName <- name(object)
    pSetDir <- paste0(objectName, '_PSet')

    tryCatch({ dir.create(file.path(filePath, pSetDir)) },
             warning=function(w)
                 message(paste0('\n', pSetDir,
                                ' directory already exists, writing .csv files there\n'))
    )

    if (!grepl(paste0('.*', objectName, '$'), filePath)) {
        filePath <- file.path(filePath, pSetDir)
    } else {
        message(paste0('\nFYI: It is not necessay to specify the ',
                       pSetDir, 'directroy in `filePath\n. We have
                       already do that for you :-)\n'))
    }

    message(paste0('Writing ', objectName, ' to csv in: \n\t', filePath, '\n'))

    message("       --> Writing molecularProfiles slot to disk\n")
    .writeMolecularProfilesToCsv(molecularProfilesSlot(object), filePath, objectName)

    # annotation and datasetType
    message("       --> Writing annotations slot to disk\n")
    .writeAnnotationToTxt(annotation(object), datasetType(object), filePath, objectName)

    if (datasetType(object) %in% c('sensitivity', 'both'))
        message("       --> Writing sensitivity slot to disk\n")
        .writeSensitivityToCsv(sensitivitySlot(object), filePath, objectName)

    if (datasetType(object) %in% c('perturbation', 'both'))
        ## TODO:: implement this for perturbation dataset!
        message("       SORRY: The method for writing the perturbation slot to .csv has not been implemented yet!\n")
        #message("   Writing perturbation slot to disk")
        #.writePerturbationToCsv(slot(object, 'perturbation'), filePath, objectName)

    .writeDFslotToCsv <- function(slotDF, filePath, objectName, fileSuffix)
        fwrite(data.table(slotDF, keep.rownames='rownames') ,
               file=file.path(filePath, paste0(objectName, fileSuffix)), compress='gzip')

    # -- drug
    message("       --> Writing drug slot to disk\n")
    .writeDFslotToCsv(slot(object, 'drug'), filePath, objectName, '@drug.csv.gz')

    # -- curation
    .writeCurationSlotToCsv <- function(curationSlot, filePath, objectName) {
        curationNames <- names(curationSlot)
        for (i in seq_along(curationSlot)) {
            .writeDFslotToCsv(curationSlot[[i]], filePath, objectName, paste0('@curation$', curationNames[[i]], '.csv.gz'))
        }
    }

    message("       --> Writing curation slot to disk\n")
    .writeCurationSlotToCsv(curation(object), filePath, objectName)

    message("       --> Writing cell slot to disk\n")
    .writeDFslotToCsv(cellInfo(object), filePath, objectName, '@cell.csv.gz')

    message('SUCCESS: Your files are in ', filePath, '\n\n')
})


#' Convert each SummarizedExperiment in a CSet into a list of data.tables, then save to disk as .csv files
#'
#' @param SElist A [`list`] of `SummarizedExperiment` objects, as returned by `molecularProfilesSlot`
#' @param filePath [`character`] Path to save the files .csv files
#' @param objectName [`character`] The name of the `CSet` object being written to disk
#'
#' @import data.table
#' @export
.writeMolecularProfilesToCsv <- function(SElist, filePath, objectName) {
    sumExperDtList <- lapply(SElist, .convertSEToDataTableList)
    sumExperNames <- names(sumExperDtList)

    # -- loop over summarized experiments
    for (i in seq_along(sumExperDtList)) {
        sumExperName <- sumExperNames[i]
        dataTableList <- sumExperDtList[[i]]
        dataTableNames <- names(dataTableList)

        # loop over each data.table, save to disk with appropraite name
        for (j in seq_along(dataTableNames)) {
            fwrite(dataTableList[[j]],
                   file=file.path(filePath,
                                  paste0(objectName, '@molecularProfiles$', # pset + slot
                                         sumExperName, '$', # summarized experiment name
                                         dataTableNames[j], '.csv.gz') # table name
                   ),
                   compress='gzip')
        }
    }
}


#' Convert the data in a SummarizedExperiment into a `list` of `data.table`s
#'
#' @section Note:
#'   Any data stored in `metadata(SummarizedExperiment)`
#'   will will be converted to an ASCII representation
#'
#' @param SE A [`SummarizedExperiment`] whos data will be converted to a `list`
#'   of `data.table`s named for their source within the object.
#'
#' @keywords internal
#' @export
.convertSEToDataTableList <- function(SE) {

    if (length(SE) == 0) return(data.table())

    # -- feature and sample annotations
    .s4DataFrameToDT <- function(DF, rownameLabel)
        data.table(as(DF, 'data.frame'), keep.rownames=rownameLabel)

    rowDataDT <- .s4DataFrameToDT(rowData(SE), '.features')

    colDataDT <- .s4DataFrameToDT(colData(SE), '.samples')

    .matrixToDT <- function(assay, rownameLabel)
        data.table(assay, keep.rownames=rownameLabel)

    # -- assay data
    assaysDtL <- mapply(FUN=.matrixToDT,
                        assay=assays(SE), rownameLabel='.features', # Arguments to function
                        SIMPLIFY=FALSE)
    names(assaysDtL) <- paste0('assay.', names(assaysDtL)) # Identifier for file name

    # -- metadata
    # Converts R code needed to recreate each object as a string; can string parse the data out of them in Python
    #   or simply use `eval(parse(text=<code as string>))`; for S4 classes, there may be issues recreating them
    #   due to different package versions using different constructor synatx
    .captureCodeAsString <- function(object) paste0(capture.output(dput(object)), collapse='')

    metadataDT <- as.data.table(lapply(metadata(SE), .captureCodeAsString))

    # -- Merge lists and return
    return(c(list('rowData'=rowDataDT, 'colData'=colDataDT, 'metadata'=metadataDT),
             assaysDtL))
}


#' Write the annotation listt in a PharmacoSet to a text file
#'
#' Each list item is one line, lines in each item are separated with `sepChar` so they
#'  can be split into to their original stucture.
#'
#' @param annotations [`list`] Annotations as returned by the `annotaiton` function
#' @param path [`character`] The path to save the output file to.
#' @param objectName [`character`] The name of the PSet the annotations are from.
#' @param sepChar [`character`] Separator to use when pasting together multiple lines of a list item.
#'
#' @return Writes to disk, does not return.
#'
#' @keywords internal
#' @export
.writeAnnotationToTxt <- function(annotations, dsType, filePath, objectName, sepChar="|||") {

    file <- file.path(filePath, paste0(objectName, '@annotations.txt'))

    # -- Date Created
    dateCreated <- annotations$dateCreated

    # -- Session Info
    sessInfo <- paste0(capture.output(annotations$sessionInfo), collapse=sepChar)

    # -- Call
    creationCall <- paste0(capture.output(annotations$call), collapse=sepChar)

    # -- Version
    version <- annotations$version

    annots <- rbind(objectName, dateCreated, sessInfo, creationCall, version, dsType)

    write.table(annots, file=file, sep="\n", row.names=FALSE, col.names=FALSE)
}


#' Convert a PSets' sensitivty slot into a long data.table and save to disk as a .csv
#'
#' @param sensSlot [`list`] PharmacoSet sensitivity slot, as returned by `sensitivitySlot`.
#'
#' @keywords internal
#' @export
.writeSensitivityToCsv <- function(sensSlot, filePath, objectName) {
    sensSlotDTs <- .sensSlotToDataTables(sensSlot)
    for (i in seq_along(sensSlotDTs))
        fwrite(sensSlotDTs[[i]],
               file=file.path(filePath,
                              paste0(objectName, '@sensitivity$', names(sensSlotDTs)[i], '.csv.gz')),
               compress='gzip')
}


#' Convert the sensitivity slot of a PSet into a list of data.tables
#'
#' Use this function to make writing data in the sensitivity slot to disk easy.
#'
#' @param sensSlot [`list`] The sensitivity slot of a PSet, as returned by `sensitivitySlot`
#'
#' @import data.table
#'
#' @keywords internal
#' @export
.sensSlotToDataTables <- function(sensSlot) {

   # -- raw
   .array3rdDimToDataTable <- function(idx, array, rownameLabel)
       data.table(array[,,idx], keep.rownames=rownameLabel)
   sensRaw <- sensSlot$raw
   sensData <- lapply(seq_len(dim(sensRaw)[3]),
                      FUN=.array3rdDimToDataTable,
                      # Arguments to function
                      array=sensRaw,
                      rownameLabel='.exp_id')
   names(sensData) <- paste0('raw.', dimnames(sensRaw)[[3]])

   # -- info, prof, n
   sensMetadata <- lapply(sensSlot[c('info', 'profiles', 'n')], data.table, keep.rownames='.rownames')

   # -- merge lists & return
   return(c(sensData, sensMetadata))
}

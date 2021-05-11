#' Dump gene annotation files to .csv for use in PharmacoDI
#'
#' @param filePath [`character`] vector indiciating path to the annotation directory
#' @param pattern [`character`] files to match using `pattern` argument in `list.files` function
#' @return NULL Writes files to `filePath`
#'
#' @import data.table
geneAnnotationToCSV <- function(filePath=file.path("../..", "data", "metadata"), pattern="Gencode.v33.annotation.*") {
    files <- list.files(filePath, pattern=paste0(pattern, '.RData'))
    fileNames <- list()
    for (file in files) {
        fileName <- load(file, verbose=TRUE)
        fileNames[[file]] <- fileName
    }

    for (names in fileNames) {
        for (name in names) {
            message(paste0("Saving ", name, " to ",
                           paste0(gsub("\\*", '', pattern), name, ".csv")))
            fwrite(data.table(get(name), keep.rownames="rownames"),
                   file=file.path(filePath, paste0(gsub("\\*", '', pattern), name, ".csv")))
        }
    }
}
#' Parse gene signature metaanalysis .rds files into a data.table.
#' 
#' @param dataDir `character(1)` Path to the top level metaanalysis directory. Folders in
#'   this directory should be named to match '*AnalyticPancan', where '*' is the molecular
#'   data type of that signature file.
#' @param saveDir `character(1)` Path to write the output .csv to.
#' 
#' @import PharmacoGx
#' @importFrom BiocParallel bplapply bpparam
#' @importFrom data.table data.table rbindlist
#' @export
convertMetaAnalysisToDT <- function(dataDir, saveDir, fileName='gene_compound.csv') {
    sigDirs <- list.files(dataDir, pattern='*AnalyticPancan', full.names=TRUE)
    sigFiles <- unlist(lapply(sigDirs, list.files, pattern='*all.rds', 
        full.names=TRUE, recursive=TRUE))
    sigList <- bplapply(sigFiles, readRDS)
    dataTableL <- bplapply(sigList, FUN=convertGeneSignatureToDT)
    metaDT <- rbindlist(dataTableL, fill=TRUE, use.names=TRUE)
    # Fix CCLE.CTRPv2
    metaDT[pSet == 'CCLE.CTRPv2', pSet := 'CTRPv2']
    fwrite(metaDT, file=file.path(saveDir, fileName))
}


if (sys.nframe() == 0 ) {
    dataDir <- '../PharmacoDI_snakemake_pipeline/rawdata/gene_signatures/metaanalysis'
    convertMetaAnalysisToDT(dataDir, saveDir=dataDir)
}
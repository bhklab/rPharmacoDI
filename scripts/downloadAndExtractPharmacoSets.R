# -- Dependencies
library(rPharmacoDI)
library(PharmacoGx)
library(data.table)
library(BiocParallel)
library(MultiAssayExperiment)

# -- Configuration
nthread <- 14
filePath <- '../PharmacoDI_snakemake_pipeline/rawdata'

# data.table config
setDTthreads(nthread)

# BiocParalle config
bp <- bpparam()
bpworkers(bp) <- nthread
bpprogressbar(bp) <- TRUE
register(bp)

# -- Download the PSets
canonicalPSetDF <- PharmacoGx::availablePSets()
pSetNames <- canonicalPSetDF$`PSet Name`
canonicalPSets <- bptry(bplapply(pSetNames, FUN=downloadPSet, saveDir=filePath, 
    timeout=1e10))
if (!all(bpok(canonicalPSets))) {
    canonicalPSets <- btry(bplapply(canonicalPSetDF$`PSet Name`, FUN=downloadPSet, 
        saveDir=filePath, BPREDO=canonicalPSets, timeout=1e20))
    if (!all(bpok(canonicalPSets))) stop("Downloading PSets failed!")
}

# -- Preprocess the PSets to use the correct molecular identifiers
procCanonicalPSets <- bplapply(canonicalPSets, FUN=stripEnsemblVersion)

## FIXME:: Remove non-UTF byte from drug metadata brand name drugs[4, 3]
##>this breaks the Snakemake pipeline everytime

# -- Extract into filePath
# This is technically bad practice, because I am using a functional looping construct
#   for it's side effects only. But it is the easiest way to parallelize this.
# TODO:: Does this break our RAM usage? No put it peaks around 80 GB, so it will break
#   if we ever lower this VMs RAM below that. Can I check that from R?
bplapply(procCanonicalPSets, FUN=writeToParquet, filePath=filePath)


## ---- Testing
if (sys.nframe() == 0) {
    object <- readRDS('rawdata/GDSC_2020(v1-8.2).rds')
    MAE <- MultiAssayExperiment(molecularProfilesSlot(object))
    colDataL <- lapply(experiments(MAE), function(x) as.data.table(colData(x)))
    colDT <- rbindlist(colDataL, fill=TRUE, use.names=TRUE, idcol='mDataType')
    colMetaDT <- colDT[, lapply(.SD, 
        function(x) paste0(unique(na.omit(x)), collapse='|')), by=rownames]
    rowDataL <- lapply(experiments(MAE), function(x) as.data.table(rowData(x)))
    rowDT <- rbindlist(rowDataL, fill=TRUE, use.names=TRUE, idcol='mDataType')
    rowMetaDT <- rowDT[, lapply(.SD, 
        function(x) paste0(unique(na.omit(x)), collapse='|')), by=rownames]
    assayL <- lapply(assays(MAE), as.data.table, keep.rownames='.feature')
    flatAssayL <- lapply(assayL, melt.data.table, id.vars='.feature', 
        variable.name='.sample', value.factor='false', variable.factor=FALSE)
    for (i in seq_along(flatAssayL)) setnames(flatAssayL[[i]], 'value', 
        names(assayL)[i])
    .merge_long_arrays <- function(x, y) 
        merge.data.table(x, y, by=c('.sample', '.feature'), all=TRUE)
    assayDT <- Reduce(.merge_long_arrays, flatAssayL)
    setnames(rowMetaDT, 'rownames', '.feature')
    setnames(colMetaDT, 'rownames', '.sample')
    MAE_DT <- merge.data.table(assayDT, colMetaDT, by='.sample')
    MAE_DT <- merge.data.table(MAE_DT, rowMetaDT, by='.feature')
}
# -- Dependencies
library(rPharmacoDI)
library(PharmacoGx)
library(data.table)
library(BiocParallel)
library(MultiAssayExperiment)
library(qs)

# -- Configuration
nthread <- 12
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
if (!file.exists(file.path('local_data', 'canonicalPSets.qs'))) {
    canonicalPSets <- bptry(bplapply(pSetNames, FUN=downloadPSet, saveDir=filePath, 
        timeout=1e10))
    canonicalPSets <- vector('list', length(pSetNames))
    # for (i in seq_along(pSetNames)) {
    #     canonicalPSets[[i]] <- downloadPSet(pSetNames[i], saveDir=filePath, 
    #         timeout=1e10)
    # }
    if (!all(bpok(canonicalPSets)))
        canonicalPSets <- btry(bplapply(canonicalPSetDF$`PSet Name`, 
            FUN=downloadPSet, saveDir=filePath, BPREDO=canonicalPSets, 
            timeout=1e20))
    if (!all(bpok(canonicalPSets))) stop("Downloading PSets failed!")
} else {
    canonicalPSets <- qread(file.path('local_data', 'canonicalPSets.qs'),
        nthread=nthread)
}

# -- Preprocess the PSets to use the correct molecular identifiers
procCanonicalPSets <- bplapply(canonicalPSets, FUN=stripEnsemblVersion)
rm(canonicalPSets); gc()
names(procCanonicalPSets) <- pSetNames

## FIXME:: Remove non-UTF byte CCLE from drug metadata brand name drugs[4, 2]
##>this breaks the Snakemake pipeline everytime
drugInfo(procCanonicalPSets[["CCLE_2015"]]) <- removeNonASCII(
        drugInfo(procCanonicalPSets[["CCLE_2015"]]))


# -- Extract into filePath
# This is technically bad practice, because I am using a functional looping 
#   construct for it's side effects only. But it is the easiest way to 
#   parallelize this.
# NOTE:: Peaks at ~55 GB RAM usage on 12 threads
bplapply(procCanonicalPSets, FUN=writeToParquet, filePath=filePath)


## ---- Testing
if (sys.nframe() == 0) {

}
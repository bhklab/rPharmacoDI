# -- Dependencies
library(rPharmacoDI)
library(PharmacoGx)
library(data.table)
library(BiocParallel)

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
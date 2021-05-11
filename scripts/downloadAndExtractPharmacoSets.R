# -- Dependencies
library(rPharmacoDI)
library(PharmaocoGx)
library(data.table)
library(BiocParallel)

# -- Configuration
nthread <- 14
filePath <- '../PharmacoDI_snakemake_pipeline/rawdata'

setDTthreads(nthread)


# -- Download the PSets
canonicalPSetDF <- PharmacoGx::availablePSets()
canonicalPSets <- bplapply(canonicalPSetDF$`PSet Name`, FUN=downloadPSet, saveDir=filePath)


# -- Extract into filePath
# This is technically bad practice, because I am using a function looping construct
#   for it's side effects only. But it is the easiest way to parallelize this.
# TODO:: Does this break our RAM usage?
bplapply(canonicalPSets, FUN=writeToCsv, filePath=filePath)
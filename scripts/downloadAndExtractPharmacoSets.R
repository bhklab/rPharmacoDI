# -- Dependencies
library(rPharmacoDI)
library(PharmacoGx)
library(data.table)
library(BiocParallel)
library(MultiAssayExperiment)
library(qs)

# -- Configuration
nthread <- 12
filePath <- "../PharmacoDI_snakemake_pipeline/rawdata"

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
if (!file.exists(file.path("local_data", "canonicalPSets.qs"))) {
    canonicalPSets <- bptry(bplapply(pSetNames, FUN=downloadPSet, saveDir=filePath, 
        timeout=1e10))
    if (!all(bpok(canonicalPSets)))
        canonicalPSets <- btry(bplapply(canonicalPSetDF$`PSet Name`, 
            FUN=downloadPSet, saveDir=filePath, BPREDO=canonicalPSets, 
            timeout=1e20))
    if (!all(bpok(canonicalPSets))) stop("Downloading PSets failed!")
    qsave(canonicalPSets, file=file.path("local_data", "canonicalPSets.qs"))
} else {
    canonicalPSets <- qread(file.path("local_data", "canonicalPSets.qs"),
        nthread=nthread)
}

# -- Preprocess the PSets to use the correct molecular identifiers
procCanonicalPSets <- bplapply(canonicalPSets, FUN=stripEnsemblVersion)
#rm(canonicalPSets); gc()
names(procCanonicalPSets) <- pSetNames

## FIXME:: Remove non-UTF byte CCLE from drug metadata brand name drugs[4, 2]
##>this breaks the Snakemake pipeline everytime
drugInfo(procCanonicalPSets[["CCLE_2015"]]) <- removeNonASCII(
        drugInfo(procCanonicalPSets[["CCLE_2015"]]))

## FIXME:: Remove this once the drugid columns are updated in the respective
##>PSets
updateDrugIDs <- function(x, old, new) {
    # Extract relevant tables as data.tables
    drugInf <- as.data.table(drugInfo(x))
    drugCur <- as.data.table(curation(x)$drug)
    sensInf <- as.data.table(sensitivityInfo(x), keep.rownames=TRUE)

    # Replace the old drug identifiers with new ones
    for (i in seq_along(old)) {
        drugInf[drugid == old[i], drugid := new[i]]
        drugCur[unique.drugid == old[i], unique.drugid := new[i]]
        sensInf[drugid == old[i], drugid := new[i]]
        sensInf[, rn := gsub(old[i], new[i], rn, fixed=TRUE)]
        colnames(sensNumber(x))[colnames(sensNumber(x)) == old[i]] <- new[i]
    }

    # Convert back to data.frames and fix the rownames
    setDF(drugInf)
    setDF(drugCur)
    setDF(sensInf)
    rownames(drugInf) <- drugInf$drugid
    rownames(drugCur) <- drugCur$unique.drugid
    rownames(sensInf) <- sensInf$rn
    rownames(sensitivityProfiles(x)) <- rownames(sensInf)
    rownames(sensitivityRaw(x)) <- rownames(sensInf)

    # Assign back to the PSet
    drugInfo(x) <- drugInf
    curation(x)$drug <- drugCur
    sensitivityInfo(x) <- sensInf[, colnames(sensInf) != 'rn']

    return(x)
}

mapToUniqueDrugId <- fread(file.path("local_data", 
    "drugid_not_in_drugs_with_ids.csv"))
for (pSet in mapToUniqueDrugId$dataset) {
    procCanonicalPSets[[pSet]] <- updateDrugIDs(
        procCanonicalPSets[[pSet]],
        old=mapToUniqueDrugId[dataset == pSet, ]$drugid,
        new=mapToUniqueDrugId[dataset == pSet, ]$unique.drugid
    )
}

# -- Extract into filePath
# This is technically bad practice, because I am using a functional looping 
#   construct for it"s side effects only. But it is the easiest way to 
#   parallelize this.
# NOTE:: Peaks at ~55 GB RAM usage on 12 threads
bplapply(procCanonicalPSets, FUN=writeToParquet, filePath=filePath)


## ---- Testing
if (sys.nframe() == 0) {

}
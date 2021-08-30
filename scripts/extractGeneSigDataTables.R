library(rPharmacoDI)
library(BiocParallel)
library(data.table)
library(qs)


# ---- 0. Script configuration

inputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/signature_files'
outputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/gene_signatures'

BP <- bpparam()
bpprogressbar(BP) <- TRUE
register(BP)

setDTthreads(14)
mDataTypes <- c('rna', 'cnv', 'mutation')

# ---- 1. Parse full gene signatures
for (i in seq_along(mDataTypes)) {
    message(paste0("Merging signatures for ", mDataTypes[i]))
    dt <- processGeneSignatureFiles(inputDir, mDataTypes[i])
    dt <- unique(dt)
    if (i == 1) {
        cols <- colnames(dt)
        fwrite(dt, file=file.path(outputDir, 
            'gene_compound_tissue_dataset.csv'))
    } else {
        setcolorder(dt, cols)
        fwrite(dt, 
            file=file.path(outputDir, 'gene_compound_tissue_dataset.csv'), 
            append=TRUE)
    }
    rm(dt); gc()
}


# ---- 2. Parse pan-cancer gene signatures
for (i in seq_along(mDataTypes)) {
    message(paste0("Merging pancancer signatures for ", mDataTypes[i]))
    dt <- processPanCancerGeneSignatureFiles(inputDir, mDataTypes[i])
    dt <- unique(dt)
    if (i == 1) {
        cols <- colnames(dt)
        fwrite(dt, file=file.path(outputDir, 'gene_compound_dataset.csv'))
    } else {
        setcolorder(dt, cols)
        fwrite(dt, file=file.path(outputDir, 'gene_compound_dataset.csv'), 
            append=TRUE)
    }
    rm(dt); gc()
}
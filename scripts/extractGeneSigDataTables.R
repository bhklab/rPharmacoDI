library(rPharmacoDI)
library(qs)
library(arrow)
library(BiocParallel)
library(doParallel)
library(data.table)

# ---- 0. Script configuration

inputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/signature_files'
outputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/gene_signatures'

BP <- bpparam()
bpprogressbar(BP) <- TRUE
register(BP)

setDTthreads(14)

# Move the files with rsync, use this when we get new signatures
#system2('bash', 'scripts/move_signatures.sh')

# ---- 1. Parse full signatures
mDataTypes <- c('rna', 'cnv', 'mutation')
for (i in seq_along(mDataTypes)) {
    dt <- processGeneSignatureFiles(inputDir, mDataTypes[i])
    if (i == 1) {
        cols <- colnames(dt)
        fwrite(dt, file=file.path(outputDir, 
            'gene_compound_tissue_dataset.csv'))
    } else {
        fwrite(dt[cols], 
            file=file.path(outputDir, 'gene_compound_tissue_dataset.csv'), 
            append=TRUE)
    }
    rm(dt); gc()
}

# ---- 2. Pancancer Results

for (i in seq_along(mDataTypes)) {
    dt <- processPanCancerGeneSignatureFiles(inputDir, mDataTypes[i])
    if (i == 1) {
        fwrite(dt, file='gene_compound_tissue_dataset.csv')
    } else {
        fwrite(dt, file='gene_compound_tissue_dataset.csv', append=TRUE)
    }
    rm(dt); gc()
}

# Remove signature files to save space
#unlink(inputDir, recursive=TRUE)
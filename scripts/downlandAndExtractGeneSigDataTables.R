library(rPharmacoDI)
library(BiocParallel)
library(data.table)

inputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/signature_files'
outputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/gene_signatures'

BP <- bpparam()
bpprogressbar(BP) <- TRUE
register(BP)
readSigToDT <- function(file) convertGeneSignatureToDT(readRDS(file))
setDTthreads(14)

# Move the files with rsync
# system2('bash scripts/move_signatures.sh')

rna_analytic_files <- list.files(file.path(inputDir, 'rna', 'rnaAnalytic'), 
    pattern='.*rds', recursive=TRUE, full.names=TRUE)
rna_analytic_dt <- rbindlist(bplapply(rna_analytic_files, readSigToDT))

rna_perm_files <- list.files(file.path(inputDir, 'rna', 'rnaPermutation'), 
    pattern='.*rds', recursive=TRUE, full.names=TRUE)
rna_perm_dt <- rbindlist(bplapply(rna_perm_files, readSigToDT))

rna_pancan_files <- list.files(file.path(inputDir, 'rna', 'rnaAnalyticPancan'),
    pattern='.*rds', recursive=TRUE, full.names=TRUE)
rna_pancan_dt <- rbindlist(bplapply(rna_pancan_files, readSigToDT))


cnv_analytic_files <- list.files(file.path(inputDir, 'cnv', 'cnvAnalytic'), 
    pattern='.*rds', recursive=TRUE, full.names=TRUE)
cnv_analytic_dt <- rbindlist(bplapply(cnv_analytic_files, readSigToDT))
fwrite(cnv_analytic_dt, file=file.path(outputDir, 'cnv_analytic.csv.gz'))

cnv_perm_files <- list.files(file.path(inputDir, 'cnv', 'cnvPermutation', 
    'pearson_perm_res'), pattern='.*rds', recursive=TRUE, full.names=TRUE)
cnv_perm_dt <- rbindlist(bplapply(cnv_perm_files, readSigToDT))

cnv_pancan_files <- list.files(file.path(inputDir, 'cnv', 'cnvAnalyticPancan'), 
    pattern='.*rds', recursive=TRUE, full.names=TRUE)
cnv_pancan_dt <- rbindlist(bplapply(cnv_pancan_files, readSigToDT))

## NOTE:: Don't parallelize this, it's done in the functions already
for (pSet in pSetNames) {
    # 
    message("Extracting sigs for ", pSet)

    # make the output directory
    saveDir <- file.path(outputDir, pSet)
    if (!dir.exists(saveDir)) dir.create(saveDir)

    # read in the data
    geneSigL <- readGeneSigsForPSet(dataDir=dataDir, pSetPattern=pSet, 
        mDataTypes=mDataTypes)

    # convert, merge and write to disk
    mergePSetGeneSigsToDT(geneSigL=geneSigL, saveDir=saveDir, 
        fileName=paste0(pSet, '_gene_sig.csv'))
}
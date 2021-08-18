library(rPharmacoDI)
library(qs)
library(arrow)
library(BiocParallel)
library(doParallel)
library(data.table)

# ---- 0. Script configuration

inputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/signature_files'
outputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/gene_signatures'

#registerDoParallel(makeCluster(16))
#BP <- DoparParam()
BP <- bpparam()
bpprogressbar(BP) <- TRUE
register(BP)

setDTthreads(14)

debug(processGeneSignatureFiles)

# Move the files with rsync, use this when we get new signatures
#system2('bash', 'scripts/move_signatures.sh')

mDataTypes <- c('rna', 'cnv', 'mutation')
for (i in seq_along(mDataTypes)) {
    dt <- processGeneSignatureFiles(inputDir, mDataType)
    if (i == 1) {
        fwrite(dt, file=file.path(outputDir, 
            'gene_compound_tissue_dataset.csv'))
    } else {
        fwrite(dt, file=file.path(outputDir, 
            'gene_compound_tissue_dataset.csv'), append=TRUE)
    }
    rm(dt); gc()
}

rna_dt <- processGeneSignatureFiles(inputDir, 'rna')
cnv_dt <- processGeneSignatureFiles(inputDir, 'cnv')
mutation_dt <- processGeneSignatureFiles(inputDir, 'mutation')

# ---- 3. Pancancer Results

## FIXME:: Sanity checks for the PanCan files?

# -- 3.1 RNA pancancer
rna_pancan_files <- list.files(file.path(inputDir, 'rna', 'rnaAnalyticPancan'),
    pattern='.*rds', recursive=TRUE, full.names=TRUE)
rna_pancan_dt <- rbindlist(bplapply(rna_pancan_files, readSigToDT))

# -- 3.2 CNV pancancer
cnv_pancan_files <- list.files(file.path(inputDir, 'cnv', 'cnvAnalyticPancan'), 
    pattern='.*rds', recursive=TRUE, full.names=TRUE)
cnv_pancan_dt <- rbindlist(bplapply(cnv_pancan_files, readSigToDT))

gene_compound_dataset <- rbindlist(list(rna_pancan_dt, cnv_pancan_dt),
    use.names=TRUE)
# Fix dataset names
gene_compound_dataset[dataset == 'CCLE.CTRPv2', dataset := 'CTRPv2']
arrow::write_parquet(gene_compound_dataset, sink=file.path(outputDir, 
    'metaanalysis', 'gene_compound_dataset.parquet'))

# Remove signature files to save space
#unlink(inputDir, recursive=TRUE)
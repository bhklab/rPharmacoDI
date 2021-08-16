library(rPharmacoDI)
library(BiocParallel)
library(data.table)
library(qs)
library(arrow)

# ---- 0. Script configuration

inputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/signature_files'
outputDir <- '../PharmacoDI_snakemake_pipeline/rawdata/gene_signatures'

BP <- bpparam()
bpprogressbar(BP) <- TRUE
register(BP)
readSigToDT <- function(file) convertGeneSignatureToDT(readGeneSig(file))
setDTthreads(14)

# Move the files with rsync, use this when we get new signatures
#system2('bash', 'scripts/move_signatures.sh')

# ---- 1. Merge RNA analytic and permutations results into a single table

# -- 1.1 Analytic results

rna_analytic_files <- list.files(file.path(inputDir, 'rna', 'rnaAnalytic'), 
    pattern='.*rds$|.*qs$', recursive=TRUE, full.names=TRUE)
rna_analytic_dt <- rbindlist(bplapply(rna_analytic_files, readSigToDT))

# Extract tissues to fix the permutation table
rna_analytic_tissues <- unique(rna_analytic_dt$tissue)

# -- 1.2 Permutation results

rna_perm_files <- list.files(file.path(inputDir, 'rna', 'rnaPermutation'), 
    pattern='.*rds$|.*qs$', recursive=TRUE, full.names=TRUE)
rna_perm_list <- bplapply(rna_perm_files, readSigToDT)
rna_perm_tissue <- gsub('^.*_|.rds$|.qs$', '', basename(rna_perm_files))
# reconstruct tissues from file names
for (i in seq_along(rna_perm_list)) {
    rna_perm_list[[i]][, tissue := rna_perm_tissue[i]]
}
rna_perm_dt <- rbindlist(rna_perm_list, use.names=TRUE, fill=TRUE)

# -- 1.3 Merge formatting

# build an fcase statement to replace tissue names
match_exprs <- paste0('tissue %like% "', make.names(rna_analytic_tissues), '"')
fcase_args <- unlist(
    Map(list, sapply(match_exprs, FUN=str2lang), rna_analytic_tissues), 
    recursive=FALSE)

# evaluate our fcase statment
## TODO: Is this slower than just using a for loop?
rna_perm_dt[, tissue := do.call('fcase', args=fcase_args)]
if (all(rna_perm_dt$tissue %in% rna_analytic_tissues)) {
    rm(rna_perm_list); gc(verbose=TRUE)
}

# 1.4 -- Join the tables, adding the appropriate suffixes
key_columns <- c('gene', 'compound', 'tissue', 'dataset', 'mDataType')
setkeyv(rna_analytic_dt, key_columns)
setkeyv(rna_perm_dt, key_columns)
rna_dt <- merge.data.table(rna_analytic_dt, rna_perm_dt, 
    suffixes=c('_analytic', '_permutation'), all=TRUE)

# Sanity check sample sizes
rna_different_n <- rna_dt[
    !is.na(n_analytic) & !is.na(n_permutation) & n_analytic != n_permutation, 
    .(gene, compound, tissue, dataset, mDataType, n_analytic, n_permutation)
    ]
if (nrow(rna_different_n) > 1) {
    warning('Sample sizes differ between analytic and permutation methods for RNA!')
    fwrite(rna_different_n, 'logs/rna_gene_sig_different_n.csv')
}

rna_dt[, c('n', 'df', 'estimate') := .(n_permutation, df_permutation, estimate_permutation)]
rna_dt[is.na(estimate), c('n', 'df', 'estimate') := .(n_analytic, df_analytic, estimate_analytic)]

# Fix CTRPv2 since it is a hybrid dataset
rna_dt[dataset == 'CCLE.CTRPv2', dataset := 'CTRPv2']

# Drop rows with NA estimates
rna_dt <- rna_dt[!is.na(estimate), ]
rna_dt[, 
    c('estimate_analytic', 'n_analytic', 'df_analytic', 
        'significant_analytic') := NULL]
rna_dt[, c('estimate_permutation', 'n_permutation', 'df_permutation') := NULL]


# Clean up so we have memory for the next table
rm(rna_perm_dt, rna_analytic_dt); gc(verbose=TRUE)


# ---- 2. Merge CNV analytic and permutation results

# -- 2.1 Analytic results

cnv_analytic_files <- list.files(file.path(inputDir, 'cnv', 'cnvAnalytic'), 
    pattern='.*rds|.*qs', recursive=TRUE, full.names=TRUE)
cnv_analytic_dt <- rbindlist(bplapply(cnv_analytic_files, readSigToDT))

# Extract tissues to fix the permutation table
cnv_analytic_tissues <- unique(cnv_analytic_dt$tissue)

# ---- 2.2 Permutation results
cnv_perm_files <- list.files(file.path(inputDir, 'cnv', 'cnvPermutation', 
    'pearson_perm_res'), pattern='.*rds$|.qs$', recursive=TRUE, full.names=TRUE)
cnv_perm_tissue <- gsub('^.*_|.rds$|.qs$', '', basename(cnv_perm_files))
cnv_perm_list <- bplapply(cnv_perm_files, readSigToDT)
# reconstruct tissues from file names
for (i in seq_along(cnv_perm_list)) {
    cnv_perm_list[[i]][, tissue := cnv_perm_tissue[i]]
}
cnv_perm_dt <- rbindlist(cnv_perm_list, use.names=TRUE, fill=TRUE)

# -- 2.3 Merge formatting

# build an fcase statement to replace tissue names
cnv_match_exprs <- paste0('tissue %like% "', make.names(cnv_analytic_tissues), '"')
cnv_fcase_args <- unlist(
    Map(list, sapply(cnv_match_exprs, FUN=str2lang), cnv_analytic_tissues), 
    recursive=FALSE)
cnv_perm_dt[, tissue := do.call('fcase', args=cnv_fcase_args)]

# -- 2.4 Join the tables, adding the appropriate suffixes
key_columns <- c('gene', 'compound', 'tissue', 'dataset', 'mDataType')
setkeyv(cnv_analytic_dt, key_columns)
setkeyv(cnv_perm_dt, key_columns)
cnv_dt <- merge.data.table(cnv_analytic_dt, cnv_perm_dt, 
    suffixes=c('_analytic', '_permutation'), all=TRUE)

# Sanity check sample sizes
cnv_different_n <- cnv_dt[
    !is.na(n_analytic) & !is.na(n_permutation) & n_analytic != n_permutation, 
    .(gene, compound, tissue, dataset, mDataType, n_analytic, n_permutation)
    ]
if (nrow(cnv_different_n) > 1) {
    warning('Sample sizes differ between analytic and permutation methods for CNV!')
    fwrite(cnv_different_n, 'logs/cnv_gene_sig_different_n.csv')
}

cnv_dt[, c('n', 'df', 'estimate') := .(n_permutation, df_permutation, estimate_permutation)]
cnv_dt[is.na(estimate), c('n', 'df', 'estimate') := .(n_analytic, df_analytic, estimate_analytic)]

# Fix CTRPv2 since it is a hybrid dataset
cnv_dt[dataset == 'CCLE.CTRPv2', dataset := 'CTRPv2']

# Drop rows with NA estimates
cnv_dt <- cn_dt[!is.na(estimate), ]
cnv_dt[, 
    c('estimate_analytic', 'n_analytic', 'df_analytic', 
        'significant_analytic') := NULL]
cnv_dt[, c('estimate_permutation', 'n_permutation', 'df_permutation') := NULL]

# Clean up memory usage
rm(cnv_perm_dt, cnv_analytic_dt); gc(verbose=TRUE)

# ---- 3. Merge RNA and CNV results and split by dataset 
#> (for pipeline compatibility)
gene_compound_tissue_dataset <- rbindlist(list(rna_dt, cnv_dt))
rm(rna_dt, cnv_dt); gc(verbose=TRUE)

for (ds in unique(gene_compound_tissue_dataset$dataset)) {
    print(ds)
    arrow::write_parquet(gene_compound_tissue_dataset[dataset == ds, ], 
        sink=file.path(outputDir, paste0(ds, '_gene_sig.parquet')))
}

rm(gene_compound_tissue_dataset); gc(v=T)

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
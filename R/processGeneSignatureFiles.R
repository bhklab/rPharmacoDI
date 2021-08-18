#' Read gene signature files for many experiments and merge them to a 
#'   `data.table`
#'
#' @param filePath `character(1)` Path to top level signature directory. This
#'   should contain a folder `<mDataType>`.
#' @param mDataType `character(1)` The molecular data type of the signatures.
#' @param ... Fall through parameters to [`BiocParallel::bplapply`]. These
#'   parameters contral internal function parallelization.
#' @param keyColumns `character` Vector of column names in the signature file
#'   to join on. For developer use, do no change.
#' @param fileFormats `character` Vector of file formats without a dot 
#'   (e.g., rds)
#' 
#' @return `data.table` The gene signatures properly formatted to a single 
#'   long format table.
#' 
#' @importFrom data.table data.table rbindlist merge.data.table `:=`
#' @importFrom BiocParallel bplapply bpparam
#' @export
processGeneSignatureFiles <- function(
    filePath,
    mDataType=c('rna', 'cnv', 'mutation'), 
    ...,
    keyColumns=c('gene', 'compound', 'tissue', 'dataset', 'mDataType'),
    fileFormats=c('rds', 'qs'),
    logDir='logs')
{
    mDataType <- match.arg(mDataType)
    format_regex_pattern <- paste0('^.*_|', paste0('.', fileFormats, '$', 
        collapse='|'))
    
    ## -- Analytic signatures
    message('Reading analytic gene signatures')
    analytic_signature_files <- list.files(
        file.path(filePath, mDataType, paste0(mDataType, 'Analytic')),
        pattern=format_regex_pattern,
        recursive=TRUE,
        full.names=TRUE)
    analytic_dt <- rbindlist(
        bplapply(analytic_signature_files, readSigToDT, ...),
        use.names=TRUE)
    print(colnames(analytic_dt))

    ## FIXME:: (A) Remove this workaround for fixing permutation tissue
    # extract tissues to fix permutation table
    analytic_tissues <- unique(analytic_dt$tissue)

    ## -- Permutation signatures
    message('Reading permutation gene signatures')
    permutation_signature_files <- list.files(
        file.path(filePath, mDataType, paste0(mDataType, 'Permutation')),
        pattern=format_regex_pattern,
        recursive=TRUE,
        full.names=TRUE)
    analytic_only <- FALSE
    if (length(permutation_signature_files) < 1) analytic_only <- TRUE
    if (!analytic_only) {
        permutation_list <- bplapply(permutation_signature_files, FUN=readSigToDT,
            ...)
        permutation_tissue <- gsub(format_regex_pattern, '',
            basename(permutation_signature_files))
        ## FIXME:: Same as (A)
        # reconstruct tissue from file names
        for (i in seq_along(permutation_list)) {
            permutation_list[[i]][, tissue := permutation_tissue[i]]
        }
        permutation_dt <- rbindlist(permutation_list, use.names=TRUE)
        rm(permutation_list); gc()

        ## -- Merge formatting
        message('Fixing permutation signature tissue column formatting')
        ## FIXME:: Same as (A)
        tissue_match_exprs <- paste0('tissue %like% "', 
            make.names(analytic_tissues), '"')
        fcase_args <- unlist(
            Map(list, sapply(tissue_match_exprs, FUN=str2lang), analytic_tissues),
            recursive=FALSE)

        # evaluation tissue matching expressions
        permutation_dt[, tissue := do.call('fcase', args=fcase_args)]
        if (!all(unique(permutation_dt$tissue) %in% unique(analytic_dt$tissue))) {
            stop('Tissue mismatches between analytic and permutation signatures!')
        }

        ## -- Join the tables and add the appropriate suffixes
        message('Joining analytic and permutation results')
        setkeyv(permutation_dt, keyColumns)
        setkeyv(analytic_dt, keyColumns)
        signature_dt <- merge.data.table(analytic_dt, permutation_dt,
            suffixes=c('_analytic', '_permutation'), all=TRUE)
        rm(permutation_dt, analytic_dt); gc()
    } else {
        # Pad with empty permutation columns if no permutation signatures are
        #>available
        signature_dt <- analytic_dt
        rm(analytic_dt); gc()
        non_key_columns <- setdiff(colnames(signature_dt), keyColumns)
        setnames(
            signature_dt, 
            old=non_key_columns, 
            new=paste0(non_key_columns, '_analytic'))
        permutation_columns <- paste0(non_key_columns, '_permutation')
        for (column in permutation_columns)
            set(signature_dt, j=column, value=NA_real_)
    }

    # Sanity check samples sizes
    message('Performing sanity checks')
    different_sample_size <- signature_dt[
        !is.na(n_analytic) & !is.na(n_permutation) & 
            n_analytic != n_permutation,
        .(gene, compound, tissue, dataset, mDataType, n_analytic, n_permutation)
    ]
    if (nrow(different_sample_size) > 1) {
        warning('Samples sizes different between analytic and permutation 
            results!')
        if (!dir.exists(logDir)) dir.create(logDir, recursive=TRUE)
        fwrite(different_sample_size, 
            file=file.path(logDir, paste0(mDataType, '_different_n.csv')))
    }

    message('Merging signature results')
    signature_dt[, 
        c('n', 'df', 'estimate') := .(n_permutation, df_permutation, 
            estimate_permutation)
        ]
    signature_dt[
        is.na(estimate), 
        c('n', 'df', 'estimate') := .(n_analytic, df_analytic, 
            estimate_analytic)
        ]

    # Check for NA estimates
    if (any(is.na(signature_dt$estimate))) {
        warning('NA estimates still in for both analytic and permuation!
            Dropping these rows!')
        na_estimate <- signature_dt[is.na(estimate), ]
        fwrite(na_estimate, 
            file=file.path(logDir, paste0(mDataType, '_na_estimate.csv')))
        rm(na_estimate); gc()
        signature_dt <- signature_dt[!is.na(estimate), ]
    }

    # Delete unwanted columns
    signature_dt[,
    c('estimate_analytic', 'n_analytic', 'df_analytic',
        'significant_analytic') := NULL
        ]
    signature_dt[, 
        c('estimate_permutation', 'n_permutation',
            'df_permutation') := NULL
        ]

    # Remove ensembl gene ids
    message('Removing Ensembl gene id versions')
    signature_dt[gene %like% '^ENS', gene := gsub('\\..*$', '', gene)]

    # Fix CTRPv2 naming
    signature_dt[dataset == 'CCLE.CTRPv2', dataset := 'CTRPv2']

    message('DONE!')
    return(signature_dt)
}

#' Read a gene signature file, then convert it to a `data.table`
#'
#' @param filePath `character(1)` The path to the signature file.
#'
#' @return `data.table` The contents of the gene signature long format.
#'
#' @export
readSigToDT <- function(filePath) convertGeneSignatureToDT(readGeneSig(filePath))

#' Read pan-cancer gene signature files and merge them into a `data.table`
#'
#' @param filePath `character(1)` Path to top level signature directory. This
#'   should contain a folder `mDataType`.
#' @param mDataType `character(1)` The molecular data type of the signatures.
#'
#' @return `data.table` The gene signatures properly formatted to a single 
#'   long format table.
#' 
#'
#' @importFrom data.table data.table rbindlist
#' @export
processPanCancerGeneSignatureFiles <- function(filePath,
    mDataType=c('rna', 'cnv', 'mutation'), fileFormats=c('rds')) 
{
    mDataType <- match.arg(mDataType)
    pancan_files <- list.files(
        file.path(filePath, mDataType, paste0(mDataType, 'AnalyticPancan')),
        pattern=paste0('.', fileFormats, '$', collapse='|'),
        recursive=TRUE,
        full.names=TRUE)
    pancan_dt <- rbindlist(pancan_files, use.names=TRUE, fill=TRUE)
    pancan_dt[gene %like% '^ENS', gene := gsub('\\..*$', '', gene)]
    pancan_dt[dataset == 'CCLE.CTRPv2', dataset := 'CTRPv2']
    return(pancan_dt)
}
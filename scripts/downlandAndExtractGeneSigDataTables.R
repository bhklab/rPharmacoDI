library(rPharmacoDI)
library(reticulate)
use_condaenv('/home/pharmdi/miniconda3/envs/data_ingestion')
di <- import('PharmacoDI')

pSetNames <- c('gCSI', 'GDSC_v1', 'GDSC_v2', 'CTRPv2', 'CCLE', 'UHNBreast', 'GRAY')
dataDir = '../PharmacoDI_snakemake_pipeline/rawdata'

# -- Download gene signatures (requires Chris' ComputeCanada password)
di$download_gene_signatures(user=Sys.getenv('UU'), password="", save_dir=dataDir)

## NOTE:: Don't parallelize this, it's done in the functions already
for (pSet in pSetNames) {

    # make the output directory
    saveDir <- file.path(dataDir, pSet)
    if (!dir.exists(saveDir)) dir.create(saveDir)

    # read in the data
    geneSigL <- readGeneSigsForPSet(dataDir=dataDir, pSetPattern=pSet)

    # convert, merge and write to disk
    mergePSetGeneSigsToDT(geneSigL=geneSigL, saveDir=saveDir, fileName=paste0(pSet, '_gene_sig.csv'))
}
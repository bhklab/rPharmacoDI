#!/bin/bash

# Move the signature files to the pipeline directory
PETR="$HOME/bhklab/psmirnov/signature_files"
TARGET="$HOME/DataIngestion/PharmacoDI_snakemake_pipeline/rawdata/signature_files"

rsync -htrP "$PETR"/cnv* "$TARGET"/cnv
rsync -htrP "$PETR"/rna* "$TARGET"/rna

#!/usr/bin/env Rscript

#######
# LOG #
#######

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, append = TRUE, type = "output")

#############
# LIBRARIES #
#############

library(data.table)
library(tidyverse)

###########
# GLOBALS #
###########

pca_covars_file <- snakemake@input[["pca_covars_file"]]
expt_covars_file <- snakemake@input[["expt_covars_file"]]

########
# MAIN #
########

pca_covars <- fread(pca_covars_file)
expt_covars <- fread(expt_covars_file)

merged_covars <- merge(pca_covars, expt_covars, by="#IID")

write_tsv(merged_covars, snakemake@output[["merged_covars"]])

# write log
sessionInfo()

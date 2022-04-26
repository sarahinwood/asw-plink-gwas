library(data.table)
library(tidyverse)
library(qqman)

## reading in all results together
adjusted_glm_res_files <- list.files("output/plink2/glm_analysis", pattern=".adjusted", recursive=T, full.names = T)

glm_results_table <- read_delim(adjusted_glm_res_files, id="path")

# sorting bonferoni adjusted sig hits
bonf_sig_glm <- data.table(subset(glm_results_table, glm_results_table$BONF<0.05))
bonf_sig_glm$covar_analysis <- ifelse(grepl("/no_covars/", bonf_sig_glm$path), "No covars", "Covars")
bonf_sig_glm$phenotype <- ifelse(grepl(".location.", bonf_sig_glm$path), "Location",
                                ifelse(grepl(".attack_status.", bonf_sig_glm$path), "Attack status", "Parasitism"))
bonf_sig_glm$input_file <- tstrsplit(bonf_sig_glm$path, "/", fixed=T, keep=c(5))


pruned_location_sig <- subset(bonf_sig_glm, bonf_sig_glm$input_file=="all_samples_pruned")
location_sig <- subset(bonf_sig_glm, bonf_sig_glm$input_file=="all_samples")
length(intersect(pruned_location_sig$ID, location_sig$ID))
# 24 variants found in both


# sorting FDR adjusted sig hits
fdr_sig_glm <- data.table(subset(glm_results_table, glm_results_table$FDR_BH<0.05))
fdr_sig_glm$covar_analysis <- ifelse(grepl("/no_covars/", fdr_sig_glm$path), "No covars", "Covars")
fdr_sig_glm$phenotype <- ifelse(grepl(".location.", fdr_sig_glm$path), "Location",
                                ifelse(grepl(".attack_status.", fdr_sig_glm$path), "Attack status", "Parasitism"))
fdr_sig_glm$input_file <- tstrsplit(fdr_sig_glm$path, "/", fixed=T, keep=c(5))


# 47 found in both - so adjustment method affects what results we detect - but is this a common way to pick assoc. variants in gwas?
length(intersect(fdr_sig_glm$ID, bonf_sig_glm$ID))

## plotting for each analysis

#for i in adjusted_glm_res_files




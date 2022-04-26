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
library(qqman)

###########
# GLOBALS #
###########

results_file <- snakemake@input[["results"]]

########
# MAIN #
########

results <- fread(results_file)

## sig results ##

sig_assoc <- sum(results$BONF < 0.05)
paste(sig_assoc, "significant assoc.s after Bonferroni correction", sep=" ")

bonferroni_threshold <- 0.05/(length(results$ID))

## plotting ##

# plotting qq
pdf(snakemake@output[["qq_plot"]])
qq(results$UNADJ)
dev.off()

# plotting manhattan
plot_threshold <- -log10(bonferroni_threshold)
# chr label
results$chr_plot <- tstrsplit(results$ID, "ASW_TRINITY_DN", keep=2)
results$chr_plot <- tstrsplit(results$chr_plot, ":", keep=1)
results$chr_plot <- gsub("_c", "", results$chr_plot)
results$chr_plot <- gsub("_g", "", results$chr_plot)
results$chr_plot <- as.numeric(results$chr_plot)
#BP
results$BP <- tstrsplit(results$ID, ":", keep=2)
results$BP <- as.numeric(results$BP)
# plot
pdf(snakemake@output[["manhattan_plot"]])
manhattan(results, chr="chr_plot", p="UNADJ", genomewideline = plot_threshold, suggestiveline = F)
dev.off()

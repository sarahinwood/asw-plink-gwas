library(data.table)
library(qqman)

nocov_res <- fread("test_output/no_covars/all_samples_pruned.assoc.logistic")
nocov_res_adj <- fread("test_output/no_covars/all_samples_pruned.assoc.logistic.adjusted")



nocov_res$bonf_adj <- p.adjust(nocov_res$P, "bonferroni")

sum(cov2_res_adj$BONF<0.05)



qq(cov2_res$P)

qq(cov2_res_adj$UNADJ)


# plotting manhattan
plot_threshold <- -log10(bonferroni_threshold)
# chr label
cov2_res$chr_plot <- tstrsplit(cov2_res$SNP, "ASW_TRINITY_DN", keep=2)
cov2_res$chr_plot <- tstrsplit(cov2_res$chr_plot, ":", keep=1)
cov2_res$chr_plot <- gsub("_c", "", cov2_res$chr_plot)
cov2_res$chr_plot <- gsub("_g", "", cov2_res$chr_plot)
cov2_res$chr_plot <- as.numeric(cov2_res$chr_plot)
#BP
cov2_res$BP <- tstrsplit(cov2_res$SNP, ":", keep=2)
cov2_res$BP <- as.numeric(cov2_res$BP)
# plot
manhattan(cov2_res, chr="chr_plot", p="UNADJ", genomewideline = plot_threshold, suggestiveline = F)


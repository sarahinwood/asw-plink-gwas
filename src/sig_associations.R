library(data.table)
library(tidyverse)
library(VennDiagram)

assoc_tables <- list.files("sig_assoc", pattern=".adjusted", recursive=T, full.names=T)

full <- read_csv(assoc_tables, id="path")
sig_assocs <- subset(full, p_Bonferroni<0.05)
#detects same parasitism assoc

##location
location_pruned <- subset(sig_assocs, path=="sig_assoc/location_all_samples_pruned.assoc.adjusted")
location_unpruned <- subset(sig_assocs, path=="sig_assoc/location_all_samples.assoc.adjusted")

intersect(location_pruned$SNP, location_unpruned$SNP)
##22 overlap - 6 unique to pruned, 23 to unpruned

#more location assoc.s than tassel
tassel_pruned_loc <- fread("../asw-gwas-rnaseq/output/pruned_location/pruned_location-only_bonferroni_annots.csv")
tassel_loc <- fread("../asw-gwas-rnaseq/output/location/location-only_bonferroni_annots.csv")


vd <- venn.diagram(x=list("plink_loc_pruned"=location_pruned$SNP, "plink_loc"=location_unpruned$SNP, "tassel_loc_pruned"=tassel_pruned_loc$Marker), filename=NULL)
grid.newpage()
grid.draw(vd)

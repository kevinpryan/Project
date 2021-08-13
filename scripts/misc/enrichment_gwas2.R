#!/usr/bin/Rscript

# This script is used to: determine whether the DBD variants overlap GWAS variants more often than is expected by chance. Also writes two tables - overlaps_genes_gwascatalog.txt which details the overlaps of the genes with the GWAS catalog and genes_any_overlap_gwas.txt which tells us whether each gene had any overlap with GWAS variants. 
# requirements: frac_rank file (all that is really required is a table with one column being variants in GTEx format and the other being the gene symboli)
# Also requires the GWAS catalog and the GTEx lookup table (https://www.ebi.ac.uk/gwas/docs/file-downloads and https://gtexportal.org/home/datasets)
library(vroom)
library(GenomicRanges, quietly = T)
library(stringr)
library(tidyr)
library(ggplot2)
library(dplyr)
# read in files
gwas <- vroom("/data/kryan/project/gtex/analysis/gwas_catalog_v1.0.2-associations_e104_r2021-07-08.tsv")
lookup <- vroom("/data/kryan/project/gtex/GTEx_Analysis_2017-06-05_v8_WholeGenomeSeq_838Indiv_Analysis_Freeze.lookup_table.txt.gz")
frac_rank_enrichment_data <- read.table("/data/kryan/project/gtex/outputs/trans_eqtl_results/frac_rank/all_targets/all_covariates/no_filter_tf/enrichment_frac_rank_covariates_AF1_qc_corrected.txt", header = T) 
# only want unique genes
frac_rank_enrichment_data <- frac_rank_enrichment_data %>% distinct(Variant_gene, .keep_all= TRUE)
# create granges object from variants of interest
var_split_eqtl <- str_split_fixed(frac_rank_enrichment_data$Variant, pattern = "_", n = 5)
chrom <- substr(var_split_eqtl[,1], start = 4, stop = nchar(var_split_eqtl[,1]))
print(chrom)
start <- as.numeric(var_split_eqtl[,2]) - 10000
print(start)
end <- as.numeric(var_split_eqtl[,2]) + 10000
print(end)
df_vars <- data.frame(chrom = chrom, start = start, end = end)
print(df_vars)
gr_vars <- as(df_vars, "GRanges")

# create granges object from significant gwas catalog variants
gwas_rm_na <- gwas[!is.na(gwas$CHR_ID) & !is.na(gwas$CHR_POS),]
gwas_rm_na_sig <- gwas_rm_na[which(gwas_rm_na$`P-VALUE` < 5e-8),]
df_gwas <- data.frame(chrom = gwas_rm_na_sig$CHR_ID, start = gwas_rm_na_sig$CHR_POS, end = gwas_rm_na_sig$CHR_POS)
gr_gwas <- as(df_gwas, "GRanges")
gr_gwas$"MAPPED_TRAIT" <- gwas_rm_na_sig$"MAPPED_TRAIT"
gr_gwas$"PUBMEDID" <- gwas_rm_na_sig$"PUBMEDID"
print(gr_gwas[1:5])
# get rid of na chromosomes or positions from gtex lookup
lookup_drop_na <- lookup %>% drop_na(chr, variant_pos)

overlaps_vars_gwas <- length(which(overlapsAny(gr_vars, gr_gwas) == TRUE))
print("overlaps_vars_gwas")
print(overlaps_vars_gwas)
print("overlaps any")
print(overlapsAny(gr_vars, gr_gwas))
var_gwas_overlaps_proportion <- overlaps_vars_gwas/length(gr_vars)
print("var_gwas_overlaps_proportion")
print(var_gwas_overlaps_proportion)
df.out <- data.frame()
for (k in 1:length(gr_vars)){
	print("gene")
	print(frac_rank_enrichment_data$Variant_gene[k])
	print("GWAS catalog entry")
	gwas_catalog_overlap <- gr_gwas[!is.na(findOverlaps(gr_gwas, gr_vars[k], select="arbitrary"))]
	print(gr_gwas[!is.na(findOverlaps(gr_gwas, gr_vars[k], select="arbitrary"))])
	gene_fordf <- rep(frac_rank_enrichment_data$Variant_gene[k], length(gwas_catalog_overlap))
	df.gene <- cbind.data.frame(gene_fordf, gwas_catalog_overlap$"MAPPED_TRAIT", gwas_catalog_overlap$PUBMEDID)
	df.out <- rbind.data.frame(df.out, df.gene)
}
colnames(df.out) <- c("Variant_gene", "MAPPED_TRAIT", "PUBMEDID")
write.table(df.out, file = "overlaps_genes_gwascatalog.txt", sep = "\t", quote = F, row.names = F)
freq_overlap_background <- c()
for (j in 1:300){
    print(paste("Iteration:", j))
    x <- sample(x = rownames(lookup_drop_na), size = 5000, replace = TRUE)
    lookup_sampled <- lookup_drop_na[x,]
    lookup_chrom_no <- substr(x = lookup_sampled$chr, start = 4, stop = nchar(lookup_sampled$chr))
    df_lookup <- data.frame(chrom = lookup_chrom_no, start = as.numeric(lookup_sampled$variant_pos)-10000, end =    as.numeric(lookup_sampled$variant_pos)+10000)
    gr_lookup <- as(df_lookup, "GRanges") 
    overlaps_proportion <- length(which(overlapsAny(gr_lookup, gr_gwas) == TRUE))/5000
    freq_overlap_background <- c(freq_overlap_background, overlaps_proportion)
}
print(freq_overlap_background)
background_frequency <- mean(freq_overlap_background)
print("background frequency")
print(background_frequency)
enrichment_fold <- var_gwas_overlaps_proportion/background_frequency
print(binom.test(x =overlaps_vars_gwas , n = length(gr_vars), p = background_frequency,conf.level = 0.95, alternative = "greater"))
print(paste("fold enrichment: ", enrichment_fold))
genes_overlapsany <- cbind.data.frame(frac_rank_enrichment_data$Variant_gene, overlapsAny(gr_vars, gr_gwas))
colnames(genes_overlapsany) <- c("Variant_gene", "overlaps_any")
write.table(genes_overlapsany, file = "genes_any_overlap_gwas.txt", sep = "\t", quote = F, row.names = F)


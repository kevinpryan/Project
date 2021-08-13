#!/usr/bin/Rscript
# This script produces a boxplot for every significant trans-eQTL seen to result from a DBD mutation
library(vroom)
library(ggplot2)
library(stringr)
library(dplyr, warn.conflicts = FALSE)
library(GenomicRanges, quietly = T)
library(tibble)
library(ggpubr)
library(rlist)
library(gridExtra)
library(grid)

# Run script by giving arguments: args[1] = "default", args[2] = "whole_blood_norm" etc, args[3] = "pass", args[4] = variant_interest, args[5] = target_interest
args = commandArgs(trailingOnly = TRUE)
# hardcoding in directory containing results
results <- read.table("/data/kryan/project/gtex/outputs/trans_eqtl_results/all_tissues/all_targets/eqtl/all_covariates/all_tissues_all_targets_eqtl_sig_results.txt", header = T)
variants_interest <- results$Variant_GTEx_ID
print(paste("variants interest:", variants_interest))
targets_interest <- results$Target_gene_name
print(paste("targets_interest:", targets_interest))

if (length(args) == 0) {
  stop("At least one argument must be supplied ('average' or 'default' )", call.=FALSE)
} else if (args[1] != "average" && args[1] != "default") {
  stop("Argument must either be 'average' or 'default'", call. = FALSE)
}

genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 1.0e-02
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"
expression_dir <- args[2]
f <- list.files(path = "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_expression_matrices/", pattern = "*.gz$")
tissues <- str_split_fixed(f,pattern = "\\.", n = 5)[,1]
print("tissues")
print(tissues)
covariates <- paste(tissues, ".v8.covariates.txt", sep = "")
psam <- "/data/kryan/project/gtex/gtex.psam"
dorothea_file <- "/data/kryan/project/gtex/analysis/dorothea_hs_with_target_info.txt"
# Read in genotype table
genotype.data <- read.table(genotypes, stringsAsFactors=FALSE)
print(paste("genotype data nrow before qc:", nrow(genotype.data)))
genotype_data_fail_VQSR <- read.table("/data/kryan/project/gtex/genotype_data_dbdp_fail_qc_29072921.txt")
genotype_data_fail_vars <- genotype_data_fail_VQSR$V3
genotype.data <- genotype.data[which(!(genotype.data$V3 %in% genotype_data_fail_vars)),]
print(paste("genotype data nrow after qc:", nrow(genotype.data)))
# Change genotype encoding from 0/0,0/1,1/1 to 0,1,2
genotype_to_numeric <- function(par){
 par.out <- par
 par.out[which(par == "0/0")] <- 0
 par.out[which(par == "0/1")] <- 1
 par.out[which(par == "1/1")] <- 2
return(as.numeric(par.out))
}

# convert to numeric - get some warnings
variants.genotypes <- genotype.data[,10:ncol(genotype.data)]
variants.numeric <- apply(variants.genotypes, 2, genotype_to_numeric)
variant.ids <- genotype.data$V3

psam <- "/data/kryan/project/gtex/gtex.psam"
dorothea_file <- "/data/kryan/project/gtex/analysis/dorothea_hs_with_target_info.txt"
# Read in genotype table
genotype.data <- read.table(genotypes, stringsAsFactors=FALSE)
# Change genotype encoding from 0/0,0/1,1/1 to 0,1,2
genotype_to_numeric <- function(par){
 par.out <- par
 par.out[which(par == "0/0")] <- 0
 par.out[which(par == "0/1")] <- 1
 par.out[which(par == "1/1")] <- 2
return(as.numeric(par.out))
}

variants.genotypes <- genotype.data[,10:ncol(genotype.data)]
variants.numeric <- apply(variants.genotypes, 2, genotype_to_numeric)
variant.ids <- genotype.data$V3
myplots <- vector('list', length(variants_interest))
x_axis = c(3,2,3,2,3,2) 
y_axis = c(2.3, 2.5,1.5, 1.25, -1.5, 2.3)
for (j in 1:length(variants_interest)){
myplots[[j]] <- local({
j <- j 
tissue <- results$tissue[j]
tissue_number_expression_covariates <- which(tissues %in% tissue)
snp_id <- results$rs_id[j]
print(paste("tissue:", f[tissue_number_expression_covariates]))
expression_file <- paste("/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_expression_matrices/", f[tissue_number_expression_covariates], sep = "")
print(paste("expression file:", expression_file))
cov_file <- paste("/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/", covariates[tissue_number_expression_covariates], sep = "")
print(paste("cov_file:", cov_file))
cov <- vroom(cov_file)
cov_t <- t(cov[-1])
colnames(cov_t) <- cov$ID

# Read in expression data
expression_data <- vroom(expression_file)
dorothea <- read.table(dorothea_file, header = T, sep = "\t", check.names = F)
Gene_full <- expression_data$"gene_id"
Gene_split <- str_split_fixed(Gene_full, pattern = "\\.", n = 2)
Gene <- Gene_split[,1]
expression_data <- expression_data[Gene %in% dorothea$ensembl_gene_id_target,]
expression_data <- expression_data %>%
                   add_column(ensembl_gene_id_target_short = str_split_fixed(expression_data$gene_id, pattern = "\\.", n = 2)[,1], .after = "gene_id")

singletons <- names(which(table(expression_data$ensembl_gene_id_target_short) == 1))
expression_data_notduplicated <- expression_data[expression_data$ensembl_gene_id_target_short %in% singletons, ]

# Read in sample ids and match up with the sample ids in the expression data file
psam.ids <- read.table(psam)
colnames(variants.numeric) <- psam.ids$V1
isect <- intersect(colnames(variants.numeric),colnames(expression_data[6,]))
variants.genotypes.matched <- variants.numeric[,isect]

# Calculate MAFs based on individuals present in data. Exclude NAs from the sum as well as the denominator of the fraction
maf <- function(par){
 par.out <- sum(par, na.rm = TRUE)/(2*length(par[!is.na(par)]))
 par.out <- pmin(par.out, 1-par.out)
 return(par.out)
}
mafs <- apply(variants.genotypes.matched, 1, maf)
print(paste("tissue:", tissue))
print(paste("sample size in this tissue:", length(isect)))
variant.ids.range <- variant.ids[mafs >= allele.frequency.desired]
print(paste("Number of variants present at desired allele frequency:", length(variant.ids.range)))
variants.range.genotypes.matched <- variants.genotypes.matched[which(mafs >= allele.frequency.desired),]
print(paste("number of rows in variants.range.genotypes.matched", nrow(variants.range.genotypes.matched)))
expression.data.samples.matched <- expression_data_notduplicated[,isect]
rownames(expression.data.samples.matched) <- expression_data_notduplicated$ensembl_gene_id_target_short

# Read in table from Barrera paper (with GTEx ids added in manually)
tf.dbdp <- read.table(barrera.variants)

# Only keep rows that are in GTEx at the desired allele frequency
tf.dbdp.variants <- tf.dbdp[tf.dbdp$gtex_var_format_tfdbdp_b38 %in% variant.ids.range,]

# Get rid of duplicated variants (the Barrera paper has separate entries for different transcript ids)
tf.dbdp.uniq <- tf.dbdp.variants[!duplicated(tf.dbdp.variants$gtex_var_format_tfdbdp_b38),]

# Read in dorothea file - created manually with added features, as of 17/06/2021

# Only keep the variants in transcription factors that have data in the dorothea database
tf.dbdp.uniq.indorothea <- tf.dbdp.uniq[tf.dbdp.uniq$Gene.symbol %in% dorothea$tf,]

i <- which(tf.dbdp.uniq.indorothea$gtex_var_format_tfdbdp_b38 == variants_interest[j])

print(paste("Number of transcription factors with GTEx variants in desired frequency range available in dorothea database:",nrow(tf.dbdp.uniq.indorothea)))

# Create GRanges object from dorothea
df <- data.frame(chrom=dorothea$target_chromosome, start=dorothea$target_start, end=dorothea$target_end)
df$chrom <- paste("chr",df$chrom, sep = "")
gr1 <- as(df, "GRanges")

# Dataframes for dorothea upregulated and downregulated targets
dorothea_upreg <- dorothea[dorothea$mor == 1,]
dorothea_downreg <- dorothea[dorothea$mor == -1,]

# Create GRanges object for filtered variants range
vars <- tf.dbdp.uniq.indorothea[,22]
vars.split <- str_split_fixed(vars, "_",5)
vars.chromosome <- vars.split[,1]
var.positions <- as.numeric(vars.split[,2])
range.start <- (var.positions - 1e06)
range.start[range.start < 1] <- 1
range.end <- as.numeric(var.positions + 1e06)

df2 <- data.frame(chrom = vars.chromosome,start=range.start, end=range.end)
gr2 <- as(df2, "GRanges")

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
#df.out <- data.frame(variant = character(), tf.interest = character(), target = character(), p_value = as.numeric())

not_all_na <- function(x) any(!is.na(x))

var <- tf.dbdp.uniq.indorothea[i,22]
print(paste("variant of interest:", var))
tf.interest <- tf.dbdp.uniq.indorothea[i,9]
tf.interest.ensembl <- tf.dbdp.uniq.indorothea$Ensembl.Gene.ID[i]
tf.genotypes <- variants.range.genotypes.matched[which(variant.ids.range == var),]
maf.var <- maf(as.numeric(tf.genotypes))
print(paste("maf:", maf.var))
# Remove targets within 1MB from consideration
dorothea.possible <- dorothea[countOverlaps(gr1,gr2[i]) == 0,]
# Filter so that variants are not on the same chromosome within 1MB of the target gene
tf.targets.interest <- dorothea.possible$ensembl_gene_id_target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
tf.targets.names <- dorothea.possible$target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
mor.targets <- dorothea.possible$mor[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
confidence.targets <- dorothea.possible$confidence[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
chrom.targets <- dorothea.possible$target_chromosome[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
targets.start <- dorothea.possible$target_start[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
targets.end <- dorothea.possible$target_end[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
target.range <-paste(targets.start,targets.end, sep = "-" )
chrom.target.range <- paste(chrom.targets,target.range, sep = ":")
chrom.target.range <- paste("chr",chrom.target.range, sep = "")
df.targets <- cbind.data.frame(tf.targets.interest, tf.targets.names, chrom.target.range, mor.targets, confidence.targets)
df.targets.intersect <- df.targets[which(df.targets$tf.targets.interest %in% rownames(expression.data.samples.matched)),]
p.values <- c()
cof <- c()
adjusted_rsquared <- c()
intersect.tf.targets.expression.data <- intersect(df.targets$tf.targets.interest, rownames(expression.data.samples.matched))
expression_data_transposed <- t(expression.data.samples.matched[intersect.tf.targets.expression.data,])
variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
tab <- cbind.data.frame(tf.genotypes, expression_data_transposed, stringsAsFactors = FALSE)
tab[,2:ncol(tab)] <- sapply(tab[,2:ncol(tab)], as.numeric)
colnames(tab) <- c("tf.genotypes", df.targets.intersect$tf.targets.names)
df.for.plotting <- data.frame(genotype = factor(tab$tf.genotypes), expression = tab[,targets_interest[j]])
title <- paste(variants_interest[j], "\n", tissue, sep = "")
label_text <- paste("TF: ", tf.interest, "\nTarget gene: ", targets_interest[j], "\nMAF: ", round(maf.var, 2), "\n n = ", length(isect), sep = "")
annotation <- data.frame(x = x_axis[j], y = y_axis[j], label = label_text)
plot_object <- ggplot(df.for.plotting, aes(genotype, expression)) +
	geom_jitter(colour="darkorange", position=position_jitter(width=0.25)) +  
	geom_boxplot(outlier.size=0, alpha=0.6, fill="steelblue") +
	geom_smooth(method = 'lm',col="darkred", aes(group=1), se=FALSE) +
	ggtitle(title) +
	theme_classic() +
	theme(plot.title = element_text(hjust = 0.5)) +
	geom_label(data = annotation, aes(x = x, y = y, label = label), color = "red3", size = 3, fontface = "bold") +	
	labs( x = NULL, y = NULL)
	print(plot_object)
	})
}
print(paste("length plots:", length(myplots)))
filename <- "eqtl_sig_results_alltissues.pdf"
print(paste("Filename:", filename))
ggar_obj <- ggarrange(plotlist = myplots, nrow = 2, ncol = 3, common.legend = TRUE)
ggar_obj_annotated <- annotate_figure(ggar_obj, left = textGrob("Target gene expression", rot = 90, vjust = 1, gp = gpar(cex = 1.6)), bottom = textGrob("Genotype", gp = gpar(cex = 1.6)))
ggsave(filename, ggar_obj_annotated, width=15, height=8.5)


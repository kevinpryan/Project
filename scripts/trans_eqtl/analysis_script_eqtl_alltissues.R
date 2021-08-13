#!/usr/bin/Rscript

# This script carried out trans-eQTL analysis on all GTEx tissues
# Script was developed with R version 4.1.0
# Files required: genotype data GTEx extracted according to instructions (see README), dbdp_nssnp_barrera_with_GTEx_varid_cadd.txt, covariates - download from GTEx, do not unzip gz files, source: https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_covariates.tar.gz
# gtex.psam - GTEx sample info, dorothea_hs_with_target_info.txt. Change file paths as necessary. GTEx expression data normalised - https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_expression_matrices.tar - again, do not gunzip files
# To run for example the default approach using the fractional rank data, filtering approach and all covariates - Rscript interaction_analysis_script_one_file.R default frac_rank pass all_covariates filtered_tf
# args: args[1] - default or average
# args[2] - full path to GTEx expression expression data (including directory name)
# args[3] - all_covariates, reduced_covariates or no_covariates
# args[4] - filtered_tf or no_filter_tf
# Edit line 250 and/or 518 to write to file of interest 

# Load in packages
library(stringr)
library(biomaRt)
library(vroom)
library(dplyr, warn.conflicts = FALSE)
library(GenomicRanges, quietly = T)
library(tibble)
print(paste( "Time", Sys.time(), sep = " "))
# Take user input (average expression or use all expression values)
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied ('average' or 'default' )", call.=FALSE)
} else if (args[1] != "average" && args[1] != "default") {
  stop("Argument must either be 'average' or 'default'", call. = FALSE)
}

print(paste("Arg 1:",args[1], "arg 2:",args[2],"arg 3:", args[3], "arg 4:", args[4], sep = " "))
print(paste("Time: ", Sys.time(), sep = ""))
start.time <- Sys.time()
# Point to file names
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

for (m in 1:length(f)){
print(paste("tissue:", tissues[m]))
print(paste("tissue number:", m))
expression_file <- paste("/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_expression_matrices/", f[m], sep = "")
print(paste("expression file:", expression_file))
cov_file <- paste("/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/", covariates[m], sep = "")
print(paste("cov_file:", cov_file))
cov <- vroom(cov_file)
cov_t <- t(cov[-1])
colnames(cov_t) <- cov$ID

if (args[3] == "reduced_covariates"){
	reduced_covariates <- c("PC1", "PC2", "PC3", "PC4", "PC5", "pcr", "platform", "sex")
	cov_t <- cov_t[,reduced_covariates]
}
# Read in expression data
expression_data <- vroom(expression_file)
print(paste("nrow expression data initially:", nrow(expression_data)))

#read in dorothea file - created manually with added features, as of 17/06/2021
dorothea <- read.table(dorothea_file, header = T, sep = "\t", check.names = F)
Gene_full <- expression_data$"gene_id"
Gene_split <- str_split_fixed(Gene_full, pattern = "\\.", n = 2)
Gene <- Gene_split[,1]
expression_data <- expression_data[Gene %in% dorothea$ensembl_gene_id_target,]
expression_data <- expression_data %>% 
		   add_column(ensembl_gene_id_target_short = str_split_fixed(expression_data$gene_id, pattern = "\\.", n = 2)[,1], .after = "gene_id")
# remove rows with duplicated gene names - we don't know the ensembl ids for the frac_rank expression data so I wanted to do the same to the whole blood data so that we can compare like with like
singletons <- names(which(table(expression_data$ensembl_gene_id_target_short) == 1))
expression_data_notduplicated <- expression_data[expression_data$ensembl_gene_id_target_short %in% singletons, ]

print(paste("dim(expression_data_notduplicated)", dim(expression_data_notduplicated)))
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
dorothea <- read.table(dorothea_file, header = T, sep = "\t", check.names = F)

# Only keep the variants in transcription factors that have data in the dorothea database
tf.dbdp.uniq.indorothea <- tf.dbdp.uniq[tf.dbdp.uniq$Gene.symbol %in% dorothea$tf,]

print(paste("Number of transcription factors with GTEx variants in desired frequency range available in dorothea database:",nrow(tf.dbdp.uniq.indorothea)))

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
df.out <- data.frame()
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

not_all_na <- function(x) any(!is.na(x))

if(args[1] == "default"){
print("entering default method")
for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	print(paste("Variant number:", i, sep = " "))
	var <- tf.dbdp.uniq.indorothea[i,22]
	tf.interest <- tf.dbdp.uniq.indorothea[i,9]
	tf.interest.ensembl <- tf.dbdp.uniq.indorothea$Ensembl.Gene.ID[i]	
	tf.genotypes <- variants.range.genotypes.matched[which(variant.ids.range == var),]
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
	if (args[4] == "filtered_tf"){
		print("entering filtering method...")
		expression.data.tf <- expression.data.samples.matched[tf.interest.ensembl,]
		if (!not_all_na(expression.data.tf)){
                next
        	}
		expression.data.tf.vec <- unlist(expression.data.tf)
		for (j in 2:ncol(tab)){
		if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
                                linear.model <- lm(tab[,j] ~ expression.data.tf.vec + as.matrix(cov_t))
                        } else if (args[3] == "no_covariates") {
                                linear.model <- lm(tab[,j] ~ expression.data.tf.vec)
                        }

                linear.model.pval <- coef(summary(linear.model))["expression.data.tf.vec", "Pr(>|t|)"]
                sign.linear.model <- sign(coef(summary(linear.model))["expression.data.tf.vec", "Estimate"])
                adjusted.rsquared.linear.model <- summary(linear.model)$adj.r.squared
                sign.mor <- sign(mor.targets[j-1])
                if (linear.model.pval < 0.05 & sign.linear.model == sign.mor){ 
			if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
                                linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes + as.matrix(cov_t))
	                 } else if (args[3] == "no_covariates"){
                                linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes)
        	         }
                	 p.values <- c(p.values, summary(linear.model.eqtl)$coefficients[2,4])
                 	 cof <- c(cof,summary(linear.model.eqtl)$coefficients[2,1])
			 adjusted_rsquared <- c(adjusted_rsquared, summary(linear.model.eqtl)$adj.r.squared)
                }  else  {
                        p.values <- c(p.values, NA)
			cof <- c(cof, NA)
			adjusted_rsquared <- c(adjusted_rsquared, NA)
			}
		}
	} else if (args[4] == "no_filter_tf"){	
	for (j in 2:ncol(tab)){
		if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
				linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes + as.matrix(cov_t))
                 } else if (args[3] == "no_covariates"){
                            	linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes)
                 }
		 p.values <- c(p.values, summary(linear.model.eqtl)$coefficients[2,4])
 		 cof <- c(cof,summary(linear.model.eqtl)$coefficients[2,1])
		}
	}
	results <- cbind(rep(var, length(intersect.tf.targets.expression.data)), rep(tf.interest, length(intersect.tf.targets.expression.data)), intersect.tf.targets.expression.data, df.targets.intersect$tf.targets.names, df.targets.intersect$chrom.target.range, p.values, cof, adjusted_rsquared, df.targets.intersect$mor.targets, df.targets.intersect$confidence.targets)
	df.out <- rbind.data.frame(df.out,results)
  }
  colnames(df.out) <- c("Variant_GTEx_ID", "TF_containing_variant", "Target_gene_id", "Target_gene_name", "Target_chromosome", "p_value", "coefficient", "adjusted_rsquared", "Mode_of_regulation","Target_confidence")
  df.out$"p_value" <- as.numeric(df.out$"p_value")
  df.out.na.remove <- na.omit(df.out)
  FDR.adjusted.pvals <- p.adjust(df.out.na.remove$"p_value", method = "BH") # adjusted using Benjamini and Hochberg method
  FDR.adjusted.results <- cbind(df.out.na.remove, FDR.adjusted.pvals)
  FDR.adjusted.results.cols.ordered <- select(FDR.adjusted.results, "Variant_GTEx_ID", "TF_containing_variant", "Target_gene_id", "Target_gene_name", "Target_chromosome", "p_value", "FDR.adjusted.pvals", "coefficient", "adjusted_rsquared", "Mode_of_regulation", "Target_confidence")
  sig.results.bh <- df.out.na.remove[FDR.adjusted.pvals <= 0.05,]
  sig.results.bh.added <- sig.results.bh %>% mutate(FDR.adjusted.pvalues = FDR.adjusted.pvals[FDR.adjusted.pvals <= 0.05]) 
  results.sorted <- FDR.adjusted.results.cols.ordered[order(FDR.adjusted.results.cols.ordered$FDR.adjusted.pvals),]
  print("writing to file...")
  write.table(results.sorted, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/all_tissues/all_targets/eqtl/", args[3], "/", tissues[m], "_", args[3], "_", args[4], "_",  args[1], "_FDR_adjusted_results.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)

} else if (args[1] == "average") {
print("entering averaging method")
	df.out.upreg <- data.frame()
	df.out.downreg <- data.frame()
	no.of.tests <- 0
	  for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	    var <- tf.dbdp.uniq.indorothea[i,22]
	    tf.interest <- tf.dbdp.uniq.indorothea[i,9]
	    tf.interest.ensembl <- tf.dbdp.uniq.indorothea$Ensembl.Gene.ID[i]
	    tf.genotypes <- variants.range.genotypes.matched[which(variant.ids.range == var),]
	    # Remove targets within 1MB from consideration
	    dorothea.possible <- dorothea[countOverlaps(gr1,gr2[i]) == 0,]
	    # Filter so that variants are not on the same chromosome within 1MB of the target gene
            tf.targets.interest <- dorothea.possible$ensembl_gene_id_target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
	    mor.targets <- dorothea$mor[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	    confidence.targets <- dorothea$confidence[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
	    chrom.targets <- dorothea.possible$target_chromosome[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
            targets.start <- dorothea.possible$target_start[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
            targets.end <- dorothea.possible$target_end[which(dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest)]
            target.range <-paste(targets.start,targets.end, sep = "-" )
            chrom.target.range <- paste(chrom.targets,target.range, sep = ":")
            chrom.target.range <- paste("chr",chrom.target.range, sep = "")
	    tf.targets.upreg <- tf.targets.interest[mor.targets == 1]
	    tf.targets.downreg <- tf.targets.interest[mor.targets == -1]
	    confidence.targets.upreg <- confidence.targets[mor.targets == 1]
	    confidence.targets.downreg <- confidence.targets[mor.targets == -1]
	    p.values.upreg <- c()
	    p.values.downreg <- c()
	    cof.upreg <- c()
	    cof.downreg <- c()
	    adjusted_rsquared_downreg <- c()
	    variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
       if (args[4] == "filtered_tf"){
                expression.data.tf <- expression.data.samples.matched[tf.interest.ensembl,]
                if (!not_all_na(expression.data.tf)){
		cof.upreg <- NA
		p.values.upreg <- NA
                next
                }
                expression.data.tf.vec <- unlist(expression.data.tf)
		intersect.tf.targets.upreg.expression.data <- intersect(tf.targets.upreg, rownames(expression.data.samples.matched))
		if (length(tf.targets.upreg) == 0 | length(intersect.tf.targets.upreg.expression.data) == 0){
				cof.upreg <- NA
				p.values.upreg <- NA
       		} else {	
			expression_data_transposed_upreg <- t(expression.data.samples.matched[intersect.tf.targets.upreg.expression.data,])
			tab_upreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_upreg, stringsAsFactors = FALSE)
			if (ncol(tab_upreg) == 2){
				if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
					linear.model.upreg <- lm(tab_upreg[,2] ~ expression.data.tf.vec + as.matrix(cov_t))
				} else if (args[3] == "no_covariates") {
					linear.model.upreg <- lm(tab_upreg[,2] ~ expression.data.tf.vec)
				}
				linear.model.pval <- coef(summary(linear.model.upreg))["expression.data.tf.vec", "Pr(>|t|)"]
				sign.linear.model.upreg <- sign(coef(summary(linear.model.upreg))["expression.data.tf.vec", "Estimate"])
				adjusted.rsquared.linear.model <- summary(linear.model.upreg)$adj.r.squared
				if (linear.model.pval < 0.05 & sign.linear.model.upreg == 1){
				        tab_rank_upreg <- rank(tab_upreg[,2])
                                        tab_rank_upreg_mean <- tab_rank_upreg
					if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
						linear.model.eqtl <- lm(tab_rank_upreg_mean ~ tab_upreg$tf.genotypes + as.matrix(cov_t))
				 	} else if (args[3] == "no_covariates"){
						linear.model.eqtl <- lm(tab_rank_upreg_mean ~ tab_upreg$tf.genotypes)
				 	}	
				 p.values.upreg <- summary(linear.model.eqtl)$coefficients[2,4]
				 cof.upreg <- summary(linear.model.eqtl)$coefficients[2,1]
				 adjusted_rsquared <- summary(linear.model.eqtl)$adj.r.squared
				} else {
					p.values.upreg <- NA
					cof.upreg <- NA
                                 	adjusted_rsquared <- c(adjusted_rsquared, NA)
				}
			} else {
				for (j in 2:ncol(tab_upreg)) {
					if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
                                        	linear.model.upreg <- lm(tab_upreg[,j] ~ expression.data.tf.vec + as.matrix(cov_t))
					} else if (args[3] == "no_covariates") {
						linear.model.upreg <- lm(tab_upreg[,j] ~ expression.data.tf.vec)
					}
					linear.model.pval <- coef(summary(linear.model.upreg))["expression.data.tf.vec", "Pr(>|t|)"]
					sign.linear.model.upreg <- sign(coef(summary(linear.model.upreg))["expression.data.tf.vec", "Estimate"])
                                	adjusted.rsquared.linear.model <- summary(linear.model.upreg)$adj.r.squared
					if (linear.model.pval >= 0.05 | sign.linear.model.upreg != 1){
						tab_upreg[,j] <- rep(NA, nrow(tab_upreg))
					}
				}
				tab_upreg <- tab_upreg %>% select(where(not_all_na))
				if (ncol(tab_upreg) > 2){
                                tab_rank_upreg <- apply(tab_upreg[,2:ncol(tab_upreg)], FUN = rank, MARGIN = 2)
                                tab_rank_upreg <- cbind(tab_upreg$tf.genotypes, tab_rank_upreg)
                                tab_rank_upreg_mean <- apply(tab_rank_upreg[,2:ncol(tab_rank_upreg)], FUN = mean, MARGIN = 1)
				} else if (ncol(tab_upreg) == 2){
					tab_rank_upreg <- rank(tab_upreg[,2])
                                        tab_rank_upreg_mean <- tab_rank_upreg
				}
				tab_rank_upreg_mean_genotypes <- cbind.data.frame(tab_upreg$tf.genotypes, tab_rank_upreg_mean)
		                colnames(tab_rank_upreg_mean_genotypes)[1] <- "tf.genotypes"
				if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
                			linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
                                } else if (args[3] == "no_covariates") {
					linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes, na.action = na.exclude)
				}
				cof.upreg <- summary(linear_model_upreg)$coefficients[2,1]
				p.values.upreg <- summary(linear_model_upreg)$coefficients[2,4]

			}
		}
                intersect.tf.targets.downreg.expression.data <- intersect(tf.targets.downreg, rownames(expression.data.samples.matched))
		if (length(tf.targets.downreg) == 0 | length(intersect.tf.targets.downreg.expression.data) == 0){
				cof.downreg <- NA
				p.values.downreg <- NA
                } else {
                        expression_data_transposed_downreg <- t(expression.data.samples.matched[intersect.tf.targets.downreg.expression.data,])
                        tab_downreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_downreg, stringsAsFactors = FALSE)
                        if (ncol(tab_downreg) == 2){
                                if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
                                        linear.model <- lm(tab_downreg[,2] ~ expression.data.tf.vec + as.matrix(cov_t))
                                } else if (args[3] == "no_covariates") {
                                        linear.model <- lm(tab_downreg[,2] ~ expression.data.tf.vec)
                                }
                                linear.model.pval <- coef(summary(linear.model))["expression.data.tf.vec", "Pr(>|t|)"]
                                sign.linear.model.downreg <- sign(coef(summary(linear.model))["expression.data.tf.vec", "Estimate"])
                                adjusted.rsquared.linear.model <- summary(linear.model)$adj.r.squared
				if (linear.model.pval < 0.05 & sign.linear.model.downreg == -1){   
					tab_rank_downreg <- rank(tab_downreg[,2])
                                        tab_rank_downreg_mean <- tab_rank_downreg
				if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
                                        linear.model.eqtl <- lm(tab_rank_downreg_mean ~ tab_downreg$tf.genotypes + as.matrix(cov_t))
                                 } else if (args[3] == "no_covariates"){
                                        linear.model.eqtl <- lm(tab_rank_downreg_mean ~ tab_downreg$tf.genotypes)
                                 }
				 p.values.downreg <- summary(linear.model.eqtl)$coefficients[2,4]
                                 cof.downreg <- summary(linear.model.eqtl)$coefficients[2,1]
                                 adjusted_rsquared_downreg <- c(adjusted_rsquared_downreg, summary(linear.model.eqtl)$adj.r.squared)
				} else {
					p.values.downreg <- NA
					cof.downreg <- NA
                                        adjusted_rsquared_downreg <- c(adjusted_rsquared_downreg, NA)
				}
			} else {
                                for (j in 2:ncol(tab_downreg)) {
                                        if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
                                                linear.model <- lm(tab_downreg[,j] ~ expression.data.tf.vec + as.matrix(cov_t))
                                        } else if (args[3] == "no_covariates") {
                                                linear.model <- lm(tab_downreg[,j] ~ expression.data.tf.vec)
                                        }
                                        linear.model.pval <- coef(summary(linear.model))["expression.data.tf.vec", "Pr(>|t|)"]
                                        sign.linear.model.downreg <- sign(coef(summary(linear.model))["expression.data.tf.vec", "Estimate"])
                                        adjusted.rsquared.linear.model <- summary(linear.model)$adj.r.squared
					if (linear.model.pval >= 0.05 | sign.linear.model.downreg != -1){
                                                tab_downreg[,j] <- rep(NA, nrow(tab_downreg))
                                        }
                                }
                                tab_downreg <- tab_downreg %>% select(where(not_all_na))
				if (ncol(tab_downreg) > 2){
                                tab_rank_downreg <- apply(tab_downreg[,2:ncol(tab_downreg)], FUN = rank, MARGIN = 2)
                                tab_rank_downreg <- cbind(tab_downreg$tf.genotypes, tab_rank_downreg)
                                tab_rank_downreg_mean <- apply(tab_rank_downreg[,2:ncol(tab_rank_downreg)], FUN = mean, MARGIN = 1)
                                } else if (ncol(tab_downreg) == 2){
                                        tab_rank_downreg <- rank(tab_downreg[,2])
                                        tab_rank_downreg_mean = tab_rank_downreg
                                } else {
					cof.downreg <- NA
					p.values.downreg <- NA
					next
				}
                                tab_rank_downreg_mean_genotypes <- cbind.data.frame(tab_downreg$tf.genotypes, tab_rank_downreg_mean)
                                colnames(tab_rank_downreg_mean_genotypes)[1] <- "tf.genotypes"
                                if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
                                        linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
                                } else if (args[3] == "no_covariates") {				                                        
					linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes, na.action = na.exclude)
				}
				cof.downreg <- summary(linear_model_downreg)$coefficients[2,1]
                                p.values.downreg <- summary(linear_model_downreg)$coefficients[2,4]
				}
		}

        } else if (args[4] == "no_filter_tf"){
	    if (length(tf.targets.upreg) == 0){
	    	cof.upreg <- NA
            	p.values.upreg <- NA
            } else {
	    	intersect.tf.targets.upreg.expression.data <- intersect(tf.targets.upreg, rownames(expression.data.samples.matched))
            	expression_data_transposed_upreg <- t(expression.data.samples.matched[intersect.tf.targets.upreg.expression.data,])
            	tab_upreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_upreg, stringsAsFactors = FALSE)
            	if (ncol(tab_upreg) == 2){
	    		tab_rank_upreg <- rank(tab_upreg[,2])
            		tab_rank_upreg_mean <- tab_rank_upreg
                } else {
	    		tab_rank_upreg <- apply(tab_upreg[,2:ncol(tab_upreg)], FUN = rank, MARGIN = 2)
	    		tab_rank_upreg <- cbind(tab_upreg$tf.genotypes, tab_rank_upreg)
	    		tab_rank_upreg_mean <- apply(tab_rank_upreg[,2:ncol(tab_rank_upreg)], FUN = mean, MARGIN = 1)
		}
   	    	tab_rank_upreg_mean_genotypes <- cbind.data.frame(tab_upreg$tf.genotypes, tab_rank_upreg_mean)
            	colnames(tab_rank_upreg_mean_genotypes)[1] <- "tf.genotypes"
		if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
			linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
		} else if (args[3] == "no_covariates"){			                        
			linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes, na.action = na.exclude)
		}
            	cof.upreg <- c(cof.upreg,summary(linear_model_upreg)$coefficients[2,1])
            	p.values.upreg <- c(p.values.upreg, summary(linear_model_upreg)$coefficients[2,4])
	    	no.of.tests <- no.of.tests + 1
		} 
            if (length(tf.targets.downreg) == 0){
	    	cof.downreg <- NA
            	p.values.downreg <- NA
	    } else {
	    	intersect.tf.targets.downreg.expression.data <- intersect(tf.targets.downreg, rownames(expression.data.samples.matched))
		if (length(intersect.tf.targets.downreg.expression.data) == 0){
		cof.downreg <- NA
                p.values.downreg <- NA
		} else {
            	expression_data_transposed_downreg <- t(expression.data.samples.matched[intersect.tf.targets.downreg.expression.data,])
            	tab_downreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_downreg, stringsAsFactors = FALSE)
		}
	    	if (ncol(tab_downreg) == 2){
	    		tab_rank_downreg <- rank(tab_downreg[,2])
	   		tab_rank_downreg_mean <- tab_rank_downreg
		} else {
	    		tab_rank_downreg <- apply(tab_downreg[,2:ncol(tab_downreg)], FUN = rank, MARGIN = 2)
	    		tab_rank_downreg <- cbind(tab_downreg$tf.genotypes, tab_rank_downreg)
            		tab_rank_downreg_mean <- apply(tab_rank_downreg[,2:ncol(tab_rank_downreg)], FUN = mean, MARGIN = 1)
	    	}
            	tab_rank_downreg_mean_genotypes <- cbind.data.frame(tab_downreg$tf.genotypes, tab_rank_downreg_mean)
            	colnames(tab_rank_downreg_mean_genotypes)[1] <- "tf.genotypes"
		if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
            		linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
        	} else if (args[3] == "no_covariates") {		                            
			linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
		}
            	colnames(tab_rank_downreg_mean_genotypes)[1] <- "tf.genotypes"
		if (args[3] == "all_covariates" | args[3] == "reduced_covariates"){
            		linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
        	} else if (args[3] == "no_covariates") {		                            
			linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
		}
		cof.downreg <- c(cof.downreg,summary(linear_model_downreg)$coefficients[2,1])
            	p.values.downreg <- c(p.values.downreg, summary(linear_model_downreg)$coefficients[2,4])
		}
            }
	    results.upreg <- cbind.data.frame(var, tf.interest, p.values.upreg, cof.upreg)
            results.downreg <- cbind.data.frame(var,tf.interest, p.values.downreg, cof.downreg)
            results.upreg$mor <- rep(1, nrow(results.upreg))
            results.downreg$mor <- rep(-1, nrow(results.downreg))
                if (ncol(results.upreg) == ncol(results.downreg)){
                        colnames(results.upreg) <- c("Variant_GTEx_ID", "TF_containing_variant", "p_value", "coefficient", "mor")
                        colnames(results.downreg) <- colnames(results.upreg)
                        df.out.upreg <- rbind.data.frame(df.out.upreg,results.upreg)
                        df.out.downreg <- rbind.data.frame(df.out.downreg, results.downreg)
                } else if (ncol(results.upreg) < ncol(results.downreg)){
                        colnames(results.downreg) <- c("Variant_GTEx_ID", "TF_containing_variant", "p_value", "coefficient", "mor")
                        df.out.downreg <- rbind.data.frame(df.out.downreg, results.downreg)
                } else if (ncol(results.upreg) > ncol(results.downreg)){
                        colnames(results.upreg) <- c("Variant_GTEx_ID", "TF_containing_variant", "p_value", "coefficient", "mor")
                        df.out.upreg <- rbind.data.frame(df.out.upreg,results.upreg)
                        }
	    }
	df.out.combined <- rbind.data.frame(df.out.upreg,df.out.downreg)    
        df.out.combined$"p_value" <- as.numeric(df.out.combined$"p_value")
        df.out.na.remove <- na.omit(df.out.combined)
        df.out.na.remove$FDR.adjusted.pvals <- p.adjust(as.numeric(df.out.na.remove$"p_value"), method = "BH")
	df.out.na.remove.cols.ordered <- select(df.out.na.remove, "Variant_GTEx_ID", "TF_containing_variant", "mor", "p_value", "FDR.adjusted.pvals", "coefficient")
	results.sorted <- df.out.na.remove.cols.ordered[order(df.out.na.remove.cols.ordered$FDR.adjusted.pvals),]
	print("writing to file...")
	write.table(results.sorted, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/all_tissues/avg_targets/eqtl/", args[3], "/", tissues[m], "_", args[3], "_", args[4], "_",  args[1], "_FDR_adjusted_results.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
	# closing bracket of else if statement
	}
 #closing bracket of for loop for files
}
end.time <- Sys.time()
print(paste("End time:", end.time, sep = ""))
print(paste("time taken:", (end.time - start.time)))


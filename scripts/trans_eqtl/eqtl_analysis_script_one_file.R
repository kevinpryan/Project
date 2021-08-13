#!/usr/bin/Rscript

# This script carried out trans-eQTL analysis on a single tissue (whole blood or fractional rank data)
# Script was developed with R version 4.1.0

# Load packages
library(stringr)
library(biomaRt)
library(vroom)
library(dplyr, warn.conflicts = FALSE)
library(GenomicRanges, quietly = T)
library(optparse)
print(paste( "Time", Sys.time(), sep = " "))
 
option_list = list(
 make_option(c("-m", "--method"), type="character", default=NULL, 
              help="either default or average", metavar="character"),
 make_option(c("-d", "--dataset"), type="character", default=NULL, 
              help="either frac_rank or whole_blood_norm", metavar="character"),
 make_option(c("-c", "--covariates"), type="character", default=NULL,
              help="either no_covariates, reduced covariates or all_covariates", metavar="character"),
 make_option(c("-f", "--filter"), type="character", default=NULL,
              help="either filtered_tf or no_filter_tf", metavar="character"),
 make_option(c("-o", "--out"), type="character", default="eqtl_out.txt",
              help="specify outfile name", metavar="character")

); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (length(opt) < 4){
  print_help(opt_parser)
  stop("At least 4 arguments must be supplied (don't have to specify outfile)", call.=FALSE)
}
print("opt:")
print(opt)
opt$method -> method
opt$dataset -> dataset
opt$covariates -> covariates
opt$filter -> filter
opt$out -> out
# Point to file names - hard coded in
genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 1.0e-02
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"
if (dataset == "whole_blood_norm") {
	expression_file <- "/data/kryan/project/gtex/sample_data/Whole_blood_v8.normalized_filtered_expression.bed"
	cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
} else if (dataset == "frac_rank") {
	expression_file <- "/data/kryan/project/gtex/analysis/frac_rank.out.all"
	cov_file <- "/data/kryan/project/gtex/sample_data/covariates_frac_rank_extracted_exp_and_geno_rownames_altered.txt"
}
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

# Read in covariates
if (dataset == "frac_rank"){
	cov <- read.table(cov_file, check.names=F)
	cov_t <- t(cov)
} else if (dataset == "whole_blood_norm") {
	cov <- vroom(cov_file)
	cov_t <- t(cov[-1])
	colnames(cov_t) <- cov$ID
}

if (covariates == "reduced_covariates"){
	reduced_covariates <- c("PC1", "PC2", "PC3", "PC4", "PC5", "pcr", "platform", "sex")
        cov_t <- cov_t[,reduced_covariates]
}

# convert to numeric - get some warnings
variants.genotypes <- genotype.data[,10:ncol(genotype.data)]
variants.numeric <- apply(variants.genotypes, 2, genotype_to_numeric)
variant.ids <- genotype.data$V3
# Read in expression data
expression_data <- vroom(expression_file)

# remove rows with duplicated gene names - we don't know the ensembl ids for the frac_rank expression data so I wanted to do the same to the whole blood data so that we can compare like with like
if (dataset == "frac_rank"){
	singletons <- names(which(table(expression_data$Gene) == 1))
	expression_data_notduplicated <- expression_data[expression_data$Gene %in% singletons, ]

} else if (dataset == "whole_blood_norm"){
	singletons <- names(which(table(expression_data$ensembl_gene_id_target_short) == 1))
	expression_data_notduplicated <- expression_data[expression_data$ensembl_gene_id_target_short %in% singletons, ]
}

# Read in sample ids and match up with the sample ids in the expression data file
psam.ids <- read.table(psam)
colnames(variants.numeric) <- psam.ids$V1
if (dataset == "frac_rank"){
	isect <- intersect(colnames(variants.numeric),colnames(expression_data[2,]))
} else if (dataset == "whole_blood_norm"){
	isect <- intersect(colnames(variants.numeric),colnames(expression_data[6,]))
}
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
if (dataset == "frac_rank"){
	rownames(expression.data.samples.matched) <- expression_data_notduplicated$Gene
} else if (dataset == "whole_blood_norm") {
	rownames(expression.data.samples.matched) <- expression_data_notduplicated$ensembl_gene_id_target_short
}

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

# Function to test whether all elements of a vector are NA - returns TRUE if all elements are NOT NA, FALSE if all elements are NA 
not_all_na <- function(x) any(!is.na(x))

if(method == "default"){
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
	if (dataset == "whole_blood_norm"){
		tf.targets.interest <- dorothea.possible$ensembl_gene_id_target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
	} else if (dataset == "frac_rank"){
		tf.targets.interest <- dorothea.possible$target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
	}
	# Gather information about transcription factor from dorothea.possible dataframe
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
	# vectors for outputs
	p.values <- c()
	cof <- c()
	adjusted_rsquared <- c()
	# get expression data for the targets of the transcription factor
	intersect.tf.targets.expression.data <- intersect(df.targets$tf.targets.interest, rownames(expression.data.samples.matched))
	# need to transpose as we want the rows to be samples and the columns to be genes
	expression_data_transposed <- t(expression.data.samples.matched[intersect.tf.targets.expression.data,])
	variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	tab <- cbind.data.frame(tf.genotypes, expression_data_transposed, stringsAsFactors = FALSE)
	tab[,2:ncol(tab)] <- sapply(tab[,2:ncol(tab)], as.numeric)
	# Give to option to filter based on linear relationship between expression of TF and expression of target gene (p < 0.05). Direction of relationship must be the same as that in dorothea to be considered
	if (filter == "filtered_tf"){
		if (dataset == "whole_blood_norm"){
			expression.data.tf <- expression.data.samples.matched[tf.interest.ensembl,]
		} else if (dataset == "frac_rank") {
			expression.data.tf <- expression.data.samples.matched[tf.interest,]
		}
		print("expression.data.tf")
		print(expression.data.tf)
		# skip this variant if the TF is not expressed
		if (!not_all_na(expression.data.tf)){
                	print(paste("No expression data for tf, skipping...", tf.interest))
                	next
        	}
		expression.data.tf.vec <- unlist(expression.data.tf)
		# Fit first linear model - exp of TF vs exp of target
		for (j in 2:ncol(tab)){
			if (covariates == "all_covariates" | covariates == "reduced_covariates"){
                                linear.model <- lm(tab[,j] ~ expression.data.tf.vec + as.matrix(cov_t))
                        } else if (covariates == "no_covariates") {
                                linear.model <- lm(tab[,j] ~ expression.data.tf.vec)
                        }
		# Extract out figures of interest
                linear.model.pval <- coef(summary(linear.model))["expression.data.tf.vec", "Pr(>|t|)"]
                sign.linear.model <- sign(coef(summary(linear.model))["expression.data.tf.vec", "Estimate"])
                adjusted.rsquared.linear.model <- summary(linear.model)$adj.r.squared
                sign.mor <- sign(mor.targets[j-1])
		# Filtering step - only fit second linear model (eQTL model) if it passes the filter. If not, skip to the next variant
                if (linear.model.pval < 0.05 & sign.linear.model == sign.mor){ 
			if (covariates == "all_covariates" | covariates == "reduced_covariates"){
                                linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes + as.matrix(cov_t))
	                 } else if (covariates == "no_covariates"){
                                linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes)
        	         }
                	 p.values <- c(p.values, summary(linear.model.eqtl)$coefficients[2,4])
                 	 cof <- c(cof,summary(linear.model.eqtl)$coefficients[2,1])
			 adjusted_rsquared <- c(adjusted_rsquared, summary(linear.model.eqtl)$adj.r.squared)
                }  else  {
                       print(paste("skipping var:", var, "target:", intersect.tf.targets.expression.data[j-1]))
                        p.values <- c(p.values, NA)
			cof <- c(cof, NA)
			adjusted_rsquared <- c(adjusted_rsquared, NA)
			}
		}
	} else if (filter == "no_filter_tf"){	
	for (j in 2:ncol(tab)){
		if (covariates == "all_covariates" | covariates == "reduced_covariates"){
				linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes + as.matrix(cov_t))
				summary(linear.model.eqtl)
                 } else if (covariates == "no_covariates"){
                            	linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes)
                 }
		 p.values <- c(p.values, summary(linear.model.eqtl)$coefficients[2,4])
 		 cof <- c(cof,summary(linear.model.eqtl)$coefficients[2,1])
		 adjusted_rsquared <- c(adjusted_rsquared, summary(linear.model.eqtl)$adj.r.squared)
		}
	}
	results <- cbind(rep(var, length(intersect.tf.targets.expression.data)), rep(tf.interest, length(intersect.tf.targets.expression.data)), intersect.tf.targets.expression.data, df.targets.intersect$tf.targets.names, df.targets.intersect$chrom.target.range, p.values, cof, adjusted_rsquared, df.targets.intersect$mor.targets, df.targets.intersect$confidence.targets)
	df.out <- rbind.data.frame(df.out,results)
  }
  colnames(df.out) <- c("Variant_GTEx_ID", "TF_containing_variant", "Target_gene_id", "Target_gene_name", "Target_chromosome", "p_value", "coefficient", "adjusted_rsquared", "Mode_of_regulation","Target_confidence")
  print(paste("Number of tests carried out", nrow(df.out)))
  sig.values.df.out <- 0.05/nrow(df.out)
  print(paste("Significance value:", sig.values.df.out))
  df.out$"p_value" <- as.numeric(df.out$"p_value")
  df.out.na.remove <- na.omit(df.out)
  sig.results <- df.out.na.remove[df.out.na.remove$"p_value" <= sig.values.df.out,]
  print("Significant results (Bonferroni)")
  print(sig.results, row.names = FALSE)
  FDR.adjusted.pvals <- p.adjust(df.out.na.remove$"p_value", method = "BH") # adjusted using Benjamini and Hochberg method
  FDR.adjusted.results <- cbind(df.out.na.remove, FDR.adjusted.pvals)
  sig.results.bh <- df.out.na.remove[FDR.adjusted.pvals <= 0.05,]
  sig.results.bh.added <- sig.results.bh %>% mutate(FDR.adjusted.pvalues = FDR.adjusted.pvals[FDR.adjusted.pvals <= 0.05]) 
  print("FDR adjusted results") 
  print(sig.results.bh.added, row.names = FALSE)
  results.sorted <- FDR.adjusted.results[order(FDR.adjusted.results$FDR.adjusted.pvals),]
  write.table(results.sorted, file = out, quote = FALSE, sep = "\t", row.names = FALSE)
  #write.table(results.sorted, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/", args[2], "/all_targets/", args[4], "/", args[5],"/filter_sign/", args[2], "_FDR_adjusted_results_AF1_",args[4], "_", args[5], "_filter_sign_and_lm_2907.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)
#write.table(results.sorted, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/", args[2], "/all_targets/", args[4], "/", args[5], "/", args[2], "_FDR_adjusted_results_AF1_",args[4], "_", args[5], "_2907.txt", sep = ""), quote = FALSE, sep = "\t", row.names = FALSE)


} else if (method == "average") {
print("entering averaging method")
	df.out.upreg <- data.frame()
	df.out.downreg <- data.frame()
	no.of.tests <- 0
	  for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	    var <- tf.dbdp.uniq.indorothea[i,22]
	    print(paste("var:", var))
	    tf.interest <- tf.dbdp.uniq.indorothea[i,9]
	    print(paste("tf.interest", tf.interest))	
	    tf.interest.ensembl <- tf.dbdp.uniq.indorothea$Ensembl.Gene.ID[i]
	    tf.genotypes <- variants.range.genotypes.matched[which(variant.ids.range == var),]
	    # Remove targets within 1MB from consideration
	    dorothea.possible <- dorothea[countOverlaps(gr1,gr2[i]) == 0,]
	    # Filter so that variants are not on the same chromosome within 1MB of the target gene
	    if (dataset == "whole_blood_norm"){
                tf.targets.interest <- dorothea.possible$ensembl_gene_id_target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
            } else if (dataset == "frac_rank"){
                tf.targets.interest <- dorothea.possible$target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
            }
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
       if (filter == "filtered_tf"){
                if (args[2] == "whole_blood_norm"){
                        expression.data.tf <- expression.data.samples.matched[tf.interest.ensembl,]
                } else if (dataset == "frac_rank") {
                        expression.data.tf <- expression.data.samples.matched[tf.interest,]
                }
                if (!not_all_na(expression.data.tf)){
                print(paste("No expression data for tf, skipping...", tf.interest))
		cof.upreg <- NA
		p.values.upreg <- NA
                next
                }
                expression.data.tf.vec <- unlist(expression.data.tf)
		if (length(tf.targets.upreg) == 0){
				cof.upreg <- NA
				p.values.upreg <- NA
				print(paste("No upregulated genes in variant:", i))
       		} else {
			intersect.tf.targets.upreg.expression.data <- intersect(tf.targets.upreg, rownames(expression.data.samples.matched))
			expression_data_transposed_upreg <- t(expression.data.samples.matched[intersect.tf.targets.upreg.expression.data,])
			tab_upreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_upreg, stringsAsFactors = FALSE)
			if (ncol(tab_upreg) == 2){
				if (covariates == "all_covariates" | covariates == "reduced_covariates"){
					linear.model.upreg <- lm(tab_upreg[,2] ~ expression.data.tf.vec + as.matrix(cov_t))
				} else if (covariates == "no_covariates") {
					linear.model.upreg <- lm(tab_upreg[,2] ~ expression.data.tf.vec)
				}
				linear.model.pval <- coef(summary(linear.model.upreg))["expression.data.tf.vec", "Pr(>|t|)"]
				sign.linear.model.upreg <- sign(coef(summary(linear.model.upreg))["expression.data.tf.vec", "Estimate"])
				adjusted.rsquared.linear.model <- summary(linear.model.upreg)$adj.r.squared
				if (linear.model.pval < 0.05 & sign.linear.model.upreg == 1){
				        tab_rank_upreg <- rank(tab_upreg[,2])
                                        tab_rank_upreg_mean <- tab_rank_upreg
					if (covariates == "all_covariates" | covariates == "reduced_covariates"){
						linear.model.eqtl <- lm(tab_rank_upreg_mean ~ tab_upreg$tf.genotypes + as.matrix(cov_t))
				 	} else if (covariates == "no_covariates"){
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
					if (covariates == "all_covariates" | covariates == "reduced_covariates"){
                                        	linear.model.upreg <- lm(tab_upreg[,j] ~ expression.data.tf.vec + as.matrix(cov_t))
					} else if (covariates == "no_covariates") {
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
				if (covariates == "all_covariates" | covariates == "reduced_covariates"){
                			linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
                                } else if (covariates == "no_covariates") {
					linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes, na.action = na.exclude)
				}
				cof.upreg <- summary(linear_model_upreg)$coefficients[2,1]
				p.values.upreg <- summary(linear_model_upreg)$coefficients[2,4]
			}
		}
		if (length(tf.targets.downreg) == 0){
				cof.downreg <- NA
				p.values.downreg <- NA
                                print(paste("No downregulated genes in variant:", i))
                } else {
                        intersect.tf.targets.downreg.expression.data <- intersect(tf.targets.downreg, rownames(expression.data.samples.matched))
                        expression_data_transposed_downreg <- t(expression.data.samples.matched[intersect.tf.targets.downreg.expression.data,])
                        tab_downreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_downreg, stringsAsFactors = FALSE)
                        if (ncol(tab_downreg) == 2){
                                if (covariates == "all_covariates" | covariates == "reduced_covariates"){
                                        linear.model <- lm(tab_downreg[,2] ~ expression.data.tf.vec + as.matrix(cov_t))
                                } else if (covariates == "no_covariates") {
                                        linear.model <- lm(tab_downreg[,2] ~ expression.data.tf.vec)
                                }
                                linear.model.pval <- coef(summary(linear.model))["expression.data.tf.vec", "Pr(>|t|)"]
                                sign.linear.model.downreg <- sign(coef(summary(linear.model))["expression.data.tf.vec", "Estimate"])
                                adjusted.rsquared.linear.model <- summary(linear.model)$adj.r.squared
				if (linear.model.pval < 0.05 & sign.linear.model.downreg == -1){   
					tab_rank_downreg <- rank(tab_downreg[,2])
                                        tab_rank_downreg_mean <- tab_rank_downreg
				if (covariates == "all_covariates" | covariates == "reduced_covariates"){
                                        linear.model.eqtl <- lm(tab_rank_downreg_mean ~ tab_downreg$tf.genotypes + as.matrix(cov_t))
                                 } else if (covariates == "no_covariates"){
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
                                        if (covariates == "all_covariates" | covariates == "reduced_covariates"){
                                                linear.model <- lm(tab_downreg[,j] ~ expression.data.tf.vec + as.matrix(cov_t))
                                        } else if (covariates == "no_covariates") {
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
                                if (covariates == "all_covariates" | covariates == "reduced_covariates"){
                                        linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
                                } else if (covariates == "no_covariates") {				                                        
					linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes, na.action = na.exclude)
				}
				cof.downreg <- summary(linear_model_downreg)$coefficients[2,1]
                                p.values.downreg <- summary(linear_model_downreg)$coefficients[2,4]
				}
		}

        } else if (filter == "no_filter_tf"){
	    if (length(tf.targets.upreg) == 0){
	    	cof.upreg <- NA
            	p.values.upreg <- NA
		print(paste("No upregulated genes in variant:", i))
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
		if (covariates == "all_covariates" | covariates == "reduced_covariates"){
			linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
		} else if (covariates == "no_covariates"){			                        
			linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes, na.action = na.exclude)
		}
            	cof.upreg <- c(cof.upreg,summary(linear_model_upreg)$coefficients[2,1])
            	p.values.upreg <- c(p.values.upreg, summary(linear_model_upreg)$coefficients[2,4])
	    	no.of.tests <- no.of.tests + 1
            	print(paste("Upreg,", i))
		} 
            if (length(tf.targets.downreg) == 0){
	    	cof.downreg <- NA
            	p.values.downreg <- NA
            print(paste("No downregulated target genes for variant:", i))
	    } else {
            	print(paste("Trying downreg:", i))
	    	intersect.tf.targets.downreg.expression.data <- intersect(tf.targets.downreg, rownames(expression.data.samples.matched))
		if (length(intersect.tf.targets.downreg.expression.data) == 0){
			print("No overlap with expression data")
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
		if (covariates == "all_covariates" | covariates == "reduced_covariates"){
            		linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
        	} else if (covariates == "no_covariates") {		                            
			linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
		}
            	colnames(tab_rank_downreg_mean_genotypes)[1] <- "tf.genotypes"
		if (covariates == "all_covariates" | covariates == "reduced_covariates"){
            		linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ tab_rank_downreg_mean_genotypes$tf.genotypes + as.matrix(cov_t), na.action = na.exclude)
        	} else if (covariates == "no_covariates") {		                            
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
                        print(paste("ncol(results.upreg):", ncol(results.upreg)))
                        print(paste("ncol(reslts.downreg):", ncol(results.downreg)))
                        print("1")
                        colnames(results.upreg) <- c("Variant_GTEx_ID", "TF_containing_variant", "p_value", "coefficient", "mor")
                        colnames(results.downreg) <- colnames(results.upreg)
                        df.out.upreg <- rbind.data.frame(df.out.upreg,results.upreg)
                        df.out.downreg <- rbind.data.frame(df.out.downreg, results.downreg)
                } else if (ncol(results.upreg) < ncol(results.downreg)){
                        print("2")
                        print(paste("ncol(results.upreg):", ncol(results.upreg)))
                        print(paste("ncol(reslts.downreg):", ncol(results.downreg)))
                        colnames(results.downreg) <- c("Variant_GTEx_ID", "TF_containing_variant", "p_value", "coefficient", "mor")
                        df.out.downreg <- rbind.data.frame(df.out.downreg, results.downreg)
                } else if (ncol(results.upreg) > ncol(results.downreg)){
                        print("3")
                        print(paste("ncol(results.upreg):", ncol(results.upreg)))
                        print(paste("ncol(reslts.downreg):", ncol(results.downreg)))
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
	print(results.sorted)
	write.table(results.sorted, file = out, quote = FALSE, sep = "\t", row.names = FALSE)
	}
	
	


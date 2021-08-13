#!/usr/bin/Rscript
# The purpose of this script is to apply the filter and then see if the targets are enriched for trans-eQTLs vs non-targets
# Script was developed with R version 4.1.0
# Files required: genotype data GTEx extracted according to instructions (see README), dbdp_nssnp_barrera_with_GTEx_varid_cadd.txt, covariates - download from GTEx, do not unzip gz files, source: https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_covariates.tar.gz
# gtex.psam - GTEx sample info, dorothea_hs_with_target_info.txt. Change file paths as necessary. GTEx expression data normalised - https://storage.googleapis.com/gtex_analysis_v8/single_tissue_qtl_data/GTEx_Analysis_v8_eQTL_expression_matrices.tar - again, do not gunzip files
# To run for example the default approach using the fractional rank data, filtering approach and all covariates - Rscript enrichment_analysis_eqtl.R default frac_rank pass all_covariates filtered_tf
# args: args[1] - default or average
# args[2] - full path to GTEx expression expression data (including directory name)
# args[3] - all_covariates, reduced_covariates or no_covariates
# args[4] - filtered_tf or no_filter_tf
# Edit line 327 to write to file of interest 
# The purpose of this script is to do our filtering and then see if this changes whether the targets are enriched for trans-eQTLs vs non-targets
# Load in packages
library(stringr)
library(biomaRt)
library(vroom)
library(dplyr, warn.conflicts = FALSE)
library(GenomicRanges, quietly = T)

# Take user input (average expression or use all expression values)
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied ('average' or 'default' )", call.=FALSE)
} else if (args[1] != "average" && args[1] != "default") {
  stop("Argument must either be 'average' or 'default'", call. = FALSE)
}

print(paste("Arg 1:",args[1], "arg 2:",args[2],"arg 3:", args[3], "arg 4:", args[4], "arg 5:", args[5], sep = " "))
start.time <- Sys.time()
print(paste("Start time:", start.time, sep = ""))

# Point to file names
genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 0.01
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"
if (args[2] == "whole_blood_tpm"){
	expression_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_2017_Whole_blood_tpm.gct"
	cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
} else if (args[2] == "whole_blood_norm") {
	expression_file <- "/data/kryan/project/gtex/sample_data/Whole_blood_v8.normalized_filtered_expression.bed"
	cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
} else if (args[2] == "frac_rank") {
	#expression_file <- "/data/kryan/project/gtex/analysis/frac_rank.out.all"
	expression_file <- "/data/kryan/project/gtex/sample_data/frac_rank_with_extra_info_v2.txt"
	cov_file <- "/data/kryan/project/gtex/sample_data/covariates_frac_rank_extracted_exp_and_geno.txt"
}
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

# Read in covariates
if (args[2] == "frac_rank"){
	cov <- read.table(cov_file, check.names=F)
	#cov <- cov[-1,]
	cov_t <- t(cov)

} else if (args[2] == "whole_blood_norm" | args[2] == "whole_blood_tpm") {
	cov <- vroom(cov_file)
	cov_t <- t(cov[-1])
	colnames(cov_t) <- cov$ID
}

if (args[4] == "reduced_covariates"){
        reduced_covariates <- c("PC1", "PC2", "PC3", "PC4", "PC5", "pcr", "platform", "sex")
        cov_t <- cov_t[,reduced_covariates]
}

# convert to numeric - get some warnings
variants.genotypes <- genotype.data[,10:ncol(genotype.data)]
variants.numeric <- apply(variants.genotypes, 2, genotype_to_numeric)
variant.ids <- genotype.data$V3
# Read in expression data
expression_data <- vroom(expression_file)
# Convert sample ids to individual ids if required ("convert_ids" as argument 3)
if (args[3] == "convert_ids"){
	split.up <- str_split_fixed(colnames(expression_data)[3:ncol(expression_data)], "-",5)
	joined.up <- paste(split.up[,1],"-",split.up[,2], sep = "")
	first_2_cols <- colnames(expression_data[1:2])
	all.cols <- c(first_2_cols,joined.up)
	colnames(expression_data) <- all.cols
}

# remove rows with duplicated gene names - we don't know the ensembl ids for the frac_rank expression data so I wanted to do the same to the whole blood data so that we can compare like with like
if (args[2] == "whole_blood_tpm"){
	singletons <- names(which(table(expression_data$Description) == 1))
	expression_data_notduplicated <- expression_data[expression_data$Description %in% singletons, ]
} else if (args[2] == "frac_rank"){
	singletons <- names(which(table(expression_data$Gene) == 1))
	expression_data_notduplicated <- expression_data[expression_data$Gene %in% singletons, ]

} else if (args[2] == "whole_blood_norm"){
	singletons <- names(which(table(expression_data$ensembl_gene_id_target_short) == 1))
	expression_data_notduplicated <- expression_data[expression_data$ensembl_gene_id_target_short %in% singletons, ]
}
# Read in sample ids and match up with the sample ids in the expression data file
psam.ids <- read.table(psam)
colnames(variants.numeric) <- psam.ids$V1
if (args[2] == "whole_blood_tpm"){
	isect <- intersect(colnames(variants.numeric),colnames(expression_data[3,]))
} else if (args[2] == "frac_rank"){
	isect <- intersect(colnames(variants.numeric),colnames(expression_data[7,]))
} else if (args[2] == "whole_blood_norm"){
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
if (args[2] == "whole_blood_tpm"){
	rownames(expression.data.samples.matched) <- expression_data_notduplicated$Description
} else if ((args[2] == "frac_rank")){
	rownames(expression.data.samples.matched) <- expression_data_notduplicated$Gene
} else if (args[2] == "whole_blood_norm") {
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

# Function for extracting p-value from linear model - don't actually need this anymore
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
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

# Create granges object for expression data
if (args[2] == "whole_blood_norm"){
	df3_norm <- data.frame(chrom = expression_data_notduplicated$`#chr`, start = expression_data_notduplicated$start, end = expression_data_notduplicated$end)
	gr3 <- as(df3_norm, "GRanges")
} else if (args[2] == "frac_rank"){
	df3_frac_rank <- data.frame(chrom = expression_data_notduplicated$chromosome, start = expression_data_notduplicated$start, end = expression_data_notduplicated$end)
	df3_frac_rank$chrom <- paste("chr", df3_frac_rank$chrom, sep = "")
	gr3 <- as(df3_frac_rank, "GRanges")
}


outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
fisher.exact.pvalues <- c()
chisq.pvalues <- c()
no.of.targets <- c()
no.of.targets.after.filtering <- c()

# function to carry out eqtl analysis on non-targets
non_targets_eqtl <- function (tab, tab_outersect.rm.na, cov_t, p.values.outersect){
		linear.model.outersect <- lm(tab ~ tab_outersect.rm.na$tf.genotypes + as.matrix(cov_t), na.action = na.omit)
                p.values.outersect <- c(p.values.outersect, summary(linear.model.outersect)$coefficients[2,4])
		return(p.values.outersect)
		}

not_all_na <- function(x) any(!is.na(x))

if(args[1] == "default"){
for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	var <- tf.dbdp.uniq.indorothea[i,22]
	print(paste("Variant number", i, ":", var, sep = " "))
	tf.interest <- tf.dbdp.uniq.indorothea[i,9]	
	tf.interest.ensembl <- tf.dbdp.uniq.indorothea$Ensembl.Gene.ID[i]
	tf.genotypes <- variants.range.genotypes.matched[which(variant.ids.range == var),]
	# Remove targets within 1MB from consideration
	dorothea.possible <- dorothea[countOverlaps(gr1,gr2[i]) == 0,]
	# Remove all genes with 1MB from expression file, only looking for transeQTLs
	expression.data.samples.matched.possible <- expression.data.samples.matched[countOverlaps(gr3,gr2[i]) == 0,]
	rownames(expression.data.samples.matched.possible) <- rownames(expression.data.samples.matched)[which((countOverlaps(gr3,gr2[i]) == 0))]
	# Filter so that variants are not on the same chromosome within 1MB of the target gene
	if (args[2] == "whole_blood_norm"){
		tf.targets.interest <- dorothea.possible$ensembl_gene_id_target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
	} else if (args[2] == "frac_rank" | args[2] == "whole_blood_tpm"){
		tf.targets.interest <- dorothea.possible$target[dorothea.possible$tf == tf.interest & dorothea.possible$target != tf.interest]
	}
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
	df.targets.intersect <- df.targets[which(df.targets$tf.targets.interest %in% rownames(expression.data.samples.matched.possible)),]
	p.values <- c()
	p.values.outersect <- c()
	intersect.tf.targets.expression.data <- intersect(df.targets$tf.targets.interest, rownames(expression.data.samples.matched.possible))
	outersect.tf.targets.expression.data <- outersect(df.targets$tf.targets.interest, rownames(expression.data.samples.matched.possible))
	print(paste("outersect.tf.targets.expression data starting length:", length(outersect.tf.targets.expression.data)))
	expression_data_transposed <- t(expression.data.samples.matched.possible[intersect.tf.targets.expression.data,])
	variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	tab <- cbind.data.frame(tf.genotypes, expression_data_transposed, stringsAsFactors = FALSE)
	tab[,2:ncol(tab)] <- sapply(tab[,2:ncol(tab)], as.numeric)
	if (args[2] == "whole_blood_norm"){
                        expression.data.tf <- expression.data.samples.matched[tf.interest.ensembl,]
                } else if (args[2] == "frac_rank" | args[2] == "whole_blood_tpm") {
                        expression.data.tf <- expression.data.samples.matched[tf.interest,]
                }
                if (!not_all_na(expression.data.tf)){
                print(paste("No expression data for tf, skipping...", tf.interest))
		fisher.exact.pvalues <- c(fisher.exact.pvalues, NA)
		no.of.targets <- c(no.of.targets, NA)
	        no.of.targets.after.filtering <- c(no.of.targets.after.filtering, NA)
                next
                }
                expression.data.tf.vec <- unlist(expression.data.tf)
                for (j in 2:ncol(tab)){
                if (args[4] == "all_covariates" | args[4] == "reduced_covariates"){
                                linear.model <- lm(tab[,j] ~ expression.data.tf.vec + as.matrix(cov_t))
                        } else if (args[4] == "no_covariates") {
                                linear.model <- lm(tab[,j] ~ expression.data.tf.vec)
                        }

                linear.model.pval <- coef(summary(linear.model))["expression.data.tf.vec", "Pr(>|t|)"]
                sign.linear.model <- sign(coef(summary(linear.model))["expression.data.tf.vec", "Estimate"])
                adjusted.rsquared.linear.model <- summary(linear.model)$adj.r.squared
                sign.mor <- sign(mor.targets[j-1])

		if (linear.model.pval < 0.05 & sign.linear.model == sign.mor){
                        if (args[4] == "all_covariates" | args[4] == "reduced_covariates"){
                                linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes + as.matrix(cov_t))
                         } else if (args[4] == "no_covariates"){
                                linear.model.eqtl <- lm(tab[,j] ~ tab$tf.genotypes)
                         }
                         p.values <- c(p.values, summary(linear.model.eqtl)$coefficients[2,4])
                }  else  {
                        p.values <- c(p.values, NA)
			outersect.tf.targets.expression.data <- c(outersect.tf.targets.expression.data,intersect.tf.targets.expression.data[j-1])
                        }
                }
	print(paste("outersect.tf.targets.expression data finishing length:", length(outersect.tf.targets.expression.data)))
	expression_data_transposed_outersect <- t(expression.data.samples.matched.possible[outersect.tf.targets.expression.data,])
	tab_outersect <- cbind.data.frame(tf.genotypes, expression_data_transposed_outersect, stringsAsFactors = FALSE)
        tab_outersect[,2:ncol(tab_outersect)] <- sapply(tab_outersect[,2:ncol(tab_outersect)], as.numeric)
	not_all_na <- function(x) any(!is.na(x))
	tab_outersect_rm_na <- tab_outersect %>% select(where(not_all_na))
	for (l in 2:ncol(tab_outersect_rm_na)){
		if (args[4] == "all_covariates" | args[4] == "reduced_covariates"){
                                linear.model.outersect <- lm(tab_outersect_rm_na[,l] ~ tab_outersect_rm_na$tf.genotypes + as.matrix(cov_t), na.action = na.omit)
				p.values.outersect <- c(p.values.outersect, summary(linear.model.outersect)$coefficients[2,4])
                 } else if (args[4] == "no_covariates"){
                                linear.model.outersect <- lm(tab_outersect_rm_na[,l] ~ tab_outersect_rm_na$tf.genotypes, na.action = na.omit)
				p.values.outersect <- c(p.values.outersect, summary(linear.model.outersect)$coefficients[2,4])
                 }
	}
	p.values.remove.na <- na.omit(p.values)
	p.values.outersect.remove.na <- na.omit(p.values.outersect)
	eqtl_targets <- length(p.values.remove.na[p.values.remove.na < 0.05])
   	eqtl_nontargets <- length(p.values.outersect.remove.na[p.values.outersect.remove.na < 0.05])
   	noneqtl_targets <- length(p.values.remove.na[p.values.remove.na >= 0.05])
   	noneqtl_nontargets <- length(p.values.outersect.remove.na[p.values.outersect.remove.na >= 0.05])
	print(paste("eqtl_targets:", eqtl_targets))
	print(paste("eqtl_nontargets:", eqtl_nontargets))
	print(paste("noneqtl_targets:", noneqtl_targets))
	print(paste("noneqtl_nontargets:", noneqtl_nontargets))
   	contingency_tab <- array(c(eqtl_targets, eqtl_nontargets, noneqtl_targets, noneqtl_nontargets), c(2,2))
	fisher.exact.pval <- fisher.test(contingency_tab, alternative="greater")$p.value
	print(paste("fisher.exact.pval:", fisher.exact.pval))
	fisher.exact.pvalues <- c(fisher.exact.pvalues, fisher.exact.pval)
	no.of.targets <- c(no.of.targets,(ncol(tab)-1))
	no.of.targets.after.filtering <- c(no.of.targets.after.filtering, length(p.values.remove.na))
  }
  FDR.adjusted.pvals.fisher <- p.adjust(fisher.exact.pvalues, method = "BH")
  results <- cbind.data.frame(tf.dbdp.uniq.indorothea[,22],tf.dbdp.uniq.indorothea[,9], fisher.exact.pvalues, FDR.adjusted.pvals.fisher, no.of.targets, no.of.targets.after.filtering)
  colnames(results) <- c("Variant", "Variant_gene", "fisher.exact.pvalues", "FDR.adjusted.pvals.fisher", "Number_of_targets", "Number_of_targets_after_filtering")
  results.sorted <- results[order(results$FDR.adjusted.pvals.fisher),]
  print(results.sorted)
  write.table(results.sorted, file = paste("/data/kryan/project/gtex/outputs/trans_eqtl_results/", args[2], "/all_targets/", args[4], "/", args[5], "/enrichment_", args[2], "_", args[4], "_filtered_norsq_AF1.txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
 end.time <- Sys.time()
 print(paste("End time:", end.time, sep = ""))
 print(paste("time taken:", (end.time - start.time)))


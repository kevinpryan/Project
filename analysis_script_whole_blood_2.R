#!/usr/bin/Rscript

# Load in packages
library(stringr)
library(dorothea)
library(biomaRt)
library(vroom)
library(dplyr, warn.conflicts = FALSE)

# Take user input (average expression or use all expression values)
args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied ('average' or 'default' )", call.=FALSE)
} else if (args[1] != "average" && args[1] != "default") {
  stop("Argument must either be 'average' or 'default'", call. = FALSE)
}

# Point to file names
genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 5.0e-02
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"
expression_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_2017_Whole_blood_tpm.gct"
psam <- "/data/kryan/project/gtex/gtex.psam"

# Read in table and extract rows with allele frequency in the desired range
genotype.data <- read.table(genotypes, stringsAsFactors=FALSE)
allele.freq <- str_split_fixed(genotype.data$V8, ";" , n=14)[,2]
allele.frequency.figure <- as.numeric(str_split_fixed(allele.freq, "=" , n=2)[,2])
variants.range <- genotype.data[allele.frequency.figure >= allele.frequency.desired & allele.frequency.figure <= 1 - allele.frequency.desired,]

# Change genotype encoding from 0/0,0/1,1/1 to 0,1,2
genotype_to_numeric <- function(par){
 par.out <- par
 par.out[which(par == "0/0")] <- 0
 par.out[which(par == "0/1")] <- 1
 par.out[which(par == "1/1")] <- 2
return(as.numeric(par.out)) 
}

# convert to numeric - get some warnings, seems to be 196 na genotype values
variants.range.genotypes <- variants.range[,10:ncol(variants.range)]
variants.range.numeric <- apply(variants.range.genotypes, 2, genotype_to_numeric)

variant.ids <- variants.range$V3
print(paste("Number of variants within desired allele frequency:",nrow(variants.range)))

# Read in table from Barrera paper (with GTEx ids added in manually)
tf.dbdp <- read.table(barrera.variants)

# Only keep rows that are in GTEx at the desired allele frequency
tf.dbdp.variants <- tf.dbdp[tf.dbdp$gtex_var_format_tfdbdp_b38 %in% variant.ids,]

# Get rid of duplicated variants (the Barrera paper has separate entries for different transcript ids)
tf.dbdp.uniq <- tf.dbdp.variants[!duplicated(tf.dbdp.variants$gtex_var_format_tfdbdp_b38),]

# Only keep the variants in transcription factors that have data in the dorothea database
tf.dbdp.uniq.indorothea <- tf.dbdp.uniq[tf.dbdp.uniq$Gene.symbol %in% dorothea_hs$tf,]

print(paste("Number of transcription factors with GTEx variants in desired frequency range available in dorothea database:",nrow(tf.dbdp.uniq.indorothea)))

# Function for extracting p-value from linear model
lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

# Read in expression_data file - takes 10-15 minues
#expression_data <- read.table(expression_file, header = TRUE, fill = TRUE,nrows = 57000, comment.char = "", check.names=FALSE)
expression_data <- vroom(expression_file)
# Fix sample ids
split.up <- str_split_fixed(colnames(expression_data)[3:ncol(expression_data)], "-",5)
joined.up <- paste(split.up[,1],"-",split.up[,2], sep = "")
first_2_cols <- colnames(expression_data[1:2])
all.cols <- c(first_2_cols,joined.up)
colnames(expression_data) <- all.cols

# Delete any rows with duplicated gene names - we don't know the ensembl ids for the expression data so I wanted to do the same to the whole blood data so that we can compare like with like
singletons <- names(which(table(expression_data$Description) == 1))
expression_data_notduplicated <- expression_data[expression_data$Description %in% singletons, ]

# Read in sample ids and match up with the sample ids in the expression data file
psam.ids <- read.table(psam)
colnames(variants.range.numeric) <- psam.ids$V1
isect <- intersect(colnames(variants.range.numeric),colnames(expression_data[3,]))
variants.range.genotypes.matched <- variants.range.numeric[,isect]
expression.data.samples.matched <- expression_data_notduplicated[,isect]
rownames(expression.data.samples.matched) <- expression_data_notduplicated$Description

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
#df.out <- data.frame(variant = character(), tf.interest = character(), target = character(), p_value = as.numeric())
df.out <- data.frame()
for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	var <- tf.dbdp.uniq.indorothea[i,22]
	tf.interest <- tf.dbdp.uniq.indorothea[i,9]	
	tf.genotypes <- variants.range.genotypes.matched[which(variants.range$V3 == var),]
	tf.targets.interest <- dorothea_hs$target[which(dorothea_hs$tf == tf.interest & dorothea_hs$target != tf.interest)]
	intersect.tf.targets.expression.data <- intersect(tf.targets.interest, rownames(expression.data.samples.matched))
	expression_data_transposed <- t(expression.data.samples.matched[intersect.tf.targets.expression.data,])
	variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	tab <- cbind.data.frame(tf.genotypes, expression_data_transposed, stringsAsFactors = FALSE)
	tab[,2:ncol(tab)] <- sapply(tab[,2:ncol(tab)], as.numeric)
	p.values <- c()
	if(args[1] == "default"){
	p.values <- c()
	for (j in 2:ncol(tab)){
 		 linear.model <- lm(tab[,j] ~ tab$tf.genotypes)
 		 p.values <- c(p.values, lmp(linear.model))
		}
	results <- cbind(rep(var, length(intersect.tf.targets.expression.data)), rep(tf.interest, length(intersect.tf.targets.expression.data)),intersect.tf.targets.expression.data, p.values )
	colnames(results) <- c("Variant GTEx ID", "TF containing variant", "Target gene", "p-value")
	}
	else {
		tab_rank <- apply(tab[,2:ncol(tab)], FUN = rank, MARGIN = 2)
		tab_rank <- cbind(tab$tf.genotypes, tab_rank)
		tab_rank_mean <- apply(tab_rank[,2:ncol(tab_rank)], FUN = mean, MARGIN = 1)
		tab_rank_mean_genotypes <- cbind.data.frame(tab$tf.genotypes, tab_rank_mean)
		colnames(tab_rank_mean_genotypes)[1] <- "tf.genotypes" 
		linear.model <- lm(tab_rank_mean_genotypes$tab_rank_mean ~ tab_rank_mean_genotypes$tf.genotypes)
		p.values <- c(p.values, lmp(linear.model))
	        results <- cbind(var, tf.interest, p.values )
		colnames(results) <- c("Variant GTEx ID", "TF containing variant", "p-value")
	}
        df.out <- rbind(df.out,results)
}
# Print the bonferroni threshold
bonferroni.cutoff <- 0.05/no.of.tests
print(paste("Number of tests carried out", nrow(df.out)))
sig.values.df.out <- 0.05/nrow(df.out)
print(paste("Significance value:", sig.values.df.out))
colnames(df.out) <- colnames(results)
#jpeg("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_trans_eqtl_analysis_hist_pvals_nocovariates.jpeg")
#hist(as.numeric(df.out$"p-value"), xlab = "Unadjusted p-values", main = "eQTL analysis p-values with Whole Blood GTEx Data")
#dev.off()
print("results table:")
#print(df.out, row.names = FALSE)
df.out$"p-value" <- as.numeric(df.out$"p-value")
df.out.na.remove <- na.omit(df.out)
sig.results <- df.out.na.remove[df.out.na.remove$"p-value" <= sig.values.df.out,]
write.table(df.out,"/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/whole_blood_transeqtl_analysis_bonf05_no_covariates_10062021_summarised_expression_allresults_AF03.txt",quote = FALSE, sep = "\t", row.names = FALSE)
print("Significant results (Bonferroni)")
#print(sig.results, row.names = FALSE)
FDR.adjusted.pvals <- p.adjust(df.out.na.remove$"p-value", method = "BH") # adjusted using Benjamini and Hochberg method

FDR.adjusted.results <- cbind(df.out.na.remove, FDR.adjusted.pvals)
print("FDR adjusted results")
#print(FDR.adjusted.results, row.names = FALSE)
sig.results.bh <- df.out.na.remove[FDR.adjusted.pvals <= 0.05,]
sig.results.bh.added <- sig.results.bh %>%
                        mutate(FDR.adjusted.pvalues = FDR.adjusted.pvals[FDR.adjusted.pvals <= 0.05])
#print(sig.results.bh.added, row.names = FALSE)
write.table(FDR.adjusted.results, "/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/all_targets/whole_blood_trans_eqtl_analysis_no_covariates_allresults_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(sig.results, "/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/all_targets/whole_blood_trans_eqtl_analysis_no_covariates_sigresults_bonf_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(sig.results.bh.added, "/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/all_targets/whole_blood_trans_eqtl_analysis_no_covariates_sigresults_bh_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE)


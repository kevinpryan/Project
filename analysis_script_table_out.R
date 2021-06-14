#!/usr/bin/Rscript

args = commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("At least one argument must be supplied ('average' or 'default' )", call.=FALSE)
} else if (args[1] != "average" && args[1] != "default") {
  stop("Argument must either be 'average' or 'default'", call. = FALSE)
}

# Load in packages
library(stringr)
library(dorothea)
library(biomaRt)
library(vroom)
library(dplyr, warn.conflicts = FALSE)
# Point to file names
genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 5.0e-02
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"
expression_file <- "/data/kryan/project/gtex/analysis/frac_rank.out.all"
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

# Read in frac_rank file - takes 10-15 minues
expression_data <- vroom(expression_file)

# Delete any rows with duplicated gene names - we don't know the ensembl ids here
singletons <- names(which(table(expression_data$Gene) == 1))
expression_data_notduplicated <- expression_data[expression_data$Gene %in% singletons, ]

# Read in sample ids and match up with the sample ids in the fractional rank file
psam.ids <- read.table(psam)
colnames(variants.range.numeric) <- psam.ids$V1
isect <- intersect(colnames(variants.range.numeric),colnames(expression_data[2,]))
variants.range.genotypes.matched <- variants.range.numeric[,isect]
expression.data.samples.matched <- expression_data_notduplicated[,isect]
rownames(expression.data.samples.matched) <- expression_data_notduplicated$Gene

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
df.out <- data.frame(variant = character(), tf.interest = character(), target = character(), p_value = as.numeric())
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

	#for (j in 2:ncol(tab)){
 	#	 linear.model <- lm(tab[,j] ~ tab$tf.genotypes)
 	#	 p.values <- c(p.values, lmp(linear.model))
	#}
	#names(p.values) <- intersect.tf.targets.frac.rank
	#results <- cbind(rep(var, length(intersect.tf.targets.expression.data)), rep(tf.interest, length(intersect.tf.targets.expression.data)), intersect.tf.targets.expression.data, p.values)
	df.out <- rbind(df.out,results)
	#colnames(df.out) <- c("Variant GTEx ID", "Transcription factor containing variant", "Target gene", "p-value")
	#rownames(df.out) <- NULL
	#outputs[[i]] <- p.values
	#no.of.tests <- no.of.tests + length(p.values) 
}
#names(outputs) <- tf.dbdp.uniq.indorothea$gtex_var_format_tfdbdp_b38
print(paste("Number of tests carried out (no.of.tests):", nrow(df.out)))

# Print the bonferroni threshold
bonferroni.cutoff <- 0.05/nrow(df.out)
print(paste("Significance value (list method):", bonferroni.cutoff))

print(paste("Number of tests carried out(df.out)", nrow(df.out)))
print(head(df.out))
sig.values.df.out <- 0.05/nrow(df.out)
print(paste("Significance value (df method):", sig.values.df.out))
colnames(df.out) <- colnames(results)

#jpeg("/data/kryan/project/gtex/outputs/trans_eqtl_results/trans_eqtl_analysis_hist_pvals_nocovariates.jpeg")
#hist(as.numeric(df.out$"p-value"), xlab = "Unadjusted p-value", main = "Histogram of p-values from eQTL analysis using frac_rank data")
#dev.off()
print("results table:")
#print(df.out, row.names = FALSE)
df.out$"p-value" <- as.numeric(df.out$"p-value")
df.out.na.remove <- na.omit(df.out)
sig.results <- df.out.na.remove[df.out.na.remove$"p-value" <= sig.values.df.out,]
print("Sig results (bonferroni)")
print(sig.results, row.names = FALSE)
#write.table(df.out.na.remove,"/data/kryan/project/gtex/outputs/trans_eqtl_results/frac_rank_results/frac_rank_trans_eqtl_analysis_bonf05_no_covariates_10062021_summarised_expression_allresults_AF10.txt",quote = FALSE, sep = "\t", row.names = FALSE)

FDR.adjusted.pvals <- p.adjust(as.numeric(df.out$"p-value"), method = "BH") # adjusted using Benjamini and Hochberg method
FDR.adjusted.results <- cbind(df.out.na.remove, FDR.adjusted.pvals)
print("FDR adjusted results")
#print(FDR.adjusted.results, row.names = FALSE)

sig.results.bh <- df.out.na.remove[FDR.adjusted.pvals <= 0.05,]
sig.results.bh.added <- sig.results.bh %>%
			mutate(FDR.adjusted.pvalues = FDR.adjusted.pvals[FDR.adjusted.pvals <= 0.05])
print("FDR adjusted significicant results")
#print(sig.results.bh.added, row.names = FALSE)
write.table(FDR.adjusted.results, "/data/kryan/project/gtex/outputs/trans_eqtl_results/frac_rank_results/all_targets/frac_rank_trans_eqtl_analysis_no_covariates_allresults_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(sig.results, "/data/kryan/project/gtex/outputs/trans_eqtl_results/frac_rank_results/all_targets/frac_rank_trans_eqtl_analysis_no_covariates_sigresults_bonf_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE)
write.table(sig.results.bh.added, "/data/kryan/project/gtex/outputs/trans_eqtl_results/frac_rank_results/all_targets/frac_rank_trans_eqtl_analysis_no_covariates_sigresults_bh_AF05_11062021.txt", quote = FALSE, sep = "\t", row.names = FALSE)

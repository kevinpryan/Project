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
dorothea_file <- "/home/rstudio/Documents/MSc/MA5105/project_data/trans_eqtl_analysis/dorothea_hs_with_target_info.txt"
cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
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

# Read in dorothea table 
dorothea <- read.table(dorothea_file, header = T, sep = "\t", check.names = F)

# Create GRanges object from dorothea
df <- data.frame(chrom=dorothea$target_chromosome, start=dorothea$target_start, end=dorothea$target_end)
gr1 <- as(df, "GRanges")

# Dataframes for dorothea upregulated and downregulated targets
dorothea_upreg <- dorothea[dorothea$mor == 1,]
dorothea_downreg <- dorothea[dorothea$mor == -1,]

# Create GRanges object for filtered variants range
vars <- tf.dbdp.uniq.indorothea[,22]
vars.chromosome <- substr(vars,4,4)
var.positions <- as.numeric(str_split_fixed(vars,"_",5)[,2])
range.start <- as.numeric(max(1,var.positions - 1000000))
range.end <- as.numeric(var.positions + 1000000)
df2 <- data.frame(chrom = vars.chromosome,start=range.start, end=range.end)
gr2 <- as(df2, "GRanges")

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
#df.out <- data.frame(variant = character(), tf.interest = character(), target = character(), p_value = as.numeric())
df.out <- data.frame()
if(args[1] == "default"){
for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	var <- tf.dbdp.uniq.indorothea[i,22]
	#var.chromosome <- substr(var,4,4)
	#var.position <- as.numeric(str_split_fixed(var,"_",5)[,2])
	#range.start <- as.numeric(max(1,var.position - 1000000))
	#range.end <- as.numeric(var.position + 1000000)
	#range <- data.frame(start = range.start, end = range.end) 
	tf.interest <- tf.dbdp.uniq.indorothea[i,9]	
	tf.genotypes <- variants.range.genotypes.matched[which(variants.range$V3 == var),]
	# Remove targets within 1MB from consideration
	dorothea.possible <- dorothea[sum(countOverlaps(gr2[i],gr1)) == 0,]
	# Filter so that variants are not on the same chromosome within 1MB of the target gene
	tf.targets.interest <- dorothea.possible$target[dorothea.possible$tf == tf.interest & 
	                                                dorothea.possible$target != tf.interest]
	mor.targets <- dorothea$mor[which(dorothea.possible$target[dorothea.possible$tf == tf.interest & 
	                                                             dorothea.possible$target != tf.interest])]
	confidence.targets <- dorothea$confidence[which(dorothea.possible$target[dorothea.possible$tf == tf.interest & 
	                                                                         dorothea.possible$target != tf.interest])]
	#tf.targets.upreg <- tf.targets.interest[mor.targets == 1]
	#tf.targets.downreg <- tf.targets.interest[mor.targets == -1]
	#confidence.targets.upreg <- confidence.targets[mor.targets == 1]
	#confidence.targets.downreg <- confidence.targets[mor.targets == -1]
	#tf.targets.df <- cbind.data.frame(tf.targets.interest, mor.targets,confidence.targets)
	#intersect.tf.targets.expression.data <- intersect(tf.targets.interest, rownames(expression.data.samples.matched))
	#tf.targets.df.intersect <- tf.targets.df[tf.targets.df$tf.targets.interest == intersect.tf.targets.expression.data,]
	#expression_data_transposed <- t(expression.data.samples.matched[intersect.tf.targets.expression.data,])
	#variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	#tab <- cbind.data.frame(tf.genotypes, expression_data_transposed, stringsAsFactors = FALSE)
	#tab[,2:ncol(tab)] <- sapply(tab[,2:ncol(tab)], as.numeric)
	#if(args[1] == "default"){
	p.values <- c()
	cof <- c()
	intersect.tf.targets.expression.data <- intersect(tf.targets.interest, rownames(expression.data.samples.matched))
	#tf.targets.df.intersect <- tf.targets.df[tf.targets.df$tf.targets.interest == intersect.tf.targets.expression.data,]
	expression_data_transposed <- t(expression.data.samples.matched[intersect.tf.targets.expression.data,])
	variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	tab <- cbind.data.frame(tf.genotypes, expression_data_transposed, stringsAsFactors = FALSE)
	tab[,2:ncol(tab)] <- sapply(tab[,2:ncol(tab)], as.numeric)
	for (j in 2:ncol(tab)){
 		 linear.model <- lm(tab[,j] ~ tab$tf.genotypes)
 		 p.values <- c(p.values, lmp(linear.model))
 		 cof <- c(cof,summary(linear.model)$coefficients[2,1])
		}
	results <- cbind(rep(var, length(intersect.tf.targets.expression.data)), 
	                 rep(tf.interest, length(intersect.tf.targets.expression.data)),
	                 intersect.tf.targets.expression.data, p.values, cof)
	colnames(results) <- c("Variant_GTEx_ID", "TF_containing_variant", "Target_gene", "p_value", "coefficient")
  }
  df.out <- rbind(df.out,results,mor.targets,confidence.targets)
  print(paste("Number of tests carried out", nrow(df.out)))
  sig.values.df.out <- 0.05/nrow(df.out)
  print(paste("Significance value:", sig.values.df.out))
  colnames(df.out) <- c(colnames(results),"Mode_of_regulation","Target_confidence")
  df.out$"p_value" <- as.numeric(df.out$"p_value")
  df.out.na.remove <- na.omit(df.out)
  sig.results <- df.out.na.remove[df.out.na.remove$"p_value" <= sig.values.df.out,]
  FDR.adjusted.pvals <- p.adjust(df.out.na.remove$"p_value", method = "BH") # adjusted using Benjamini and Hochberg method
  FDR.adjusted.results <- cbind(df.out.na.remove, FDR.adjusted.pvals)
  
} else {
	  for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	    var <- tf.dbdp.uniq.indorothea[i,22]
	    tf.interest <- tf.dbdp.uniq.indorothea[i,9]	
	    tf.genotypes <- variants.range.genotypes.matched[which(variants.range$V3 == var),]
	    # Remove targets within 1MB from consideration
	    dorothea.possible <- dorothea[sum(countOverlaps(gr2[i],gr1)) == 0,]
	    # Filter so that variants are not on the same chromosome within 1MB of the target gene
	    tf.targets.interest <- dorothea.possible$target[dorothea.possible$tf == tf.interest & 
	                                                      dorothea.possible$target != tf.interest]
	    mor.targets <- dorothea$mor[which(dorothea.possible$target[dorothea.possible$tf == tf.interest & 
	                                                                 dorothea.possible$target != tf.interest])]
	    confidence.targets <- dorothea$confidence[which(dorothea.possible$target[dorothea.possible$tf == tf.interest & 
	                                                                               dorothea.possible$target != tf.interest])]
	    tf.targets.upreg <- tf.targets.interest[mor.targets == 1]
	    tf.targets.downreg <- tf.targets.interest[mor.targets == -1]
	    confidence.targets.upreg <- confidence.targets[mor.targets == 1]
	    confidence.targets.downreg <- confidence.targets[mor.targets == -1]
	    p.values.upreg <- c()
	    p.values.downreg <- c()
	    cof.upreg <- c()
	    cof.downreg <- c()
	    intersect.tf.targets.upreg.expression.data <- intersect(tf.targets.upreg, rownames(expression.data.samples.matched))
	    intersect.tf.targets.downreg.expression.data <- intersect(tf.targets.downreg, rownames(expression.data.samples.matched))
	    expression_data_transposed_upreg <- t(expression.data.samples.matched[intersect.tf.targets.upreg.expression.data,])
	    expression_data_transposed_downreg <- t(expression.data.samples.matched[intersect.tf.targets.downreg.expression.data,])
	    variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	    tab_upreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_upreg, stringsAsFactors = FALSE)
	    tab_downreg <- cbind.data.frame(tf.genotypes, expression_data_transposed_downreg, stringsAsFactors = FALSE)
		  tab_rank_upreg <- apply(tab_upreg[,2:ncol(tab_rank_upreg)], FUN = rank, MARGIN = 2)
		  tab_rank_downreg <- apply(tab_downreg[,2:ncol(tab_rank_downreg)], FUN = rank, MARGIN = 2)
		  tab_rank_upreg <- cbind(tab$tf.genotypes, tab_rank_upreg)
		  tab_rank_downreg <- cbind(tab$tf.genotypes, tab_rank_downreg)
		  tab_rank_upreg_mean <- apply(tab_rank_upreg[,2:ncol(tab_rank_upreg)], FUN = mean, MARGIN = 1)
		  tab_rank_downreg_mean <- apply(tab_rank_downreg[,2:ncol(tab_rank_downreg)], FUN = mean, MARGIN = 1)
		  tab_rank_upreg_mean_genotypes <- cbind.data.frame(tab_upreg$tf.genotypes, tab_rank_upreg_mean)
		  tab_rank_downreg_mean_genotypes <- cbind.data.frame(tab_downreg$tf.genotypes, tab_rank_downreg_mean)
		  colnames(tab_rank_upreg_mean_genotypes)[1] <- "tf.genotypes" 
		  colnames(tab_rank_downreg_mean_genotypes)[1] <- "tf.genotypes" 
		  #tab_rank_mean_genotypes <- cbind.data.frame(tab$tf.genotypes, tab_rank_mean)
		  #colnames(tab_rank_mean_genotypes)[1] <- "tf.genotypes" 
		  #linear.model <- lm(tab_rank_mean_genotypes$tab_rank_mean ~ tab_rank_mean_genotypes$tf.genotypes)
		  linear_model_upreg <- lm(tab_rank_upreg_mean_genotypes$tab_rank_upreg_mean ~ tab_rank_upreg_mean_genotypes$tf.genotypes)
		  linear_model_downreg <- lm(tab_rank_downreg_mean_genotypes$tab_rank_downreg_mean ~ 
		                               tab_rank_downreg_mean_genotypes$tf.genotypes)
		  p.values <- c(p.values, lmp(linear.model))
		  cof.upreg <- c(cof.upreg,summary(linear_model_upreg)$coefficients[2,1])
		  cof.downreg <- c(cof.downreg,summary(linear_model_downreg)$coefficients[2,1])
	    results <- cbind(var, tf.interest, p.values.upreg, p.values.upreg, cof.upreg, cof.downreg)
		  colnames(results) <- c("Variant_GTEx_ID", "TF_containing_variant", "p_value_upreg", 
		                       "p_value_downreg", "coefficient_upreg", "coefficient_downreg")
	}
        df.out <- rbind(df.out,results)
                        #,mor.targets,confidence.targets)
        # Print the bonferroni threshold
        no.of.tests <- 2*nrow(df.out)
        print(paste("Number of tests carried out", no.of.tests))
        sig.values.df.out <- 0.05/no.of.tests
        print(paste("Significance value:", bonferroni.cutoff))
        #colnames(df.out) <- c(colnames(results),"Mode_of_regulation","Target_confidence")
        #df.out$"p_value" <- as.numeric(df.out$"p_value")
        df.out.na.remove <- na.omit(df.out)
        sig.results <- df.out.na.remove[df.out.na.remove$"p_value" <= sig.values.df.out,]
        all.pvals <- c(df.out.na.remove$"p_value_upreg",df.out.na.remove$"p_value_downreg")
        FDR.adjusted.pvals <- p.adjust(all.pvals, method = "BH")
        FDR.adjusted.pvals.upreg <- FDR.adjusted.pvals[1:length(FDR.adjusted.pvals)/2]
        FDR.adjusted.pvals.downreg <- FDR.adjusted.pvals[length(FDR.adjusted.pvals)/2:length(FDR.adjusted.pvals)]
        FDR.adjusted.results <- cbind(df.out.na.remove, FDR.adjusted.pvals.upreg, FDR.adjusted.pvals.downreg)
        print("FDR adjusted results")
        print(FDR.adjusted.results, row.names = FALSE)
}


# Print the bonferroni threshold
#bonferroni.cutoff <- 0.05/no.of.tests
#print(paste("Number of tests carried out", nrow(df.out)))
#sig.values.df.out <- 0.05/nrow(df.out)
#print(paste("Significance value:", sig.values.df.out))
#colnames(df.out) <- c(colnames(results),"Mode_of_regulation","Target_confidence")
#jpeg("/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_trans_eqtl_analysis_hist_pvals_nocovariates.jpeg")
#hist(as.numeric(df.out$"p-value"), xlab = "Unadjusted p-values", main = "eQTL analysis p-values with Whole Blood GTEx Data")
#dev.off()
print("results table:")
#print(df.out, row.names = FALSE)
#df.out$"p-value" <- as.numeric(df.out$"p_value")
#df.out.na.remove <- na.omit(df.out)
#sig.results <- df.out.na.remove[df.out.na.remove$"p_value" <= sig.values.df.out,]
#write.table(df.out,"/data/kryan/project/gtex/outputs/trans_eqtl_results/whole_blood_results/whole_blood_transeqtl_analysis_bonf05_no_covariates_10062021_summarised_expression_allresults_AF03.txt",quote = FALSE, sep = "\t", row.names = FALSE)
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


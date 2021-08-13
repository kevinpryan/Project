#!/usr/bin/Rscript
# script to create scatterplot coloured by genotype. Script contains lots of superfluous code as an interaction analysis script was used as a template
library(vroom)
library(ggplot2)
library(stringr)
library(dplyr, warn.conflicts = FALSE)
library(GenomicRanges, quietly = T)
library(ggpubr)
library(rlist)
library(gridExtra)
library(grid)
args = commandArgs(trailingOnly = TRUE)
variants_interest <- args[4]
print(paste("variants interest:", variants_interest))
targets_interest <- c("ZNF14", "ZSCAN12", "ZNF461",  "MTERF1",  "ZNF879",  "ZNF285" )

print(paste("targets_interest:", targets_interest))

if (length(args) == 0) {
  stop("At least one argument must be supplied ('average' or 'default' )", call.=FALSE)
} else if (args[1] != "average" && args[1] != "default") {
  stop("Argument must either be 'average' or 'default'", call. = FALSE)
}

genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 1.0e-02
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"

if (args[2] == "whole_blood_tpm"){
        expression_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_2017_Whole_blood_tpm.gct"
        cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
} else if (args[2] == "whole_blood_norm") {
        expression_file <- "/data/kryan/project/gtex/sample_data/Whole_blood_v8.normalized_filtered_expression.bed"
        cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
} else if (args[2] == "frac_rank") {
        expression_file <- "/data/kryan/project/gtex/analysis/frac_rank.out.all"
        cov_file <- "/data/kryan/project/gtex/sample_data/covariates_frac_rank_extracted_exp_and_geno_rownames_altered.txt"
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
        isect <- intersect(colnames(variants.numeric),colnames(expression_data[2,]))
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

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
#df.out <- data.frame(variant = character(), tf.interest = character(), target = character(), p_value = as.numeric())
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

#df.out <- data.frame(variant = character(), tf.interest = character(), target = character(), p_value = as.numeric())

not_all_na <- function(x) any(!is.na(x))
#plots <- list()
myplots <- vector('list', length(targets_interest))
i <- which(tf.dbdp.uniq.indorothea$gtex_var_format_tfdbdp_b38 == variants_interest)
var <- tf.dbdp.uniq.indorothea[i,22]
tf.interest <- tf.dbdp.uniq.indorothea[i,9]
print(paste("tf.interest:", tf.interest))
tf.interest.ensembl <- tf.dbdp.uniq.indorothea$Ensembl.Gene.ID[i]
tf.genotypes <- variants.range.genotypes.matched[which(variant.ids.range == var),]
# Remove targets within 1MB from consideration
dorothea.possible <- dorothea[countOverlaps(gr1,gr2[i]) == 0,]
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
print("tf.targets.interest")
print(tf.targets.interest)
print("tf.targets.names")
print(tf.targets.names)
print("chrom.target.range")
print(chrom.target.range)
print("mor.targets")
print(mor.targets)
print("confidence.targets")
print(confidence.targets)
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
if (args[2] == "whole_blood_norm"){
                        expression.data.tf <- expression.data.samples.matched[tf.interest.ensembl,]
                } else if (args[2] == "frac_rank" | args[2] == "whole_blood_tpm") {
                        expression.data.tf <- expression.data.samples.matched[tf.interest,]
                }
expression.data.tf.vec <- unlist(expression.data.tf)

for (j in 1:length(targets_interest)){
	myplots[[j]] <- local({
	j <- j
	df.for.plotting <- data.frame(genotype = factor(tab$tf.genotypes), expression_target_gene = tab[,targets_interest[j]], expression_transcription_factor = expression.data.tf.vec)
	title <- paste(targets_interest[j])
	plot_object <- ggplot(df.for.plotting, aes(x=expression_transcription_factor, y=expression_target_gene, color=genotype, shape=genotype)) +
  	geom_point(size = 0.5) + 
  	geom_smooth(method=lm, se=FALSE, fullrange=TRUE, size = 0.5)+
  	ggtitle(title) +
  	theme_bw() +
  	theme(plot.title = element_text(hjust = 0.5, size = rel(1.5)), 
	      panel.grid.major = element_blank(), 
	      panel.grid.minor = element_blank(), 
	      plot.margin = unit(c(1,1,1,3), "lines"),
	      legend.key.size = unit(3, 'cm'),
	      legend.title = element_text(size=rel(1.5)),
	      legend.position = "right") +
	labs(x = NULL, y = NULL)
	print(plot_object)
	})
}

print(paste("length plots:", length(myplots)))
filename <- paste("scatterplot_expression", var[1], sep = "_")
end <- ".pdf"
filename_complete <- paste(filename, end, sep = "")
print(paste("Filename:", filename_complete))
ggar_obj <- ggarrange(plotlist = myplots, nrow = 2, ncol = 3, common.legend = TRUE)#, legend = "bottom")
ggar_obj_annotated <- annotate_figure(ggar_obj, left = textGrob("Target gene expression (mean fractional rank)", rot = 90, vjust = 1, gp = gpar(cex = 1.6)), bottom = textGrob("ZNF502 expression (mean fractional rank) ", gp = gpar(cex = 1.6)))
ggsave(filename_complete, ggar_obj_annotated, width=15, height=8.5)




#!/usr/bin/Rscript

# Version of analysis script but with allele frequency of 1%
# Load in packages
library(stringr)
library(dorothea)
library(biomaRt)
# Point to file names
genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 1.0e-02
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"
frac_rank_file <- "/data/kryan/project/gtex/analysis/frac_rank.out.all"
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
# Read in Barrera table with gtex variant ids (build 38)

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
frac_rank <- read.table(frac_rank_file, header = TRUE, fill = TRUE,nrows = 55914, comment.char = "", check.names=FALSE)

# Delete any rows with duplicated gene names - we don't know the ensembl ids here
singletons <- names(which(table(frac_rank$Gene) == 1))
frac_rank_notduplicated <- frac_rank[frac_rank$Gene %in% singletons, ]

# Read in sample ids and match up with the sample ids in the fractional rank file
psam.ids <- read.table(psam)
colnames(variants.range.numeric) <- psam.ids$V1
isect <- intersect(colnames(variants.range.numeric),colnames(frac_rank[2,]))
variants.range.genotypes.matched <- variants.range.numeric[,isect]
frac.rank.samples.matched <- frac_rank_notduplicated[,isect]
rownames(frac.rank.samples.matched) <- frac_rank_notduplicated$Gene

# For each transcription factor variant, get the targets for the transcription factor and fit a linear model for the expression of the targets wrt increasing numbers of the minor allele
outputs <- list()
no.of.tests <- 0
for (i in 1:nrow(tf.dbdp.uniq.indorothea)){
	var <- tf.dbdp.uniq.indorothea[i,22]
	tf.interest <- tf.dbdp.uniq.indorothea[i,9]	
	tf.genotypes <- variants.range.genotypes.matched[which(variants.range$V3 == var),]
	tf.targets.interest <- dorothea_hs$target[which(dorothea_hs$tf == tf.interest & dorothea_hs$target != tf.interest)]
	intersect.tf.targets.frac.rank <- intersect(tf.targets.interest, rownames(frac.rank.samples.matched))
	frac_rank_transposed <- t(frac.rank.samples.matched[intersect.tf.targets.frac.rank,])
	variants.range.genotypes.matched.vector <- unlist(variants.range.genotypes.matched)
	tab <- cbind.data.frame(tf.genotypes, frac_rank_transposed, stringsAsFactors = FALSE)
	#rownames(tab) <- isect
	tab[,2:ncol(tab)] <- sapply(tab[,2:ncol(tab)], as.numeric)
	p.values <- c()
	for (j in 2:ncol(tab)){
 		 linear.model <- lm(tab[,j] ~ tab$tf.genotypes)
 		 p.values <- c(p.values, lmp(linear.model))
	}
	names(p.values) <- intersect.tf.targets.frac.rank
	outputs[[i]] <- p.values
	no.of.tests <- no.of.tests + length(p.values) 
}
names(outputs) <- tf.dbdp.uniq.indorothea$gtex_var_format_tfdbdp_b38
print(paste("Number of tests carried out:", no.of.tests))

# Print the bonferroni threshold
bonferroni.cutoff <- 0.05/no.of.tests
#sink("trans_eqtl_output.txt")
#print(outputs)
#sink()
#dump(outputs, file = "trans_eqtl_output_dump_05.txt")
# dput(outputs, "trans_eqtl_outputs_dput_05.txt")
#dput(outputs, "test_remove_same_gene_qtl.txt")
dput(outputs, file = "/data/kryan/project/gtex/outputs/2805/test_outputs_1percent_1036.txt")
# read in dput outputs and look for significant ones @ bonferroni threshold

for (i in 1:length(outputs)){
	print(names(outputs[i]))
	print(outputs[[i]][which(outputs[[i]] <= bonferroni.cutoff)])
}


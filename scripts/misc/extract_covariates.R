#!/usr/bin/Rscript
# Purpose of this script is to extract covariates of interest from the GTEx samples to add to the fractional rank data. We want to extract the first 5 genotyping PCs,
# PCR conditions (pcr free or not), sex, sequencing platform for the GTEx samples in the frac_rank data. To do this, we start with the 
# whole blood covariates and extract the relevant covariates. Then we loop through each GTEx tissue and extract the covariates for any sample ID in the frac_rank
# data which has not yet been encountered. To finish we add the covariates to the expression PCs which have already been determined using PCAtools
# 
library(stringr)
library(biomaRt)
library(vroom)
library(dplyr, warn.conflicts = FALSE)
library(GenomicRanges, quietly = T)
genotypes <- "/data/kryan/project/gtex/genotype_data_dbdp.vcf"
allele.frequency.desired <- 5.0e-02
barrera.variants <- "/data/kryan/project/gtex/analysis/dbdp_nssnp_barrera_with_GTEx_varid.txt"
#expression_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_2017_Whole_blood_tpm.gct"
#expression_file <- "/data/kryan/project/gtex/sample_data/Whole_blood_v8.normalized_filtered_expression.bed"
expression_file <- "/data/kryan/project/gtex/analysis/frac_rank.out.all"
psam <- "/data/kryan/project/gtex/gtex.psam"
dorothea_file <- "/data/kryan/project/gtex/analysis/dorothea_hs_with_target_info.txt"
cov_file <- "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/Whole_Blood.v8.covariates.txt"
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
cov <- read.table(cov_file, header = T, row.names = 1, check.names = F)

# convert to numeric - get some warnings
variants.genotypes <- genotype.data[,10:ncol(genotype.data)]
variants.numeric <- apply(variants.genotypes, 2, genotype_to_numeric)
variant.ids <- genotype.data$V3
# Read in expression data
expression_data <- vroom(expression_file)
singletons <- names(which(table(expression_data$Gene) == 1))
expression_data_notduplicated <- expression_data[expression_data$Gene %in% singletons, ]
psam.ids <- read.table(psam)
colnames(variants.numeric) <- psam.ids$V1
isect <- intersect(colnames(variants.numeric),colnames(expression_data[2,]))
variants.genotypes.matched <- variants.numeric[,isect]
expression_file_wb <- "/data/kryan/project/gtex/sample_data/Whole_blood_v8.normalized_filtered_expression.bed"
expression_data_wb <- vroom(expression_file_wb)
singletons_wb <- names(which(table(expression_data_wb$ensembl_gene_id_target_short) == 1))
expression_data_notduplicated_wb <- expression_data_wb[expression_data_wb$ensembl_gene_id_target_short %in% singletons_wb, ]
isect_wb <- intersect(colnames(variants.numeric),colnames(expression_data_wb[6,]))
isect.x <- intersect(isect_wb, isect)
cov_frac_rank <- cov[,isect.x]
x <- c(1:5)
y <- c((nrow(cov)-2):nrow(cov))
cov_frac_rank <- cov_frac_rank[c(x,y),]
outersect <- function(x, y) {
  sort(c(setdiff(x, y),
         setdiff(y, x)))
}
out <- outersect(isect, isect_wb)
f <- list.files(path = "/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/", pattern = "^[^W]*.txt")
for (i in 1:length(f)){
#for (i in c(1,2)) {
	if (length(out) != 0){
		print(paste("Trying file:", f[i]))
		tab <- read.table(paste("/data/kryan/project/gtex/sample_data/GTEx_Analysis_v8_eQTL_covariates/",f[i], sep = ""), header = T, row.names = 1, check.names = F)
		isect.x <- c(isect.x,intersect(out, colnames(tab)))
		new.samples <- intersect(out, colnames(tab))
		out <- outersect(out, new.samples)
		if (length(new.samples) > 1){
			cov_frac_rank <- cbind(cov_frac_rank, tab[c(x,y),new.samples])
		} else if (length(new.samples) == 1) {
			cov_frac_rank <- cbind(cov_frac_rank, tab[c(x,y),new.samples])
			colnames(cov_frac_rank)[ncol(cov_frac_rank)] <- new.samples
		} 
		print(paste("new samples:",new.samples)) 
		print(paste("Length of out:", length(out)))	
	} else {
		break
	}
}
cov_frac_rank_out <- cov_frac_rank[,isect]
print("writing to file")
#write.table(cov_frac_rank_out, file = "/data/kryan/project/gtex/sample_data/covariates_frac_rank_extracted7.txt", quote = F, sep = "\t")
# read in expression principal components
pcs <- read.table("/data/kryan/project/gtex/sample_data/principal_components_frac_rank_pcatools.txt", check.names=F)
pcs_t <- t(pcs)
pcs_t_ordered <- pcs_t[,isect]
pcs_t_ordered_reduced <- pcs_t_ordered[1:60,]
rownames(cov_frac_rank_out)[1:5] <- c("PC1_geno", "PC2_geno", "PC3_geno", "PC4_geno", "PC5_geno")
#colnames(cov_frac_rank_out)[1] <- "ID"
#colnames(pcs_t_ordered_reduced)[1] <- "ID"
covariates_frac_rank <- rbind.data.frame(cov_frac_rank_out, pcs_t_ordered_reduced)
write.table(covariates_frac_rank, file = "/data/kryan/project/gtex/sample_data/covariates_frac_rank_extracted_exp_and_geno.txt", quote = F, sep = "\t")


	





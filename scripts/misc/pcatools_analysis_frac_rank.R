#!/usr/bin/Rscript
# script to run PCA on the frac_rank data. Also finds elbow point, carries out Horn's parallel analysis (find number of covariates which explain x% of the variance), and produces a scree plot.
library(vroom)
library(PCAtools)
frac_rank_file <- "/data/kryan/project/gtex/analysis/frac_rank.out.all"
frac_rank <- vroom(frac_rank_file)
frac_rank <- frac_rank[,2:(ncol(frac_rank)-1)]
frac_rank_transpose <- t(frac_rank)
remove_zero_var <- frac_rank_transpose[ , which(apply(frac_rank_transpose, 2, var) != 0)]
p <- pca(na.omit(remove_zero_var))
elbow <- findElbowPoint(p$variance)
# takes a while to run
horn <- parallelPCA(na.omit(remove_zero_var))
which(cumsum(p$variance) > 80)[1] # 303
which(cumsum(p$variance) > 70)[1] # 137
which(cumsum(p$variance) > 60)[1] # 49           
which(cumsum(p$variance) > 62)[1] # 60
pdf("scree_advanced_wb.pdf")
screeplot(p.wb,
    components = getComponents(p.wb, 1:20),
    vline = c(horn.wb$n, elbow.wb)) +
    geom_label(aes(x = horn.wb$n + 1, y = 50,
    label = 'Horn\'s', vjust = -1, size = 8)) +
    geom_label(aes(x = elbow.wb + 1, y = 50,
    label = 'Elbow method', vjust = -1, size = 8))
dev.off()
pcs <- getLoadings(p)
write.table(pcs, file = "/data/kryan/project/gtex/sample_data/principal_components_frac_rank_pcatools.txt", col.names = T, row.names = T, quote = F, sep = "\t")
		   
           


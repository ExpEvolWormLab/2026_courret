############################################################
## permGWAS - preparation & analyses 
## C. elegans – CeMEE panel
############################################################


## ===============================
## PREPARATION
## ===============================

library(AGHmatrix)
library(data.table)
library(vcfR)
library(ggplot2)
library(dplyr)

# Load the phenotypes
pheno_var <- fread("data/gemma_pheno_withLines_var.txt")
pheno_mean <- fread("data/gemma_pheno_withLines.txt")
pheno_CV <- fread("data/gemma_pheno_withLines_CV.txt")

str(pheno_var)

#Modif fichier pheno
pheno_var <- pheno_var[, .(line = V2, phenotype_value= V1)]
pheno_mean <- pheno_mean[, .(line = V2, phenotype_value= V1)]
pheno_CV <- pheno_CV[, .(line = V2, phenotype_value= V1)]


#Modif fichier geno
#geno_raw <- fread("data/genotype_raw_GT.txt")

geno <- fread("data/genotype_012.txt")

sum(is.na(geno))

str(geno)
length(unique(geno$CHROM))
unique(geno$CHROM)

geno[, rs := paste(CHROM, POS, sep = "_")]   # crée la colonne rs
setcolorder(geno, c("rs", setdiff(names(geno), c("CHROM","POS","rs"))))  # place rs en premier
geno[, c("CHROM","POS") := NULL]  # supprime CHROM et POS

geno_matrix <- as.matrix(geno[, -1])     # tout sauf la colonne 'rs'
rownames(geno_matrix) <- geno$rs         # mettre les rs comme noms de lignes

# Transposer
geno_t <- t(geno_matrix)

str(geno_t)

sum(is.na(geno_t))

geno_t_df <- data.frame(rs = rownames(geno_t), geno_t, check.names = FALSE)

str(geno_t_df)

sum(is.na(geno_t))

#Kinship


Kinship <- as.matrix(read.table("plink_pruning/pruning01/CeMEE.pruned.0.1.rel"))
fam <- read.table("plink_pruning/pruning01/CeMEE.pruned.0.1.fam")
ids <- fam$V2  # IID
ids <- sub("_.*", "", ids)

rownames(Kinship) <- ids
colnames(Kinship) <- ids

str(Kinship)

# modif des pop label (en commun)
ids_pheno <- unique(pheno_mean$line)
ids_geno  <- rownames(geno_t_df) 
ids_kin   <- rownames(Kinship)

str(geno_t_df)

common_ids <- Reduce(intersect, list(ids_pheno, ids_geno, ids_kin)) 

pheno_mean_perm <- pheno_mean[line %in% common_ids, ]
geno_perm <- geno_t_df[rownames(geno_t_df) %in% common_ids, ]
kinship_perm <- Kinship[common_ids, common_ids]

str(kinship_perm)
str(geno_perm)

#RETIRER NA DE GENO
sum(is.na(geno_perm))
geno_clean <- geno_perm
geno_clean <- geno_clean[, colSums(is.na(geno_clean)) == 0] # enlever SNPs avec au moins un NA
sum(is.na(geno_clean))  # doit être 0
dim(geno_clean)

#En csv

fwrite(pheno_mean, "data/pheno_mean_permGWAS.csv", sep = ",")
fwrite(pheno_var,  "data/pheno_var_permGWAS.csv",  sep = ",")
fwrite(pheno_CV,   "data/pheno_CV_permGWAS.csv",   sep = ",")


write.table(geno_clean,
            file = "data/geno_permGWAS.csv",
            sep = ",",
            quote = FALSE,
            row.names = F,
            col.names = T)

geno_perm= as.data.frame(geno_perm)
str(geno_perm)


write.table(geno_perm,
            file = "data/geno_grid.txt",
            sep = "\t",        
            quote = FALSE,
            row.names = TRUE,  
            col.names = TRUE)


kinship_dt <- as.data.table(kinship_perm, keep.rownames = "line")

fwrite(kinship_dt, "data/kinship_permGWAS.csv", sep = ",")


pheno_mean <- fread("data/pheno_mean_permGWAS.csv", sep=",", header=TRUE, encoding="UTF-8")
fwrite(pheno_mean, "data/pheno_mean_permGWAS.csv", sep=",")


## ===============================
## ANALYSES
## ===============================


#MEAN

permGWAS_mean <- fread("output_permGWAS/p_values_phenotype_value_mean.csv")
str(permGWAS_mean)

permGWAS_mean <- as.data.table(permGWAS_mean)

res <- permGWAS_mean

# ordre correct des chromosomes (important pour C. elegans)
chr_levels <- c("I", "II", "III", "IV", "V", "X", "MtDNA")

# transformer en facteur ordonné
res$CHR <- factor(res$CHR, levels = chr_levels)

# trier
setorder(res, CHR, POS)

res[, BP_cum := POS]
offset <- 0

for(chr in levels(res$CHR)){
  idx <- which(res$CHR == chr)
  res[idx, BP_cum := POS + offset]
  offset <- max(res$BP_cum[idx])
}

axis_df <- res[, .(
  center = (min(BP_cum) + max(BP_cum)) / 2
), by = CHR]

res$log10p <- -log10(res$p_value)

ggplot(res, aes(x = BP_cum, y = log10p, color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  # couleurs chromosomes
  scale_color_manual(values = rep(c("grey70", "grey40"), length(chr_levels))) +
  
  # -------------------------
# SEUILS (avec légende)
# -------------------------

geom_hline(aes(yintercept = -log10(1.581517750955237e-07),
               linetype = "Bonferroni 5%"),
           color = "red", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(3.163035501910473e-08),
                 linetype = "Bonferroni 1%"),
             color = "darkred", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(2.7600882732707975e-06),
                 linetype = "Permutation 5%"),
             color = "blue", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(8.807478080315957e-07),
                 linetype = "Permutation 1%"),
             color = "darkblue", linewidth = 0.8) +
  
  scale_linetype_manual(
    name = "Significance thresholds",
    values = c(
      "Bonferroni 5%" = "dashed",
      "Bonferroni 1%" = "dotted",
      "Permutation 5%" = "twodash",
      "Permutation 1%" = "solid"
    )
  ) +
  
  theme_classic() +
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](p)),
    color = "Chromosome"
  ) +
  
  theme(
    legend.position = "right"
  )



#VARIANCE

permGWAS_var <- fread("output_permGWAS/p_values_phenotype_value_var.csv")
str(permGWAS_var)

library(ggplot2)
library(dplyr)
library(data.table)

permGWAS_var <- as.data.table(permGWAS_var)

res <- permGWAS_var

# ordre correct des chromosomes (important pour C. elegans)
chr_levels <- c("I", "II", "III", "IV", "V", "X", "MtDNA")

# transformer en facteur ordonné
res$CHR <- factor(res$CHR, levels = chr_levels)

# trier
setorder(res, CHR, POS)

res[, BP_cum := POS]
offset <- 0

for(chr in levels(res$CHR)){
  idx <- which(res$CHR == chr)
  res[idx, BP_cum := POS + offset]
  offset <- max(res$BP_cum[idx])
}

axis_df <- res[, .(
  center = (min(BP_cum) + max(BP_cum)) / 2
), by = CHR]

res$log10p <- -log10(res$p_value)

ggplot(res, aes(x = BP_cum, y = log10p, color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  # couleurs chromosomes
  scale_color_manual(values = rep(c("grey70", "grey40"), length(chr_levels))) +
  
  # -------------------------
# SEUILS (avec légende)
# -------------------------

geom_hline(aes(yintercept = -log10(1.581517750955237e-07),
               linetype = "Bonferroni 5%"),
           color = "red", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(3.163035501910473e-08),
                 linetype = "Bonferroni 1%"),
             color = "darkred", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(9.503197388682205e-08),
                 linetype = "Permutation 5%"),
             color = "blue", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(2.006834383080637e-08),
                 linetype = "Permutation 1%"),
             color = "darkblue", linewidth = 0.8) +
  
  scale_linetype_manual(
    name = "Significance thresholds",
    values = c(
      "Bonferroni 5%" = "dashed",
      "Bonferroni 1%" = "dotted",
      "Permutation 5%" = "twodash",
      "Permutation 1%" = "solid"
    )
  ) +
  
  theme_classic() +
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](p)),
    color = "Chromosome"
  ) +
  
  theme(
    legend.position = "right"
  )


#CV

permGWAS_CV <- fread("output_permGWAS/p_values_phenotype_value_CV.csv")
str(permGWAS_CV)

library(ggplot2)
library(dplyr)
library(data.table)

permGWAS_mean <- as.data.table(permGWAS_CV)

res <- permGWAS_CV

# ordre correct des chromosomes (important pour C. elegans)
chr_levels <- c("I", "II", "III", "IV", "V", "X", "MtDNA")

# transformer en facteur ordonné
res$CHR <- factor(res$CHR, levels = chr_levels)

# trier
setorder(res, CHR, POS)

res[, BP_cum := POS]
offset <- 0

for(chr in levels(res$CHR)){
  idx <- which(res$CHR == chr)
  res[idx, BP_cum := POS + offset]
  offset <- max(res$BP_cum[idx])
}

axis_df <- res[, .(
  center = (min(BP_cum) + max(BP_cum)) / 2
), by = CHR]

res$log10p <- -log10(res$p_value)

ggplot(res, aes(x = BP_cum, y = log10p, color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  # couleurs chromosomes
  scale_color_manual(values = rep(c("grey70", "grey40"), length(chr_levels))) +
  
  # -------------------------
# SEUILS (avec légende)
# -------------------------

geom_hline(aes(yintercept = -log10(1.581517750955237e-07),
               linetype = "Bonferroni 5%"),
           color = "red", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(3.163035501910473e-08),
                 linetype = "Bonferroni 1%"),
             color = "darkred", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(3.9945021064056265e-07),
                 linetype = "Permutation 5%"),
             color = "blue", linewidth = 0.8) +
  
  geom_hline(aes(yintercept = -log10(9.758057779994083e-08),
                 linetype = "Permutation 1%"),
             color = "darkblue", linewidth = 0.8) +
  
  scale_linetype_manual(
    name = "Significance thresholds",
    values = c(
      "Bonferroni 5%" = "dashed",
      "Bonferroni 1%" = "dotted",
      "Permutation 5%" = "twodash",
      "Permutation 1%" = "solid"
    )
  ) +
  
  theme_classic() +
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](p)),
    color = "Chromosome"
  ) +
  
  theme(
    legend.position = "right"
  )


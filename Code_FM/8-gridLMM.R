############################################################
## DATA PREPARATION – GridLMM
## C. elegans – CeMEE panel
############################################################


############################################################
## STEP 1 — PHENOTYPE CLEAN (avec correction isotypes)
############################################################

rm(list = ls())
library(data.table)

## ===============================
## 1. LOAD DATA
## ===============================

phen_data <- fread("data/phen_data.txt")
corr_isotypes <- fread("data/correspondance_isotypes.csv")
removed_isotypes <- fread("data/removed_isotypes_Cemee_2025.txt", header=F)

## ===============================
## 5. CORRECTION ISOTYPES
## ===============================

# clean noms
phen_data[, pop_label := gsub("CeMee", "", pop_label)]
corr_isotypes[, Line := gsub("CeMee", "", Line)]
corr_isotypes[, Representant_line := gsub("CeMee", "", Representant_line)]
removed_isotypes[, V1 := gsub("CeMee", "", V1)]

# merge correspondance
phen_data <- merge(
  phen_data,
  corr_isotypes[, .(Line, Representant_line)],
  by.x = "pop_label",
  by.y = "Line",
  all.x = TRUE
)

# remplacer par isotype représentatif
phen_data[, pop_label := ifelse(
  is.na(Representant_line),
  pop_label,
  Representant_line
)]

phen_data[, Representant_line := NULL]

# enlever lignes à exclure
phen_data <- phen_data[!pop_label %in% removed_isotypes$V1]

## ===============================
## 6. FORMAT FINAL
## ===============================

phen_clean <- phen_data[, .(
  pop_label,
  area.F,
  area.F.var,
  logD,
  rel_humidity,
  temperature,
  date_str,
  session_id, 
  env_label
)]

############################################################
## STEP 2 — EXTRACT GENOTYPE
############################################################

#Sur ordi ENS, transfo en 0,1 ou 2
geno <- fread("data/genotype_012.txt")

geno[, SNP_ID := paste(CHROM, POS, sep = "_")] # créer ID SNP
geno_num <- geno[, !c("CHROM", "POS", "SNP_ID"), with = FALSE] # garder uniquement SNP numériques
geno_mat <- as.matrix(geno_num) # convertir en matrice
geno_mat <- t(geno_mat) # transpose
colnames(geno_mat) <- geno$SNP_ID # ajouter noms SNP
rownames(geno_mat) <- colnames(geno_num) # ajouter noms individus


############################################################
## STEP 3 — KINSHIP (VCF + PLINK + vcfR)
############################################################

library(vcfR)
library(data.table)
library(AGHmatrix)

## ===============================
## 1. Creer matrice de genotype pruné 
## ===============================

# créer dossier
system("mkdir -p plink_pruning/pruning_grid")
#pruning de 0.1
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 0.1 --out plink_pruning/pruning_grid/CeMEE.pruned.0.1 --allow-extra-chr --double-id")
# extraire matrice de genotype pruné
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --extract plink_pruning/pruning_grid/CeMEE.pruned.0.1.prune.in --recode A --out plink_pruning/pruning_grid/genotype_pruned --allow-extra-chr --double-id")

geno_pruned <- fread("plink_pruning/pruning_grid/genotype_pruned.raw")
str(geno_pruned)

############################################################
## STEP 4 — Bonne taille de fichier
############################################################

# modif des pop label (en commun)
ids_pheno <- unique(phen_clean$pop_label)
ids_geno  <- rownames(geno_mat)
ids_pruned <- unique(geno_pruned$FID)

common_ids <- Reduce(intersect, list(ids_pheno, ids_geno, ids_pruned))

phen_clean_sub <- phen_clean[phen_clean$pop_label %in% common_ids, ]
geno_snp_sub <- geno_mat[common_ids, ]
pruned_sub <- geno_pruned[geno_pruned$FID %in% common_ids, ]

# Pour le phenotype
fwrite(phen_clean_sub, "data/phenotype_grid.txt", sep = "\t")

# Pour le genotype
write.table(geno_mat,
            "data/geno_grid.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

# Pour la kinship
write.table(pruned_sub,
            file = "data/pruned_grid.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)



#exemple Francois 
gwas.NGM = GridLMM_GWAS(formula = meanTrait ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                        test_formula = ~1, data = all_BLUEs_genotyped_NGM, X = genotype_matrix,
                        X_ID = 'pop_label', relmat = list(pop_label = kinship_matrix) ,method = 'REML')


############################################################
## STEP 1 — PHENOTYPE CLEAN (avec correction isotypes)
############################################################

rm(list = ls())

devtools::install_github("deruncie/GridLMM")
library(GridLMM)
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

## ===============================
## 7. SAVE PHENOTYPE
## ===============================

fwrite(phen_clean, "data/phenotype_grid.txt", sep = "\t")



############################################################
## STEP 2 — EXTRACTION GENOTYPE DEPUIS VCF (BRUT)
############################################################

#Sur ordi ENS, transfo en 0,1 ou 2
geno <- fread("data/genotype_012.txt")

# créer ID SNP
geno[, SNP_ID := paste(CHROM, POS, sep = "_")]

# garder uniquement SNP numériques
geno_num <- geno[, !c("CHROM", "POS", "SNP_ID"), with = FALSE]

# convertir en matrice
geno_mat <- as.matrix(geno_num)

# transpose (IMPORTANT : sur la matrice, pas sur geno)
geno_mat <- t(geno_mat)

# ajouter noms SNP (colonnes)
colnames(geno_mat) <- geno$SNP_ID

# ajouter noms individus (lignes)
rownames(geno_mat) <- colnames(geno_num)


############################################################
## STEP 3 — KINSHIP (VCF + PLINK + vcfR)
############################################################

library(vcfR)
library(data.table)
library(AGHmatrix)

## ===============================
## 1. PLINK LD pruning (IDENTIQUE)
## ===============================

# créer dossier
system("mkdir -p plink_pruning/pruning_grid")
#pruning de 0.1
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 0.1 --out plink_pruning/pruning_grid/CeMEE.pruned.0.1 --allow-extra-chr --double-id")
# extraire matrice de genotype pruné
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --extract plink_pruning/pruning_grid/CeMEE.pruned.0.1.prune.in --recode A --out plink_pruning/pruning_grid/genotype_pruned --allow-extra-chr --double-id")

geno_pruned <- fread("plink_pruning/pruning_grid/genotype_pruned.raw")
str(geno_pruned)
# enlever colonnes meta PLINK
geno_pruned <- geno_pruned[, -(1:6)]
genotype_pruned_matrix <- as.matrix(geno_pruned)

str(genotype_pruned_matrix)

colnames(genotype_pruned_matrix) <- sapply(
  strsplit(colnames(genotype_pruned_matrix), "_"),
  function(x) paste(x[1], x[2], sep = "_")
)

## ===============================
## 3. Bonne taille de fichier 
## ===============================

# modif des pop label (en commun)
ids_pheno <- unique(phen_clean$pop_label)
ids_geno  <- rownames(geno_mat)
ids_pruned <- unique(geno_pruned$FID)

common_ids <- Reduce(intersect, list(ids_pheno, ids_geno))

phen_clean_sub <- phen_clean[phen_clean$pop_label %in% common_ids, ]
geno_snp_sub <- geno_mat[common_ids, ]
pruned_sub <- geno_pruned[geno_pruned$FID %in% common_ids, ]

str(geno_snp_sub)

fwrite(phen_clean_sub, "data/phenotype_grid.txt", sep = "\t")

# Pour le genotype
write.table(geno_mat,
            "data/geno_grid.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)

str(geno_snp_sub)
str(phen_clean_sub)
str(pruned_sub)

# Pour la kinship
write.table(pruned_sub,
            file = "data/pruned_grid.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


# creer matrice de genotype pruné 

# créer dossier
system("mkdir -p plink_pruning/pruning_grid")
#pruning de 0.1
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 0.1 --out plink_pruning/pruning_grid/CeMEE.pruned.0.1 --allow-extra-chr --double-id")
# extraire matrice de genotype pruné
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --extract plink_pruning/pruning_grid/CeMEE.pruned.0.1.prune.in --recode A --out plink_pruning/pruning_grid/genotype_pruned --allow-extra-chr --double-id")

geno_pruned <- fread("plink_pruning/pruning_grid/genotype_pruned.raw")
# enlever colonnes meta PLINK
geno_pruned <- geno_pruned[, -(1:6)]
genotype_pruned_matrix <- as.matrix(geno_pruned)

str(genotype_pruned_matrix)

colnames(genotype_pruned_matrix) <- sapply(
  strsplit(colnames(genotype_pruned_matrix), "_"),
  function(x) paste(x[1], x[2], sep = "_")
)






















## ===============================
## 0. Extraire IDs communs
## ===============================

ids_pheno  <- unique(phen_clean_sub$pop_label)
ids_geno   <- rownames(geno_snp_sub)
ids_kin    <- pruned_sub$FID

common_ids <- Reduce(intersect, list(ids_pheno, ids_geno, ids_kin))

cat("Common IDs:", length(common_ids), "\n")


## ===============================
## 1. Reordonner PHENOTYPE
## ===============================

phen_clean_sub <- phen_clean_sub[
  match(common_ids, phen_clean_sub$pop_label)
]

phen_clean_sub <- phen_clean_sub[
  phen_clean_sub$pop_label %in% common_ids
]


## ===============================
## 2. Reordonner GENOTYPE
## ===============================

geno_snp_sub <- geno_snp_sub[common_ids, , drop = FALSE]


## ===============================
## 3. Reordonner KINSHIP
## ===============================

pruned_sub <- pruned_sub[
  match(common_ids, pruned_sub$FID)
]


## ===============================
## 4. CHECK FINAL (OBLIGATOIRE)
## ===============================

stopifnot(all(phen_clean_sub$pop_label == rownames(geno_snp_sub)))
stopifnot(all(pruned_sub$FID == rownames(geno_snp_sub)))
stopifnot(nrow(phen_clean_sub) == nrow(geno_snp_sub))


## ===============================
## 5. EXPORT PROPRE
## ===============================

fwrite(phen_clean_sub,
       "data/phenotype_grid.txt",
       sep = "\t")

write.table(geno_snp_sub,
            file = "data/geno_grid.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)

write.table(pruned_sub,
            file = "data/pruned_grid.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)


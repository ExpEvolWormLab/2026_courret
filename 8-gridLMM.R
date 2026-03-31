install.packages("gridLMM")  
library(GridLMM)
library(data.table)
devtools::install_github('deruncie/GridLMM')


# Installer la version GitHub de GridLMM
devtools::install_github("deruncie/GridLMM", force = TRUE)


#exemple Francois 
gwas.NGM = GridLMM_GWAS(formula = meanTrait ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                        test_formula = ~1, data = all_BLUEs_genotyped_NGM, X = genotype_matrix,
                        X_ID = 'pop_label', relmat = list(pop_label = kinship_matrix) ,method = 'REML')




############################################################
## STEP 1 — PHENOTYPE CLEAN (avec correction isotypes)
############################################################

rm(list = ls())

library(data.table)

## ===============================
## 1. LOAD DATA
## ===============================

phen_data <- fread("data/all_phenotypes_from_movies_fixed.txt")
corr_isotypes <- fread("data/correspondance_isotypes.csv")
removed_isotypes <- fread("data/removed_isotypes_Cemee_2025.txt", header=F)

## ===============================
## 2. VARIABLES EXPÉRIMENTALES
## ===============================

# session ID
phen_data[, session_id := substr(date_str, 1, 4)]
phen_data[session_id == "2014", session_id := "2013"]
phen_data[substr(block_id, 1, 1) == "F", session_id := "F"]
phen_data[substr(block_id, 1, 2) == "SB", session_id := "SB"]

# covariables numériques
for(i in c("temperature","rel_humidity","logD")){
  phen_data[, i := as.numeric(get(i))]
}

## ===============================
## 3. FILTRAGE
## ===============================

phen_data <- phen_data[!is.na(pop_label)]
phen_data <- phen_data[!block_id %in% paste0("B", 200:400)]
phen_data <- phen_data[!pop_label %in% c("A00", "CX12311L0", "LSJ2L0")]

phen_data <- phen_data[
  !is.na(temperature) &
    !is.na(rel_humidity) &
    !is.na(logD)
]

## ===============================
## 4. PHÉNOTYPES D'INTÉRÊT
## ===============================

phen_data <- phen_data[
  !is.na(area.F) &
    !is.na(area.F.var)
]

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
  session_id
)]

## ===============================
## 7. SAVE PHENOTYPE
## ===============================

fwrite(phen_clean, "output/phenotype_clean.txt", sep = "\t")

str(phen_clean)

cat("✅ Phénotype clean sauvegardé\n")
cat("Nb lignes:", nrow(phen_clean), "\n")
cat("Nb individus:", length(unique(phen_clean$pop_label)), "\n")



############################################################
## STEP 2 — EXTRACTION GENOTYPE DEPUIS VCF (BRUT)
############################################################

library(vcfR)

## ===============================
## 1. LOAD VCF
## ===============================

vcf <- read.vcfR("data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz")

## ===============================
## 2. EXTRAIRE GT (FORMAT BRUT)
## ===============================

gt <- extract.gt(vcf, element = "GT")
# garde "0/0", "0/1", etc.

## ===============================
## 3. FORMAT TABLEAU
## ===============================

geno_dt <- as.data.table(gt)

# ajouter info SNP
geno_dt[, CHROM := getCHROM(vcf)]
geno_dt[, POS := getPOS(vcf)]

# mettre CHROM/POS au début
setcolorder(geno_dt, c("CHROM", "POS", setdiff(colnames(geno_dt), c("CHROM","POS"))))

## ===============================
## 4. SAVE GENOTYPE
## ===============================

fwrite(geno_dt, "output/genotype_raw_GT.txt", sep = "\t")

cat("✅ Génotype brut sauvegardé\n")
cat("Nb SNPs:", nrow(geno_dt), "\n")
cat("Nb individus:", ncol(geno_dt) - 2, "\n")

dim(geno)


geno <- fread("data/genotype_012.txt")
str(geno)


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
system("mkdir -p plink_pruning/pruning01")
#pruning de 0.1
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 0.1 --out plink_pruning/pruning01/CeMEE.pruned.0.1 --allow-extra-chr --double-id")
# Extraire SNPs prunés
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --extract plink_pruning/pruning01/CeMEE.pruned.0.1.prune.in --recode vcf --out plink_pruning/pruning01/CeMEE.pruned.0.1 --allow-extra-chr --double-id")
# compression
system("bgzip plink_pruning/pruning01/CeMEE.pruned.0.1.vcf")

## ===============================
## 2. Lire VCF pruné
## ===============================

genotype_pruned <- read.vcfR("plink_pruning/pruning01/CeMEE.pruned.0.1.vcf.gz")

## ===============================
## 3. Convertir en matrice
## ===============================

genotype_pruned_matrix <- apply(genotype_pruned@gt[, -1], 1, function(x){
  as.numeric(substring(x, 1, 1))
})

# convertir en dosage 0/1/2
genotype_pruned_matrix <- 2 * genotype_pruned_matrix

# mettre individus en lignes (IMPORTANT)
genotype_pruned_matrix <- t(genotype_pruned_matrix)

## ===============================
## 4. IDs individus
## ===============================

strainID <- tstrsplit(colnames(genotype_pruned@gt)[-1], '_')[[1]]

genotype_pruned_matrix <- t(genotype_pruned_matrix)

rownames(genotype_pruned_matrix) <- strainID

## ===============================
## 5. Charger phénotype
## ===============================

ids_pheno <- unique(phen_clean$pop_label)

## ===============================
## 6. Subset commun
## ===============================

common_ids <- intersect(strainID, ids_pheno)

geno_sub <- genotype_pruned_matrix[common_ids, , drop = FALSE]

## ===============================
## 7. Reorder (CRUCIAL 🔥)
## ===============================

ids_order <- unique(phen_clean$pop_label)
ids_order <- ids_order[ids_order %in% common_ids]

geno_sub <- geno_sub[ids_order, ]

## ===============================
## 8. Kinship (VanRaden)
## ===============================

GRM <- Gmatrix(geno_sub, method = "VanRaden")

## ===============================
## 9. Vérification
## ===============================

stopifnot(all(rownames(GRM) == ids_order))

## ===============================
## 10. Export
## ===============================

write.table(GRM,
            file = "data/kinship_grid.txt",
            sep = '\t',
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)







# Supprimer CHROM et POS
geno_snp <- geno[, !c("CHROM", "POS"), with = FALSE]


gwas.NGM <- GridLMM_GWAS(
  formula = area.F ~ 1 + logD + rel_humidity + temperature + (1|pop_label) + (1|session_id),
  reduced_formula = area.F ~ 1 + logD + rel_humidity + temperature + (1|pop_label) + (1|session_id),
  test_formula = ~1,
  data = phen_clean,
  X = geno_snp,
  X_ID = 'pop_label',
  relmat = list(pop_label = GRM_sub),
  method = 'REML',
  verbose = TRUE 
)





# IDs communs
common_ids <- intersect(ids_pheno, ids_kinship)

# Filtrer le phénotype
phen_clean_sub <- phen_clean[phen_clean$pop_label %in% common_ids, ]

# Réordonner GRM_sub pour que l'ordre corresponde à phen_clean_sub
GRM_sub <- GRM[common_ids, common_ids]

# Si ton X (geno_snp) doit être aligné sur pop_label :
geno_snp_sub <- geno_snp[, common_ids ]  # ou cols selon l'orientation de ton tableau

# Vérification
all(phen_clean_sub$pop_label %in% rownames(GRM_sub)) # doit être TRUE

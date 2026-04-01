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
  session_id
)]

## ===============================
## 7. SAVE PHENOTYPE
## ===============================

fwrite(phen_clean, "data/phenotype_grid.txt", sep = "\t")

str(phen_clean)



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

dim(geno)

#Sur ordi ENS, transfo en 0,1 ou 2
geno <- fread("data/genotype_012.txt")
str(geno)
geno_snp <- geno[, -(1:2)]
geno_snp <- t(geno_snp) #fichier final

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
## 3. 
## ===============================

system("/opt/homebrew/bin/plink --vcf plink_pruning/pruning01/CeMEE.pruned.0.1.vcf.gz --make-bed --out plink_pruning/pruning01/CeMEE.pruned.0.1 --allow-extra-chr --double-id")

system("/opt/homebrew/bin/plink --bfile plink_pruning/pruning01/CeMEE.pruned.0.1 --make-rel square --out plink_pruning/pruning01/CeMEE.pruned.0.1 --allow-extra-chr") 


GRM <- as.matrix(read.table("plink_pruning/pruning01/CeMEE.pruned.0.1.rel"))

fam <- read.table("plink_pruning/pruning01/CeMEE.pruned.0.1.fam")

ids <- fam$V2  # IID

ids <- sub("_.*", "", ids)

rownames(GRM) <- ids
colnames(GRM) <- ids

str(GRM)

ids_pheno <- unique(phen_clean$pop_label)

common_ids <- intersect(ids, ids_pheno) # 

phen_clean_sub <- phen_clean[phen_clean$pop_label %in% common_ids, ]
GRM_sub <- GRM[common_ids, common_ids]


stopifnot(all(rownames(GRM_sub) == colnames(GRM_sub)))

relmat = list(pop_label = GRM_sub)


# modif des pop label (en commun)
ids_pheno <- unique(phen_clean$pop_label)
ids_geno  <- rownames(geno_snp)
ids_kin   <- rownames(GRM_sub)

common_ids <- Reduce(intersect, list(ids_pheno, ids_geno, ids_kin))

phen_clean_sub <- phen_clean[phen_clean$pop_label %in% common_ids, ]
geno_snp_sub <- geno_snp[common_ids, ]
GRM_sub2 <- GRM_sub[common_ids, common_ids]



fwrite(phen_clean_sub, "data/phenotype_grid.txt", sep = "\t")

# Pour le genotype
write.table(geno_snp_sub,
            file = "data/geno_grid.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)     # met un NA dans la première case du header

# Pour la kinship
write.table(GRM_sub2,
            file = "data/kinship_grid.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)


str(phen_clean_sub)


######## R ENS #########

#Modele MEAN
gwas.NGM <- GridLMM_GWAS(
  formula = area.F ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
  reduced_formula = area.F ~ 1 + logD + rel_humidity + temperature,
  test_formula = ~1,
  data = phen_clean_sub,
  X = geno_snp_sub,
  X_ID = 'pop_label',
  relmat = list(pop_label = GRM_sub2),
  method = 'REML'
)
#Modele VAR
gwas.NGM <- GridLMM_GWAS(
  formula = area.F.var ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
  reduced_formula = area.F ~ 1 + logD + rel_humidity + temperature,
  test_formula = ~1,
  data = phen_clean_sub,
  X = geno_snp_sub,
  X_ID = 'pop_label',
  relmat = list(pop_label = GRM_sub2),
  method = 'REML'
)
#calcul du CV
#Modele CV
gwas.NGM <- GridLMM_GWAS(
  formula = area.F ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
  reduced_formula = area.F ~ 1 + logD + rel_humidity + temperature,
  test_formula = ~1,
  data = phen_clean_sub,
  X = geno_snp_sub,
  X_ID = 'pop_label',
  relmat = list(pop_label = GRM_sub2),
  method = 'REML'
)

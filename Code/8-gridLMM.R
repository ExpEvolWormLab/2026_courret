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

plot(phen_clean$area.F, phen_clean$area.F.var)

ggplot(phen_clean, aes(x = area.F, y = area.F.var)) +
  
  geom_point(
    alpha = 0.6,
    size = 2,
    color = "steelblue"
  ) +
  
  scale_y_log10() +
  
  theme_classic() +
  
  labs(
    x = "Mean",
    y = "log(variance)"
  )
cor.test(phen_clean$area.F, phen_clean$area.F.var)

############################################################
## STEP 2 — EXTRACT GENOTYPE
############################################################

# 1. lecture VCF
vcf <- read.vcfR(
  "data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz"
)

# 2. extraire génotypes
gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

# 3. fonction de conversion robuste
gt_to_012 <- function(x) {
  ifelse(x %in% c("0/0", "0|0"), 0,
         ifelse(x %in% c("0/1", "1/0", "0|1", "1|0"), 1,
                ifelse(x %in% c("1/1", "1|1"), 2, NA)))
}

# 4. appliquer conversion
geno_012 <- apply(gt, 2, gt_to_012)

# 5. reconstruire data.frame
geno_012 <- data.frame(
  CHROM = vcf@fix[, "CHROM"],
  POS   = as.numeric(vcf@fix[, "POS"]),
  geno_012,
  check.names = FALSE
)

str(geno_012)

saveRDS(geno_012, file = "data/geno_012_clean.rds")
geno <- readRDS("data/geno_012_clean.rds")

#6. mise en forme 
geno$SNP_ID <- paste(geno$CHROM, geno$POS, sep = "_")
geno_num <- geno[, !(names(geno) %in% c("CHROM", "POS", "SNP_ID"))]
geno_mat <- as.matrix(geno_num)
geno_mat <- t(geno_mat)
colnames(geno_mat) <- geno$SNP_ID
rownames(geno_mat) <- colnames(geno_num)


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

common_ids <- Reduce(intersect, list(ids_pheno, ids_geno))

phen_clean_sub <- phen_clean[phen_clean$pop_label %in% common_ids, ]
geno_snp_sub <- geno_mat[common_ids, ]
pruned_sub <- geno_pruned[geno_pruned$FID %in% common_ids, ]

# Pour le phenotype
fwrite(phen_clean_sub, "data/phenotype_grid.txt", sep = "\t")

# Pour le genotype
write.table(geno_snp_sub,
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



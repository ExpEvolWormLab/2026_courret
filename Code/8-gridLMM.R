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

#gt <- extract.gt(vcf)
#sum(is.na(gt))

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
str(geno_snp)
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

install.packages("devtools")
library("devtools")
install_github("jdstorey/qvalue")
library("qvalue")


#MEAN

gwas.mean <- readRDS("output_gridLMM/gwas_mean.rds")

str(gwas.mean)

hist(pvals, breaks = 50, col = "grey80",
     main = "Distribution des p-values",
     xlab = "p-values") #hyp validé

pvals <- gwas.mean$results$p_value_REML

valid <- !is.na(pvals)

qobj <- qvalue(p = pvals[valid])

gwas.mean$results$qvalue <- NA
gwas.mean$results$qvalue[valid] <- qobj$qvalues

gwas.mean$results$significant <- gwas.mean$results$qvalue < 0.05

summary(qobj)
plot(qobj)

max(qobj$qvalues[pvals[valid] <= 0.01])






#comparer beta 

gemma_mean= fread("output_GEMMA/gwas_results.assoc.txt")
str(gemma_mean)

# Compare the beta
plot(gwas.mean$results$beta.5,gemma_mean$beta,pch=10)
#Compare the pvalues
plot(-log10(gwas.mean$results$p_value_REML),-log10(gemma_mean$p_wald),pch=10)




#VARIANCE 












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
















#code francois 
gwas.NGM = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                        test_formula = ~1, reduced_formula = ~0,
                        data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                        X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

test.glm <- lmer(y~1+logD + rel_humidity + temperature + (1|date_str)  + (1|session_id),data=all_BLUEs_genotyped_NGM)
all_BLUEs_genotyped_NGM$y_corrected_residuals <- resid(test.glm)

gwas.NGM.2 = GridLMM_GWAS(formula = y_corrected_residuals ~ 1 + (1|pop_label),
                          test_formula = ~1, reduced_formula = ~0,
                          data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                          X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gwas.NGM.3 = GridLMM_GWAS(formula = y_corrected_residuals ~ 1 + (1|pop_label),
                          test_formula = ~1, reduced_formula = ~1,
                          data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                          X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gwas.NGM.4 = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                          test_formula = ~ 1 + logD + rel_humidity + temperature,
                          reduced_formula = ~0,
                          data = all_BLUEs_genotyped_NGM, X = genotype_matrix[,is.NA==0],
                          X_ID = 'pop_label', relmat = list(pop_label = GRM_sub), mc.cores = 2 ,method = 'REML')

gemma_T7=(read.table("gemma_R02/gwas_results_BLUE_R02_T7.assoc.txt",h=TRUE))

# Compare the beta
plot(gwas.NGM$results$beta.5,gemma_T7$beta,pch=16)
#Compare the pvalues
plot(-log10(gwas.NGM$results$p_value_REML),-log10(gemma_T7$p_wald),pch=16)

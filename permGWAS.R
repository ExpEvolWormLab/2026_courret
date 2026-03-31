## ===============================
## 0. Packages
## ===============================

library(AGHmatrix)
library(data.table)
library(vcfR)

# Load the phenotypes
pheno_mean <- read.table("data/pheno_mean.txt", header = T)
pheno_var <- read.table("data/pheno_var.txt", header = T)
pheno_CV <- read.table("data/pheno_CV.txt", header = T)

# Load the isotype information
removed_isotypes <- read.table("data/removed_isotypes_Cemee_2025.txt")
corr_isotypes <- read.csv("data/correspondance_isotypes.csv")

## ===============================
## 1. Isogenic lines 
## ===============================

### Remove the additional name info
removed_isotypes$V1 <- gsub('CeMee', '', removed_isotypes$V1)
removed_isotypes$V1 <- gsub('_sorted', '', removed_isotypes$V1)
removed_isotypes=unique(removed_isotypes)
corr_isotypes$Line <- gsub('CeMee', '', corr_isotypes$Line)
corr_isotypes$Line <- gsub('_sorted', '', corr_isotypes$Line)
corr_isotypes$Representant_line <- gsub('CeMee', '', corr_isotypes$Representant_line)
corr_isotypes$Representant_line <- gsub('_sorted', '', corr_isotypes$Representant_line)
corr_isotypes <- corr_isotypes[,c(1,3)]
corr_isotypes=subset(corr_isotypes,Line!=Representant_line)
table(removed_isotypes$V1%in%corr_isotypes$Representant_line)
replaced_iso <- subset(corr_isotypes,Line%in%removed_isotypes$V1)

## ===============================
## 1.1 Cemee data 
## ===============================

vcf <- read.vcfR("data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz")
vcf@fix
snp_positions <- as.data.frame(getFIX(vcf))

snp_positions=as.data.frame(snp_positions)
snp_positions$CHROM=factor(snp_positions$CHROM,levels=c("I","II","III","IV","V","X","MtDNA"))
snp_positions$POS=as.numeric(as.character(snp_positions$POS))

## Order the genotypes
genotype_info <- extract.gt(vcf, element = "GT")
genotype_info <- genotype_info[order(snp_positions$CHROM,snp_positions$POS),]
snp_positions <- snp_positions[order(snp_positions$CHROM,snp_positions$POS),]

## Order the genotype - keep shared lines
genotyped_lines <- colnames(genotype_info)
all_lines <- sort(unique(pheno_CV$pop_label)) # pour la moyenne / var / CV
phenotyped_lines = rep(TRUE,length(genotyped_lines))
gemma_pheno <- NULL
### Now we will create the phenotype file for GEMMA
k=0
for(i in genotyped_lines){
  k=k+1
  ## Is there multiple isotypes ?
  temp_iso <- subset(replaced_iso,i==Line | i==Representant_line)
  all_lines_ID <- unique(c(i,temp_iso$Line,temp_iso$Representant_line))
  # Phenotyped once ?
  if(sum(pheno_CV$pop_label%in%all_lines_ID)==0){ # pour la moyenne / var / CV
    phenotyped_lines[k]=FALSE
  }else{
    temp_pheno <- subset(pheno_CV,pop_label%in%all_lines_ID) # pour la moyenne / var / CV
    ## Now extract all phenotypic data for these lines
    if(nrow(temp_pheno)==0){
      gemma_pheno=rbind(gemma_pheno,rep(NA,1))
    }else{
      gemma_pheno <- rbind(gemma_pheno,mean(temp_pheno$emmean,na.rm=TRUE))
    }
  }
}

table(phenotyped_lines) # 966 TRUE
dim(gemma_pheno) # 966 OK 

genotype_info_kept <- genotype_info[,phenotyped_lines]


#preparation fichier permGWAS

# transformer gemma_pheno en vecteur propre
pheno_vec <- as.numeric(gemma_pheno)


## ===============================
## 4. TRANSPOSE (clé pour permGWAS)
## ===============================

geno_t <- t(geno_num)

## ===============================
## 5. Ajouter noms propres
## ===============================

# noms SNPs = CHROM_POS
colnames(geno_t) <- paste0(snp_positions$CHROM, "_", snp_positions$POS)

# noms individus
rownames(geno_t) <- genotyped_lines_kept

## ===============================
## 6. Construire dataframe final
## ===============================

# Génotype
geno_df <- as.data.frame(geno_t)
geno_df$id <- rownames(geno_df)
geno_df <- geno_df[, c("id", setdiff(colnames(geno_df), "id"))]

# Phénotype
pheno_df <- data.frame(
  id = genotyped_lines_kept,
  phenotype = pheno_vec
)

## ===============================
## 7. Vérification CRUCIALE
## ===============================

if (!all(pheno_df$id == geno_df$id)) {
  stop("❌ ERREUR : mismatch géno / phéno")
}

## ===============================
## 8. Export
## ===============================

write.csv(geno_df, "data/genotype_permGWAS.csv", row.names = FALSE)
write.csv(pheno_df, "data/phenotype_permGWAS_CV.csv", row.names = FALSE)






#################

## ===============================
## 1. Charger les données
## ===============================

# Génotypes déjà convertis (0/1/2)
geno_raw <- read.csv("data/genotype_for_permGWAS.csv", header = FALSE)

str(geno_raw)

# Phénotype
pheno_vec <- as.numeric(read.table("data/gemma_pheno_log_var.txt")[,1])

## ===============================
## 2. TRANSPOSE (clé pour permGWAS)
## ===============================

geno_t <- t(geno_raw)

## ===============================
## 3. Ajouter noms propres
## ===============================

# ⚠️ À adapter si besoin selon tes objets déjà existants
colnames(geno_t) <- paste0(snp_positions$CHROM, "_", snp_positions$POS)
rownames(geno_t) <- genotyped_lines_kept

## ===============================
## 4. Construire dataframe final
## ===============================

# Génotype
geno_df <- as.data.frame(geno_t)
geno_df$id <- rownames(geno_df)
geno_df <- geno_df[, c("id", setdiff(colnames(geno_df), "id"))]

# Phénotype
pheno_df <- data.frame(
  id = genotyped_lines_kept,
  phenotype = pheno_vec
)

## ===============================
## 5. Vérification CRUCIALE
## ===============================

if (!all(pheno_df$id == geno_df$id)) {
  stop("❌ ERREUR : mismatch géno / phéno")
}

## ===============================
## 6. Export
## ===============================

write.csv(geno_df, "data/genotype_permGWAS.csv", row.names = FALSE)
write.csv(pheno_df, "data/phenotype_permGWAS_logvar.csv", row.names = FALSE)




############################################################
## ISOGNIC LINES + DATA FOR GEMMA – area.F / area.F.var
## C. elegans – CeMEE panel
############################################################

## ===============================
## 0. Packages
## ===============================

library(AGHmatrix)
library(data.table)

# Load the phenotypes
pheno_mean <- read.table("data/pheno_mean.txt", header = T)
pheno_var <- read.table("data/pheno_var.txt", header = T)
pheno_CV <- read.table("data/pheno_CV.txt", header = T)

str(pheno_CV)

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


### We will use the bimbam format which is simple for GEMMA
## First load the genotypes
## Load Cemee data

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

bimbam <- data.frame(id=snp_positions$ID,major=snp_positions$REF,minor=snp_positions$ALT)
bimbam <- cbind(bimbam,genotype_info_kept)
write.table(bimbam,"data/genotypes_for_gemma_CV.txt",row.names=FALSE,col.names=FALSE,quote=FALSE) # table geno
write.table(gemma_pheno,file="data/gemma_pheno_CV.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE) #table pheno

gemma_pheno_withlines=data.frame(gemma_pheno)
gemma_pheno_withlines$strain=as.character(genotyped_lines[phenotyped_lines])
write.table(gemma_pheno_withlines,file="data/gemma_pheno_withLines_CV.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)

## ===============================
## 1.2 Distribution of data
## ===============================

pheno_variso <- fread("data/gemma_pheno_withLines_var.txt")
pheno_meaniso <- fread("data/gemma_pheno_withLines.txt")
pheno_CViso <- fread("data/gemma_pheno_withLines_CV.txt")

plot(pheno_meaniso$V1,pheno_variso$V1)
plot(pheno_meaniso$V1, pheno_CViso$V1)
plot(pheno_variso$V1, pheno_CViso$V1)

ggplot(pheno_meaniso, aes(x = V1)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.002, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "Estimated mean", y = "Density",
       title = "Distribution of emmean with density curve")

ggplot(pheno_variso, aes(x = log(V1))) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "Estimated variance", y = "Density",
       title = "Distribution of emmean with density curve")

cor(pheno_meaniso$V1,pheno_variso$V1) #cor de 0.14
cor(pheno_meaniso$V1, pheno_CViso$V1) # cor de -0.47
cor(pheno_variso$V1, pheno_CViso$V1) # cor de 0.0.75

log_emmeaniso_var=log(pheno_variso$V1)
plot(pheno_meaniso$V1, log_emmeaniso_var)

log_emmeaniso_CV=log(pheno_CViso$V1)
plot(pheno_variso$V1, log_emmeaniso_CV)

plot(log_emmeaniso_var, log_emmeaniso_CV)
cor(log_emmeaniso_var, log_emmeaniso_CV) #cor de 0.73

#CV
ggplot(pheno_CViso, aes(x = V1)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.2, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "CV", y = "Density",
       title = "Distribution of CV with density curve")

ggplot(pheno_CViso, aes(x = log(V1))) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "log(CV)", y = "Density",
       title = "Distribution of CV with density curve")

#log-tranformation for GEMMA (CV and variance)

pheno_variso_gemma <- fread("data/gemma_pheno_var.txt")
pheno_CViso_gemma <- fread("data/gemma_pheno_CV.txt")

str(pheno_CViso_gemma)
str(pheno_variso_gemma)

pheno_CViso_gemma$log_V1 <- log(pheno_CViso_gemma$V1)
pheno_variso_gemma$log_V1 <- log(pheno_variso_gemma$V1)

write.table(pheno_CViso_gemma$log_V1, "gemma_pheno_log_CV.txt",sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(pheno_variso_gemma$log_V1, "gemma_pheno_log_var.txt",sep='\t', row.names = FALSE, col.names = FALSE, quote = FALSE)

## ===============================
## 2. Data for permGWAS
## ===============================

#pas de log-transformation
#juste genotype et phenotype




## ===============================
## 3. Data for GEMMA (kinship matrix + pruning)
## ===============================

# créer dossier
system("mkdir -p plink_pruning/pruning01")
#pruning de 0.1
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 0.1 --out plink_pruning/pruning01/CeMEE.pruned.0.1 --allow-extra-chr --double-id")
# Extraire SNPs prunés
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --extract plink_pruning/pruning01/CeMEE.pruned.0.1.prune.in --recode vcf --out plink_pruning/pruning01/CeMEE.pruned.0.1 --allow-extra-chr --double-id")
# compression
system("bgzip plink_pruning/pruning01/CeMEE.pruned.0.1.vcf")

genotype_pruned <- read.vcfR("plink_pruning/pruning01/CeMEE.pruned.0.1.vcf.gz")
dim(genotype_pruned@gt)
genotype_pruned_matrix <- apply(genotype_pruned@gt[,2:1002],1,function(x){
  as.numeric(substring(x,1,1))
})
str(genotype_pruned_matrix)
genotype_pruned_matrix=2*genotype_pruned_matrix
strainID <- tstrsplit(colnames(genotype_pruned@gt)[2:1002],'_')[[1]]
GRM <- Gmatrix((genotype_pruned_matrix), method="VanRaden")
dim(GRM)
## Re-order
GRM_sub = GRM[strainID %in% colnames(genotype_info_kept),strainID %in% colnames(genotype_info_kept)]
strainID_sub = strainID[strainID %in% colnames(genotype_info_kept)]
table(strainID_sub==colnames(genotype_info_kept)) # Same order

head(bimbam)
head(gemma_pheno)

write.table(GRM_sub,file="data/kinship.txt",sep='\t',quote=FALSE,row.names=FALSE,col.names=FALSE)



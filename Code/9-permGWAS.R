## ===============================
## 0. Packages
## ===============================

library(AGHmatrix)
library(data.table)
library(vcfR)

# Load the phenotypes
pheno_var <- fread("data/gemma_pheno_withLines_var.txt")
pheno_mean <- fread("data/gemma_pheno_withLines.txt")
pheno_CV <- fread("data/gemma_pheno_withLines_CV.txt")

str(pheno_var)

#Modif fichier pheno
pheno_var <- pheno_var[, .(line = V2, variance = V1)]
pheno_mean <- pheno_mean[, .(line = V2, mean = V1)]
pheno_CV <- pheno_CV[, .(line = V2, CV = V1)]

#Modif fichier geno

geno <- fread("data/genotype_012.txt")

str(geno)

geno[, rs := paste(CHROM, POS, sep = "_")]   # crée la colonne rs
setcolorder(geno, c("rs", setdiff(names(geno), c("CHROM","POS","rs"))))  # place rs en premier
geno[, c("CHROM","POS") := NULL]  # supprime CHROM et POS

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
ids_geno  <- colnames(geno)
ids_kin   <- rownames(Kinship)

str(geno)

common_ids <- Reduce(intersect, list(ids_pheno, ids_geno, ids_kin))

pheno_mean_perm <- pheno_mean[line %in% common_ids, ]
geno_perm <- geno[, c("rs", intersect(common_ids, colnames(geno))), with = FALSE]
kinship_perm <- Kinship[common_ids, common_ids]


#En csv

fwrite(pheno_mean, "data/pheno_mean_permGWAS.csv", sep = ",")
fwrite(pheno_var,  "data/pheno_var_permGWAS.csv",  sep = ",")
fwrite(pheno_CV,   "data/pheno_CV_permGWAS.csv",   sep = ",")

write.table(geno_perm,
            file = "data/geno_permGWAS.csv",
            sep = ",",
            quote = FALSE,
            row.names = FALSE,
            col.names = TRUE)

write.table(kinship_perm,
            file = "data/kinship_permGWAS.csv",
            sep = ",",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA)




#####################






fwrite(pheno_mean, "data/pheno_mean_permGWAS.txt", sep = "\t")
fwrite(pheno_var, "data/pheno_var_permGWAS.txt", sep = "\t")
fwrite(pheno_CV, "data/pheno_CV_permGWAS.txt", sep = "\t")

str(geno_perm)
str(kinship_perm)


# Pour le genotype
write.table(geno_perm,
            file = "data/geno_permGWAS.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE) 

# Pour la kinship
write.table(kinship_perm,
            file = "data/kinship_permGWAS.txt",
            sep = "\t",
            quote = FALSE,
            row.names = TRUE)






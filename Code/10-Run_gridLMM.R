#preparation resultats gridLMM

library(ggplot2)
library(rutilstimflutre)
library(GridLMM)
library(data.table)
library(matrixStats)
library(lme4)
library(AGHmatrix)
library(qqman)

all_BLUEs_genotyped=fread("phenotype_grid.txt")
all_BLUEs_genotyped_NGM <- all_BLUEs_genotyped[env_label == "NGM"]

genotype <- read.table(
  "geno_grid.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

genotype_matrix <- as.matrix(genotype)
geno_pruned <- fread("pruned_grid.txt")

# enlever colonnes meta PLINK
geno_pruned_ids <- geno_pruned$IID # garder les IDs
geno_pruned <- geno_pruned[, -(1:6)] # enlever colonnes PLINK
genotype_pruned_matrix <- as.matrix(geno_pruned)
rownames(genotype_pruned_matrix) <- geno_pruned_ids #remettre les rownames
colnames(genotype_pruned_matrix) <- sapply(
  strsplit(colnames(genotype_pruned_matrix), "_"),
  function(x) paste(x[1], x[2], sep = "_")
)

# Je définis mes traits - les noms des colonnes sont dans la table des phénotypes
all_BLUEs_genotyped_NGM[, CV_areaF := sqrt(area.F.var) / area.F]
all_BLUEs_genotyped_NGM[, `:=`(date_str = as.factor(date_str), logD = scale(logD), temperature = scale(temperature), rel_humidity = scale(rel_humidity))] #mettre en facteur et scale 

vect_P_traits=c("area.F","area.F.var", "CV_areaF")

# Je récupére les noms des chromosomes dans la table des génotypes
CHR_genotype_matrix = tstrsplit(colnames(genotype_matrix),"_")[[1]]
# Je récupére les noms des chromosomes dans la table des génotypes prunés
CHR_genotype_matrix_pruned = tstrsplit(colnames(genotype_pruned_matrix),"_")[[1]]

## Boucle for sur chaque trait
for(traitID in c("area.F","area.F.var","CV_areaF")){
  print(traitID)
  
  #Je crée une colonne temporaire "y" avec le trait sélectionné
  all_BLUEs_genotyped_NGM$y = all_BLUEs_genotyped_NGM[[traitID]]
  
  # Je scale le trait pour le GWAS
  all_BLUEs_genotyped_NGM$y=scale(all_BLUEs_genotyped_NGM$y)
  
  #AJOUT : si on analyse la variance, on ajoute la moyenne comme covariable
  if(traitID == "area.F.var"){
    all_BLUEs_genotyped_NGM$areaF_scaled = scale(all_BLUEs_genotyped_NGM$area.F)
  }
  
  # Je boucle sur les chromosomes
  for(chrID in c("I","II","III","IV","V","X","MtDNA")){
    
    # Je ne garde que la sous matrice des génotypes pour mon CHR d'interet
    genotype_matrix_CHR = genotype_matrix[,CHR_genotype_matrix==chrID]
    
    # A l'inverse, je ne garde que la sous matrice des génotypes prunés qui exclue ce chromosome
    CHR_genotype_matrix_pruned_CHR=genotype_pruned_matrix[,!CHR_genotype_matrix_pruned==chrID]
    
    # Je recalcule la matrice d'apparentement avec cette nouvelle sous matrice
    GRM_CHR <- Gmatrix(CHR_genotype_matrix_pruned_CHR, method="VanRaden") 
    
    # vecteur pour exclure les NA
    is.NA=colSums(is.na(genotype_matrix_CHR))
    
    #AJOUT : formule conditionnelle selon le trait
    if(traitID == "area.F.var"){
      formula_model <- y ~ 1 + areaF_scaled + logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id)
    } else {
      formula_model <- y ~ 1 + logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id)
    }
    
    # Je fais tourner mon modèle
    gwas.NGM.CHR = GridLMM_GWAS(
      formula = formula_model, # utilisation de la formule conditionnelle
      test_formula = ~1,
      reduced_formula = ~0,
      data = all_BLUEs_genotyped_NGM,
      X = genotype_matrix_CHR[,is.NA==0],
      X_ID = 'pop_label',
      relmat = list(pop_label = GRM_CHR),
      mc.cores = 10,
      method = 'REML'
    )
    
    # Je sauvegarde ma table
    NGM_results <- (gwas.NGM.CHR$results)
    write.table(
      NGM_results,
      paste0("GridLMM_results_mixed_perCHR/gwas_results_",traitID,"_CHR_",chrID,"_NGM.txt"),
      quote=FALSE,
      row.names=FALSE
    )
    
  }
  
}
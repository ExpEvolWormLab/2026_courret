rm(list=ls());gc()

library(vcfR)
library(ggplot2)
library(rutilstimflutre)
library(ASRgenomics)
library(GridLMM)
library(data.table)
library(matrixStats)
library(performance)
library(lme4)
library(AGHmatrix)
library(qqman)

## Je charge les 2 matrices de génotypes + la table des phénotypes avec les cov. environnementales
load('RData/Preparatory_GLMM.RData')

# Je définis mes traits - les noms des colonnes sont dans la table des péhnotypes
vect_P_traits=c("SF",'SB','FS','FB','BS','BF',"area.F","beta")

# Je récupére les noms des chromosomes dans la table des génotypes
CHR_genotype_matrix = tstrsplit(colnames(genotype_matrix),"_")[[1]]
# Je récupére les noms des chromosomes dans la table des génotypes prunés
CHR_genotype_matrix_pruned = tstrsplit(colnames(genotype_pruned_matrix),"_")[[1]]


## Boucle for sur chaque trait
for(traitID in c("beta","SF",'SB','FS','FB','BS','BF',"area.F")){
  print(traitID)
  #Je crée une colonne temporaire "y" avec le trait sélectionné
  all_BLUEs_genotyped_NGM$y=all_BLUEs_genotyped_NGM[,traitID]
  
  # Je scale le trait pour le GWAS
  all_BLUEs_genotyped_NGM$y= scale(all_BLUEs_genotyped_NGM$y)

  
  # Je boucle sur les chromosomes
  for(chrID in c("I","II","III","IV","V","X")){
    # Je ne garde que la sous matrice des génotypes pour mon CHR d'interet
    genotype_matrix_CHR = genotype_matrix[,CHR_genotype_matrix==chrID]
    # A l'inverse, je ne garde que la sous matrice des génotypes prunés qui exclue ce chromosome
    CHR_genotype_matrix_pruned_CHR=genotype_pruned_matrix[,!CHR_genotype_matrix_pruned==chrID]
    # Je recalcule la matrice d'apparentement avec cette nouvelle sous matrice
    GRM_CHR <- Gmatrix(CHR_genotype_matrix_pruned_CHR, method="VanRaden") 
    # vecteur pour exclure les NA
    is.NA=colSums(is.na(genotype_matrix_CHR))
    # Je fais tourner mon modèle
    gwas.NGM.CHR = GridLMM_GWAS(formula = y ~ 1 +logD + rel_humidity + temperature + (1|date_str) + (1|pop_label) + (1|session_id),
                            test_formula = ~1, reduced_formula = ~0,
                            data = all_BLUEs_genotyped_NGM, X = genotype_matrix_CHR[,is.NA==0],
                            X_ID = 'pop_label', relmat = list(pop_label = GRM_CHR), mc.cores = 10 ,method = 'REML')
    # Je sauvegarde ma table
    NGM_results <- (gwas.NGM.CHR$results)
    write.table(NGM_results,paste0("GridLMM_results_mixed_perCHR/gwas_results_",traitID,"_CHR_",chrID,"_NGM.txt"),quote=FALSE,row.names=FALSE)
    
  }
  
}

## La commande pour faire tourner un script depuis le terminal
# nohup Rscript ./10-Run_gridLMM.R > output_grid.txt &


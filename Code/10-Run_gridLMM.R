#preparation resultats gridLMM

library(ggplot2)
library(rutilstimflutre)
library(GridLMM)
library(data.table)
library(matrixStats)
library(lme4)
library(AGHmatrix)
library(qqman)

all_BLUEs_genotyped=fread("data/phenotype_grid.txt")
all_BLUEs_genotyped_NGM <- all_BLUEs_genotyped[env_label == "NGM"]

plot(density(all_BLUEs_genotyped_NGM$area.F.var))
plot(density(log(all_BLUEs_genotyped_NGM$area.F.var)))
length(unique(all_BLUEs_genotyped_NGM$log_area.F.var))
table(table(all_BLUEs_genotyped_NGM$pop_label))

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

# Je définis mes traits
all_BLUEs_genotyped_NGM[, CV_areaF := sqrt(area.F.var) / area.F]

# AJOUT LOG TRANSFORMATION
all_BLUEs_genotyped_NGM[, log_area.F.var := log(area.F.var)]
all_BLUEs_genotyped_NGM[, log_CV_areaF := log(CV_areaF)]

all_BLUEs_genotyped_NGM[, `:=`(
  date_str = as.factor(date_str),
  logD = scale(logD),
  temperature = scale(temperature),
  rel_humidity = scale(rel_humidity)
)]

# Traits analysés
vect_P_traits=c("area.F","log_area.F.var", "log_CV_areaF")

# Je récupére les noms des chromosomes dans la table des génotypes
CHR_genotype_matrix = tstrsplit(colnames(genotype_matrix),"_")[[1]]

# Je récupére les noms des chromosomes dans la table des génotypes prunés
CHR_genotype_matrix_pruned = tstrsplit(colnames(genotype_pruned_matrix),"_")[[1]]

## Boucle for sur chaque trait
for(traitID in vect_P_traits){
  
  print(traitID)
  
  # Je crée une colonne temporaire "y" avec le trait sélectionné
  all_BLUEs_genotyped_NGM$y = all_BLUEs_genotyped_NGM[[traitID]]
  
  # Je scale le trait pour le GWAS
  all_BLUEs_genotyped_NGM$y = scale(all_BLUEs_genotyped_NGM$y)
  
  # AJOUT : si on analyse la variance, on ajoute la moyenne comme covariable
  if(traitID == "log_area.F.var"){
    
    all_BLUEs_genotyped_NGM$areaF_scaled = scale(all_BLUEs_genotyped_NGM$area.F)
    
  }
  
  # Je boucle sur les chromosomes
  for(chrID in c("I","II","III","IV","V","X","MtDNA")){
    
    # Je ne garde que la sous matrice des génotypes pour mon CHR d'interet
    genotype_matrix_CHR = genotype_matrix[,CHR_genotype_matrix==chrID]
    
    # A l'inverse, je ne garde que la sous matrice des génotypes prunés qui exclue ce chromosome
    CHR_genotype_matrix_pruned_CHR = genotype_pruned_matrix[,!CHR_genotype_matrix_pruned==chrID]
    
    # Je recalcule la matrice d'apparentement
    GRM_CHR <- Gmatrix(CHR_genotype_matrix_pruned_CHR, method="VanRaden") 
    
    # vecteur pour exclure les NA
    is.NA = colSums(is.na(genotype_matrix_CHR))
    
    # Formule conditionnelle selon le trait
    if(traitID == "log_area.F.var"){
      
      formula_model <- y ~ 1 + areaF_scaled + logD + rel_humidity + temperature +
        (1|date_str) + (1|pop_label) + (1|session_id)
      
    } else {
      
      formula_model <- y ~ 1 + logD + rel_humidity + temperature +
        (1|date_str) + (1|pop_label) + (1|session_id)
    }
    
    # GWAS
    gwas.NGM.CHR = GridLMM_GWAS(
      formula = formula_model,
      test_formula = ~1,
      reduced_formula = ~0,
      data = all_BLUEs_genotyped_NGM,
      X = genotype_matrix_CHR[,is.NA==0],
      X_ID = 'pop_label',
      relmat = list(pop_label = GRM_CHR),
      mc.cores = 15,
      method = 'REML'
    )
    
    # Sauvegarde
    NGM_results <- (gwas.NGM.CHR$results)
    
    write.table(
      NGM_results,
      paste0(
        "GridLMM_results_mixed_perCHR/gwas_results_",
        traitID,
        "_CHR_",
        chrID,
        "_NGM.txt"
      ),
      quote = FALSE,
      row.names = FALSE
    )
    
  }
}
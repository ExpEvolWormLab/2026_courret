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
all_BLUEs_genotyped_NaCl <- all_BLUEs_genotyped[env_label == "NaCl"]

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
rownames(genotype_pruned_matrix) <- geno_pruned_ids # remettre les rownames
colnames(genotype_pruned_matrix) <- sapply(
  strsplit(colnames(genotype_pruned_matrix), "_"),
  function(x) paste(x[1], x[2], sep = "_")
)

# Définition des traits
all_BLUEs_genotyped_NaCl[, CV_areaF := sqrt(area.F.var) / area.F]

# AJOUT LOG TRANSFORMATION
all_BLUEs_genotyped_NaCl[, log_area.F.var := log(area.F.var)]
all_BLUEs_genotyped_NaCl[, log_CV_areaF := log(CV_areaF)]

all_BLUEs_genotyped_NaCl[, `:=`(
  date_str = as.factor(date_str),
  logD = scale(logD),
  temperature = scale(temperature),
  rel_humidity = scale(rel_humidity)
)]

# Traits analysés
vect_P_traits = c("area.F", "log_area.F.var", "log_CV_areaF")

# Chromosomes
CHR_genotype_matrix = tstrsplit(colnames(genotype_matrix), "_")[[1]]
CHR_genotype_matrix_pruned = tstrsplit(colnames(genotype_pruned_matrix), "_")[[1]]

## Boucle for sur chaque trait
for(traitID in vect_P_traits){
  
  print(traitID)
  
  # Création variable y
  all_BLUEs_genotyped_NaCl$y = all_BLUEs_genotyped_NaCl[[traitID]]
  
  # Scale du trait
  all_BLUEs_genotyped_NaCl$y = scale(all_BLUEs_genotyped_NaCl$y)
  
  # Si analyse de variance → ajouter moyenne
  if(traitID == "log_area.F.var"){
    
    all_BLUEs_genotyped_NaCl$areaF_scaled = scale(all_BLUEs_genotyped_NaCl$area.F)
  }
  
  # Boucle chromosomes
  for(chrID in c("I","II","III","IV","V","X","MtDNA")){
    
    # SNPs du chromosome
    genotype_matrix_CHR = genotype_matrix[,CHR_genotype_matrix == chrID]
    
    # SNPs prunés sans le chromosome courant
    CHR_genotype_matrix_pruned_CHR =
      genotype_pruned_matrix[,!CHR_genotype_matrix_pruned == chrID]
    
    # Matrice apparentement
    GRM_CHR <- Gmatrix(
      CHR_genotype_matrix_pruned_CHR,
      method = "VanRaden"
    ) 
    
    # Exclusion NA
    is.NA = colSums(is.na(genotype_matrix_CHR))
    
    # Formule conditionnelle
    if(traitID == "log_area.F.var"){
      
      formula_model <- y ~ 1 +
        areaF_scaled +
        logD +
        rel_humidity +
        temperature +
        (1|date_str) +
        (1|pop_label) +
        (1|session_id)
      
    } else {
      
      formula_model <- y ~ 1 +
        logD +
        rel_humidity +
        temperature +
        (1|date_str) +
        (1|pop_label) +
        (1|session_id)
    }
    
    # GWAS
    gwas.NaCl.CHR = GridLMM_GWAS(
      formula = formula_model,
      test_formula = ~1,
      reduced_formula = ~0,
      data = all_BLUEs_genotyped_NaCl,
      X = genotype_matrix_CHR[,is.NA == 0],
      X_ID = 'pop_label',
      relmat = list(pop_label = GRM_CHR),
      mc.cores = 10,
      method = 'REML'
    )
    
    # Sauvegarde résultats
    NaCl_results <- (gwas.NaCl.CHR$results)
    
    write.table(
      NaCl_results,
      paste0(
        "GridLMM_results_mixed_perCHR_NaCl/gwas_results_",
        traitID,
        "_CHR_",
        chrID,
        "_NaCl.txt"
      ),
      quote = FALSE,
      row.names = FALSE
    )
    
  }
}
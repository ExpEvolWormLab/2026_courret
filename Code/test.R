library(AGHmatrix)
library(qqman)

all_BLUEs_genotyped=fread("data/phenotype_grid.txt")
all_BLUEs_genotyped_NGM <- all_BLUEs_genotyped[env_label == "NGM"]

genotype <- read.table(
  "data/geno_grid.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

genotype_matrix <- as.matrix(genotype)
geno_pruned <- fread("data/pruned_grid.txt")

# enlever colonnes meta PLINK
geno_pruned_ids <- geno_pruned$IID # garder les IDs
geno_pruned <- geno_pruned[, -(1:6)] # enlever colonnes PLINK
genotype_pruned_matrix <- as.matrix(geno_pruned)
rownames(genotype_pruned_matrix) <- geno_pruned_ids # remettre les rownames

colnames(genotype_pruned_matrix) <- sapply(
  strsplit(colnames(genotype_pruned_matrix), "_"),
  function(x) paste(x[1], x[2], sep = "_")
)

# scale et format
all_BLUEs_genotyped_NGM[, `:=`(
  date_str = as.factor(date_str),
  logD = scale(logD),
  temperature = scale(temperature),
  rel_humidity = scale(rel_humidity)
)]

############################################################
## TRAIT = area.F uniquement
############################################################

traitID = "area.F"

all_BLUEs_genotyped_NGM$y = all_BLUEs_genotyped_NGM[[traitID]]

# scale du trait
all_BLUEs_genotyped_NGM$y = scale(all_BLUEs_genotyped_NGM$y)

############################################################
## CHROMOSOME I uniquement
############################################################

# noms chromosomes
CHR_genotype_matrix = tstrsplit(colnames(genotype_matrix),"_")[[1]]

# chromosomes pruned
CHR_genotype_matrix_pruned = tstrsplit(
  colnames(genotype_pruned_matrix),
  "_"
)[[1]]

chrID = "I"

# garder uniquement chromosome I
genotype_matrix_CHR = genotype_matrix[
  ,
  CHR_genotype_matrix == chrID
]

# enlever chromosome I pour kinship LOCO
CHR_genotype_matrix_pruned_CHR = genotype_pruned_matrix[
  ,
  !CHR_genotype_matrix_pruned == chrID
]

# matrice apparentement
GRM_CHR <- Gmatrix(
  CHR_genotype_matrix_pruned_CHR,
  method = "VanRaden"
)

# enlever SNP avec NA
is.NA = colSums(is.na(genotype_matrix_CHR))

############################################################
## GWAS
############################################################

gwas.NGM.CHR = GridLMM_GWAS(
  
  formula = y ~ 1 + 
    logD + 
    rel_humidity + 
    temperature + 
    (1|date_str) + 
    (1|pop_label) + 
    (1|session_id),
  
  test_formula = ~1,
  reduced_formula = ~0,
  
  data = all_BLUEs_genotyped_NGM,
  
  X = genotype_matrix_CHR[,is.NA==0],
  X_ID = 'pop_label',
  
  relmat = list(pop_label = GRM_CHR),
  
  mc.cores = 3,
  method = 'REML'
)

############################################################
## SAVE RESULTS
############################################################

NGM_results <- gwas.NGM.CHR$results

str(NGM_results)

write.table(
  NGM_results,
  "gwas_results_areaF_CHR_I_NGM.txt",
  quote = FALSE,
  row.names = FALSE
)


library(data.table)
library(ggplot2)
library(dplyr)

############################################################
## LOAD RESULTS
############################################################

GWAS.mean <- fread("gwas_results_areaF_CHR_I_NGM.txt")

############################################################
## PREPARE DATA
############################################################

GWAS.mean <- GWAS.mean %>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP  = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  arrange(CHR, BP)

############################################################
## BH CORRECTION
############################################################

GWAS.mean$p.adjust_BH <- p.adjust(
  GWAS.mean$p_value_REML,
  method = "BH"
)

GWAS.mean$log10p_adj <- -log10(GWAS.mean$p.adjust_BH)

############################################################
## MANHATTAN PLOT — adjusted p-values
############################################################

ggplot(
  GWAS.mean,
  aes(x = BP, y = log10p_adj)
) +
  
  geom_point(
    size = 0.5,
    alpha = 0.7,
    color = "grey40"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_hline(
    yintercept = -log10(0.1),
    color = "orange",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_hline(
    yintercept = -log10(0.2),
    color = "blue",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  labs(
    x = "Position chromosome I",
    y = expression(-log[10](adjusted~p))
  ) +
  
  theme_classic()

############################################################
## MANHATTAN PLOT — raw p-values
############################################################

ggplot(
  GWAS.mean,
  aes(x = BP, y = -log10(p_value_REML))
) +
  
  geom_point(
    size = 0.5,
    alpha = 0.7,
    color = "grey40"
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_hline(
    yintercept = -log10(0.1),
    color = "orange",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  geom_hline(
    yintercept = -log10(0.2),
    color = "blue",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  labs(
    x = "Position chromosome I",
    y = expression(-log[10](p~value))
  ) +
  
  theme_classic()

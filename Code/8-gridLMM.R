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

## 1.2 Outliers de variance

cutoff <- 3

# =========================================================
# NGM
# =========================================================

mean_pheno.varA <- mean(phen_data$area.F.var[phen_data$env_label == "NGM"], na.rm = TRUE)
sd_pheno.varA   <- sd(phen_data$area.F.var[phen_data$env_label == "NGM"], na.rm = TRUE)

z_scores.varA <- (phen_data$area.F.var - mean_pheno.varA) / sd_pheno.varA

table(abs(z_scores.varA[phen_data$env_label == "NGM"]) <= cutoff)

phen_data$exper_data_id[phen_data$env_label == "NGM"][abs(z_scores.varA[phen_data$env_label == "NGM"]) > cutoff]

dim(phen_data)
phen_data <- phen_data[!(phen_data$env_label == "NGM" & abs(z_scores.varA) > cutoff), ]
dim(phen_data)

# =========================================================
# NaCl
# =========================================================

mean_pheno.varA <- mean(phen_data$area.F.var[phen_data$env_label == "NaCl"], na.rm = TRUE)
sd_pheno.varA   <- sd(phen_data$area.F.var[phen_data$env_label == "NaCl"], na.rm = TRUE)

z_scores.varA <- (phen_data$area.F.var - mean_pheno.varA) / sd_pheno.varA

table(abs(z_scores.varA[phen_data$env_label == "NaCl"]) <= cutoff)

phen_data$exper_data_id[phen_data$env_label == "NaCl"][abs(z_scores.varA[phen_data$env_label == "NaCl"]) > cutoff]

dim(phen_data)
phen_data <- phen_data[!(phen_data$env_label == "NaCl" & abs(z_scores.varA) > cutoff), ]
dim(phen_data)

str(phen_data)
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

common_ids <- Reduce(intersect, list(ids_pheno, ids_geno, ids_pruned))

phen_clean_sub <- phen_clean[phen_clean$pop_label %in% common_ids, ]
geno_snp_sub <- geno_mat[common_ids, ]
pruned_sub <- geno_pruned[geno_pruned$FID %in% common_ids, ]

# Pour le phenotype
fwrite(phen_clean_sub, "data/phenotype_grid.txt", sep = "\t")

phen_clean_sub <- fread("data/phenotype_grid.txt")

str(phen_clean_sub)
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


#distrubution des donnees 
str(phen_clean_sub)

library(ggplot2)
library(patchwork)

# Garder uniquement NGM
phen_NGM <- subset(phen_clean, env_label == "NGM")

# Calcul du coefficient de variation
phen_NGM$CV_areaF <- sqrt(phen_NGM$area.F.var) / phen_NGM$area.F

# Distribution variance
p1 <- ggplot(phen_NGM, aes(x = area.F.var)) +
  geom_histogram(bins = 400, fill = "#A6CEE3", color = "black") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of variance (NGM)",
    x = "Variance of area.F",
    y = "Count"
  )

# Distribution coefficient de variation
p2 <- ggplot(phen_NGM, aes(x = log(CV_areaF))) +
  geom_histogram(bins = 4000, fill = "#B2DF8A", color = "black") +
  theme_classic(base_size = 14) +
  labs(
    title = "Distribution of coefficient of variation (NGM)",
    x = "CV of area.F",
    y = "Count"
  )

# Affichage
p1
p2

#correlation mean-var (+env.)

str(phen_clean_sub)
plot(phen_clean_sub$area.F, phen_clean_sub$area.F.var)

table(table(phen_clean_sub$pop_label))

ggplot(phen_clean_sub, aes(x = area.F, y = area.F.var)) +
  
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
cor.test(phen_clean_sub$area.F, phen_clean_sub$area.F.var)


library(ggplot2)
library(dplyr)

# ---------------------------
# Préparation des données
# ---------------------------

plot_data <- phen_clean_sub %>%
  dplyr::select(area.F, area.F.var, env_label) %>%
  na.omit()

# ---------------------------
# Corrélations par environnement
# ---------------------------

cor_NGM <- cor.test(
  plot_data$area.F[plot_data$env_label == "NGM"],
  plot_data$area.F.var[plot_data$env_label == "NGM"],
  method = "pearson"
)

cor_NaCl <- cor.test(
  plot_data$area.F[plot_data$env_label == "NaCl"],
  plot_data$area.F.var[plot_data$env_label == "NaCl"],
  method = "pearson"
)

# Affichage résultats
cor_NGM
cor_NaCl

# ---------------------------
# Labels pour le plot
# ---------------------------

label_NGM <- paste0(
  "NGM: r = ",
  round(cor_NGM$estimate, 3)
)

label_NaCl <- paste0(
  "NaCl: r = ",
  round(cor_NaCl$estimate, 3)
)

# ---------------------------
# Plot
# ---------------------------

ggplot(plot_data,
       aes(x = area.F,
           y = area.F.var,
           color = env_label)) +
  
  geom_point(alpha = 0.6, size = 2) +
  
  annotate(
    "text",
    x = max(plot_data$area.F) * 0.7,
    y = max(plot_data$area.F.var) * 0.95,
    label = label_NGM,
    color = "black",
    hjust = 0
  ) +
  
  annotate(
    "text",
    x = max(plot_data$area.F) * 0.7,
    y = max(plot_data$area.F.var) * 0.80,
    label = label_NaCl,
    color = "black",
    hjust = 0
  ) +
  
  labs(
    x = "Mean",
    y = "Variance",
    color = "Environment"
  ) +
  
  theme_classic()



#distribution de la variance et de la moyenne 

library(ggplot2)
library(dplyr)
library(patchwork)

# =========================================================
# Préparation des données
# =========================================================

phen_clean_sub <- phen_clean_sub %>%
  mutate(env_label = factor(env_label, levels = c("NGM", "NaCl")))

# =========================================================
# Plot distribution variance
# =========================================================

p_var <- ggplot(phen_clean_sub, aes(x = log(area.F.var), fill = env_label)) +
  
  geom_density(alpha = 0.45, linewidth = 1) +
  
  scale_fill_manual(values = c("NGM" = "#6BAED6",
                               "NaCl" = "#FC9272")) +
  
  labs(
    x = "log(variance of body size)",
    y = "Density",
    fill = "Environment"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    legend.position = "top",
    axis.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  )

# =========================================================
# Plot distribution moyenne
# =========================================================

p_mean <- ggplot(phen_clean_sub, aes(x = area.F, fill = env_label)) +
  
  geom_density(alpha = 0.45, linewidth = 1) +
  
  scale_fill_manual(values = c("NGM" = "#6BAED6",
                               "NaCl" = "#FC9272")) +
  
  labs(
    x = "Mean body size",
    y = "Density",
    fill = "Environment"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    legend.position = "none",
    axis.title = element_text(face = "bold")
  )

# =========================================================
# Figure finale
# =========================================================

final_plot <- p_var / p_mean +
  plot_layout(heights = c(3, 1))

final_plot

mean_NaCl<- fread("data/pheno_mean_NaCl.txt")
var_NaCl<- fread("data/pheno_var_NaCl.txt")

str(var_NaCl)








library(data.table)
library(qqman)

chr_I   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_I_NGM.txt")
chr_II  <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_II_NGM.txt")
chr_III <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_III_NGM.txt")
chr_IV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_IV_NGM.txt")
chr_V   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_V_NGM.txt")
chr_X   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_X_NGM.txt")
chr_MtDNA <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_MtDNA_NGM.txt")

chr_I_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F.var_CHR_I_NGM.txt")
chr_II_var  <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F.var_CHR_II_NGM.txt")
chr_III_var <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F.var_CHR_III_NGM.txt")
chr_IV_var  <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F.var_CHR_IV_NGM.txt")
chr_V_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F.var_CHR_V_NGM.txt")
chr_X_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F.var_CHR_X_NGM.txt")
chr_MtDNA_var <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F.var_CHR_MtDNA_NGM.txt")

chr_I_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_CV_areaF_CHR_I_NGM.txt")
chr_II_CV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_CV_areaF_CHR_II_NGM.txt")
chr_III_CV <- fread("GridLMM_results_mixed_perCHR/gwas_results_CV_areaF_CHR_III_NGM.txt")
chr_IV_CV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_CV_areaF_CHR_IV_NGM.txt")
chr_V_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_CV_areaF_CHR_V_NGM.txt")
chr_X_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_CV_areaF_CHR_X_NGM.txt")
chr_MtDNA_CV <- fread("GridLMM_results_mixed_perCHR/gwas_results_CV_areaF_CHR_MtDNA_NGM.txt")

GWAS.mean <- rbind(chr_I, chr_II, chr_III, chr_IV, chr_V, chr_X, chr_MtDNA)
GWAS.var <- rbind(chr_I_var, chr_II_var, chr_III_var, chr_IV_var, chr_V_var, chr_X_var, chr_MtDNA_var)
GWAS.cv <- rbind(chr_I_CV, chr_II_CV, chr_III_CV, chr_IV_CV, chr_V_CV, chr_X_CV, chr_MtDNA_CV)

### Mean ###

# extraire chromosome et position
GWAS.mean[, c("CHR", "BP") := tstrsplit(X_ID, "_", fixed = TRUE)]
GWAS.mean[, BP := as.numeric(BP)] # convertir BP en numérique
GWAS.mean[, P := p_value_REML]
GWAS.mean[, CHR_num := as.numeric(factor(CHR, levels = c("I","II","III","IV","V","X","MtDNA")))]``
GWAS.mean[, SNP := X_ID]

manhattan(GWAS.mean,
          chr = "CHR_num",
          bp = "BP",
          p = "P",
          snp = "SNP",
          col = c("#A6CEE3", "#1F78B4"),
          suggestiveline = -log10(0.2),
          genomewideline = -log10(0.05 / nrow(GWAS.mean)),
          main = "GWAS mean trait")



GWAS.mean[, qval := p.adjust(p_value_REML, method = "BH")]
GWAS.mean[, log10p := -log10(p_value_REML)]
GWAS.mean[, log10q := -log10(qval)]

chr_order <- c("I", "II", "III", "IV", "V", "X", "MtDNA")
GWAS.mean[, CHR := factor(CHR, levels = chr_order)]

setorder(GWAS.mean, CHR, BP)

GWAS.mean[, BP_cum := BP]

offset <- 0
gap <- 1e6  # espace entre chromosomes (à ajuster)

chr_centers <- c()  # pour les labels

for(chr in chr_order){
  idx <- GWAS.mean$CHR == chr
  
  GWAS.mean[idx, BP_cum := BP + offset]
  
  # centre du chromosome pour label
  chr_centers[chr] <- mean(GWAS.mean$BP_cum[idx], na.rm = TRUE)
  
  # update offset avec espace
  offset <- max(GWAS.mean$BP_cum[idx], na.rm = TRUE) + gap
}

#Manhattan plot
ggplot(GWAS.mean, aes(x = BP_cum, y = log10q, color = CHR)) +
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(values = rep(c("#A6CEE3", "#1F78B4"), length(chr_order))) +
  
  # ligne seuil
  geom_hline(yintercept = -log10(0.2),
             linetype = "dashed",
             color = "red") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "yellow") +
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed",
             color = "orange") +
  
  # axes avec chromosomes
  scale_x_continuous(label = names(chr_centers),
                     breaks = chr_centers) +
  
  labs(x = "Chromosomes", y = expression(-log[10](p))) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

#effect



### Variance ###

GWAS.var[, c("CHR", "BP") := tstrsplit(X_ID, "_", fixed = TRUE)]
GWAS.var[, BP := as.numeric(BP)]

GWAS.var[, qval := p.adjust(p_value_REML, method = "BH")]
GWAS.var[, log10p := -log10(p_value_REML)]
GWAS.var[, log10q := -log10(qval)]

chr_order <- c("I", "II", "III", "IV", "V", "X", "MtDNA")
GWAS.var[, CHR := factor(CHR, levels = chr_order)]

setorder(GWAS.var, CHR, BP)

GWAS.var[, BP_cum := BP]

offset <- 0
gap <- 1e6

chr_centers_var <- numeric(length(chr_order))
names(chr_centers_var) <- chr_order

for(chr in chr_order){
  idx <- GWAS.var$CHR == chr
  
  GWAS.var[idx, BP_cum := BP + offset]
  
  chr_centers_var[chr] <- mean(GWAS.var$BP_cum[idx], na.rm = TRUE)
  
  offset <- max(GWAS.var$BP_cum[idx], na.rm = TRUE) + gap
}

ggplot(GWAS.var, aes(x = BP_cum, y = log10q, color = CHR)) +
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(values = rep(c("#A6CEE3", "#1F78B4"), length(chr_order))) +
  
  geom_hline(yintercept = -log10(0.2),
             linetype = "dashed",
             color = "red") +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed",
             color = "yellow") +
  geom_hline(yintercept = -log10(0.1),
             linetype = "dashed",
             color = "orange") +
  
  scale_x_continuous(
    breaks = chr_centers_var,
    labels = names(chr_centers_var)
  ) +
  
  labs(x = "Chromosomes", y = expression(-log[10](p))) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )


### Coef of variation ###

GWAS.cv[, c("CHR", "BP") := tstrsplit(X_ID, "_", fixed = TRUE)]
GWAS.cv[, BP := as.numeric(BP)]

GWAS.cv[, qval := p.adjust(p_value_REML, method = "BH")]
GWAS.cv[, log10p := -log10(p_value_REML)]
GWAS.cv[, log10q := -log10(qval)]

GWAS.cv[, CHR := factor(CHR, levels = chr_order)]

setorder(GWAS.cv, CHR, BP)

GWAS.cv[, BP_cum := BP]

offset <- 0
gap <- 1e6

chr_centers_cv <- numeric(length(chr_order))
names(chr_centers_cv) <- chr_order

for(chr in chr_order){
  idx <- GWAS.cv$CHR == chr
  
  GWAS.cv[idx, BP_cum := BP + offset]
  
  chr_centers_cv[chr] <- mean(GWAS.cv$BP_cum[idx], na.rm = TRUE)
  
  offset <- max(GWAS.cv$BP_cum[idx], na.rm = TRUE) + gap
}

ggplot(GWAS.cv, aes(x = BP_cum, y = log10p, color = CHR)) +
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(values = rep(c("#A6CEE3", "#1F78B4"), length(chr_order))) +
  
  geom_hline(yintercept = -log10(0.2),
             linetype = "dashed",
             color = "red") +
  
  scale_x_continuous(
    breaks = chr_centers_cv,
    labels = names(chr_centers_cv)
  ) +
  
  labs(x = "Chromosomes", y = expression(-log[10](p))) +
  
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )







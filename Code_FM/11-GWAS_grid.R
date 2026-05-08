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
str(GWAS.mean)

library(ggplot2)
library(dplyr)

# ---------------------------
# Préparation des données
# ---------------------------

GWAS.mean <- GWAS.mean %>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP  = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  arrange(CHR, BP)

# ---------------------------
# Position cumulée
# ---------------------------

chr_sizes <- GWAS.mean %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP), .groups = "drop") %>%
  mutate(offset = cumsum(lag(max_BP, default = 0)))

GWAS.mean <- GWAS.mean %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

# ---------------------------
# Axis
# ---------------------------

axis_df <- GWAS.mean %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum), .groups = "drop")

# ---------------------------
# p.adjust uniquement
# ---------------------------

GWAS.mean$p.adjust_BH <- p.adjust(GWAS.mean$p_value_REML, method = "BH")
GWAS.mean$log10p_adj <- -log10(GWAS.mean$p.adjust_BH)

names(GWAS.mean)

summary(GWAS.mean$p.adjust_BH)
length(unique(GWAS.mean$p.adjust_BH))

summary(GWAS.mean$log10p_adj)
summary(-log10(GWAS.mean$p.adjust_BH))

tapply(GWAS.mean$p_value_REML, GWAS.mean$CHR, length)
tapply(GWAS.mean$p.adjust_BH, GWAS.mean$CHR, range)

# ---------------------------
# Manhattan plot
# ---------------------------

ggplot(GWAS.mean, aes(x = BP_cum, y = log10p_adj, color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(values = rep(c("grey70", "grey40"), length(unique(GWAS.mean$CHR)))) +
  
  geom_hline(yintercept = -log10(0.05),
             color = "red", linetype = "dashed", linewidth = 1)+

  geom_hline(yintercept = -log10(0.1),
           color = "orange", linetype = "dashed", linewidth = 1)+

  geom_hline(yintercept = -log10(0.2),
           color = "blue", linetype = "dashed", linewidth = 1)+
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](adjusted~p)),
    color = "Chromosome"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right"
  )

#p-value RAW
ggplot(GWAS.mean, aes(x = BP_cum, y = -log10(p_value_REML), color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(values = rep(c("grey70", "grey40"), length(unique(GWAS.mean$CHR)))) +
  
  geom_hline(yintercept = -log10(0.05),
             color = "red", linetype = "dashed", linewidth = 1) +
  
  geom_hline(yintercept = -log10(0.1),
             color = "orange", linetype = "dashed", linewidth = 1) +
  
  geom_hline(yintercept = -log10(0.2),
             color = "blue", linetype = "dashed", linewidth = 1) +
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](p~value)),
    color = "Chromosome"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "right"
  )



## par chromosome 
add_bh <- function(df) {
  df$p_adj_BH <- p.adjust(df$p_value_REML, method = "BH")
  df$log10p_adj <- -log10(df$p_adj_BH)
  return(df)
}
chr_I     <- add_bh(chr_I)
chr_II    <- add_bh(chr_II)
chr_III   <- add_bh(chr_III)
chr_IV    <- add_bh(chr_IV)
chr_V     <- add_bh(chr_V)
chr_X     <- add_bh(chr_X)
chr_MtDNA <- add_bh(chr_MtDNA)

res <- chr_I
ggplot(res, aes(x = X_ID, y = log10p_adj)) +
  
  geom_point(aes(color = p_adj_BH < 0.05),
             size = 0.6, alpha = 0.7) +
  
  scale_color_manual(values = c("grey70", "red")) +
  
  geom_hline(yintercept = -log10(0.05),
             color = "red", linetype = "dashed", linewidth = 1) +
  
  geom_hline(yintercept = -log10(0.1),
             color = "orange", linetype = "dashed", linewidth = 1) +
  
  geom_hline(yintercept = -log10(0.2),
             color = "blue", linetype = "dashed", linewidth = 1) +
  
  theme_classic() +
  
  labs(
    x = "Position (bp)",
    y = expression(-log[10](BH~adjusted~p)),
    color = "Significant (BH < 0.05)"
  )





#var 

# ---------------------------
# 🧬 Préparation des données
# ---------------------------

GWAS.var <- GWAS.var %>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP  = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  arrange(CHR, BP)

# ---------------------------
# 📍 Position cumulée
# ---------------------------

chr_sizes <- GWAS.var %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP), .groups = "drop") %>%
  mutate(offset = cumsum(lag(max_BP, default = 0)))

GWAS.var <- GWAS.var %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

# ---------------------------
# 🧭 Axis
# ---------------------------

axis_df <- GWAS.var %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum), .groups = "drop")

# ---------------------------
# 🔥 p.adjust uniquement
# ---------------------------

GWAS.var$p_adj <- p.adjust(GWAS.var$p_value_REML, method = "BH")
GWAS.var$log10p_adj <- -log10(GWAS.var$p_adj)

# ---------------------------
# 🎨 Manhattan plot
# ---------------------------

ggplot(GWAS.var, aes(x = BP_cum, y = log10p_adj, color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(values = rep(c("grey70", "grey40"), length(unique(GWAS.var$CHR)))) +
  
  geom_hline(yintercept = -log10(0.05),
             color = "red", linetype = "dashed", linewidth = 1)+
  
  geom_hline(yintercept = -log10(0.1),
             color = "orange", linetype = "dashed", linewidth = 1)+
  
  geom_hline(yintercept = -log10(0.2),
             color = "blue", linetype = "dashed", linewidth = 1)+
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](adjusted~p)),
    color = "Chromosome"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right"
  )

# CV


# ---------------------------
# 🧬 Préparation des données
# ---------------------------

GWAS.cv <- GWAS.cv %>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP  = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  arrange(CHR, BP)

# ---------------------------
# 📍 Position cumulée
# ---------------------------

chr_sizes <- GWAS.cv %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP), .groups = "drop") %>%
  mutate(offset = cumsum(lag(max_BP, default = 0)))

GWAS.cv <- GWAS.cv %>%
  left_join(chr_sizes, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

# ---------------------------
# 🧭 Axis
# ---------------------------

axis_df <- GWAS.cv %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum), .groups = "drop")

# ---------------------------
# 🔥 p.adjust uniquement
# ---------------------------

GWAS.cv$p_adj <- p.adjust(GWAS.cv$p_value_REML, method = "BH")
GWAS.cv$log10p_adj <- -log10(GWAS.cv$p_adj)

# ---------------------------
# 🎨 Manhattan plot
# ---------------------------

ggplot(GWAS.cv, aes(x = BP_cum, y = log10p_adj, color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(values = rep(c("grey70", "grey40"), length(unique(GWAS.cv$CHR)))) +
  
  geom_hline(yintercept = -log10(0.05),
             color = "red", linetype = "dashed", linewidth = 1)+
  
  geom_hline(yintercept = -log10(0.1),
             color = "orange", linetype = "dashed", linewidth = 1)+
  
  geom_hline(yintercept = -log10(0.2),
             color = "blue", linetype = "dashed", linewidth = 1)+
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](adjusted~p)),
    color = "Chromosome"
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right"
  )



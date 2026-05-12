library(data.table)
library(qqman)

rm(list = ls())

chr_I   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_I_NGM.txt")
chr_II  <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_II_NGM.txt")
chr_III <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_III_NGM.txt")
chr_IV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_IV_NGM.txt")
chr_V   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_V_NGM.txt")
chr_X   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_X_NGM.txt")
chr_MtDNA <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_MtDNA_NGM.txt")

chr_I_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_I_NGM.txt")
chr_II_var  <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_II_NGM.txt")
chr_III_var <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_III_NGM.txt")
chr_IV_var  <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_IV_NGM.txt")
chr_V_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_V_NGM.txt")
chr_X_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_X_NGM.txt")
chr_MtDNA_var <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_MtDNA_NGM.txt")

chr_I_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_I_NGM.txt")
chr_II_CV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_II_NGM.txt")
chr_III_CV <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_III_NGM.txt")
chr_IV_CV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_IV_NGM.txt")
chr_V_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_V_NGM.txt")
chr_X_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_X_NGM.txt")
chr_MtDNA_CV <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_MtDNA_NGM.txt")

GWAS.mean <- rbind(chr_I, chr_II, chr_III, chr_IV, chr_V, chr_X, chr_MtDNA)
GWAS.var <- rbind(chr_I_var, chr_II_var, chr_III_var, chr_IV_var, chr_V_var, chr_X_var, chr_MtDNA_var)
GWAS.cv <- rbind(chr_I_CV, chr_II_CV, chr_III_CV, chr_IV_CV, chr_V_CV, chr_X_CV, chr_MtDNA_CV)


# =========================================================
# MEAN
# =========================================================


str(GWAS.mean)

library(ggplot2)
library(dplyr)

### Mean ###
str(GWAS.mean)

library(ggplot2)
library(dplyr)

# ---------------------------
# Préparation des données
# ---------------------------

chr_order <- c("I", "II", "III", "IV", "V", "X", "MtDNA")

GWAS.mean <- GWAS.mean %>%
  
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP  = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  
  mutate(
    CHR = factor(CHR, levels = chr_order)
  ) %>%
  
  arrange(CHR, BP)

# ---------------------------
# Espacement chromosomes
# ---------------------------

spacing <- 5e5

chr_sizes <- GWAS.mean %>%
  
  group_by(CHR) %>%
  
  summarise(
    max_BP = max(BP),
    .groups = "drop"
  ) %>%
  
  mutate(
    chr_len = max_BP + spacing,
    offset = cumsum(lag(chr_len, default = 0))
  )

# ---------------------------
# Position cumulée
# ---------------------------

GWAS.mean <- GWAS.mean %>%
  
  left_join(chr_sizes, by = "CHR") %>%
  
  mutate(
    BP_cum = BP + offset
  )

# ---------------------------
# Axis
# ---------------------------

axis_df <- GWAS.mean %>%
  
  group_by(CHR) %>%
  
  summarise(
    center = mean(BP_cum),
    .groups = "drop"
  )

# ---------------------------
# p.adjust
# ---------------------------

GWAS.mean$p.adjust_BH <- p.adjust(
  GWAS.mean$p_value_REML,
  method = "BH"
)

GWAS.mean$log10p_adj <- -log10(
  GWAS.mean$p.adjust_BH
)

# ---------------------------
# Manhattan plot
# ---------------------------

ggplot(
  GWAS.mean,
  aes(
    x = BP_cum,
    y = log10p_adj,
    color = CHR
  )
) +
  
  geom_point(
    size = 0.5,
    alpha = 0.7
  ) +
  
  scale_color_manual(
    values = rep(
      c("grey70", "grey40"),
      length(unique(GWAS.mean$CHR))
    )
  ) +
  
  geom_hline(
    yintercept = -log10(0.05),
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  
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
    legend.position = "right",
    axis.text.x = element_text(size = 12),
    axis.title = element_text(size = 14)
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



# SNP significatif ou non
GWAS.mean$signif <- GWAS.mean$p.adjust_BH < 0.05

# ---------------------------
# Manhattan plot des beta
# ---------------------------

ggplot(GWAS.mean,
       aes(x = BP_cum,
           y = beta.5,
           color = signif)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "FALSE" = "grey70",
      "TRUE"  = "red"
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(beta[1]),
    color = "FDR < 0.05"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "right"
  )


# =========================================================
# Variance
# =========================================================


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
  
  #geom_hline(yintercept = -log10(0.1),
  #          color = "orange", linetype = "dashed", linewidth = 1)+
  
  #geom_hline(yintercept = -log10(0.2),
  #           color = "blue", linetype = "dashed", linewidth = 1)+
  
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


# SNP significatif ou non
GWAS.var$signif <- GWAS.var$p_adj < 0.05

# ---------------------------
# Manhattan plot des beta
# ---------------------------

ggplot(GWAS.var,
       aes(x = BP_cum,
           y = beta.6,
           color = signif)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "FALSE" = "grey70",
      "TRUE"  = "red"
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(beta[1]),
    color = "FDR < 0.05"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "right"
  )





# =========================================================
# CV
# =========================================================


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
  
  #geom_hline(yintercept = -log10(0.1),
  #           color = "orange", linetype = "dashed", linewidth = 1)+
  
  #geom_hline(yintercept = -log10(0.2),
  #           color = "blue", linetype = "dashed", linewidth = 1)+
  
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

# SNP significatif ou non
GWAS.cv$signif <- GWAS.cv$p_adj < 0.05

# ---------------------------
# Manhattan plot des beta
# ---------------------------

ggplot(GWAS.cv,
       aes(x = BP_cum,
           y = beta.5,
           color = signif)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(
    values = c(
      "FALSE" = "grey70",
      "TRUE"  = "red"
    )
  ) +
  
  geom_hline(
    yintercept = 0,
    color = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  
  scale_x_continuous(
    breaks = axis_df$center,
    labels = axis_df$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(beta[1]),
    color = "FDR < 0.05"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "right"
  )












# =========================================================
# SNP significatifs communs Mean vs Variance
# =========================================================
install.packages("VennDiagram")
library(VennDiagram)
library(grid)

# ---------------------------
# SNP significatifs
# ---------------------------

sig_mean <- GWAS.mean %>%
  filter(p.adjust_BH < 0.05) %>%
  pull(X_ID)

sig_var <- GWAS.var %>%
  filter(p_adj < 0.05) %>%
  pull(X_ID)

# Vérification
length(sig_mean)
length(sig_var)

# SNP communs
common_snps <- intersect(sig_mean, sig_var)

length(common_snps)

# Affichage SNP communs
common_snps

# =========================================================
# Venn diagram
# =========================================================

venn.plot <- draw.pairwise.venn(
  
  area1 = length(unique(sig_mean)),
  area2 = length(unique(sig_var)),
  cross.area = length(unique(common_snps)),
  
  category = c("Mean", "Variance"),
  
  fill = c("skyblue", "orange"),
  
  alpha = c(0.5, 0.5),
  
  cex = 2,
  cat.cex = 2,
  
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05)
)

grid.draw(venn.plot)

graphics.off()

library(VennDiagram)
library(grid)

# SNP significatifs mean
SNP_mean <- GWAS.mean$X_ID[GWAS.mean$p.adjust_BH < 0.05]

# SNP significatifs variance
SNP_var <- GWAS.var$X_ID[GWAS.var$p_adj < 0.05]

# Nombre en commun
length(intersect(SNP_mean, SNP_var))

# Venn diagram
grid.newpage()

draw.pairwise.venn(
  
  area1 = length(unique(SNP_mean)),
  area2 = length(unique(SNP_var)),
  cross.area = length(intersect(SNP_mean, SNP_var)),
  
  category = c("Mean", "Variance"),
  
  fill = c("skyblue", "tomato"),
  alpha = c(0.5, 0.5),
  
  cex = 2,
  cat.cex = 2,
  
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  
  scaled = FALSE
)








library(dplyr)
library(ggplot2)

# ---------------------------
# Préparation mean
# ---------------------------

mean_df <- GWAS.mean%>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP  = as.numeric(sub(".*_", "", X_ID)),
    trait = "mean"
  )

mean_df$p_adj <- p.adjust(mean_df$p_value_REML, method = "BH")
mean_df$logP <- -log10(mean_df$p_adj)

# ---------------------------
# Préparation variance
# ---------------------------

var_df <- GWAS.var %>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP  = as.numeric(sub(".*_", "", X_ID)),
    trait = "variance"
  )

var_df$p_adj <- p.adjust(var_df$p_value_REML, method = "BH")
var_df$logP <- -log10(var_df$p_adj)

# ---------------------------
# Combine
# ---------------------------

plot_df <- bind_rows(
  mean_df %>% mutate(y = logP),
  var_df %>% mutate(y = -logP)   # miroir
)

# ---------------------------
# Ordre chromosomes
# ---------------------------

chr_order <- c("I","II","III","IV","V","X","MtDNA")

plot_df$CHR <- factor(plot_df$CHR, levels = chr_order)

plot_df <- plot_df %>%
  arrange(CHR, BP) %>%
  group_by(CHR) %>%
  mutate(BP_cum = BP + cumsum(lag(max(BP, na.rm = TRUE), default = 0))) %>%
  ungroup()

# centres pour axes
axis_df <- plot_df %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum), .groups = "drop")

# ---------------------------
# Plot miroir
# ---------------------------

ggplot(plot_df, aes(x = BP_cum, y = y, color = trait)) +
  
  geom_hline(yintercept = 0, color = "black") +
  
  geom_point(size = 0.4, alpha = 0.7) +
  
  scale_color_manual(values = c("mean" = "steelblue", "variance" = "tomato")) +
  
  scale_x_continuous(breaks = axis_df$center, labels = axis_df$CHR) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](adjusted~p) ~ "(mean en haut / variance en bas)"),
    color = "Trait"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "right"
  )


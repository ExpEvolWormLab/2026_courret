library(data.table)
library(qqman)
library(ggplot2)
library(dplyr)

# =========================================================
# LOAD GWAS RESULTS NaCl
# =========================================================

chr_I.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_I_NaCl.txt")
chr_II.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_II_NaCl.txt")
chr_III.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_III_NaCl.txt")
chr_IV.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_IV_NaCl.txt")
chr_V.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_V_NaCl.txt")
chr_X.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_X_NaCl.txt")
chr_MtDNA.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_MtDNA_NaCl.txt")

chr_I_var.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_I_NaCl.txt")
chr_II_var.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_II_NaCl.txt")
chr_III_var.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_III_NaCl.txt")
chr_IV_var.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_IV_NaCl.txt")
chr_V_var.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_V_NaCl.txt")
chr_X_var.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_X_NaCl.txt")
chr_MtDNA_var.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_MtDNA_NaCl.txt")

chr_I_CV.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_I_NaCl.txt")
chr_II_CV.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_II_NaCl.txt")
chr_III_CV.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_III_NaCl.txt")
chr_IV_CV.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_IV_NaCl.txt")
chr_V_CV.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_V_NaCl.txt")
chr_X_CV.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_X_NaCl.txt")
chr_MtDNA_CV.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_MtDNA_NaCl.txt")

# =========================================================
# MERGE RESULTS
# =========================================================

GWAS.mean.NaCl <- rbind(
  chr_I.NaCl,
  chr_II.NaCl,
  chr_III.NaCl,
  chr_IV.NaCl,
  chr_V.NaCl,
  chr_X.NaCl,
  chr_MtDNA.NaCl
)

GWAS.var.NaCl <- rbind(
  chr_I_var.NaCl,
  chr_II_var.NaCl,
  chr_III_var.NaCl,
  chr_IV_var.NaCl,
  chr_V_var.NaCl,
  chr_X_var.NaCl,
  chr_MtDNA_var.NaCl
)

GWAS.cv.NaCl <- rbind(
  chr_I_CV.NaCl,
  chr_II_CV.NaCl,
  chr_III_CV.NaCl,
  chr_IV_CV.NaCl,
  chr_V_CV.NaCl,
  chr_X_CV.NaCl,
  chr_MtDNA_CV.NaCl
)

# =========================================================
# MEAN
# =========================================================

GWAS.mean.NaCl <- GWAS.mean.NaCl %>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  arrange(CHR, BP)

chr_sizes.mean.NaCl <- GWAS.mean.NaCl %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP), .groups = "drop") %>%
  mutate(offset = cumsum(lag(max_BP, default = 0)))

GWAS.mean.NaCl <- GWAS.mean.NaCl %>%
  left_join(chr_sizes.mean.NaCl, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

axis_df.mean.NaCl <- GWAS.mean.NaCl %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum), .groups = "drop")

GWAS.mean.NaCl$p.adjust_BH <- p.adjust(
  GWAS.mean.NaCl$p_value_REML,
  method = "BH"
)

GWAS.mean.NaCl$log10p_adj <- -log10(GWAS.mean.NaCl$p.adjust_BH)

ggplot(GWAS.mean.NaCl,
       aes(x = BP_cum,
           y = log10p_adj,
           color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(
    values = rep(c("grey70", "grey40"),
                 length(unique(GWAS.mean.NaCl$CHR)))
  ) +
  
  geom_hline(yintercept = -log10(0.05),
             color = "red",
             linetype = "dashed",
             linewidth = 1) +
  
  #geom_hline(yintercept = -log10(0.1),
  #           color = "orange",
  #           linetype = "dashed",
  #           linewidth = 1) +
  
  #geom_hline(yintercept = -log10(0.2),
  #           color = "blue",
  #           linetype = "dashed",
  #           linewidth = 1) +
  
  scale_x_continuous(
    breaks = axis_df.mean.NaCl$center,
    labels = axis_df.mean.NaCl$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](adjusted~p)),
    color = "Chromosome"
  ) +
  
  theme_classic()

# SNP significatif ou non
GWAS.mean.NaCl$signif <- GWAS.mean.NaCl$p.adjust_BH < 0.05

# ---------------------------
# Manhattan plot des beta
# ---------------------------

ggplot(GWAS.mean.NaCl,
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
    y = expression(beta),
    color = "FDR < 0.05"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "right"
  )



# =========================================================
# VARIANCE
# =========================================================

GWAS.var.NaCl <- GWAS.var.NaCl %>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  arrange(CHR, BP)

chr_sizes.var.NaCl <- GWAS.var.NaCl %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP), .groups = "drop") %>%
  mutate(offset = cumsum(lag(max_BP, default = 0)))

GWAS.var.NaCl <- GWAS.var.NaCl %>%
  left_join(chr_sizes.var.NaCl, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

axis_df.var.NaCl <- GWAS.var.NaCl %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum), .groups = "drop")

GWAS.var.NaCl$p.adjust_BH <- p.adjust(
  GWAS.var.NaCl$p_value_REML,
  method = "BH"
)

GWAS.var.NaCl$log10p_adj <- -log10(GWAS.var.NaCl$p.adjust_BH)

ggplot(GWAS.var.NaCl,
       aes(x = BP_cum,
           y = log10p_adj,
           color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(
    values = rep(c("grey70", "grey40"),
                 length(unique(GWAS.var.NaCl$CHR)))
  ) +
  
  geom_hline(yintercept = -log10(0.05),
             color = "red",
             linetype = "dashed",
             linewidth = 1) +
  
  #geom_hline(yintercept = -log10(0.1),
  #          color = "orange",
  #           linetype = "dashed",
  #           linewidth = 1) +
  
  #geom_hline(yintercept = -log10(0.2),
  #           color = "blue",
  #           linetype = "dashed",
  #           linewidth = 1) +
  
  scale_x_continuous(
    breaks = axis_df.var.NaCl$center,
    labels = axis_df.var.NaCl$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](adjusted~p)),
    color = "Chromosome"
  ) +
  
  theme_classic()

# SNP significatif ou non
GWAS.var.NaCl$signif <- GWAS.var.NaCl$p.adjust_BH < 0.05

# ---------------------------
# Manhattan plot des beta
# ---------------------------

ggplot(GWAS.var.NaCl,
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
    y = expression(beta),
    color = "FDR < 0.05"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "right"
  )

# =========================================================
# CV
# =========================================================

GWAS.cv.NaCl <- GWAS.cv.NaCl %>%
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  arrange(CHR, BP)

chr_sizes.cv.NaCl <- GWAS.cv.NaCl %>%
  group_by(CHR) %>%
  summarise(max_BP = max(BP), .groups = "drop") %>%
  mutate(offset = cumsum(lag(max_BP, default = 0)))

GWAS.cv.NaCl <- GWAS.cv.NaCl %>%
  left_join(chr_sizes.cv.NaCl, by = "CHR") %>%
  mutate(BP_cum = BP + offset)

axis_df.cv.NaCl <- GWAS.cv.NaCl %>%
  group_by(CHR) %>%
  summarise(center = mean(BP_cum), .groups = "drop")

GWAS.cv.NaCl$p.adjust_BH <- p.adjust(
  GWAS.cv.NaCl$p_value_REML,
  method = "BH"
)

GWAS.cv.NaCl$log10p_adj <- -log10(GWAS.cv.NaCl$p.adjust_BH)

ggplot(GWAS.cv.NaCl,
       aes(x = BP_cum,
           y = log10p_adj,
           color = CHR)) +
  
  geom_point(size = 0.5, alpha = 0.7) +
  
  scale_color_manual(
    values = rep(c("grey70", "grey40"),
                 length(unique(GWAS.cv.NaCl$CHR)))
  ) +
  
  geom_hline(yintercept = -log10(0.05),
             color = "red",
             linetype = "dashed",
             linewidth = 1) +
  
  #geom_hline(yintercept = -log10(0.1),
  #           color = "orange",
  #           linetype = "dashed",
  #           linewidth = 1) +
  
  #geom_hline(yintercept = -log10(0.2),
  #           color = "blue",
  #           linetype = "dashed",
  #           linewidth = 1) +
  
  scale_x_continuous(
    breaks = axis_df.cv.NaCl$center,
    labels = axis_df.cv.NaCl$CHR
  ) +
  
  labs(
    x = "Chromosome",
    y = expression(-log[10](adjusted~p)),
    color = "Chromosome"
  ) +
  
  theme_classic()



# SNP significatif ou non
GWAS.cv.NaCl$signif <- GWAS.cv.NaCl$p.adjust_BH < 0.05

# ---------------------------
# Manhattan plot des beta
# ---------------------------

ggplot(GWAS.cv.NaCl,
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
    y = expression(beta),
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

sig_mean <- GWAS.mean.NaCl %>%
  filter(p.adjust_BH < 0.05) %>%
  pull(X_ID)

sig_var <- GWAS.var.NaCl %>%
  filter(p.adjust_BH < 0.05) %>%
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
SNP_mean <- GWAS.var.NaCl$X_ID[GWAS.var.NaCl$p.adjust_BH < 0.05]

# SNP significatifs variance
SNP_var <- GWAS.cv.NaCl$X_ID[GWAS.cv.NaCl$p.adjust_BH < 0.05]

# Nombre en commun
length(intersect(SNP_mean, SNP_var))

# Venn diagram
grid.newpage()

draw.pairwise.venn(
  
  area1 = length(unique(SNP_mean)),
  area2 = length(unique(SNP_var)),
  cross.area = length(intersect(SNP_mean, SNP_var)),
  
  category = c("variance", "CV"),
  
  fill = c("skyblue", "tomato"),
  alpha = c(0.5, 0.5),
  
  cex = 2,
  cat.cex = 2,
  
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  
  scaled = FALSE
)




library(dplyr)
library(VennDiagram)
library(grid)

# =========================================================
# 1. FONCTION PROPRE DE CLEAN SNPs
# =========================================================

get_sig_snps <- function(df, pcol){
  df %>%
    filter(!is.na(X_ID)) %>%
    mutate(X_ID = as.character(X_ID)) %>%
    filter(.data[[pcol]] < 0.05) %>%
    pull(X_ID) %>%
    unique()
}

# =========================================================
# 2. EXTRACTION DES 4 SETS
# =========================================================

NGM_mean   <- get_sig_snps(GWAS.mean, "p.adjust_BH")
NGM_var    <- get_sig_snps(GWAS.var, "p_adj")

NaCl_mean  <- get_sig_snps(GWAS.mean.NaCl, "p.adjust_BH")
NaCl_var   <- get_sig_snps(GWAS.var.NaCl, "p.adjust_BH")

# =========================================================
# 3. CHECK RAPIDE (IMPORTANT DEBUG)
# =========================================================

length(NGM_mean)
length(NGM_var)
length(NaCl_mean)
length(NaCl_var)

# =========================================================
# 4. INTERSECTIONS (SAFE VERSION)
# =========================================================

n12 <- length(intersect(NGM_mean, NGM_var))
n13 <- length(intersect(NGM_mean, NaCl_mean))
n14 <- length(intersect(NGM_mean, NaCl_var))
n23 <- length(intersect(NGM_var, NaCl_mean))
n24 <- length(intersect(NGM_var, NaCl_var))
n34 <- length(intersect(NaCl_mean, NaCl_var))

n123 <- length(Reduce(intersect, list(NGM_mean, NGM_var, NaCl_mean)))
n124 <- length(Reduce(intersect, list(NGM_mean, NGM_var, NaCl_var)))
n134 <- length(Reduce(intersect, list(NGM_mean, NaCl_mean, NaCl_var)))
n234 <- length(Reduce(intersect, list(NGM_var, NaCl_mean, NaCl_var)))

n1234 <- length(Reduce(intersect, list(NGM_mean, NGM_var, NaCl_mean, NaCl_var)))

# =========================================================
# 5. VENN 4 SETS
# =========================================================

grid.newpage()

draw.quad.venn(
  
  area1 = length(NGM_mean),
  area2 = length(NGM_var),
  area3 = length(NaCl_mean),
  area4 = length(NaCl_var),
  
  n12 = n12,
  n13 = n13,
  n14 = n14,
  n23 = n23,
  n24 = n24,
  n34 = n34,
  
  n123 = n123,
  n124 = n124,
  n134 = n134,
  n234 = n234,
  
  n1234 = n1234,
  
  category = c("NGM Mean", "NGM Var", "NaCl Mean", "NaCl Var"),
  
  fill = c("steelblue", "skyblue", "tomato", "orange"),
  alpha = 0.5,
  
  cex = 1.2,
  cat.cex = 1.2,
  
  cat.pos = c(-20, 20, 140, -140),
  cat.dist = c(0.05, 0.05, 0.05, 0.05)
)

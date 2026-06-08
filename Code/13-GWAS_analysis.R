library(data.table)
library(qqman)
library(ggplot2)
library(dplyr)


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









# =========================================================
# EXTRACTION DES SNPs SIGNIFICATIFS - VARIANCE NaCl
# =========================================================

# SNPs significatifs
sig_var_table <- GWAS.var.NaCl %>%
  
  filter(p.adjust_BH < 0.05) %>%
  
  mutate(
    CHR = sub("_.*", "", X_ID),
    BP  = as.numeric(sub(".*_", "", X_ID))
  ) %>%
  
  arrange(CHR, BP) %>%
  
  select(
    SNP            = X_ID,
    Chromosome     = CHR,
    Position_BP    = BP,
    beta           = beta.6,
    p_value        = p_value_REML,
    p_adjusted_BH  = p.adjust_BH,
    minus_log10_p  = log10p_adj
  )

# =========================================================
# APERCU
# =========================================================

head(sig_var_table)

# Nombre total de SNPs significatifs
nrow(sig_var_table)

# =========================================================
# RESUME PAR CHROMOSOME
# =========================================================

summary_chr <- sig_var_table %>%
  
  group_by(Chromosome) %>%
  
  summarise(
    
    n_SNPs = n(),
    
    min_position = min(Position_BP),
    
    max_position = max(Position_BP),
    
    mean_beta = mean(beta, na.rm = TRUE),
    
    .groups = "drop"
  )








##### Miroir  ######


library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)
library(cowplot)

# =========================================================
# CHROMOSOME ORDER
# =========================================================

chr_order <- c("I", "II", "III", "IV", "V", "X", "MtDNA")

# =========================================================
# COLORS
# =========================================================

col_mean <- "#5B6CFF"
col_var  <- "#E68613"
col_cv   <- "#1B9E77"

chr_cols <- c(
  "I"     = "grey20",
  "II"    = "grey55",
  "III"   = "grey20",
  "IV"    = "grey55",
  "V"     = "grey20",
  "X"     = "grey55",
  "MtDNA" = "grey20"
)

trait_cols <- c(
  "Mean"     = col_mean,
  "Variance" = col_var,
  "CV"       = col_cv
)

# =========================================================
# THEME
# =========================================================

theme_manhattan <- theme_classic(base_size = 14) +
  
  theme(
    
    axis.title = element_text(
      face = "bold",
      size = 15
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 10,
      angle = 45,
      hjust = 1
    ),
    
    axis.text.y = element_text(
      size = 11
    ),
    
    plot.title = element_text(
      face = "bold",
      size = 15,
      hjust = 0.5
    ),
    
    legend.position = "none",
    
    plot.margin = margin(8,8,8,8)
    
  )

# =========================================================
# PREP FUNCTION
# =========================================================

prepare_gwas <- function(df){
  
  df <- df %>%
    
    mutate(
      CHR = sub("_.*", "", X_ID),
      BP  = as.numeric(sub(".*_", "", X_ID))
    ) %>%
    
    filter(CHR %in% chr_order) %>%
    
    mutate(
      CHR = factor(CHR,
                   levels = chr_order)
    ) %>%
    
    arrange(CHR, BP)
  
  chr_sizes <- df %>%
    
    group_by(CHR) %>%
    
    summarise(
      max_BP = max(BP),
      .groups = "drop"
    ) %>%
    
    mutate(
      offset = cumsum(
        lag(max_BP,
            default = 0)
      )
    )
  
  df <- df %>%
    
    left_join(chr_sizes,
              by = "CHR") %>%
    
    mutate(
      BP_cum = BP + offset
    )
  
  axis_df <- df %>%
    
    group_by(CHR) %>%
    
    summarise(
      center = mean(BP_cum),
      .groups = "drop"
    )
  
  df$p.adjust_BH <- p.adjust(
    df$p_value_REML,
    method = "BH"
  )
  
  df$log10p <- -log10(df$p.adjust_BH)
  
  df$signif <- df$p.adjust_BH < 0.05
  
  list(
    data = df,
    axis = axis_df
  )
}

# =========================================================
# PREP DATA
# =========================================================

mean_NGM  <- prepare_gwas(GWAS.mean.pruned)
var_NGM   <- prepare_gwas(GWAS.var.pruned)
cv_NGM    <- prepare_gwas(GWAS.cv.pruned)

mean_NaCl <- prepare_gwas(GWAS.mean.NaCl.pruned)
var_NaCl  <- prepare_gwas(GWAS.var.NaCl.pruned)
cv_NaCl   <- prepare_gwas(GWAS.cv.NaCl.pruned)

# =========================================================
# MIRROR FUNCTION
# =========================================================

make_mirror_plot <- function(
    top_data,
    bottom_data,
    axis_df,
    top_col,
    bottom_col,
    title = NULL,
    show_y = FALSE
){
  
  bottom_data <- bottom_data %>%
    mutate(log10p = -log10p)
  
  ymax <- max(top_data$log10p, na.rm = TRUE)
  ymin <- min(bottom_data$log10p, na.rm = TRUE)
  
  ggplot() +
    
    # =====================================================
  # SNPs
  # =====================================================
  
  geom_point(
    data = top_data,
    aes(
      x = BP_cum,
      y = log10p,
      colour = CHR
    ),
    size = 0.14,
    alpha = 0.70
  ) +
    
    geom_point(
      data = bottom_data,
      aes(
        x = BP_cum,
        y = log10p,
        colour = CHR
      ),
      size = 0.14,
      alpha = 0.70
    ) +
    
    # =====================================================
  # SIGNIFICANT SNPs
  # =====================================================
  
  geom_point(
    data = top_data %>%
      filter(signif),
    aes(
      x = BP_cum,
      y = log10p
    ),
    colour = top_col,
    size = 0.35
  ) +
    
    geom_point(
      data = bottom_data %>%
        filter(signif),
      aes(
        x = BP_cum,
        y = log10p
      ),
      colour = bottom_col,
      size = 0.35
    ) +
    
    # =====================================================
  # CENTRAL LINE AT 0
  # =====================================================
  
  geom_hline(
    yintercept = 0,
    colour = "black",
    linewidth = 0.9,
    alpha = 1
  ) +
    
    # =====================================================
  # FDR THRESHOLD
  # =====================================================
  
  geom_hline(
    yintercept = c(
      -log10(0.05),
      log10(0.05)
    ),
    colour = "firebrick",
    linetype = "dashed",
    linewidth = 0.7
  ) +
    
    # =====================================================
  # AXES
  # =====================================================
  
  scale_colour_manual(
    values = chr_cols
  ) +
    
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$CHR
    ) +
    
    coord_cartesian(
      ylim = c(ymin * 1.05,
               ymax * 1.05),
      expand = FALSE
    ) +
    
    labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    
    theme_manhattan +
    
    theme(
      panel.border = element_blank()
    ) +
    
    {
      if(!show_y)
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()
        )
    }
}

# =========================================================
# PLOTS
# =========================================================

p1 <- make_mirror_plot(
  mean_NGM$data,
  var_NGM$data,
  mean_NGM$axis,
  col_mean,
  col_var,
  "Mean vs Variance",
  show_y = TRUE
)

p2 <- make_mirror_plot(
  mean_NGM$data,
  cv_NGM$data,
  mean_NGM$axis,
  col_mean,
  col_cv,
  "Mean vs CV",
  show_y = FALSE
)

p3 <- make_mirror_plot(
  var_NGM$data,
  cv_NGM$data,
  var_NGM$axis,
  col_var,
  col_cv,
  "Variance vs CV",
  show_y = FALSE
)

p4 <- make_mirror_plot(
  mean_NaCl$data,
  var_NaCl$data,
  mean_NaCl$axis,
  col_mean,
  col_var,
  NULL,
  show_y = TRUE
)

p5 <- make_mirror_plot(
  mean_NaCl$data,
  cv_NaCl$data,
  mean_NaCl$axis,
  col_mean,
  col_cv,
  NULL,
  show_y = FALSE
)

p6 <- make_mirror_plot(
  var_NaCl$data,
  cv_NaCl$data,
  var_NaCl$axis,
  col_var,
  col_cv,
  NULL,
  show_y = FALSE
)

# =========================================================
# LEFT LABELS
# =========================================================

env_top <- ggdraw() +
  draw_label(
    "NGM",
    angle = 90,
    fontface = "bold",
    size = 18
  )

env_bottom <- ggdraw() +
  draw_label(
    "NaCl",
    angle = 90,
    fontface = "bold",
    size = 18
  )

yleft_top <- ggdraw() +
  draw_label(
    expression(-log[10](p-value)),
    angle = 90,
    fontface = "bold",
    size = 14
  )

yleft_bottom <- ggdraw() +
  draw_label(
    expression(-log[10](p-value)),
    angle = 90,
    fontface = "bold",
    size = 14
  )

# =========================================================
# LEGEND
# =========================================================

legend_df <- data.frame(
  Trait = c("Mean", "Variance", "CV"),
  x = 1:3,
  y = 1
)

legend_plot <- ggplot(
  legend_df,
  aes(
    x = x,
    y = y,
    colour = Trait
  )
) +
  
  geom_point(size = 5) +
  
  scale_colour_manual(
    values = trait_cols
  ) +
  
  theme_void(base_size = 14) +
  
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(
      face = "bold",
      size = 12
    )
  )

legend <- get_legend(legend_plot)

# =========================================================
# FINAL LAYOUT
# =========================================================

top_row <- plot_grid(
  
  env_top,
  
  yleft_top,
  
  p1,
  
  p2,
  
  p3,
  
  nrow = 1,
  
  rel_widths = c(0.04, 0.05, 1, 1, 1)
  
)

bottom_row <- plot_grid(
  
  env_bottom,
  
  yleft_bottom,
  
  p4,
  
  p5,
  
  p6,
  
  nrow = 1,
  
  rel_widths = c(0.04, 0.05, 1, 1, 1)
  
)

final_plot <- plot_grid(
  
  top_row,
  
  legend,
  
  bottom_row,
  
  ncol = 1,
  
  rel_heights = c(1, 0.08, 1)
  
)

final_plot



######### Cluster #########



library(dplyr)
library(VennDiagram)
library(grid)

# =========================================================
# FUNCTION
# =========================================================

get_sig_snps <- function(df, pcol = "p.adjust_BH"){
  
  df %>%
    
    filter(!is.na(X_ID)) %>%
    
    mutate(
      X_ID = as.character(X_ID)
    ) %>%
    
    filter(.data[[pcol]] < 0.05) %>%
    
    pull(X_ID) %>%
    
    unique()
}

# =========================================================
# SNP SETS
# =========================================================

# -------------------------
# MEAN
# -------------------------

NGM_mean  <- get_sig_snps(mean_NGM$data)
NaCl_mean <- get_sig_snps(mean_NaCl$data)

# -------------------------
# VAR
# -------------------------

NGM_var   <- get_sig_snps(var_NGM$data)
NaCl_var  <- get_sig_snps(var_NaCl$data)

# -------------------------
# CV
# -------------------------

NGM_cv    <- get_sig_snps(cv_NGM$data)
NaCl_cv   <- get_sig_snps(cv_NaCl$data)

# =========================================================
# =========================================================
# VENN 1 : MEAN vs VAR
# =========================================================
# =========================================================

# =========================================================
# INTERSECTIONS
# =========================================================

n12 <- length(intersect(NGM_mean, NGM_var))
n13 <- length(intersect(NGM_mean, NaCl_mean))
n14 <- length(intersect(NGM_mean, NaCl_var))
n23 <- length(intersect(NGM_var, NaCl_mean))
n24 <- length(intersect(NGM_var, NaCl_var))
n34 <- length(intersect(NaCl_mean, NaCl_var))

n123  <- length(Reduce(intersect,
                       list(NGM_mean, NGM_var, NaCl_mean)))

n124  <- length(Reduce(intersect,
                       list(NGM_mean, NGM_var, NaCl_var)))

n134  <- length(Reduce(intersect,
                       list(NGM_mean, NaCl_mean, NaCl_var)))

n234  <- length(Reduce(intersect,
                       list(NGM_var, NaCl_mean, NaCl_var)))

n1234 <- length(Reduce(intersect,
                       list(NGM_mean,
                            NGM_var,
                            NaCl_mean,
                            NaCl_var)))

# =========================================================
# PLOT
# =========================================================

grid.newpage()

v=draw.quad.venn(
  
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
  
  category = c(
    "Mean low-salt",
    "Variance low-salt",
    "Mean high-salt",
    "Variance high-salt"
  ),
  
  fill = c(
    "#E69F00",
    "#56B4E9",
    "#F0E442",
    "#5A9BD4"
  ),
  
  alpha = 0.5,
  
  cex = 1.4,
  cat.cex = 0.9,
  
  fontface = "bold",
  cat.fontface = "bold",
  
  cat.pos = c(-20, 20, 140, -140),
  cat.dist = c(0.05, 0.05, 0.05, 0.05),
  
  lwd = 2,
  col = "black"
)



# =========================================================
# =========================================================
# VENN 2 : MEAN vs CV
# =========================================================
# =========================================================

# =========================================================
# INTERSECTIONS
# =========================================================

n12 <- length(intersect(NGM_mean, NGM_cv))
n13 <- length(intersect(NGM_mean, NaCl_mean))
n14 <- length(intersect(NGM_mean, NaCl_cv))
n23 <- length(intersect(NGM_cv, NaCl_mean))
n24 <- length(intersect(NGM_cv, NaCl_cv))
n34 <- length(intersect(NaCl_mean, NaCl_cv))

n123  <- length(Reduce(intersect,
                       list(NGM_mean, NGM_cv, NaCl_mean)))

n124  <- length(Reduce(intersect,
                       list(NGM_mean, NGM_cv, NaCl_cv)))

n134  <- length(Reduce(intersect,
                       list(NGM_mean, NaCl_mean, NaCl_cv)))

n234  <- length(Reduce(intersect,
                       list(NGM_cv, NaCl_mean, NaCl_cv)))

n1234 <- length(Reduce(intersect,
                       list(NGM_mean,
                            NGM_cv,
                            NaCl_mean,
                            NaCl_cv)))

# =========================================================
# PLOT
# =========================================================

grid.newpage()

cv=draw.quad.venn(
  
  area1 = length(NGM_mean),
  area2 = length(NGM_cv),
  area3 = length(NaCl_mean),
  area4 = length(NaCl_cv),
  
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
  
  category = c(
    "Mean low-salt",
    "CV low-salt",
    "Mean high-salt",
    "CV high-salt"
  ),
  
  fill = c(
    "#E69F00",
    "#009E73",
    "#F0E442",
    "#D55E00"
  ),
  
  alpha = 0.5,
  
  cex = 1.4,
  cat.cex = 0.9,
  
  fontface = "bold",
  cat.fontface = "bold",
  
  cat.pos = c(-20, 20, 140, -140),
  cat.dist = c(0.05, 0.05, 0.05, 0.05),
  
  lwd = 2,
  col = "black"
)



(v | cv)

shared_CV_sigma2 <- Reduce(
  intersect,
  list(NGM_cv,
       NGM_var,
       NaCl_cv,
       NaCl_var)
)

length(shared_CV_sigma2)

length(intersect(NGM_cv, NGM_var))
length(intersect(NaCl_cv, NaCl_var))




############################### manhattan definitif ######################


library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

############################################################
##################### PREP FUNCTION ########################
############################################################

prepare_gwas <- function(df, trait_name, strip_color) {
  
  # enlever MtDNA
  df <- df %>%
    filter(!grepl("^MtDNA", X_ID))
  
  # ordre des chromosomes
  chr_order <- c("I", "II", "III", "IV", "V", "X")
  
  # extraction chromosome + position
  df <- df %>%
    mutate(
      CHR = sub("_.*", "", X_ID),
      BP  = as.numeric(sub(".*_", "", X_ID))
    ) %>%
    filter(CHR %in% chr_order) %>%
    mutate(CHR = factor(CHR, levels = chr_order)) %>%
    arrange(CHR, BP)
  
  # correction BH
  df$p.adjust_BH <- p.adjust(df$p_value_REML, method = "BH")
  
  # -log10(p)
  df$log10p <- -log10(df$p.adjust_BH)
  
  # AJOUT IMPORTANT
  df$trait <- trait_name
  
  # cumulative genomic position
  chr_sizes <- df %>%
    group_by(CHR) %>%
    summarise(max_BP = max(BP), .groups = "drop") %>%
    mutate(offset = cumsum(lag(max_BP, default = 0)))
  
  df <- df %>%
    left_join(chr_sizes, by = "CHR") %>%
    mutate(BP_cum = BP + offset)
  
  # centres des chromosomes
  axis_df <- df %>%
    group_by(CHR) %>%
    summarise(center = mean(BP_cum), .groups = "drop")
  
  # SNP significatifs
  sig_snps <- df %>%
    filter(p.adjust_BH < 0.05)
  
  list(
    df = df,
    axis_df = axis_df,
    sig_snps = sig_snps,
    color = strip_color
  )
}


make_plot <- function(obj, show_x = FALSE) {
  
  p <- ggplot(obj$df,
              aes(x = BP_cum,
                  y = log10p)) +
    
    # points
    geom_point(aes(color = CHR),
               size = 0.5,
               alpha = 0.5) +
    
    # alternating chromosome colors
    scale_color_manual(
      values = rep(c("grey20", "grey70"), 3)
    ) +
    
    # SNP significatifs
    geom_point(
      data = obj$sig_snps,
      color = "red",
      size = 1.5
    ) +
    
    # seuil significatif
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dotted",
      color = "grey40",
      linewidth = 0.5
    ) +
    
    # espace entre chromosomes
    scale_x_continuous(
      breaks = obj$axis_df$center,
      labels = obj$axis_df$CHR,
      expand = expansion(mult = c(0.01, 0.03))
    ) +
    
    scale_y_continuous(
      limits = c(0, ymax + 0.5),
      expand = c(0, 0)
    ) +
    
    labs(
      x = if(show_x) "Genomic position (Mb)" else NULL,
      y = expression(-log[10](italic(p)))
    ) +
    
    facet_grid(
      rows = vars(trait),
      switch = "y"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      strip.position = "right",
      
      strip.switch.pad.grid = unit(0.1, "cm"),
      
      legend.position = "none",
      
      axis.title.x = if(show_x)
        element_text(size = 16)
      else
        element_blank(),
      
      axis.title.y = element_text(size = 16),
      
      axis.text.x = if(show_x)
        element_text(size = 12, color = "black")
      else
        element_blank(),
      
      axis.ticks.x = if(show_x)
        element_line()
      else
        element_blank(),
      
      axis.text.y = element_text(size = 11),
      
      # bordure panel
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      # BOITES COLOREES A DROITE
      strip.placement = "outside",
      
      strip.background.y = element_rect(
        fill = obj$color,
        color = "black"
      ),
      
      strip.text.y.right = element_text(
        angle = -90,
        size = 12,
        face = "plain"
      ),
      
      # marges
      plot.margin = margin(2, 2, 2, 2)
    )
  
  return(p)
}

############################################################
###################### PREP DATA ###########################
############################################################

mean_low  <- prepare_gwas(GWAS.mean.pruned,
                          "Mean low-salt",
                          "#E69F00")

var_low   <- prepare_gwas(GWAS.var.pruned,
                          "ln(\u03C3²) low-salt",
                          "#56B4E9")

cv_low    <- prepare_gwas(GWAS.cv.pruned,
                          "ln(CV) low-salt",
                          "#009E73")

mean_high <- prepare_gwas(GWAS.mean.NaCl.pruned,
                          "Mean high-salt",
                          "#F0E442")

var_high  <- prepare_gwas(GWAS.var.NaCl.pruned,
                          "ln(\u03C3²) high-salt",
                          "#5A9BD4")

cv_high   <- prepare_gwas(GWAS.cv.NaCl.pruned,
                          "ln(CV) high-salt",
                          "#D55E00")

all_gwas <- list(
  mean_low,
  var_low,
  cv_low,
  mean_high,
  var_high,
  cv_high
)

############################################################
##################### GLOBAL LIMIT #########################
############################################################

ymax <- max(
  sapply(all_gwas, function(x) max(x$df$log10p, na.rm = TRUE))
)

############################################################
###################### PLOT FUNCTION #######################
############################################################

make_plot <- function(obj, show_x = FALSE) {
  
  p <- ggplot(obj$df,
              aes(x = BP_cum,
                  y = log10p)) +
    
    # alternating chromosome colors
    geom_point(aes(color = CHR),
               size = 0.5,
               alpha = 0.5) +
    
    scale_color_manual(
      values = rep(c("grey20", "grey70"), 3)
    ) +
    
    # significant SNPs
    geom_point(
      data = obj$sig_snps,
      color = "red",
      size = 1.5
    ) +
    
    # significance threshold
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dotted",
      color = "grey40",
      linewidth = 0.5
    ) +
    
    scale_x_continuous(
      breaks = obj$axis_df$center,
      labels = obj$axis_df$CHR,
      expand = c(0.01, 0.01)
    ) +
    
    scale_y_continuous(
      limits = c(0, ymax + 0.5),
      expand = c(0, 0)
    ) +
    
    labs(
      
      x = if(show_x) "Genomic position" else NULL,
      
      y = NULL
      
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      
      legend.position = "none",
      
      axis.title.x = if(show_x)
        
        element_text(size = 16, face = "bold")
      
      else
        
        element_blank(),
      
      axis.title.y = element_text(
        
        size = 16,
        
        face = "bold"
        
      ),
      
      axis.text.x = if(show_x)
        element_text(size = 12, color = "black")
      else
        element_blank(),
      
      axis.ticks.x = if(show_x)
        element_line()
      else
        element_blank(),
      
      axis.text.y = element_text(size = 11),
      
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      strip.background = element_rect(
        fill = obj$color,
        color = "black"
      ),
      
      strip.text = element_text(
        angle = -90,
        face = "bold",
        size = 12
      ),
      
      plot.margin = margin(2, 2, 2, 2)
    ) +
    
    facet_grid(
      rows = vars(trait),
      switch = "y"
    )
  
  return(p)
}

############################################################
########################## PLOTS ###########################
############################################################

p1 <- make_plot(mean_low)
p2 <- make_plot(var_low)
p3 <- make_plot(cv_low)
p4 <- make_plot(mean_high)
p5 <- make_plot(var_high)
p6 <- make_plot(cv_high, show_x = TRUE)

############################################################
###################### FINAL FIGURE ########################
############################################################

final_plot <- (
  p1 /
    p2 /
    p3 /
    p4 /
    p5 /
    p6
) +
  
  plot_layout(heights = rep(1, 6)) &
  
  theme(
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16)
  )

final_plot



final_plot +
  plot_annotation(
    theme = theme(
      plot.margin = margin(5, 5, 5, 5)
    )
  )

wrap_elements(full = 
                grid::textGrob(
                  expression(-log[10](italic(p))),
                  rot = 90,
                  gp = grid::gpar(fontsize = 16, fontface = "bold")
                )
) | final_plot

library(cowplot)

ggdraw(final_plot) +
  
  draw_label(
    
    expression(-log[10](italic(p))),
    
    x = 0.015,
    
    y = 0.5,
    
    angle = 90,
    
    fontface = "bold",
    
    size = 16
    
  )









############################################################
###################### PLOT FUNCTION #######################
############################################################

make_plot <- function(obj,
                      show_x = FALSE,
                      show_y = FALSE) {
  
  p <- ggplot(obj$df,
              aes(x = BP_cum,
                  y = log10p)) +
    
    # SNPs
    geom_point(aes(color = CHR),
               size = 0.5,
               alpha = 0.5) +
    
    # alternating chromosome colors
    scale_color_manual(
      values = rep(c("grey20", "grey70"), 3)
    ) +
    
    # significant SNPs
    geom_point(
      data = obj$sig_snps,
      color = "red",
      size = 1.2
    ) +
    
    # significance threshold
    geom_hline(
      yintercept = -log10(0.05),
      linetype = "dotted",
      color = "grey40",
      linewidth = 0.5
    ) +
    
    # chromosome labels
    scale_x_continuous(
      breaks = obj$axis_df$center,
      labels = obj$axis_df$CHR,
      expand = c(0.01, 0.01)
    ) +
    
    # same y scale everywhere
    scale_y_continuous(
      limits = c(0, ymax + 0.5),
      expand = c(0, 0)
    ) +
    
    # labels
    labs(
      x = if(show_x)
        "Genomic position"
      else
        NULL,
      
      y = if(show_y)
        "-log10(p)"
      else
        NULL
    ) +
    
    # one facet per trait
    facet_grid(
      rows = vars(trait),
      switch = "y"
    ) +
    
    theme_classic(base_size = 14) +
    
    theme(
      
      legend.position = "none",
      
      # X axis title
      axis.title.x = if(show_x)
        element_text(size = 16,
                     face = "bold")
      else
        element_blank(),
      
      # Y axis title
      axis.title.y = if(show_y)
        element_text(size = 16,
                     face = "bold")
      else
        element_blank(),
      
      # x labels only at bottom
      axis.text.x = if(show_x)
        element_text(size = 12,
                     color = "black")
      else
        element_blank(),
      
      axis.ticks.x = if(show_x)
        element_line()
      else
        element_blank(),
      
      # y tick labels
      axis.text.y = element_text(
        size = 11,
        color = "black"
      ),
      
      # panel border
      panel.border = element_rect(
        color = "black",
        fill = NA,
        linewidth = 0.8
      ),
      
      # colored labels at right
      strip.background = element_rect(
        fill = obj$color,
        color = "black"
      ),
      
      strip.text = element_text(
        angle = -90,
        face = "bold",
        size = 12
      ),
      
      # margins
      plot.margin = margin(2, 2, 2, 2)
    )
  
  return(p)
}

############################################################
########################## PLOTS ###########################
############################################################

p1 <- make_plot(mean_low)

p2 <- make_plot(var_low)

p3 <- make_plot(cv_low,
                show_y = TRUE)

p4 <- make_plot(mean_high)

p5 <- make_plot(var_high)

p6 <- make_plot(cv_high,
                show_x = TRUE)

############################################################
###################### FINAL FIGURE ########################
############################################################

final_plot <- (
  p1 /
    p2 /
    p3 /
    p4 /
    p5 /
    p6
) +
  
  plot_layout(heights = rep(1, 6))

final_plot





############################ Effets Volcano #######################

library(data.table)
library(dplyr)
library(ggplot2)
library(patchwork)

###########################################################
# COULEURS
###########################################################

col_mean_low  <- "#E69F00"
col_var_low   <- "#56B4E9"
col_cv_low    <- "#009E73"

col_mean_high <- "#F0E442"
col_var_high  <- "#5A9BD4"
col_cv_high   <- "#D55E00"

###########################################################
# PREPARATION
###########################################################

mean_low <- GWAS.mean.pruned %>%
  mutate(
    Trait = "Mean",
    beta_SNP = beta.5,
    p.adjust_BH = p.adjust(p_value_REML, method = "BH")
  )

var_low <- GWAS.var.pruned %>%
  mutate(
    Trait = "ln(σ²)",
    beta_SNP = beta.6,
    p.adjust_BH = p.adjust(p_value_REML, method = "BH")
  )

cv_low <- GWAS.cv.pruned %>%
  mutate(
    Trait = "ln(CV)",
    beta_SNP = beta.5,
    p.adjust_BH = p.adjust(p_value_REML, method = "BH")
  )

volcano_low <- bind_rows(
  mean_low,
  var_low,
  cv_low
)

###########################################################

mean_high <- GWAS.mean.NaCl.pruned %>%
  mutate(
    Trait = "Mean",
    beta_SNP = beta.5,
    p.adjust_BH = p.adjust(p_value_REML, method = "BH")
  )

var_high <- GWAS.var.NaCl.pruned %>%
  mutate(
    Trait = "ln(σ²)",
    beta_SNP = beta.6,
    p.adjust_BH = p.adjust(p_value_REML, method = "BH")
  )

cv_high <- GWAS.cv.NaCl.pruned %>%
  mutate(
    Trait = "ln(CV)",
    beta_SNP = beta.5,
    p.adjust_BH = p.adjust(p_value_REML, method = "BH")
  )

volcano_high <- bind_rows(
  mean_high,
  var_high,
  cv_high
)

###########################################################
# LOG10(P)
###########################################################

volcano_low <- volcano_low %>%
  mutate(
    log10p = -log10(p.adjust_BH),
    Significant = p.adjust_BH < 0.05
  )

volcano_high <- volcano_high %>%
  mutate(
    log10p = -log10(p.adjust_BH),
    Significant = p.adjust_BH < 0.05
  )

###########################################################
# VOLCANO FUNCTION
###########################################################

make_volcano <- function(df,
                         title,
                         colours){
  
  threshold <- -log10(0.05)
  
  ggplot(df,
         aes(x = beta_SNP,
             y = log10p,
             colour = Trait)) +
    
    geom_point(
      data = subset(df, !Significant),
      alpha = 0.25,
      size = 0.7
    ) +
    
    geom_point(
      data = subset(df, Significant),
      alpha = 1,
      size = 2.2
    ) +
    
    geom_hline(
      yintercept = threshold,
      linetype = "dashed",
      colour = "grey40"
    ) +
    
    scale_colour_manual(
      values = colours
    ) +
    
    guides(
      colour = guide_legend(
        override.aes = list(
          size = 5,
          alpha = 1
        )
      )
    ) +
    
    labs(
      x = expression(bold("SNP effect size (" * beta * ")")),
      y = expression(bold(-log[10](italic(p)))),
      title = title,
      colour = NULL
    ) +
    
    theme_classic(base_size = 16) +
    
    theme(
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      ),
      
      legend.position = "top",
      
      legend.text = element_text(
        face = "bold",
        size = 14
      ),
      
      axis.title = element_text(
        face = "bold"
      )
    )
}

###########################################################
# PLOTS
###########################################################

p_low <- make_volcano(
  volcano_low,
  "Low-salt",
  c(
    "Mean"   = col_mean_low,
    "ln(σ²)" = col_var_low,
    "ln(CV)" = col_cv_low
  )
)

p_high <- make_volcano(
  volcano_high,
  "High-salt",
  c(
    "Mean"   = col_mean_high,
    "ln(σ²)" = col_var_high,
    "ln(CV)" = col_cv_high
  )
)

###########################################################
# FIGURE FINALE
###########################################################

p_low | p_high







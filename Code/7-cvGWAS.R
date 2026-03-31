## ===============================
## cvGWAS (log)
## ===============================

pheno_CViso <- fread("data/gemma_pheno_withLines_CV.txt")
SNP_perm_CV =  fread("output/perm_log_CV_pwald_matrix.txt")
SNP_CV = fread("output/gwas_log_CV_results.assoc.txt")


str(pheno_CViso)
str(SNP_perm_CV)
str(SNP_CV)

SNP_CV[, chr := sub("_.*", "", rs)]
SNP_CV[, ps  := as.numeric(sub(".*_", "", rs))]
str(SNP_CV)

all(SNP_CV$rs == SNP_perm_CV$V1) #correspondnace des SNP


perm_matrix <- as.matrix(SNP_perm_CV[, -1, with = FALSE])

perm_threshold <- apply(perm_matrix, 1, function(x) {
  quantile(x, 0.05, na.rm = TRUE)
})

SNP_CV[, threshold_perm := perm_threshold]
SNP_CV[, significant_perm := p_wald < threshold_perm]


SNP_CV[, chr := factor(chr, levels = c("I","II","III","IV","V","X", "MtDNA"))]
SNP_CV[, logp := -log10(p_wald)]


SNP_CV <- SNP_CV[order(chr, ps)] #ordonner SNP
SNP_CV[, pos_cum := ps + cumsum(c(0, head(tapply(ps, chr, max), -1)))[chr]]


#manhattan plot 
ggplot(SNP_CV, aes(x = pos_cum, y = logp)) +
  
  geom_point(aes(color = chr), alpha = 0.7, size = 0.5) +
  
  geom_point(
    data = SNP_CV[significant_perm == TRUE],
    color = "#D62728",
    size = 1.2
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Genomic position",
    y = expression(-log[10](italic(p))),
    title = "Genome-wide association analysis"
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    
    panel.grid = element_blank()
  )


# Couleurs
chrom_colors <- c("#B2DF8A", "#CAB2D6")  # alternance vert/violet pastel
sig_color <- "#E31A1C"  # rouge visible pour les significatifs

# Ajouter un espace entre chromosomes
SNP_CV <- SNP_CV %>%
  arrange(chr, pos_cum) %>%
  group_by(chr) %>%
  mutate(pos_cum_chr = pos_cum + (as.numeric(chr)-1) * 1e6) %>%
  ungroup()

chrom_labels <- SNP_CV %>%
  group_by(chr) %>%
  summarise(
    pos_label = (min(pos_cum_chr) + max(pos_cum_chr)) / 2
  ) %>%
  ungroup()

ggplot() +
  # SNP non significatifs
  geom_point(
    data = SNP_CV %>% filter(!significant_perm),
    aes(x = pos_cum_chr, y = logp, color = factor(chr)),
    size = 0.6,
    alpha = 0.3
  ) +
  # SNP significatifs
  geom_point(
    data = SNP_CV %>% filter(significant_perm),
    aes(x = pos_cum_chr, y = logp),
    color = sig_color,
    size = 0.6,
    alpha = 0.6
  ) +
  scale_color_manual(values = rep(chrom_colors, length.out = length(unique(SNP_CV$chr)))) +
  scale_x_continuous(
    breaks = chrom_labels$pos_label,
    labels = chrom_labels$chr
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "GWAS coefficient of variation"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    panel.grid = element_blank()
  )




ggplot(SNP_CV, aes(x = pos_cum, y = beta)) +
  
  # SNP non significatifs (fond)
  geom_point(
    data = SNP_CV[significant_perm == FALSE],
    color = "grey80",
    size = 0.5,
    alpha = 0.6
  ) +
  
  # SNP significatifs (par-dessus)
  geom_point(
    data = SNP_CV[significant_perm == TRUE],
    aes(color = chr),
    size = 1.2,
    alpha = 0.9
  ) +
  
  # ligne effet nul
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Genomic position",
    y = "SNP effect size (beta)",
    title = "Genetic effects on body size CV"
  ) +
  
  scale_color_manual(
    values = rep(c("#F28E2B", "#D55E00"), length.out = length(levels(SNP_CV$chr)))
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    
    panel.grid = element_blank()
  )




####### pour chr 1

# Filtrer uniquement le chromosome 1
SNP_chr1 <- SNP_CV[chr == "II"]

ggplot(SNP_chr1, aes(x = pos_cum, y = logp)) +
  
  # Tous les SNPs
  geom_point(aes(color = chr), alpha = 0.7, size = 0.5) +
  
  # SNPs significatifs
  geom_point(
    data = SNP_chr1[significant_perm == TRUE],
    color = "#D62728",
    size = 1.2
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Genomic position",
    y = expression(-log[10](italic(p))),
    title = "GWAS - Chromosome 1"
  ) +
  
  theme(
    legend.position = "none",
    
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    
    panel.grid = element_blank()
  )



ggplot(SNP_chr1, aes(x = pos_cum, y = beta)) +
  
  # SNP non significatifs
  geom_point(
    data = SNP_chr1[significant_perm == FALSE],
    color = "grey80",
    size = 0.5,
    alpha = 0.6
  ) +
  
  # SNP significatifs
  geom_point(
    data = SNP_chr1[significant_perm == TRUE],
    color = "#D55E00",  # une seule couleur suffit ici
    size = 1.2,
    alpha = 0.9
  ) +
  
  # ligne effet nul
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Position on chromosome II",
    y = "SNP effect size (beta)",
    title = "Genetic effects on body size variance (Chr II)"
  ) +
  
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    panel.grid = element_blank()
  )







################################################# calcul avec LD


library(data.table)

# S'assurer que c'est un data.table
setDT(SNP_CV)

# Filtrer chr II + effet positif
lead_snp_chrII_pos <- SNP_CV[chr == "II" & beta > 0][order(p_wald)][1]

lead_snp_chrII_pos

### Plink step ##


ld <- read.table("/Users/edouardcourret/Documents/GitHub/2026_courret/data/LD/ld_chrII.ld", header = TRUE)
str(ld)


library(dplyr)

SNP_CV_ld <- SNP_CV %>%
  left_join(ld, by = c("rs" = "SNP_B"))

library(dplyr)

SNP_CV_ld <- SNP_CV %>%
  left_join(ld %>% select(SNP_B, R2), by = c("rs" = "SNP_B"))

library(ggplot2)
library(dplyr)
library(scales)  # pour alpha()

# Définir la couleur des SNP significatifs
sig_color_palette <- c("low" = "#FEE0D2", "high" = "#D55E00")  # du pastel au rouge intense

ggplot() +
  # SNP non significatifs sur Chr II
  geom_point(
    data = SNP_CV_ld %>% filter(!significant_perm, chr == "II"),
    aes(x = pos_cum, y = beta),
    color = "grey80",
    size = 0.5,
    alpha = 0.6
  ) +
  # SNP significatifs sur Chr II, gradient selon LD
  geom_point(
    data = SNP_CV_ld %>% filter(significant_perm, chr == "II"),
    aes(x = pos_cum, y = beta, color = R2),
    size = 1.2,
    alpha = 0.9
  ) +
  # Ligne effet nul
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 0.5) +
  # Gradient de couleur rouge pour LD
  scale_color_gradient(low = "#FEE0D2", high = "#D55E00", na.value = "#FEE0D2") +
  
  theme_classic(base_size = 15) +
  labs(
    x = "Position on Chromosome II",
    y = "SNP effect size (beta)",
    color = expression(LD~(r^2)),
    title = "Genetic effects on body size variance — Chr II"
  ) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    panel.grid = element_blank()
  )



library(ggplot2)
library(dplyr)
library(scales)

# Définir la couleur des SNP significatifs (gradient LD)
sig_color_palette <- c("low" = "#FEE0D2", "high" = "#D55E00")  # pastel → rouge

ggplot() +
  # SNP non significatifs sur Chr II
  geom_point(
    data = SNP_CV_ld %>% filter(!significant_perm, chr == "II"),
    aes(x = pos_cum, y = -log10(p_wald)),
    color = "grey80",
    size = 0.5,
    alpha = 0.6
  ) +
  # SNP significatifs sur Chr II, gradient selon LD
  geom_point(
    data = SNP_CV_ld %>% filter(significant_perm, chr == "II"),
    aes(x = pos_cum, y = -log10(p_wald), color = R2),
    size = 1.2,
    alpha = 0.9
  ) +
  # Gradient de couleur rouge pour LD
  scale_color_gradient(low = "#FEE0D2", high = "#D55E00", na.value = "#FEE0D2") +
  
  theme_classic(base_size = 15) +
  labs(
    x = "Position on Chromosome II",
    y = expression(-log[10](italic(p))),
    color = expression(LD~(r^2)),
    title = "GWAS — Chr II with LD gradient"
  ) +
  theme(
    legend.position = "right",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    axis.line = element_line(linewidth = 0.6),
    axis.ticks = element_line(linewidth = 0.6),
    panel.grid = element_blank()
  )




#volcano plots 

ggplot(SNP_CV, aes(x = beta, y = logp)) +
  
  geom_point(
    data = SNP_CV[significant_perm == FALSE],
    color = "grey80",
    size = 0.6,
    alpha = 0.6
  ) +
  
  geom_point(
    data = SNP_CV[significant_perm == TRUE],
    color = "#6A3D9A",
    size = 1.2,
    alpha = 0.9
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Effect size (beta)",
    y = expression(-log[10](italic(p))),
    title = "Volcano plot — CV GWAS"
  )

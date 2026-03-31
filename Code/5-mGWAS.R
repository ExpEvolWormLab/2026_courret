## ===============================
## mGWAS
## ===============================

pheno_meaniso <- fread("data/gemma_pheno_withLines.txt")
SNP_perm_mean =  fread("output/perm_pwald_matrix.txt")
SNP_mean = fread("output/gwas_results.assoc.txt")


str(pheno_meaniso)
str(SNP_perm_mean)
str(SNP_mean)

SNP_mean[, chr := sub("_.*", "", rs)]
SNP_mean[, ps  := as.numeric(sub(".*_", "", rs))]
str(SNP_mean)

all(SNP_mean$rs == SNP_perm_mean$V1) #correspondnace des SNP


perm_matrix <- as.matrix(SNP_perm_mean[, -1, with = FALSE])

perm_threshold <- apply(perm_matrix, 1, function(x) {
  quantile(x, 0.05, na.rm = TRUE)
})

SNP_mean[, threshold_perm := perm_threshold]
SNP_mean[, significant_perm := p_wald < threshold_perm]


SNP_mean[, chr := factor(chr, levels = c("I","II","III","IV","V","X", "MtDNA"))]
SNP_mean[, logp := -log10(p_wald)]


SNP_mean <- SNP_mean[order(chr, ps)] #ordonner SNP
SNP_mean[, pos_cum := ps + cumsum(c(0, head(tapply(ps, chr, max), -1)))[chr]]


#manhattan plot 
ggplot(SNP_mean, aes(x = pos_cum, y = logp)) +
  
  geom_point(aes(color = chr), alpha = 0.7, size = 0.5) +
  
  geom_point(
    data = SNP_mean[significant_perm == TRUE],
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
SNP_mean <- SNP_mean %>%
  arrange(chr, pos_cum) %>%
  group_by(chr) %>%
  mutate(pos_cum_chr = pos_cum + (as.numeric(chr)-1) * 1e6) %>%
  ungroup()

chrom_labels <- SNP_mean %>%
  group_by(chr) %>%
  summarise(
    pos_label = (min(pos_cum_chr) + max(pos_cum_chr)) / 2
  ) %>%
  ungroup()

ggplot() +
  # SNP non significatifs
  geom_point(
    data = SNP_mean %>% filter(!significant_perm),
    aes(x = pos_cum_chr, y = logp, color = factor(chr)),
    size = 0.6,
    alpha = 0.3
  ) +
  # SNP significatifs
  geom_point(
    data = SNP_mean %>% filter(significant_perm),
    aes(x = pos_cum_chr, y = logp),
    color = sig_color,
    size = 0.6,
    alpha = 0.6
  ) +
  scale_color_manual(values = rep(chrom_colors, length.out = length(unique(SNP_mean$chr)))) +
  scale_x_continuous(
    breaks = chrom_labels$pos_label,
    labels = chrom_labels$chr
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "GWAS mean"
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



ggplot(SNP_mean, aes(x = pos_cum, y = beta)) +
  
  # SNP non significatifs (fond)
  geom_point(
    data = SNP_mean[significant_perm == FALSE],
    color = "grey80",
    size = 0.5,
    alpha = 0.6
  ) +
  
  # SNP significatifs (par-dessus)
  geom_point(
    data = SNP_mean[significant_perm == TRUE],
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
    title = "Genetic effects on body size mean"
  ) +
  
  scale_color_manual(
    values = rep(c("#F28E2B", "#D55E00"), length.out = length(levels(SNP_mean$chr)))
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
SNP_chr1 <- SNP_mean[chr == "I"]

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




#volcano plots 

library(ggplot2)

ggplot(SNP_mean, aes(x = beta, y = logp)) +
  
  # SNP non significatifs
  geom_point(
    data = SNP_mean[significant_perm == FALSE],
    color = "grey80",
    size = 0.6,
    alpha = 0.6
  ) +
  
  # SNP significatifs
  geom_point(
    data = SNP_mean[significant_perm == TRUE],
    color = "#1F78B4",
    size = 1.2,
    alpha = 0.9
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Effect size (beta)",
    y = expression(-log[10](italic(p))),
    title = "Volcano plot — Mean GWAS"
  ) +
  
  theme(
    legend.position = "none",
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    plot.title = element_text(face = "bold", size = 15, hjust = 0)
  )





#################### manhattan miroir #############################

library(data.table)
library(dplyr)

# SNP significatifs uniquement
mean_sig <- SNP_mean[significant_perm == TRUE]
var_sig  <- SNP_var[significant_perm == TRUE]

# On garde les colonnes utiles
mean_sig <- mean_sig[, .(rs, chr, pos_cum, logp)]
var_sig  <- var_sig[, .(rs, chr, pos_cum, logp)]

# Mean au-dessus (positif)
mean_sig[, logp_plot := logp]

# Variance en dessous (négatif)
var_sig[, logp_plot := -logp]

# Ajouter un label de type
mean_sig[, trait := "Mean"]
var_sig[, trait := "Variance"]

# Combiner
SNP_mirror <- rbind(mean_sig, var_sig)


library(ggplot2)

ggplot(SNP_mirror, aes(x = pos_cum, y = logp_plot)) +
  
  # Mean (haut)
  geom_point(
    data = SNP_mirror %>% filter(trait == "Mean"),
    color = "#1F78B4",   # bleu propre
    size = 1,
    alpha = 0.8
  ) +
  
  # Variance (bas)
  geom_point(
    data = SNP_mirror %>% filter(trait == "Variance"),
    color = "#E31A1C",   # rouge propre
    size = 1,
    alpha = 0.8
  ) +
  
  # Ligne centrale (0)
  geom_hline(yintercept = 0, color = "black", linewidth = 0.6) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Genomic position",
    y = expression(paste("Mean  ", -log[10](italic(p)), "   |   Variance")),
    title = "GWAS comparison: Mean vs Variance"
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




# SNP significatifs dans chaque analyse
mean_sig <- SNP_mean[significant_perm == TRUE]
var_sig  <- SNP_var[significant_perm == TRUE]

# Trouver les SNP en commun
common_snps <- intersect(mean_sig$rs, var_sig$rs)

# Garder seulement ceux-là
SNP_common <- SNP_mean[rs %in% common_snps]

# Ajouter logp (au cas où)
SNP_common[, logp := -log10(p_wald)]


library(ggplot2)

ggplot(SNP_common, aes(x = pos_cum, y = logp)) +
  
  geom_point(
    color = "#6A3D9A",  # violet propre (différent de tes autres figures)
    size = 1.2,
    alpha = 0.9
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Genomic position",
    y = expression(-log[10](italic(p))),
    title = "Shared genetic architecture: Mean & Variance"
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



ggplot() +
  # Tous les SNP en fond
  geom_point(
    data = SNP_mean,
    aes(x = pos_cum, y = logp),
    color = "grey85",
    size = 0.4,
    alpha = 0.5
  ) +
  
  # SNP communs (highlight)
  geom_point(
    data = SNP_common,
    aes(x = pos_cum, y = logp),
    color = "#6A3D9A",
    size = 1.3,
    alpha = 1
  ) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = "Genomic position",
    y = expression(-log[10](italic(p))),
    title = "Shared SNPs between Mean and Variance GWAS"
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

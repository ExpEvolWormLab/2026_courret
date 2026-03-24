## ===============================
## vGWAS (log)
## ===============================

library(ggplot2)

pheno_variso <- fread("data/gemma_pheno_withLines_var.txt")
SNP_perm_var =  fread("output/perm_log_var_pwald_matrix.txt")
SNP_var = fread("output/gwas_log_var_results.assoc.txt")


str(pheno_variso)
str(SNP_perm_var)
str(SNP_var)

SNP_var[, chr := sub("_.*", "", rs)]
SNP_var[, ps  := as.numeric(sub(".*_", "", rs))]
str(SNP_var)

all(SNP_var$rs == SNP_perm_var$V1) #correspondnace des SNP


perm_matrix <- as.matrix(SNP_perm_var[, -1, with = FALSE])

perm_threshold <- apply(perm_matrix, 1, function(x) {
  quantile(x, 0.05, na.rm = TRUE)
})

SNP_var[, threshold_perm := perm_threshold]
SNP_var[, significant_perm := p_wald < threshold_perm]


SNP_var[, chr := factor(chr, levels = c("I","II","III","IV","V","X", "MtDNA"))]
SNP_var[, logp := -log10(p_wald)]


SNP_var <- SNP_var[order(chr, ps)] #ordonner SNP
SNP_var[, pos_cum := ps + cumsum(c(0, head(tapply(ps, chr, max), -1)))[chr]]


#manhattan plot 
ggplot(SNP_var, aes(x = pos_cum, y = logp)) +
  
  geom_point(aes(color = chr), alpha = 0.7, size = 0.5) +
  
  geom_point(
    data = SNP_var[significant_perm == TRUE],
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


library(ggplot2)
library(dplyr)

# Couleurs
chrom_colors <- c("#B2DF8A", "#CAB2D6")  # alternance vert/violet pastel
sig_color <- "#E31A1C"  # rouge visible pour les significatifs

# Ajouter un espace entre chromosomes
SNP_var <- SNP_var %>%
  arrange(chr, pos_cum) %>%
  group_by(chr) %>%
  mutate(pos_cum_chr = pos_cum + (as.numeric(chr)-1) * 1e6) %>%
  ungroup()

ggplot() +
  # SNP non significatifs
  geom_point(
    data = SNP_var %>% filter(!significant_perm),
    aes(x = pos_cum_chr, y = logp, color = factor(chr)),
    size = 0.6,
    alpha = 0.4
  ) +
  # SNP significatifs
  geom_point(
    data = SNP_var %>% filter(significant_perm),
    aes(x = pos_cum_chr, y = logp),
    color = sig_color,
    size = 0.6,
    alpha = 0.6
  ) +
  scale_color_manual(values = rep(chrom_colors, length.out = length(unique(SNP_var$chr)))) +
  theme_classic(base_size = 15) +
  labs(
    x = "Genomic position",
    y = expression(-log[10](italic(p))),
    title = "GWAS Manhattan Plot — Significant SNPs Highlighted"
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



chrom_labels <- SNP_var %>%
  group_by(chr) %>%
  summarise(
    pos_label = (min(pos_cum_chr) + max(pos_cum_chr)) / 2
  ) %>%
  ungroup()

ggplot() +
  # SNP non significatifs
  geom_point(
    data = SNP_var %>% filter(!significant_perm),
    aes(x = pos_cum_chr, y = logp, color = factor(chr)),
    size = 0.6,
    alpha = 0.3
  ) +
  # SNP significatifs
  geom_point(
    data = SNP_var %>% filter(significant_perm),
    aes(x = pos_cum_chr, y = logp),
    color = sig_color,
    size = 0.6,
    alpha = 0.6
  ) +
  scale_color_manual(values = rep(chrom_colors, length.out = length(unique(SNP_var$chr)))) +
  scale_x_continuous(
    breaks = chrom_labels$pos_label,
    labels = chrom_labels$chr
  ) +
  theme_classic(base_size = 15) +
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "GWAS variance"
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




ggplot(SNP_var, aes(x = pos_cum, y = beta)) +
  
  # SNP non significatifs (fond)
  geom_point(
    data = SNP_var[significant_perm == FALSE],
    color = "grey80",
    size = 0.5,
    alpha = 0.6
  ) +
  
  # SNP significatifs (par-dessus)
  geom_point(
    data = SNP_var[significant_perm == TRUE],
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
    title = "Genetic effects on body size variance"
  ) +
  
  scale_color_manual(
    values = rep(c("#F28E2B", "#D55E00"), length.out = length(levels(SNP_var$chr)))
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
SNP_chr1 <- SNP_var[chr == "X"]

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
























































##########################################################################################

gemma_var_results <- read.table("output/gwas_var_results.assoc.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)

str(gemma_var_results)
dim(gemma_var_results)

gemma_var_results <- gemma_var_results %>%
  mutate(
    chr = str_extract(rs, "^[^_]+"),   # tout avant le _
    ps  = as.numeric(str_extract(rs, "(?<=_).*"))  # tout après le _
  )

gemma_var_results$chr <- factor(
  gemma_var_results$chr,
  levels = c("I","II","III","IV","V","X", "MtDNA")
)

plot(gemma_var_results$beta,-log10(gemma_var_results$p_wald), pch=16, cex=0.5)

observed <- sort(gemma_var_results$p_wald)
expected <- ppoints(length(observed))

qqplot(-log10(expected),
       -log10(observed),
       pch = 16, cex = 0.6,
       xlab = "Expected -log10(p)",
       ylab = "Observed -log10(p)")
abline(0,1,col="red")


ggplot(gemma_var_results, aes(x = af)) +
  geom_histogram(bins = 100) +
  theme_classic()

ggplot(gemma_var_results,
       aes(x = af, y = beta)) +
  geom_point(alpha = 0.5, size = 0.8) +
  theme_classic() +
  labs(
    x = "Allele frequency",
    y = "Effect size (beta)"
  )


ggplot(gemma_var_results, aes(x = ps, y = -log10(p_wald))) +
  geom_point(aes(color = chr), size = 0.6, alpha = 0.8) +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  theme_classic(base_size = 14) +
  labs(
    x = "Position",
    y = expression(-log[10](p))
  ) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 45, hjust = 1, size=8)) +
  geom_hline(yintercept = bonf_threshold, linetype = "dashed",color = "red")


# Calcul du seuil Bonferroni
bonf_threshold <- 0.05 / nrow(gemma_var_results)

ggplot(gemma_var_results, aes(x = ps, y = -log10(p_wald))) +
  # points colorés selon p-value brute
  geom_point(aes(color = p_wald < 0.05), size = 0.6, alpha = 0.8) +
  # facettes par chromosome
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  # couleurs
  scale_color_manual(values = c("TRUE" = "red",
                                "FALSE" = "grey70")) +
  # ligne seuil Bonferroni
  geom_hline(yintercept = -log10(bonf_threshold), linetype = "dashed", color = "red") +
  theme_classic(base_size = 14) +
  labs(
    x = "Position",
    y = expression(-log[10](p))
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )



# Calcul des p-values ajustées BH
gemma_var_results$p_adj_BH <- p.adjust(gemma_var_results$p_wald,
                                       method = "BH")

# Seuil FDR 5%
fdr_threshold <- 0.05

ggplot(gemma_var_results, aes(x = ps, y = -log10(p_wald))) +
  
  # Points colorés selon FDR BH
  geom_point(aes(color = p_adj_BH < fdr_threshold),
             size = 0.6, alpha = 0.8) +
  
  # Facettes par chromosome
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  
  # Couleurs
  scale_color_manual(values = c("TRUE" = "red",
                                "FALSE" = "grey70")) +
  
  theme_classic(base_size = 14) +
  
  labs(
    x = "Position",
    y = expression(-log[10](p))
  ) +
  
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8)
  )
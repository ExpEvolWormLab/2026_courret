## ===============================
## cvGWAS (log)
## ===============================

pheno_CViso <- fread("data/gemma_pheno_withLines_CV.txt")
SNP_perm_CV =  fread("output/perm_pwald_CV_matrix.txt")
SNP_CV = fread("output/gwas_CV_results.assoc.txt")


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
    title = "Genetic effects on body size variance"
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
SNP_chr1 <- SNP_CV[chr == "I"]

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
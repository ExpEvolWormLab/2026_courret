## ===============================
## GWAS
## ===============================


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
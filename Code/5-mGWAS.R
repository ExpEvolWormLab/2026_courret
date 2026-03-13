## ===============================
## GWAS
## ===============================


gemma_results <- read.table("output/gwas_results.assoc.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)

str(gemma_results)
dim(gemma_results)

gemma_results <- gemma_results %>%
  mutate(
    chr = str_extract(rs, "^[^_]+"),   # tout avant le _
    ps  = as.numeric(str_extract(rs, "(?<=_).*"))  # tout après le _
  )

gemma_results$chr <- factor(
  gemma_results$chr,
  levels = c("I","II","III","IV","V","X", "MtDNA")
)


ggplot(gemma_results, aes(x = ps, y = -log10(p_wald))) +
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

ggplot(gemma_results, aes(x = ps, y = -log10(p_wald))) +
  geom_point(aes(color = p_wald < 0.05), size = 0.6, alpha = 0.8) +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
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

# Calcul du seuil Bonferroni
bonf_threshold <- 0.05 / nrow(gemma_results)

ggplot(gemma_results, aes(x = ps, y = -log10(p_wald))) +
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



observed <- sort(gemma_results$p_wald)
expected <- ppoints(length(observed))

qqplot(-log10(expected),
       -log10(observed),
       pch = 16, cex = 0.6,
       xlab = "Expected -log10(p)",
       ylab = "Observed -log10(p)")
abline(0,1,col="red")

hist(gemma_results$beta, breaks=50)

plot(gemma_results$beta,-log10(gemma_results$p_wald), pch=16, cex=0.5)

ggplot(gemma_results,
       aes(x = af, y = beta)) +
  geom_point(alpha = 0.5, size = 0.8) +
  theme_classic() +
  labs(
    x = "Allele frequency",
    y = "Effect size (beta)"
  )

ggplot(gemma_results, aes(x = af)) +
  geom_histogram(bins = 100) +
  theme_classic()





# Calcul des p-values ajustées BH
gemma_results$p_adj_BH <- p.adjust(gemma_results$p_wald,
                                   method = "BH")

# Seuil FDR 5%
fdr_threshold <- 0.05

ggplot(gemma_results, aes(x = ps, y = -log10(p_wald))) +
  
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



Meff <- length(readLines("plink_pruning/pruning01/CeMEE.pruned.0.1.prune.in"))
eigen_threshold <- 0.05 / Meff #bonferroni mais avec SNP prunées


ggplot(gemma_results, aes(x = ps, y = -log10(p_wald))) +
  geom_point(aes(color = p_wald < eigen_threshold),
             size = 0.6, alpha = 0.8) +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("TRUE"="red", "FALSE"="grey70")) +
  geom_hline(yintercept = -log10(eigen_threshold),
             linetype = "dashed",
             color = "blue") +
  theme_classic()

ggplot(gemma_var_results, aes(x = ps, y = -log10(p_wald))) +
  geom_point(aes(color = p_wald < eigen_threshold),
             size = 0.6, alpha = 0.8) +
  facet_wrap(~ chr, scales = "free_x", nrow = 1) +
  scale_color_manual(values = c("TRUE"="red", "FALSE"="grey70")) +
  geom_hline(yintercept = -log10(eigen_threshold),
             linetype = "dashed",
             color = "blue") +
  theme_classic()




###################### permutations 

gemma_results <- read.table("output/gwas_results.assoc.txt", header = TRUE, sep = "", stringsAsFactors = FALSE)

real_p <- gemma_results$p_wald
snp_ids <- gemma_results$rs
n_snp <- length(real_p)

n_perm <- 100
perm_matrix <- matrix(NA, nrow = n_snp, ncol = n_perm)
colnames(perm_matrix) <- paste0("perm_", 1:n_perm)
rownames(perm_matrix) <- snp_ids

























perm_threshold_snp <- apply(perm_matrix, 1, quantile, probs = 0.05)

significant_snp <- real_p < perm_threshold_snp

final_table <- data.frame(
  SNP = snp_ids,
  p_real = real_p,
  perm_threshold_5perc = perm_threshold_snp,
  significant = significant_snp
)

write.table(final_table,
            "permutation_results_100.txt",
            sep="\t",
            row.names=FALSE,
            quote=FALSE)
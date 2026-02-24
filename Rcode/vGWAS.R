############################################################
## GWAS – area.forward
## C. elegans
############################################################

## ===============================
## 0. Packages
## ===============================

library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(vcfR)
library(broom)
library(qqman)
library(car)
library(performance)
library(forcats)

## ===============================
## 1.2 Visualisation brute
## ===============================

#Variance 
ggplot(phen_data_step1, aes(x = area.F, y = log(area.F.var))) +
  geom_point(alpha = 0.3) +
  theme_classic() 
ggplot(phen_data_step1, aes(x = area.F, y = area.F.var)) +
  geom_point(alpha = 0.3) +
  theme_classic() 

ggplot(phen_data_step1, aes(x = area.F.mean_corr, y = area.F.var)) +
  geom_point(alpha = 0.3) +
  theme_classic() 

#logD
ggplot(phen_data_step1, aes(x = logD, y = area.F.var)) +
  geom_point(alpha = 0.6) +
  theme_classic() 

#Temperature 
ggplot(phen_data_step1, aes(x = temperature, y = area.F.var)) +
  geom_point(alpha = 0.6) +
  theme_classic() 

#Relative humidity
ggplot(phen_data_step1, aes(x = rel_humidity, y = area.F.var)) +
  geom_point(alpha = 0.6) +
  theme_classic() 

#pop
ggplot(phen_data_step1, aes(x = pop_label, y = area.F.var)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "area.F par population",
       x = "Population",
       y = "area.F.var") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Block
ggplot(phen_data_step1, aes(x = block_id, y = area.F.var)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "area.F par block",
       x = "Block",
       y = "area.F.var") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Session
ggplot(phen_data_step1,
       aes(x = session_id, y = area.F.var, color = session_id)) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_classic() 

## ===============================
## 1.3 Modèle environnemental
## ===============================

#1 modele mean

var_model <- lmer(
  area.F.var ~
    area.F +
    logD +
    temperature +
    rel_humidity +
    pop_label +
    (1 | session_id/block_id),
  data = phen_data_step1,
  REML = TRUE
)

summary(var_model)
Anova(var_model)

check_model(var_model)

#2 modele corrected mean 

phen_data_step1$area.F.mean_corr <- predict(mean_model)

var_model_corr <- lmer(
  area.F.var ~
    area.F.mean_corr +
    logD +
    temperature +
    rel_humidity +
    pop_label +
    (1 | session_id/block_id),
  data = phen_data_step1,
  REML = TRUE
)

summary(var_model_corr)
Anova(var_model_corr)

check_model(var_model_corr)

phen_data_step1$area.F.var_corr2 <- predict(var_model_corr) # mean corr
phen_data_step1$area.F.var_corr = predict(var_model) #mean non corr
plot(phen_data_step1$area.F.var_corr, phen_data_step1$area.F.var_corr2)

plot(phen_data_step1$area.F.var, phen_data_step1$area.F.var_corr)

sd(phen_data_step1$area.F.var_corr)
sd(phen_data_step1$area.F.var)

abline(v=mean(phen_data_step1$area.F.var)+3*mean(phen_data_step1$area.F.var))



#3 coefficient of variation (CV)

phen_data_step1$CV_areaF <- (sqrt(phen_data_step1$area.F.var_corr) / phen_data_step1$area.F.mean_corr)*100

ggplot(phen_data_step1, aes(x = area.F.mean_corr, y = sqrt(area.F.var_corr))) +
  geom_point(alpha = 0.3) +
  theme_classic() 

## ===============================
## 1.4 Extraction du phénotype final
## ===============================

## Résidus = variance intrinsèque
phen_data_step1$areaF_var_resid <- resid(var_model)

## Diagnostic post-correction
ggplot(phen_data_step1, aes(x = areaF_var_resid)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "Residual moy distribution")

ggplot(phen_data_step1,
       aes(x = area.F, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = logD, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = temperature, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = rel_humidity, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = pop_label, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = block_id, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = session_id, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

## ===============================
## 2. EMMEANS
## ===============================

pheno_line_var <- emmeans(var_model, "pop_label")
pheno_var <- as.data.frame(pheno_line_var)
dim(pheno_var)

write.csv(pheno_var,
          file = "pheno_var.csv",
          row.names = FALSE)

write.table(pheno_var,
            file = "pheno_var.txt",
            row.names = FALSE,
            col.names = T,
            quote = FALSE)

pheno_var <- read.table("data/pheno_var.txt")

pheno_var <- pheno_var %>%
  mutate(
    population = ifelse(
      substr(pop_label, 1, 3) %in% c("SMR", "SLR"),
      substr(pop_label, 1, 3),
      substr(pop_label, 1, 2)
    ),
    population2 = str_extract(pop_label, ".*(?=L)")
  )


# Réordonner les populations par emmean
pheno_var$pop_label <- fct_reorder(pheno_var$pop_label, pheno_var$emmean)

# Plot
ggplot(pheno_var, aes(x = pop_label, y = emmean, color = population)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Population",
    y = "Estimated var",
    title = "Trait var with 95% confidence intervals by population"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 2)  # taille modifiée pour lisibilité
  )

ggplot(pheno_var, aes(x = log(emmean))) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "Estimated var", y = "Density",
       title = "Distribution of emmean with density curve")

ggplot(pheno_var, aes(x = emmean)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.000002, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "Estimated var", y = "Density",
       title = "Distribution of emmean with density curve")

qqnorm(pheno_var$emmean)
qqline(pheno_var$emmean, col = "red", lwd = 2)

shapiro.test(pheno_var$emmean)



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


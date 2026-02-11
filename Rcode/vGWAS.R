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

# points extremes 
#phen_data_highvar <- phen_data_step1 %>%
#  filter(area.F.var > 3e-4)

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

phen_data_step1$area.F.mean_corr <- predict(moy_model)

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

#3 coefficient of variation (CV) NON car 2fois mean

phen_data_step1$area.F.var_corr <- predict(var_model_corr)

pred_var <- predict(var_model_corr)
phen_data_step1$CV_areaF <- (sqrt(pred_var) / phen_data_step1$area.F.mean_corr)*100

ggplot(phen_data_step1, aes(x = area.F, y = CV_areaF)) +
  geom_point(alpha = 0.3) +
  theme_classic() 

## ===============================
## 1.4 Extraction du phénotype final
## ===============================

## Résidus = variance intrinsèque
phen_data_step1$areaF_var_resid <- resid(var_model_corr)

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

pheno_line_var <- emmeans(var_model_corr, "pop_label")
pheno_var <- as.data.frame(pheno_line_var)



## 2.1 Lecture du VCF
vcf <- read.vcfR("/Users/edouardcourret/Documents/GitHub/2026_courret/data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz")

str(vcf@gt)
gt_df <- as.data.frame(vcf@gt)

#GEMMA
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
library(emmeans)
library(stringr)

## ===============================
## 1. Définition du phénotype
## ===============================

## 1.1 Filtrage des données
phen_data_step1 <- phen_data %>%
  filter(
    env_label == "NGM",
  ) %>%
  mutate(
    session_id = factor(session_id),
    block_id   = factor(block_id),
    pop_label  = factor(pop_label)
  )

#standard deviation method
mean_pheno.varA <- mean(phen_data_step1$area.F.var, na.rm = TRUE)
sd_pheno.varA <- sd(phen_data_step1$area.F.var, na.rm = TRUE)
cutoff <- 3
z_scores.varA <- (phen_data_step1$area.F.var - mean_pheno.varA) / sd_pheno.varA
table(abs(z_scores.varA) <= cutoff) #6 avec 3x SD
table(abs(z_scores.varA) >= -cutoff) #0 avec -3x SD (par symetrie)
phen_data_step1$exper_data_id[abs(z_scores.varA) > cutoff] # retire 6 outliers de var trop haut 
dim(phen_data_step1)
phen_data_step1 <- phen_data_step1[abs(z_scores.varA) <= cutoff, ]
dim(phen_data_step1)

## ===============================
## 1.2 Visualisation brute
## ===============================

#logD
ggplot(phen_data_step1, aes(x = logD, y = area.F)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(title = "area.F ~ logD",
       x = "logD",
       y = "area.F")

#Temperature 
ggplot(phen_data_step1, aes(x = temperature, y = area.F)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(title = "area.F ~ temperature",
       x = "Temperature",
       y = "area.F")

#Relative humidity
ggplot(phen_data_step1, aes(x = rel_humidity, y = area.F)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(title = "area.F ~ rel_humidity",
       x = "Relative humidity",
       y = "area.F")

#pop
ggplot(phen_data_step1, aes(x = pop_label, y = area.F)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "area.F par population",
       x = "Population",
       y = "area.F") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Block
ggplot(phen_data_step1, aes(x = block_id, y = area.F)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "area.F par block",
       x = "Block",
       y = "area.F") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Session
ggplot(phen_data_step1,
  aes(x = session_id, y = area.F, color = session_id)) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_classic() +
  labs(title = "area.F par session",
       x = "Session",
       y = "area.F",
       color = "Session") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

## ===============================
## 1.3 Model
## ===============================

mean_model <- lmer(
  area.F ~
    logD +
    temperature +
    rel_humidity +
    pop_label +
    (1 | session_id/block_id),
  data = phen_data_step1,
  REML = TRUE
)

summary(mean_model)
Anova(mean_model)

check_model(mean_model)

## ===============================
## 1.4 Extraction du phénotype final
## ===============================

## Résidus 
phen_data_step1$areaF_moy_resid <- resid(mean_model)

## Diagnostic post-correction
ggplot(phen_data_step1, aes(x = areaF_moy_resid)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "Residual moy distribution")

ggplot(phen_data_step1,
       aes(x = logD, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = temperature, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = rel_humidity, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = pop_label, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = block_id, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_step1,
       aes(x = session_id, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

## ===============================
## 2. EMMEANS
## ===============================

pheno_line_mean <- emmeans(mean_model, "pop_label")
pheno_mean <- as.data.frame(pheno_line_mean)

str(pheno_mean)
dim(pheno_mean)

str(phen_data_step1)

write.csv(pheno_mean,
          file = "pheno_mean.csv",
          row.names = FALSE)

write.table(pheno_mean,
            file = "pheno_mean.txt",
            row.names = FALSE,
            col.names = T,
            quote = FALSE)


pheno_mean <- read.table("data/pheno_mean.txt")

dim(pheno_mean)

pheno_mean <- pheno_mean %>%
  mutate(
    population = ifelse(
      substr(pop_label, 1, 3) %in% c("SMR", "SLR"),
      substr(pop_label, 1, 3),
      substr(pop_label, 1, 2)
    ),
    population2 = str_extract(pop_label, ".*(?=L)")
  )

#plot de la taille 
ggplot(pheno_mean, aes(x = pop_label, y = emmean, color = population)) +
  geom_point(size = 1, alpha = 0.8) +           
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL),    
                width = 0.1, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "lines",
    y = "Estimated mean",
    title = "Trait means with 95% confidence intervals by population"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size=2),
  )




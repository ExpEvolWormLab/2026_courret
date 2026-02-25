############################################################
## MEAN ESTIMATE – area.forward
## C. elegans
############################################################

## ===============================
## 0. Packages
## ===============================

phen_data <- fread("data/phen_data.txt", sep = "\t", header = TRUE)

library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)
library(broom)
library(qqman)
library(car)
library(performance)
library(emmeans)
library(stringr)

## ===============================
## 1. Effect of envrionmental covariate 
## ===============================

## 1.1 Data filtering
phen_data_NGM <- phen_data %>%
  filter(
    env_label == "NGM",
  ) %>%
  mutate(
    session_id = factor(session_id),
    block_id   = factor(block_id),
    pop_label  = factor(pop_label)
  )

## 1.2 Outliers 
mean_pheno.varA <- mean(phen_data_NGM$area.F.var, na.rm = TRUE)
sd_pheno.varA <- sd(phen_data_NGM$area.F.var, na.rm = TRUE)
cutoff <- 3
z_scores.varA <- (phen_data_NGM$area.F.var - mean_pheno.varA) / sd_pheno.varA
table(abs(z_scores.varA) <= cutoff) #6 avec 3x SD
table(abs(z_scores.varA) >= -cutoff) #0 avec -3x SD (par symetrie)
phen_data_NGM$exper_data_id[abs(z_scores.varA) > cutoff] # retire 6 outliers de var trop haut 
dim(phen_data_NGM)
phen_data_NGM <- phen_data_NGM[abs(z_scores.varA) <= cutoff, ]
dim(phen_data_NGM)

## ===============================
## 1.2 RAW data
## ===============================

#logD
ggplot(phen_data_NGM, aes(x = logD, y = area.F)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(title = "area.F ~ logD",
       x = "logD",
       y = "area.F")

#Temperature 
ggplot(phen_data_NGM, aes(x = temperature, y = area.F)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(title = "area.F ~ temperature",
       x = "Temperature",
       y = "area.F")

#Relative humidity
ggplot(phen_data_NGM, aes(x = rel_humidity, y = area.F)) +
  geom_point(alpha = 0.6) +
  theme_classic() +
  labs(title = "area.F ~ rel_humidity",
       x = "Relative humidity",
       y = "area.F")

#pop
ggplot(phen_data_NGM, aes(x = pop_label, y = area.F)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "area.F par population",
       x = "Population",
       y = "area.F") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Block
ggplot(phen_data_NGM, aes(x = block_id, y = area.F)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "area.F par block",
       x = "Block",
       y = "area.F")

#Session
ggplot(phen_data_NGM,
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
## 1.3 Mean model
## ===============================

mean_model <- lmer(
  area.F ~
    logD +
    temperature +
    rel_humidity +
    pop_label +
    (1 | session_id/block_id),
  data = phen_data_NGM,
  REML = TRUE
)

summary(mean_model)
Anova(mean_model)

check_model(mean_model)

## ===============================
## 1.4  Residuals
## ===============================

## Résidus 
phen_data_NGM$areaF_moy_resid <- resid(mean_model)

## Diagnostic post-correction
ggplot(phen_data_NGM, aes(x = areaF_moy_resid)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "Residual moy distribution")

ggplot(phen_data_NGM,
       aes(x = logD, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = temperature, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = rel_humidity, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = pop_label, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = block_id, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = session_id, y = areaF_moy_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

## ===============================
## 2. EMMEANS
## ===============================

pheno_line_mean <- emmeans(mean_model, "pop_label")
pheno_mean <- as.data.frame(pheno_line_mean)

pheno_mean <- pheno_mean %>%
  mutate(
    population = ifelse(
      substr(pop_label, 1, 3) %in% c("SMR", "SLR"),
      substr(pop_label, 1, 3),
      substr(pop_label, 1, 2)
    ),
    population2 = str_extract(pop_label, ".*(?=L)")
  )

write.table(pheno_mean, file = "pheno_mean.txt", row.names = FALSE, col.names = T, quote = FALSE)

## ===============================
## 2.1  Distribution of mean
## ===============================

pheno_mean <- fread("data/pheno_mean.txt", header = TRUE)

# plot of line/ populations
pheno_mean$pop_label <- fct_reorder(pheno_mean$pop_label, pheno_mean$emmean) # Réordonner les populations par emmean

ggplot(pheno_mean, aes(x = pop_label, y = emmean, color = population)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Population",
    y = "Estimated mean",
    title = "Trait mean with 95% confidence intervals by population"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 2)  
  )

# Histogram of mean frequency
ggplot(pheno_mean, aes(x = emmean)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.002, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "Estimated mean", y = "Density",
       title = "Distribution of emmean with density curve")

# Normality of mean
qqnorm(pheno_mean$emmean)
qqline(pheno_mean$emmean, col = "red", lwd = 2)

shapiro.test(pheno_mean$emmean)

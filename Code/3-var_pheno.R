############################################################
## VARIANCE ESTIMATE – area.forward
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
library(vcfR)
library(broom)
library(qqman)
library(car)
library(performance)
library(forcats)

## ===============================
## 1. Effect of environmental covariate
## ===============================

#Data filtering
phen_data_NGM <- phen_data %>%
  filter(
    env_label == "NaCl",
  ) %>%
  mutate(
    session_id = factor(session_id),
    block_id   = factor(block_id),
    pop_label  = factor(pop_label)
  )

## ===============================
## 1.2 RAW data
## ===============================

#Variance 
ggplot(phen_data_NGM, aes(x = area.F, y = log(area.F.var))) +
  geom_point(alpha = 0.3) +
  theme_classic() 
ggplot(phen_data_NGM, aes(x = area.F, y = area.F.var)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM, aes(x = area.F.var)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.000002, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "Estimated var", y = "Density",
       title = "Distribution of emmean with density curve")



highlight_pops <- c("LR3L107", "LR3L84", 
                    "LR3L105", "GA150L22", 
                    "CA150L23", "CA1100L77")
ggplot() +
  # Tous les points en gris
  geom_point(data = phen_data_NGM,
             aes(x = area.F, y = area.F.var),
             color = "grey80",
             alpha = 0.4) +
  
  # Les 9 points par-dessus
  geom_point(data = subset(phen_data_NGM, 
                           pop_label %in% highlight_pops),
             aes(x = area.F, y = area.F.var, color = pop_label),
             size = 3) +
  
  theme_classic()



#logD
ggplot(phen_data_NGM, aes(x = logD, y = area.F.var)) +
  geom_point(alpha = 0.6) +
  theme_classic() 

#Temperature 
ggplot(phen_data_NGM, aes(x = temperature, y = area.F.var)) +
  geom_point(alpha = 0.6) +
  theme_classic() 

#Relative humidity
ggplot(phen_data_NGM, aes(x = rel_humidity, y = area.F.var)) +
  geom_point(alpha = 0.6) +
  theme_classic() 

#pop
ggplot(phen_data_NGM, aes(x = pop_label, y = area.F.var)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "area.F par population",
       x = "Population",
       y = "area.F.var") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Block
ggplot(phen_data_NGM, aes(x = block_id, y = area.F.var)) +
  geom_jitter(width = 0.2, alpha = 0.5) +
  theme_classic() +
  labs(title = "area.F par block",
       x = "Block",
       y = "area.F.var") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

#Session
ggplot(phen_data_NGM,
       aes(x = session_id, y = area.F.var, color = session_id)) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_classic() 


ggplot() +
  # Tous les points en gris
  geom_jitter(data = phen_data_NGM,
              aes(x = session_id, y = area.F.var),
              width = 0.2,
              color = "grey80",
              alpha = 0.5) +
  
  # Les 9 points au-dessus, colorés par session
  geom_jitter(data = subset(phen_data_NGM, 
                            pop_label %in% highlight_pops),
              aes(x = session_id, 
                  y = area.F.var,
                  color = pop_label),
              width = 0.2,
              size = 3) +
  
  theme_classic()

head(
  phen_data_NGM[order(-phen_data_NGM$area.F.var), 
                c("pop_label", "area.F.var")],
  15
)

head(
  phen_data_NGM[order(-phen_data_NGM$area.F.var), 
                c("pop_label", "area.F.var", "exper_data_id", "block_id")],
  15
)


head(
  phen_data_NGM[
    phen_data_NGM$pop_label %in% highlight_pops,
    c("pop_label", "exper_data_id")
  ]
)

head(
  phen_data_NGM[
    phen_data_NGM$pop_label %in% highlight_pops,
    c("pop_label", "exper_data_id")
  ],
  9
)
## ===============================
## 1.3 Variance models
## ===============================

#1  mean model

var_model <- lmer(
  area.F.var ~
    area.F +
    logD +
    temperature +
    rel_humidity +
    pop_label +
    (1 | session_id/block_id),
  data = phen_data_NGM,
  REML = TRUE
)

summary(var_model)
Anova(var_model)

check_model(var_model)

#2  Corrected mean model

phen_data_NGM$area.F.mean_corr <- predict(mean_model)

var_model_corr <- lmer(
  area.F.var ~
    area.F.mean_corr +
    logD +
    temperature +
    rel_humidity +
    pop_label +
    (1 | session_id/block_id),
  data = phen_data_NGM,
  REML = TRUE
)

summary(var_model_corr)
Anova(var_model_corr)

check_model(var_model_corr)

# graph distribution

phen_data_NGM$area.F.var_corr = predict(var_model) #mean non corr
#phen_data_NGM$area.F.var_corr2 <- predict(var_model_corr) # mean corr

plot(phen_data_NGM$area.F.var_corr, phen_data_NGM$area.F.var_corr2) #pareil donc je garde avec mean non corrigé

plot(phen_data_NGM$area.F.var, phen_data_NGM$area.F.var_corr)

#3  Coefficient of variation (CV)

phen_data_NGM$CV_areaF <- (sqrt(phen_data_NGM$area.F.var) / phen_data_NGM$area.F)*100

CV_model <- lmer(
  CV_areaF ~
    logD + #without mean
    temperature +
    rel_humidity +
    pop_label +
    (1 | session_id/block_id),
  data = phen_data_NGM,
  REML = TRUE
)

summary(CV_model)
Anova(CV_model)

check_model(CV_model)

phen_data_NGM$CV=predict(CV_model)

plot(phen_data_NGM$area.F, phen_data_NGM$area.F.mean_corr)

plot(phen_data_NGM$CV_areaF, phen_data_NGM$area.F.var_corr)
plot(phen_data_NGM$CV_areaF, phen_data_NGM$area.F.mean_corr)

## ===============================
## 1.4 Residuals
## ===============================

## Résidus = variance intrinsèque
phen_data_NGM$areaF_var_resid <- resid(var_model)

## Diagnostic post-correction
ggplot(phen_data_NGM, aes(x = areaF_var_resid)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "Residual moy distribution")

ggplot(phen_data_NGM,
       aes(x = area.F, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = logD, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = temperature, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = rel_humidity, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = pop_label, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = block_id, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

ggplot(phen_data_NGM,
       aes(x = session_id, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic()

## ===============================
## 2. EMMEANS - variance
## ===============================

pheno_line_var <- emmeans(var_model, "pop_label")
pheno_var <- as.data.frame(pheno_line_var)

pheno_var <- pheno_var %>%
  mutate(
    population = ifelse(
      substr(pop_label, 1, 3) %in% c("SMR", "SLR"),
      substr(pop_label, 1, 3),
      substr(pop_label, 1, 2)
    ),
    population2 = str_extract(pop_label, ".*(?=L)")
  )

write.table(pheno_var, file = "pheno_var_NaCl.txt", row.names = FALSE, col.names = T, quote = FALSE)

## ===============================
## 2.1  Distribution of variance
## ===============================

pheno_var <- fread("data/pheno_var.txt", header = TRUE)

# plot of line/ populations
pheno_var$pop_label <- fct_reorder(pheno_var$pop_label, pheno_var$emmean) # Réordonner les populations par emmean
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
    axis.text.x = element_text(angle = 45, hjust = 1, size = 2)  
  )

# Histogram of variance frequency
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

# Normality of variance
qqnorm(pheno_var$emmean)
qqline(pheno_var$emmean, col = "red", lwd = 2)

shapiro.test(pheno_var$emmean)

## ===============================
## 2.2  Correlation between mean-variance
## ===============================

log_emmean_var=log(pheno_var$emmean)
plot(pheno_mean$emmean, log_emmean_var) #log de var

plot(pheno_mean$emmean, pheno_var$emmean) #var

cor(pheno_mean$emmean, pheno_var$emmean) # corr de 0.15


## ===============================
## 3. EMMEAN - CV
## ===============================

pheno_line_CV <- emmeans(CV_model, "pop_label")
pheno_CV <- as.data.frame(pheno_line_CV)

pheno_CV <- pheno_CV %>%
  mutate(
    population = ifelse(
      substr(pop_label, 1, 3) %in% c("SMR", "SLR"),
      substr(pop_label, 1, 3),
      substr(pop_label, 1, 2)
    ),
    population2 = str_extract(pop_label, ".*(?=L)")
  )

write.table(pheno_CV, file = "data/pheno_CV.txt", row.names = FALSE, col.names = T, quote = FALSE)

## ===============================
## 3.1  Distribution of CV
## ===============================

#histogram CV
ggplot(pheno_CV, aes(x = emmean)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.2, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "CV", y = "Density",
       title = "Distribution of CV with density curve")

# histogram log(CV)
ggplot(pheno_CV, aes(x = log(emmean))) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.02, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "log(CV)", y = "Density",
       title = "Distribution of CV with density curve")


# plot of line/ populations
pheno_CV$pop_label <- fct_reorder(pheno_CV$pop_label, pheno_CV$emmean) # Réordonner les populations par emmean
ggplot(pheno_CV, aes(x = pop_label, y = emmean, color=population)) +
  geom_point(size = 1, alpha = 0.8) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, alpha = 0.5) +
  theme_classic() +
  labs(
    x = "Population",
    y = "CV",
    title = "Trait var with 95% confidence intervals by population"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 2)  
  )

## ===============================
## 3.2  Correlation with CV
## ===============================

plot(pheno_mean$emmean, pheno_CV$emmean)
cor(pheno_mean$emmean, pheno_CV$emmean) # corr de -0.46

plot(pheno_var$emmean, pheno_CV$emmean)
cor(pheno_var$emmean, pheno_CV$emmean) # corr de 0.74

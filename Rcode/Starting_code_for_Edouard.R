rm(list=ls());gc()
library(data.table)
library(matrixStats)
##### Load the phenotypic data

phen_data <- read.table("/Users/edouardcourret/Documents/GitHub/2026_courret/data/all_phenotypes_from_movies.txt",sep='\t',h=TRUE)

# ## Create the session ID 
phen_data$session_id <- substring(phen_data$date_str,1,4)
phen_data$session_id[phen_data$session_id=="2014"] <- "2013"
phen_data$session_id[substring(phen_data$block_id,1,1)=="F"] <- "F"
phen_data$session_id[substring(phen_data$block_id,1,2)=="SB"] <- "SB"

## Set the env. covariates as numeric
for(i in c("temperature","rel_humidity","logD")){
  phen_data[,i]=as.numeric(phen_data[,i])
}

phen_data=subset(phen_data,!is.na(phen_data$pop_label))  # Mostly MA lines and populations
# Remove additonal populations
phen_data=subset(phen_data,!phen_data$block_id%in%paste0("B",300:350))
phen_data=subset(phen_data,!phen_data$pop_label=="A00")
phen_data=subset(phen_data,!phen_data$pop_label=="CX12311L0")
phen_data=subset(phen_data,!phen_data$pop_label=="LSJ2L0")

## Create the population ID
phen_data$population = tstrsplit(phen_data$pop_label,"L")[[1]]
phen_data$population[phen_data$population==""]=paste0("L",tstrsplit(phen_data$pop_label,"L")[[2]])[phen_data$population==""]
phen_data$population[phen_data$population=="S"]="SL3"
#phen_data$pop_label[phen_data$population=="LSJ2"]
phen_data$population2=phen_data$population
phen_data$population2=substring(phen_data$population,1,2)
phen_data$population2[phen_data$population=="EEV1401"]="A6"

phen_data$date_str = as.factor(as.character(phen_data$date_str))

### Remove NA values
phen_data=subset(phen_data,!is.na(phen_data$temperature))
phen_data=subset(phen_data,!is.na(phen_data$rel_humidity))
phen_data=subset(phen_data,!is.na(phen_data$logD))

for(i in c("temperature","rel_humidity","logD")){
  phen_data[,i]= (phen_data[,i]-mean(phen_data[,i],na.rm=TRUE))/sd(phen_data[,i])
}

vect_P_traits = c("area.F","area.F.var")

# Remove NA values for phenotypes
phen_data=subset(phen_data,!is.na(area.F) & !is.na(area.F.var))
dim(phen_data) # 5281 movies

## Investigate the data
plot(phen_data$area.F,phen_data$area.F.var,log="y")
## Some extreme variance, due to density ?
plot(phen_data$logD,phen_data$area.F.var,log="y") # Not only

## Remove outliers of both density and variance
# Standard deviation method
mean_pheno.varA <- mean(phen_data$area.F.var, na.rm = TRUE)
mean_logD <- mean(phen_data$logD, na.rm = TRUE)
sd_pheno.varA <- sd(phen_data$area.F.var, na.rm = TRUE)
sd_logD <- sd(phen_data$logD, na.rm = TRUE)
# Chose a cutoff
cutoff <- 3
z_scores.varA <- (phen_data$area.F.var - mean_pheno.varA) / sd_pheno.varA
z_scores.logD <- (phen_data$logD - mean_logD) / sd_logD

table(abs(z_scores.varA) <= cutoff) # 15 values
table(abs(z_scores.logD) <= cutoff) # 111 values

plot(phen_data$logD,phen_data$area.F.var,log="y",pch=16,col=c("black","red")[as.numeric(abs(z_scores.varA) <= cutoff)+1]) # Not only
plot(phen_data$logD,phen_data$area.F.var,log="y",pch=16,col=c("black","red")[as.numeric(abs(z_scores.logD) <= cutoff & abs(z_scores.varA) <= cutoff)+1]) # Not only

#### Ok great
phen_data_filt <- subset(phen_data,abs(z_scores.logD) <= cutoff & abs(z_scores.varA) <= cutoff)
dim(phen_data_filt)

phen_data_NGM = subset(phen_data_filt,env_label=='NGM')
phen_data_NaCl = subset(phen_data_filt,env_label=='NaCl')

plot(phen_data_filt$area.F,phen_data_filt$area.F.var,log="y",pch=16)
plot(phen_data_filt$logD,phen_data_filt$area.F,log="y",pch=16)

### Finally we scale the phenotypic data

phen_data_filt$area.F=scale(phen_data_filt$area.F)
phen_data_filt$area.F.var=scale(phen_data_filt$area.F.var)



######## Analyse Edouard ###########

str(phen_data_filt)

############################################
## STEP 1 — DEFINE PHENOTYPE (variance trait)
############################################

library(lme4)
library(dplyr)
library(ggplot2)
library(car)

## -------------------------------
## 1. Filtrage des données
## -------------------------------

phen_data_step1 <- phen_data_filt %>%
  filter(
    env_label == "NGM",        # environnement unique
    QC_label == "ok"           # qualité OK
  ) %>%
  mutate(
    session_id = factor(session_id),
    block_id   = factor(block_id),
    pop_label  = factor(pop_label)
  )

## -------------------------------
## 2. Vérification brute 
## -------------------------------

ggplot(phen_data_step1, aes(x = area.F, y = area.F.var)) +
  geom_point(alpha = 0.3) +
  theme_classic() +
  labs(
    x = "Mean area.forward (raw)",
    y = "Variance area.forward (raw)",
    title = "Raw mean–variance relationship"
  )

## -------------------------------
## 3. Modèle environnemental sur la VARIANCE
## -------------------------------

## Y = variance du trait
## On corrige :
## - densité
## - température
## - humidité
## - session (caméra / setup)
## - block (date)

var_model <- lmer(
  area.F.var ~ 
    area.F +
    logD +
    temperature +
    rel_humidity +
    (1 | block_id / session_id),
  data = phen_data_step1,
  REML = TRUE
)

summary(var_model)

Anova(var_model)

## -------------------------------
## 4. Extraction du phénotype corrigé
## -------------------------------

## Résidus conditionnels = variance intrinsèque
phen_data_step1$areaF_var_resid <- resid(var_model)

## Centrage / standardisation 
phen_data_step1$areaF_var_resid_scaled <- scale(
  phen_data_step1$areaF_var_resid
)

## -------------------------------
## 5. Diagnostic post-correction
## -------------------------------

## a) Résidus ~ 0 ?
ggplot(phen_data_step1, aes(x = areaF_var_resid)) +
  geom_histogram(bins = 50) +
  theme_classic() +
  labs(title = "Distribution of variance residuals")

## b) Plus de structure session ?
ggplot(phen_data_step1,
       aes(x = session_id, y = areaF_var_resid)) +
  geom_boxplot(outlier.alpha = 0.2) +
  theme_classic() +
  theme(axis.text.x = element_blank()) +
  labs(title = "Residual variance by session")

## c) Mean–variance APRÈS correction
ggplot(phen_data_step1,
       aes(x = area.F, y = areaF_var_resid)) +
  geom_point(alpha = 0.3) +
  theme_classic() +
  labs(
    x = "Mean area.forward",
    y = "Corrected variance (residuals)",
    title = "Mean–variance relationship after correction"
  )

## -------------------------------
## 6. Objet final pour la suite
## -------------------------------

phenotype_final <- phen_data_step1 %>%
  select(
    pop_label,
    block_id,
    session_id,
    areaF_var_resid_scaled
  )





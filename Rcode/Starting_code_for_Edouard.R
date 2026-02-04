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















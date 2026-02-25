############################################################
## PHENOTYPE DATA PREPARATION – area.F / area.F.var
## C. elegans – CeMEE panel
############################################################

rm(list = ls())

## ===============================
## 0. Packages
## ===============================

library(data.table)
library(matrixStats)

## ===============================
## 1. Load phenotypic data
## ===============================

phen_data <- fread("data/all_phenotypes_from_movies_fixed.txt", sep = "\t", header = TRUE)

## ===============================
## 2. Create experimental variables
## ===============================

### 2.1 Session ID
phen_data[, session_id := substr(date_str, 1, 4)]
phen_data[session_id == "2014", session_id := "2013"]
phen_data[substr(block_id, 1, 1) == "F", session_id := "F"]
phen_data[substr(block_id, 1, 2) == "SB", session_id := "SB"]

### 2.2 Convert environmental covariates to numeric
for(i in c("temperature","rel_humidity","logD")){
  phen_data[,i]=as.numeric(phen_data[,i])
}

## ===============================
## 3. Data filtering
## ===============================

### 3.1 Remove unwanted populations
phen_data <- phen_data[!is.na(pop_label)]
phen_data <- phen_data[!block_id %in% paste0("B", 200:400)]
phen_data <- phen_data[!pop_label %in% c("A00", "CX12311L0", "LSJ2L0")]

### 3.2 Remove missing environmental data
phen_data <- phen_data[
  !is.na(temperature) &
    !is.na(rel_humidity) &
    !is.na(logD)
]

## ===============================
## 4. Define population identifiers
## ===============================

### 4.1 Main population ID
phen_data[, population := tstrsplit(pop_label, "L")[[1]]]

empty_idx <- phen_data$population == ""
phen_data[empty_idx,
          population := paste0("L", tstrsplit(pop_label, "L")[[2]][empty_idx])]

phen_data[population == "S", population := "SL3"]

### 4.2 Short population grouping (first two letters)
phen_data[, population2 := substr(population, 1, 2)]
phen_data[population == "EEV1401", population2 := "A6"]

## ===============================
## 5. Variable formatting
## ===============================

phen_data[, date_str := as.factor(as.character(date_str))]

## ===============================
## 6. Define phenotypes of interest
## ===============================

vect_P_traits <- c("area.F", "area.F.var")

phen_data <- phen_data[
  !is.na(area.F) &
    !is.na(area.F.var)
]

dim(phen_data) #4999 movies 

## ===============================
## 7. Final table
## ===============================

write.table(phen_data, file = "data/phen_data.txt", sep = "\t", quote = FALSE, row.names = FALSE)

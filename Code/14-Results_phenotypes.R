############################################################
## PHENOTYPE DATA PREPARATION & ANALYSES
## C. elegans – CeMEE panel
############################################################


## ===============================
## I. PREPARATION & MODELS
## ===============================

## ===============================
## 1. LOAD DATA
## ===============================

library(data.table)
library(lme4)
library(car)
library(dplyr)
library(ggplot2)
library(patchwork)
library(scales)


rm(list = ls())
phen_data <- fread("data/phen_data.txt")

cutoff <- 3

# =========================================================
# NGM
# =========================================================

mean_pheno.varA <- mean(phen_data$area.F.var[phen_data$env_label == "NGM"], na.rm = TRUE)
sd_pheno.varA   <- sd(phen_data$area.F.var[phen_data$env_label == "NGM"], na.rm = TRUE)

z_scores.varA <- (phen_data$area.F.var - mean_pheno.varA) / sd_pheno.varA

table(abs(z_scores.varA[phen_data$env_label == "NGM"]) <= cutoff)

phen_data$exper_data_id[phen_data$env_label == "NGM"][abs(z_scores.varA[phen_data$env_label == "NGM"]) > cutoff]

dim(phen_data)
phen_data <- phen_data[!(phen_data$env_label == "NGM" & abs(z_scores.varA) > cutoff), ]
dim(phen_data)

# =========================================================
# NaCl
# =========================================================

mean_pheno.varA <- mean(phen_data$area.F.var[phen_data$env_label == "NaCl"], na.rm = TRUE)
sd_pheno.varA   <- sd(phen_data$area.F.var[phen_data$env_label == "NaCl"], na.rm = TRUE)

z_scores.varA <- (phen_data$area.F.var - mean_pheno.varA) / sd_pheno.varA

table(abs(z_scores.varA[phen_data$env_label == "NaCl"]) <= cutoff)

phen_data$exper_data_id[phen_data$env_label == "NaCl"][abs(z_scores.varA[phen_data$env_label == "NaCl"]) > cutoff]

dim(phen_data)
phen_data <- phen_data[!(phen_data$env_label == "NaCl" & abs(z_scores.varA) > cutoff), ]
dim(phen_data)

str(phen_data)
corr_isotypes <- fread("data/correspondance_isotypes.csv")
removed_isotypes <- fread("data/removed_isotypes_Cemee_2025.txt", header=F)

## ===============================
## 5. CORRECTION ISOTYPES
## ===============================

# clean noms
phen_data[, pop_label := gsub("CeMee", "", pop_label)]
corr_isotypes[, Line := gsub("CeMee", "", Line)]
corr_isotypes[, Representant_line := gsub("CeMee", "", Representant_line)]
removed_isotypes[, V1 := gsub("CeMee", "", V1)]

# merge correspondance
phen_data <- merge(
  phen_data,
  corr_isotypes[, .(Line, Representant_line)],
  by.x = "pop_label",
  by.y = "Line",
  all.x = TRUE
)

# remplacer par isotype représentatif
phen_data[, pop_label := ifelse(
  is.na(Representant_line),
  pop_label,
  Representant_line
)]

phen_data[, Representant_line := NULL]

# enlever lignes à exclure
phen_data <- phen_data[!pop_label %in% removed_isotypes$V1]

## ===============================
## 6. FORMAT FINAL
## ===============================

phen_clean <- phen_data[, .(
  pop_label,
  area.F,
  area.F.var,
  logD,
  rel_humidity,
  temperature,
  date_str,
  session_id, 
  env_label, 
  block_id
)]

str(phen_clean)

mean_model <- lmer(
  area.F ~ logD + temperature + rel_humidity + pop_label + env_label +
    (1 | session_id/block_id),
  data = phen_clean,
  REML = TRUE
)

Anova(mean_model)

var_model <- lmer(
  log(area.F.var) ~ area.F + logD + temperature + rel_humidity + pop_label + env_label +
    (1 | session_id/block_id),
  data = phen_clean,
  REML = TRUE
)

Anova(var_model)

CV_model <- lmer(
  log((sqrt(area.F.var) / area.F) * 100) ~
    logD + temperature + rel_humidity + pop_label + env_label +
    (1 | session_id/block_id),
  data = phen_clean,
  REML = TRUE
)

Anova(CV_model)


library(lme4)
library(dplyr)

phen_clean <- as.data.frame(phen_clean)

# =========================================================
# SPLIT
# =========================================================

phen_NGM  <- phen_clean %>% filter(env_label == "NGM")
phen_NaCl <- phen_clean %>% filter(env_label == "NaCl")

# Correlation area.F ~ area.F.var
cor_NGM <- cor.test(
  phen_NGM$area.F,
  phen_NGM$area.F.var,
  method = "pearson"
)

cor_NGM

# Correlation area.F ~ area.F.var
cor_NaCl <- cor.test(
  phen_NaCl$area.F,
  phen_NaCl$area.F.var,
  method = "pearson"
)

cor_NaCl

# =========================================================
# MODELS - NGM
# =========================================================

mean_model_NGM <- lmer(
  area.F ~ logD + temperature + rel_humidity + pop_label +
    (1 | session_id/block_id),
  data = phen_NGM,
  REML = TRUE
)

Anova(mean_model_NGM)

var_model_NGM <- lmer(
  log(area.F.var) ~ area.F + logD + temperature + rel_humidity + pop_label +
    (1 | session_id/block_id),
  data = phen_NGM,
  REML = TRUE
)

Anova(var_model_NGM)

CV_model_NGM <- lmer(
  log((sqrt(area.F.var) / area.F) * 100) ~
    logD + temperature + rel_humidity + pop_label +
    (1 | session_id/block_id),
  data = phen_NGM,
  REML = TRUE
)

Anova(CV_model_NGM)

# =========================================================
# MODELS - NaCl
# =========================================================

mean_model_NaCl <- lmer(
  area.F ~ logD + temperature + rel_humidity + pop_label +
    (1 | session_id/block_id),
  data = phen_NaCl,
  REML = TRUE
)

Anova(mean_model_NaCl)

var_model_NaCl <- lmer(
  log(area.F.var) ~ area.F + logD + temperature + rel_humidity + pop_label +
    (1 | session_id/block_id),
  data = phen_NaCl,
  REML = TRUE
)

Anova(var_model_NaCl)

CV_model_NaCl <- lmer(
  log((sqrt(area.F.var) / area.F) * 100) ~
    logD + temperature + rel_humidity + pop_label +
    (1 | session_id/block_id),
  data = phen_NaCl,
  REML = TRUE
)

Anova(CV_model_NaCl)


library(emmeans)

# =========================================================
# NGM
# =========================================================

pheno_line_mean_NGM <- emmeans(mean_model_NGM, "pop_label")
pheno_mean_NGM <- as.data.frame(pheno_line_mean_NGM)

pheno_line_var_NGM <- emmeans(var_model_NGM, "pop_label")
pheno_var_NGM <- as.data.frame(pheno_line_var_NGM)

pheno_line_CV_NGM <- emmeans(CV_model_NGM, "pop_label")
pheno_CV_NGM <- as.data.frame(pheno_line_CV_NGM)

# =========================================================
# NaCl
# =========================================================

pheno_line_mean_NaCl <- emmeans(mean_model_NaCl, "pop_label")
pheno_mean_NaCl <- as.data.frame(pheno_line_mean_NaCl)

pheno_line_var_NaCl <- emmeans(var_model_NaCl, "pop_label")
pheno_var_NaCl <- as.data.frame(pheno_line_var_NaCl)

pheno_line_CV_NaCl <- emmeans(CV_model_NaCl, "pop_label")
pheno_CV_NaCl <- as.data.frame(pheno_line_CV_NaCl)

# =========================================================
# SAVE TABLES - NGM
# =========================================================

fwrite(
  pheno_mean_NGM,
  file = "data/pheno_mean_NGM.txt",
  sep = "\t"
)

fwrite(
  pheno_var_NGM,
  file = "data/pheno_var_NGM.txt",
  sep = "\t"
)

fwrite(
  pheno_CV_NGM,
  file = "data/pheno_CV_NGM.txt",
  sep = "\t"
)

# =========================================================
# SAVE TABLES - NaCl
# =========================================================

fwrite(
  pheno_mean_NaCl,
  file = "data/pheno_mean_NaCl.txt",
  sep = "\t"
)

fwrite(
  pheno_var_NaCl,
  file = "data/pheno_var_NaCl.txt",
  sep = "\t"
)

fwrite(
  pheno_CV_NaCl,
  file = "data/pheno_CV_NaCl.txt",
  sep = "\t"
)

str(pheno_CV_NaCl)

table_NGM <- pheno_mean_NGM %>%
  select(pop_label, mean = emmean) %>%
  left_join(pheno_var_NGM %>% select(pop_label, var = emmean), by = "pop_label") %>%
  left_join(pheno_CV_NGM %>% select(pop_label, CV = emmean), by = "pop_label")


table_NaCl <- pheno_mean_NaCl %>%
  select(pop_label, mean = emmean) %>%
  left_join(pheno_var_NaCl %>% select(pop_label, var = emmean), by = "pop_label") %>%
  left_join(pheno_CV_NaCl %>% select(pop_label, CV = emmean), by = "pop_label")

fwrite(table_NaCl, "data/table_NaCl.txt")
fwrite(table_NGM, "data/table_NGM.txt")


## ===============================
## II. ANALYSES
## ===============================


BLUE_NGM <- fread("data/table_NGM.txt")
BLUE_NaCl <- fread("data/table_NaCl.txt")

str(BLUE_NaCl)
str(BLUE_NGM)

# =========================================================
# DATA
# =========================================================

BLUE_NGM  <- as.data.frame(BLUE_NGM)
BLUE_NaCl <- as.data.frame(BLUE_NaCl)

# =========================================================
# COLORS
# =========================================================

col_NGM  <- "#4C72B0"
col_NaCl <- "#C44E52"

# =========================================================
# THEME
# =========================================================

theme_science <- theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 18, face = "bold"),
    axis.text  = element_text(size = 14, colour = "black"),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none",
    plot.margin = margin(12,12,12,12)
  )

# =========================================================
# CORRELATIONS (var / CV + mean / var + mean / CV)
# =========================================================

# -------------------------
# NGM
# -------------------------

cor_var_CV_NGM   <- cor.test(BLUE_NGM$var,  BLUE_NGM$CV)
cor_mean_var_NGM <- cor.test(BLUE_NGM$mean, BLUE_NGM$var)
cor_mean_CV_NGM  <- cor.test(BLUE_NGM$mean, BLUE_NGM$CV)

# -------------------------
# NaCl
# -------------------------

cor_var_CV_NaCl   <- cor.test(BLUE_NaCl$var,  BLUE_NaCl$CV)
cor_mean_var_NaCl <- cor.test(BLUE_NaCl$mean, BLUE_NaCl$var)
cor_mean_CV_NaCl  <- cor.test(BLUE_NaCl$mean, BLUE_NaCl$CV)

# =========================================================
# DISTRIBUTION OF VARIANCE
# =========================================================

BLUE_NGM_sorted <- BLUE_NGM %>%
  arrange(var) %>%
  mutate(rank = 1:n())

BLUE_NaCl_sorted <- BLUE_NaCl %>%
  arrange(var) %>%
  mutate(rank = 1:n())

# -------------------------
# NGM
# -------------------------

p_var_NGM <- ggplot(BLUE_NGM_sorted,
                    aes(x = rank,
                        y = var)) +
  geom_line(linewidth = 1.5,
            colour = col_NGM) +
  labs(title = "NGM",
       x = "Strains ranked by variance",
       y = "Predicted variance") +
  theme_science

# -------------------------
# NaCl
# -------------------------

p_var_NaCl <- ggplot(BLUE_NaCl_sorted,
                     aes(x = rank,
                         y = var)) +
  geom_line(linewidth = 1.5,
            colour = col_NaCl) +
  labs(title = "NaCl",
       x = "Strains ranked by variance",
       y = "Predicted variance") +
  theme_science

# =========================================================
# FUNCTION FOR CORRELATION PLOTS
# =========================================================

make_cor_plot <- function(data,
                          xvar,
                          yvar,
                          color,
                          title){
  
  cor_res <- cor.test(data[[xvar]],
                      data[[yvar]])
  
  r_val <- round(cor_res$estimate, 3)
  p_val <- signif(cor_res$p.value, 3)
  
  ggplot(data,
         aes_string(x = xvar,
                    y = yvar)) +
    
    geom_point(size = 2.3,
               alpha = 0.65,
               colour = color) +
    
    geom_smooth(method = "lm",
                se = TRUE,
                linewidth = 1.2,
                colour = "black") +
    
    annotate("text",
             x = Inf,
             y = Inf,
             label = paste0("r = ", r_val,
                            "\np = ", p_val),
             hjust = 1.1,
             vjust = 1.5,
             size = 5.2,
             fontface = "bold") +
    
    labs(title = title,
         x = xvar,
         y = yvar) +
    
    theme_science
}

# =========================================================
# NGM PLOTS
# =========================================================

p1 <- make_cor_plot(BLUE_NGM, "var",  "CV",   col_NGM, "NGM: Variance vs CV")
p2 <- make_cor_plot(BLUE_NGM, "mean", "var",  col_NGM, "NGM: Mean vs Variance")
p3 <- make_cor_plot(BLUE_NGM, "mean", "CV",   col_NGM, "NGM: Mean vs CV")

# =========================================================
# NaCl PLOTS
# =========================================================

p4 <- make_cor_plot(BLUE_NaCl, "var",  "CV",  col_NaCl, "NaCl: Variance vs CV")
p5 <- make_cor_plot(BLUE_NaCl, "mean", "var", col_NaCl, "NaCl: Mean vs Variance")
p6 <- make_cor_plot(BLUE_NaCl, "mean", "CV",  col_NaCl, "NaCl: Mean vs CV")

# =========================================================
# DISPLAY
# =========================================================

# variance ranking
p_var_NGM
p_var_NaCl

# correlations
(p1 | p2 | p3)
(p4 | p5 | p6)

# cross comparisons
(p1 | p4) # var vs CV
(p2 | p5) # mean vs var
(p3 | p6) # mean vs CV

# =========================================================
# PRINT CORRELATION RESULTS
# =========================================================

cor_var_CV_NGM
cor_mean_var_NGM
cor_mean_CV_NGM

cor_var_CV_NaCl
cor_mean_var_NaCl
cor_mean_CV_NaCl


# =========================================================
# MERGE DATA
# =========================================================

BLUE_NGM$env  <- "NGM"
BLUE_NaCl$env <- "NaCl"

str(BLUE_NaCl)
str(BLUE_NGM)

all_data <- bind_rows(BLUE_NGM, BLUE_NaCl)
str(all_data)

# =========================================================
# FUNCTION FOR COMBINED ENV PLOTS
# =========================================================

make_combined_plot <- function(data,
                               xvar,
                               yvar,
                               title){
  
  # correlations
  cor_NGM <- cor.test(
    data[[xvar]][data$env == "NGM"],
    data[[yvar]][data$env == "NGM"]
  )
  
  cor_NaCl <- cor.test(
    data[[xvar]][data$env == "NaCl"],
    data[[yvar]][data$env == "NaCl"]
  )
  
  # labels
  label_NGM <- paste0(
    "NGM: r = ",
    round(cor_NGM$estimate, 3),
    "\np = ",
    signif(cor_NGM$p.value, 3)
  )
  
  label_NaCl <- paste0(
    "NaCl: r = ",
    round(cor_NaCl$estimate, 3),
    "\np = ",
    signif(cor_NaCl$p.value, 3)
  )
  
  # plot
  ggplot(data,
         aes_string(x = xvar,
                    y = yvar,
                    colour = "env")) +
    
    # points
    geom_point(size = 0.5,
               alpha = 0.5) +
    
    # regression lines
    geom_smooth(method = "lm",
                se = TRUE,
                linewidth = 1.2) +
    
    # colors
    scale_colour_manual(values = c(
      "NGM"  = col_NGM,
      "NaCl" = col_NaCl
    )) +
    
    # annotations
    annotate("text",
             x = Inf,
             y = Inf,
             label = label_NGM,
             hjust = 1.05,
             vjust = 1.5,
             colour = col_NGM,
             size = 4.7,
             fontface = "bold") +
    
    annotate("text",
             x = Inf,
             y = Inf,
             label = label_NaCl,
             hjust = 1.05,
             vjust = 4.2,
             colour = col_NaCl,
             size = 4.7,
             fontface = "bold") +
    
    labs(title = title,
         x = xvar,
         y = yvar,
         colour = "Environment") +
    
    theme_science +
    
    theme(
      legend.position = "top"
    )
}

# =========================================================
# COMBINED PLOTS
# =========================================================

# Mean vs Variance
p_combined_mean_var <- make_combined_plot(
  all_data,
  "mean",
  "var",
  "Mean vs Variance"
)

# Variance vs CV
p_combined_var_CV <- make_combined_plot(
  all_data,
  "var",
  "CV",
  "Variance vs CV"
)

# Mean vs CV
p_combined_mean_CV <- make_combined_plot(
  all_data,
  "mean",
  "CV",
  "Mean vs CV"
)

# =========================================================
# DISPLAY
# =========================================================

(p_combined_mean_var |
   p_combined_var_CV |
   p_combined_mean_CV)






# =========================================================
# FUNCTION FOR COMBINED ENV PLOTS
# with panel labels a/b/c
# =========================================================

make_combined_plot <- function(data,
                               xvar,
                               yvar,
                               xlab,
                               ylab,
                               panel_label = NULL,
                               show_legend = FALSE){
  
  # panel label position
  x_min <- min(data[[xvar]], na.rm = TRUE)
  x_range <- diff(range(data[[xvar]], na.rm = TRUE))
  
  y_max <- max(data[[yvar]], na.rm = TRUE)
  y_range <- diff(range(data[[yvar]], na.rm = TRUE))
  
  ggplot(data,
         aes_string(x = xvar,
                    y = yvar,
                    colour = "env")) +
    
    geom_point(size = 0.75, alpha = 0.4) +
    
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 1.4) +
    
    # panel label (a, b, c)
    annotate("text",
             x = x_min - 0.02 * x_range,
             y = y_max + 0.02 * y_range,
             label = panel_label,
             fontface = "bold",
             size = 6,
             hjust = 0,
             vjust = 1) +
    
    scale_colour_manual(
      
      values = c(
        
        "NGM"  = col_NGM,
        
        "NaCl" = col_NaCl
        
      ),
      
      labels = c(
        
        "NGM"  = "Low-salt",
        
        "NaCl" = "High-salt"
        
      )
      
    ) +
    
    labs(
      x = xlab,
      y = ylab,
      colour = NULL
    ) +
    
    theme_classic(base_size = 16) +
    
    theme(
      legend.position = ifelse(show_legend, "top", "none"),
      legend.justification = "center",
      legend.direction = "horizontal",
      
      legend.text = element_text(face = "bold", size = 13),
      
      axis.title = element_text(face = "bold", size = 16),
      axis.text = element_text(size = 13),
      
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(linewidth = 0.7),
      
      plot.margin = margin(10, 10, 10, 10)
    )
}


p_combined_var_CV <- make_combined_plot(
  all_data,
  "var",
  "CV",
  xlab = expression(bold("ln(" * sigma^2 * ")")),
  ylab = expression(bold("ln(CV)")),
  panel_label = "a)"
)

p_combined_mean_var <- make_combined_plot(
  all_data,
  "mean",
  "var",
  xlab = expression(bold("Mean")),
  ylab = expression(bold("ln(" * sigma^2 * ")")),
  panel_label = "b)",
  show_legend = TRUE
)

p_combined_mean_CV <- make_combined_plot(
  all_data,
  "mean",
  "CV",
  xlab = expression(bold("Mean")),
  ylab = expression(bold("ln(CV)")),
  panel_label = "c)"
)

(p_combined_var_CV |
    p_combined_mean_var |
    p_combined_mean_CV)



# =========================================================
# DISTRIBUTION
# =========================================================

BLUE_NGM  <- as.data.frame(BLUE_NGM)
BLUE_NaCl <- as.data.frame(BLUE_NaCl)

BLUE_NGM$env  <- "NGM"
BLUE_NaCl$env <- "NaCl"

all_data <- bind_rows(BLUE_NGM, BLUE_NaCl)

# =========================================================
# THEME (Nature / Science style)
# =========================================================

theme_science <- theme_classic(base_size = 16) +
  theme(
    axis.title = element_text(size = 16, face = "bold"),
    axis.text  = element_text(size = 13),
    plot.title = element_text(size = 17, face = "bold", hjust = 0.5),
    legend.title = element_blank(),
    legend.text = element_text(size = 13)
  )

# =========================================================
# COLORS
# =========================================================

cols <- c("NGM" = "#4C72B0",
          "NaCl" = "#C44E52")

# =========================================================
# DENSITY PLOT FUNCTION
# =========================================================

plot_density <- function(data, variable, title){
  
  ggplot(data, aes_string(x = variable, fill = "env", colour = "env")) +
    
    geom_density(alpha = 0.35, linewidth = 1.1) +
    
    scale_fill_manual(values = cols) +
    scale_colour_manual(values = cols) +
    
    labs(title = title,
         x = variable,
         y = "Density") +
    
    theme_science
}

# =========================================================
# PLOTS
# =========================================================

p_mean <- plot_density(all_data,
                       "mean",
                       "Distribution of Mean")

p_var <- plot_density(all_data,
                      "var",
                      "Distribution of Variance")

p_cv <- plot_density(all_data,
                     "CV",
                     "Distribution of CV")

# =========================================================
# DISPLAY
# =========================================================

(p_mean / p_var / p_cv)




# =========================================================
# BOXPLOT 
# =========================================================

all_data$env <- factor(all_data$env,
                       levels = c("NGM", "NaCl"))

plot_box_env <- function(data, variable, panel_label = NULL, ylab = NULL){
  
  y_max <- max(data[[variable]], na.rm = TRUE)
  y_range <- diff(range(data[[variable]], na.rm = TRUE))
  y_pos <- y_max + 0.08 * y_range
  
  ggplot(data, aes(x = env, y = .data[[variable]])) +
    
    # boxplot
    geom_boxplot(aes(fill = env),
                 alpha = 0.65,
                 width = 0.6,
                 outlier.shape = NA) +
    
    # points
    geom_jitter(aes(color = env),
                width = 0.15,
                size = 1.4,
                alpha = 0.35) +
    
    # couleurs + labels propres
    scale_fill_manual(
      values = cols,
      labels = c(
        "NGM"  = "Low-salt",
        "NaCl" = "High-salt"
      )
    ) +
    
    scale_color_manual(
      values = cols,
      labels = c(
        "NGM"  = "Low-salt",
        "NaCl" = "High-salt"
      )
    ) +
    
    scale_x_discrete(labels = c(
      "NGM"  = "Low-salt",
      "NaCl" = "High-salt"
    )) +
    
    # barre de signification
    annotate("segment",
             x = 1,
             xend = 2,
             y = y_pos,
             yend = y_pos,
             linewidth = 0.9) +
    
    annotate("text",
             x = 1.5,
             y = y_pos + 0.03 * y_range,
             label = "***",
             size = 6,
             fontface = "bold") +
    
    # panel label
    annotate("text",
             x = 0.75,
             y = y_max + 0.12 * y_range,
             label = panel_label,
             fontface = "bold",
             size = 6) +
    
    labs(
      x = "Environment",
      y = ylab
    ) +
    
    theme_science +
    
    theme(
      legend.position = "none",
      axis.title = element_text(face = "bold")
    )
}

p_mean_box <- plot_box_env(all_data,
                           "mean",
                           panel_label = "a)",
                           ylab = "Mean")

p_var_box <- plot_box_env(all_data,
                          "var",
                          panel_label = "b)",
                          ylab = expression(bold(ln(sigma^2))))

p_cv_box <- plot_box_env(all_data,
                         "CV",
                         panel_label = "c)",
                         ylab = expression(bold(ln(CV))))

(p_mean_box | p_var_box | p_cv_box)



# =========================================================
# CORRELATIONS
# =========================================================


BLUE_NGM_clean <- BLUE_NGM %>%
  select(pop_label, mean, var, CV) %>%
  rename(
    mean_NGM = mean,
    var_NGM  = var,
    CV_NGM   = CV
  )

BLUE_NaCl_clean <- BLUE_NaCl %>%
  select(pop_label, mean, var, CV) %>%
  rename(
    mean_NaCl = mean,
    var_NaCl  = var,
    CV_NaCl   = CV
  )

BLUE_merge <- inner_join(BLUE_NGM_clean,
                         BLUE_NaCl_clean,
                         by = "pop_label")


BLUE_merge_var <- BLUE_merge %>%
  mutate(env = "comparison") %>%
  select(pop_label, var_NGM, var_NaCl)

ggplot(BLUE_merge, aes(x = var_NGM, y = var_NaCl)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  theme_classic(base_size = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8)+
  labs(
    x = expression("ln(" * sigma[NGM]^2 * ")"),
    y = expression("ln(" * sigma[NaCl]^2 * ")")
  )

p=ggplot(BLUE_merge, aes(x = var_NaCl, y = var_NGM)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
  theme_classic(base_size = 15) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.8)+
  labs(
    x = expression("ln(" * sigma[NGM]^2 * ")"),
    y = expression("ln(" * sigma[NaCl]^2 * ")")
  )

cor_res <- cor.test(BLUE_merge$var_NGM,
                    BLUE_merge$var_NaCl)

cor_res


cor_cv <- cor.test(BLUE_merge$CV_NGM,
                   BLUE_merge$CV_NaCl)
cor_cv

l=ggplot(BLUE_merge, aes(x = CV_NGM, y = CV_NaCl)) +
  
  geom_point(alpha = 0.5, size = 1) +
  
  geom_smooth(method = "lm",
              se = FALSE,
              linewidth = 1.2) +
  
  geom_abline(slope = 1,
              intercept = 0,
              linetype = "dashed",
              linewidth = 0.8) +
  
  theme_classic(base_size = 15) +
  
  labs(
    x = expression(ln(CV)[NGM]),
    y = expression(ln(CV)[NaCl])
  )

(p | l)

make_env_correlation_plot <- function(data,
                                      xvar,
                                      yvar,
                                      xlab,
                                      ylab,
                                      panel_label = NULL){
  
  x_min <- min(data[[xvar]], na.rm = TRUE)
  x_range <- diff(range(data[[xvar]], na.rm = TRUE))
  
  y_max <- max(data[[yvar]], na.rm = TRUE)
  y_range <- diff(range(data[[yvar]], na.rm = TRUE))
  
  ggplot(data,
         aes_string(x = xvar,
                    y = yvar)) +
    
    geom_point(size = 0.75,
               alpha = 0.4) +
    
    geom_smooth(method = "lm",
                se = FALSE,
                linewidth = 1.4,
                colour = "grey30")+
    
    annotate("text",
             x = x_min - 0.02 * x_range,
             y = y_max + 0.02 * y_range,
             label = panel_label,
             fontface = "bold",
             size = 6,
             hjust = 0,
             vjust = 1) +
    
    labs(
      x = xlab,
      y = ylab
    ) +
    
    theme_classic(base_size = 16) +
    
    theme(
      axis.title = element_text(face = "bold",
                                size = 20),
      axis.text = element_text(size = 15),
      axis.line = element_line(linewidth = 0.7),
      axis.ticks = element_line(linewidth = 0.7),
      plot.margin = margin(10,10,10,10)
    )
}

p_var <- make_env_correlation_plot(
  BLUE_merge,
  "var_NGM",
  "var_NaCl",
  xlab = expression(bold(ln(sigma^2)["Low-salt"])),
  ylab = expression(bold(ln(sigma^2)["High-salt"])),
  panel_label = "a)"
)

p_cv <- make_env_correlation_plot(
  BLUE_merge,
  "CV_NGM",
  "CV_NaCl",
  xlab = expression(bold(ln(CV)["Low-salt"])),
  ylab = expression(bold(ln(CV)["High-salt"])),
  panel_label = "b)"
)

p_var | p_cv


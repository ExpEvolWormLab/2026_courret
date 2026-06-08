library(data.table)
library(stringr)
library(dplyr)

# =========================================================
# FUNCTION TO EXTRACT H2
# =========================================================

extract_h2 <- function(file, env_name){
  
  txt <- readLines(file)
  
  trait_names <- c(
    "area.F",
    "log_area.F.var",
    "log_CV_areaF"
  )
  
  results <- list()
  
  current_trait <- NA
  chr_count <- 0
  
  for(i in seq_along(txt)){
    
    line <- txt[i]
    
    # -----------------------------------------------------
    # DETECT TRAIT
    # -----------------------------------------------------
    
    for(tr in trait_names){
      
      if(grepl(paste0('"', tr, '"'), line, fixed = TRUE)){
        
        current_trait <- tr
        chr_count <- 0
      }
    }
    
    # -----------------------------------------------------
    # DETECT H2
    # -----------------------------------------------------
    
    if(grepl("GridLMM_posterior h2s:", line)){
      
      chr_count <- chr_count + 1
      
      # line containing values
      value_line <- txt[i + 2]
      
      # remove [1,]
      value_line <- gsub("\\[1,\\]", "", value_line)
      
      # extract numeric values
      vals <- as.numeric(
        str_extract_all(value_line, "[0-9]+\\.?[0-9]*")[[1]]
      )
      
      # first value = pop_label
      h2 <- vals[1]
      
      results[[length(results) + 1]] <- data.table(
        Environment = env_name,
        Trait = current_trait,
        Chromosome = chr_count,
        H2 = h2
      )
    }
  }
  
  rbindlist(results)
}

# =========================================================
# LOAD FILES
# =========================================================

h2_NGM <- extract_h2(
  "data/output_grid.txt",
  "NGM"
)

h2_NaCl <- extract_h2(
  "data/output_grid_NaCl.txt",
  "NaCl"
)

# combine
h2_all <- rbind(h2_NGM, h2_NaCl)

print(h2_all)

# =========================================================
# SUMMARY
# =========================================================

h2_summary <- h2_all %>%
  group_by(Environment, Trait) %>%
  summarise(
    mean_H2 = mean(H2),
    sd_H2 = sd(H2),
    min_H2 = min(H2),
    max_H2 = max(H2),
    n_chr = n()
  )

print(h2_summary)




h2_summary <- h2_all %>%
  group_by(Environment, Trait) %>%
  summarise(
    mean_H2 = mean(H2),
    sd_H2 = sd(H2),
    n = n(),
    se = sd_H2 / sqrt(n),
    ci_low = mean_H2 - 1.96 * se,
    ci_high = mean_H2 + 1.96 * se,
    .groups = "drop"
  )

str(h2_summary)



values = c(
  "NGM"  = col_NGM,
  "NaCl" = col_NaCl
)

h2_summary$Trait <- as.character(h2_summary$Trait)

h2_summary$Trait[h2_summary$Trait == "area.F"] <- "Mean"
h2_summary$Trait[h2_summary$Trait == "log_area.F.var"] <- "ln(sigma2)"
h2_summary$Trait[h2_summary$Trait == "log_CV_areaF"] <- "ln(CV)"

ggplot(h2_summary,
       aes(x = Trait,
           y = mean_H2,
           color = Environment,
           group = Environment)) +
  
  geom_line(linewidth = 1, alpha = 0.8) +
  
  geom_point(alpha = 0.5, size = 2.5,
             position = position_dodge(0)) +
  
  geom_errorbar(aes(ymin = ci_low,
                    ymax = ci_high),
                width = 0.12,
                linewidth = 0.8,
                position = position_dodge(0)) +
  
  scale_color_manual(values = c(
    "NGM" = "#4C72B0",
    "NaCl" = "#C44E52"
  )) +
  
  labs(
    x = NULL,
    y = expression(bold(h^2))
  ) +
  
  scale_x_discrete(labels = function(x) parse(text = x)) +
  
  theme_classic(base_size = 15) +
  
  theme(
    legend.position = "top",
    legend.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    axis.title.y = element_text(face = "bold")
  )




ggplot(h2_summary,
       aes(x = Trait,
           y = mean_H2,
           colour = Environment,
           group = Environment)) +
  
  geom_line(linewidth = 1.4,
            alpha = 0.9) +
  
  geom_point(size = 2.5,
             alpha = 0.8) +
  
  geom_errorbar(aes(ymin = ci_low,
                    ymax = ci_high),
                width = 0.08,
                linewidth = 0.8) +
  
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
    x = NULL,
    y = expression(bold(h^2)),
    colour = NULL
  ) +
  
  scale_x_discrete(labels = c(
    "Mean" = "Mean",
    "ln(CV)" = "ln(CV)",
    "ln(sigma2)" = expression(ln(sigma^2))
  )) +
  
  theme_classic(base_size = 16) +
  
  theme(
    
    legend.position = "top",
    legend.justification = "center",
    legend.direction = "horizontal",
    
    legend.text = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.title = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.text.y = element_text(
      size = 13
    ),
    
    axis.line = element_line(
      linewidth = 0.7
    ),
    
    axis.ticks = element_line(
      linewidth = 0.7
    ),
    
    plot.margin = margin(10, 10, 10, 10)
  )


###########################

h2_summary$Trait <- factor(
  h2_summary$Trait,
  levels = c("Mean", "ln(sigma2)", "ln(CV)")
)

ggplot(h2_summary,
       aes(x = Trait,
           y = mean_H2,
           colour = Environment,
           group = Environment)) +
  
  geom_line(linewidth = 1.4, alpha = 0.9) +
  
  geom_point(size = 2.6, alpha = 0.85) +
  
  geom_errorbar(aes(ymin = ci_low,
                    ymax = ci_high),
                width = 0.08,
                linewidth = 0.8) +
  
  scale_colour_manual(
    values = c(
      "NGM"  = "#4C72B0",
      "NaCl" = "#C44E52"
    ),
    labels = c(
      "NGM"  = "Low-salt",
      "NaCl" = "High-salt"
    )
  ) +
  
  scale_x_discrete(labels = c(
    "Mean" = expression(bold("Mean")),
    "ln(CV)" = expression(bold(ln(CV))),
    "ln(sigma2)" = expression(bold(ln(sigma^2)))
  ))+
  
  labs(
    x = NULL,
    y = expression(bold(H^2)),
    colour = NULL
  ) +
  
  scale_y_continuous(
    
    limits = c(0, NA)
    
  ) +
  
  theme_classic(base_size = 16) +
  
  theme(
    legend.position = "top",
    
    legend.text = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.title.y = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text.x = element_text(
      face = "bold",
      size = 16
    ),
    
    axis.text.y = element_text(size = 13),
    
    axis.line = element_line(linewidth = 0.7),
    axis.ticks = element_line(linewidth = 0.7)
  )

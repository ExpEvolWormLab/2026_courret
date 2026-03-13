# retirer valeurs extremes (6 points de la queue histogram des BLUE var)

library(dplyr)

phen_data_NGM_corr <- phen_data_NGM %>%
  filter(
    !(pop_label == "CA150L23") &
      !(pop_label == "LR3L84") &
      !(pop_label == "LR3L105") &
      !(pop_label == "LR3L107") &
      !(pop_label == "GA150L22") &
      !(pop_label == "CA1100L77" & block_id == "F033")
  )

phen_data_NGM_corr %>%
  filter(pop_label %in% c(
    "CA150L23",
    "LR3L84",
    "LR3L105",
    "LR3L107",
    "GA150L22",
    "CA1100L77"
  )) %>%
  select(pop_label, block_id) #les points ont bien été enelevées

#distribution
ggplot(phen_data_NGM_corr,
       aes(x = session_id, y = area.F.var, color = session_id)) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  theme_classic()

# modele 
var_model_corr_queue <- lmer(
  area.F.var ~
    area.F +
    logD +
    temperature +
    rel_humidity +
    pop_label +
    (1 | session_id/block_id),
  data = phen_data_NGM_corr,
  REML = TRUE
)

summary(var_model_corr_queue)
Anova(var_model_corr_queue)

check_model(var_model_corr_queue)

#emmeans

pheno_line_var_corr <- emmeans(var_model_corr_queue, "pop_label")
pheno_var_corr <- as.data.frame(pheno_line_var_corr)

pheno_var_corr <- pheno_var_corr %>%
  mutate(
    population = ifelse(
      substr(pop_label, 1, 3) %in% c("SMR", "SLR"),
      substr(pop_label, 1, 3),
      substr(pop_label, 1, 2)
    ),
    population2 = str_extract(pop_label, ".*(?=L)")
  )

write.table(pheno_var_corr, file = "pheno_var_corr_queue.txt", row.names = FALSE, col.names = T, quote = FALSE)


df_plot <- merge(
  pheno_var[, c("pop_label","emmean")],
  pheno_var_corr[, c("pop_label","emmean")],
  by = "pop_label",
  suffixes = c("_raw","_corr")
)

plot(df_plot$emmean_raw, df_plot$emmean_corr)
abline(0,1,col="red")

cor(df_plot$emmean_raw, df_plot$emmean_corr)

ggplot(pheno_var_corr, aes(x = emmean)) +
  geom_histogram(aes(y = ..density..), 
                 binwidth = 0.000002, fill = "skyblue", color = "black", alpha = 0.6) +
  geom_density(color = "red", size = 1) +
  theme_classic() +
  labs(x = "Estimated var", y = "Density",
       title = "Distribution of emmean with density curve")

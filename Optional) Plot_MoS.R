setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data")

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Chargement et nettoyage des données
data <- read.csv2("ACP_Clustering_DATA.csv", sep = ";", dec = ".", check.names = TRUE)

# Définition des groupes de variables (noms générés par check.names=TRUE)
vars_L0 <- c("Mean_MoS.ML.HS...L0.", "Mean_MoS.ML.Stance...L0.", 
             "Mean_MoS.AP.HS...L0.", "Mean_MoS.AP.Stance...L0.")

vars_mm <- c("Mean_MoS.ML.HS..mm.", "Mean_MoS.ML.Stance..mm.", 
             "Mean_MoS.AP.HS..mm.", "Mean_MoS.AP.Stance..mm.")

# 2. Fonction de création de graphique
create_mos_plot_with_points <- function(df, variables, title_label, unit_label) {
  
  # 1. Préparation des données individuelles (pour les points)
  df_raw <- df %>%
    dplyr::select(AgeGroup, Surface, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "Variable", values_to = "Value") %>%
    mutate(
      Value = as.numeric(as.character(Value)),
      AgeGroup = factor(AgeGroup, 
                        levels = c("JeunesEnfants", "Enfants", "Adolescents", "Adultes"),
                        labels = c("Young Children", "Children", "Adolescents", "Adults")),
      Surface = factor(Surface, 
                       levels = c("Plat", "Medium", "High"),
                       labels = c("Even", "Medium", "High")),
      Variable = gsub("Mean_MoS.", "", Variable),
      Variable = gsub("\\.", " ", Variable)
    ) %>%
    filter(!is.na(Value))
  
  # 2. Résumé statistique (pour les barres et erreurs)
  df_summary <- df_raw %>%
    group_by(AgeGroup, Surface, Variable) %>%
    summarise(
      mean_val = mean(Value, na.rm = TRUE),
      sd_val = sd(Value, na.rm = TRUE),
      .groups = 'drop'
    )
  
  # 3. Plot
  ggplot() +
    # REMPLACEMENT : Boxplots (Distribution)
    # On utilise df_raw car le boxplot a besoin de toutes les données, pas juste du résumé
    geom_boxplot(data = df_raw, aes(x = AgeGroup, y = Value, fill = Surface),
                 position = position_dodge(0.85), width = 0.7, 
                 alpha = 0.7, outlier.shape = NA, color = "black") +
    
    # Points individuels (On garde geom_jitter mais on s'assure qu'il utilise df_raw)
    geom_jitter(data = df_raw, aes(x = AgeGroup, y = Value, group = Surface),
                position = position_dodge(0.85), 
                size = 1, alpha = 0.4, color = "black") +
    
    facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
    # Note : Vérifiez que vos couleurs correspondent aux labels (Level, Medium, Uneven)
    scale_fill_manual(values = c("Even" = "blue", "Medium" = "green", "High" = "red")) + 
    labs(
      title = paste("Margin Of Stability:", title_label),
      y = paste("Value", unit_label), # On enlève "Mean ± SD" car c'est une distribution
      x = "Age Group",
      fill = "Surface Condition"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 11, face = "italic"), 
      strip.text = element_text(face = "bold", size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

# 3. Génération des deux figures
fig_L0 <- create_mos_plot_with_points(data, vars_L0, "Normalized data (%L0)", "(%L0)")
ggsave("MoS_Normalise_L0.pdf", width = 10, height = 8)
ggsave("MoS_Normalise_L0.png", plot = fig_L0, width = 10, height = 7, dpi = 300)

fig_mm <- create_mos_plot_with_points(data, vars_mm, "Raw data (mm)", "(mm)")
ggsave("MoS_Raw_mm.pdf", width = 10, height = 8)
ggsave("MoS_Raw_mm.png", plot = fig_mm, width = 10, height = 7, dpi = 300)

# Affichage (dans RStudio, elles apparaîtront l'une après l'autre)
print(fig_L0)
print(fig_mm)
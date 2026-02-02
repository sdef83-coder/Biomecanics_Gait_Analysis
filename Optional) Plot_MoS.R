setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data")

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Chargement
data <- read.csv2("ACP_Clustering_DATA.csv", sep = ";", dec = ".", check.names = TRUE)

# Définition des groupes de variables (noms générés par check.names=TRUE)
vars_L0 <- c("Mean_MoS.ML.HS...L0.", "Mean_MoS.ML.Stance...L0.",
             "Mean_MoS.AP.HS...L0.", "Mean_MoS.AP.Stance...L0.")

vars_mm <- c("Mean_MoS.ML.HS..mm.", "Mean_MoS.ML.Stance..mm.",
             "Mean_MoS.AP.HS..mm.", "Mean_MoS.AP.Stance..mm.")

# 2. Fonction : Surface en abscisse
create_mos_plot_with_points_surfaceX <- function(df, variables, title_label, unit_label) {
  
  df_raw <- df %>%
    dplyr::select(AgeGroup, Surface, all_of(variables)) %>%
    pivot_longer(cols = all_of(variables), names_to = "Variable", values_to = "Value") %>%
    mutate(
      Value = as.numeric(as.character(Value)),
      AgeGroup = factor(
        AgeGroup,
        levels = c("JeunesEnfants", "Enfants", "Adolescents", "Adultes"),
        labels = c("Young Children", "Children", "Adolescents", "Adults")
      ),
      Surface = factor(
        Surface,
        levels = c("Plat", "Medium", "High"),
        labels = c("Even", "Medium", "High")
      ),
      Variable = gsub("Mean_MoS\\.", "", Variable),
      Variable = gsub("\\.", " ", Variable)
    ) %>%
    filter(!is.na(Value))
  
  ggplot(df_raw, aes(x = Surface, y = Value, fill = AgeGroup)) +
    
    # Boxplots (distribution)
    geom_boxplot(
      position = position_dodge(0.85),
      width = 0.7,
      alpha = 0.7,
      outlier.shape = NA,
      color = "black"
    ) +
    
    # Points individuels (même style que tes autres figures)
    geom_jitter(
      aes(group = AgeGroup),
      position = position_dodge(0.85),
      size = 1,
      alpha = 0.4,
      color = "black"
    ) +
    
    facet_wrap(~ Variable, scales = "free_y", ncol = 2) +
    
    # Couleurs par âge (à ajuster si tu as une palette standard)
    scale_fill_manual(values = c(
      "Young Children" = "blue",
      "Children"       = "orange",
      "Adolescents"    = "green",
      "Adults"         = "purple"
    )) +
    
    labs(
      title = paste("Margin Of Stability:", title_label, "across surfaces and age groups"),
      y = paste("Value", unit_label),
      x = "Surface",
      fill = "Age Group"
    ) +
    
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      strip.text = element_text(face = "bold", size = 10),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
}

# 3. Génération des deux figures
fig_L0 <- create_mos_plot_with_points_surfaceX(data, vars_L0, "Normalized data (%L0)", "(%L0)")
fig_mm <- create_mos_plot_with_points_surfaceX(data, vars_mm, "Raw data (mm)", "(mm)")

# 4. Affichage
print(fig_L0)
print(fig_mm)

# 5. Sauvegardes (explicites)
ggsave("MoS_Normalise_L0.pdf", plot = fig_L0, width = 10, height = 8)
ggsave("MoS_Normalise_L0.png", plot = fig_L0, width = 10, height = 7, dpi = 300)

ggsave("MoS_Raw_mm.pdf", plot = fig_mm, width = 10, height = 8)
ggsave("MoS_Raw_mm.png", plot = fig_mm, width = 10, height = 7, dpi = 300)

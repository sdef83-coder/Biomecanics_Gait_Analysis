setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data")

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Chargement des données
data <- read.csv2("ACP_Clustering_DATA.csv", sep = ";", dec = ".", check.names = TRUE)

# 2. Fonction GVI en version Boxplot (Surface en X)
create_gvi_boxplot_surfaceX <- function(df) {
  
  # Préparation des données
  df_plot <- df %>%
    dplyr::select(AgeGroup, Surface, Mean_GVI..ua.) %>%
    mutate(
      Value = as.numeric(as.character(Mean_GVI..ua.)),
      AgeGroup = factor(
        AgeGroup,
        levels = c("JeunesEnfants", "Enfants", "Adolescents", "Adultes"),
        labels = c("Young Children", "Children", "Adolescents", "Adults")
      ),
      Surface = factor(
        Surface,
        levels = c("Plat", "Medium", "High"),
        labels = c("Even", "Medium", "High")
      )
    ) %>%
    filter(!is.na(Value))
  
  # Création du graphique
  ggplot(df_plot, aes(x = Surface, y = Value, fill = AgeGroup)) +
    
    # Ligne de référence à 100
    geom_hline(yintercept = 100, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    
    # BOXPLOT
    geom_boxplot(
      position = position_dodge(0.85),
      width = 0.7,
      alpha = 0.7,
      outlier.shape = NA,
      color = "black"
    ) +
    
    # Points individuels
    geom_jitter(
      aes(group = AgeGroup),
      position = position_dodge(0.85),
      size = 1.2,
      alpha = 0.4,
      color = "black"
    ) +
    
    # Couleurs par AgeGroup
    scale_fill_manual(values = c(
      "Young Children" = "blue",
      "Children"       = "orange",
      "Adolescents"    = "green",
      "Adults"         = "purple"
    )) +
    
    # Labels
    labs(
      title = "Gait Variability Index (GVI) across surfaces and age groups",
      subtitle = "Boxplots show median and quartiles. Dashed line = Adult Level reference.",
      y = "GVI Score (ua)",
      x = "Surface",
      fill = "Age Group"
    ) +
    
    # Thème
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
}

# 3. Génération et affichage
fig_gvi_box <- create_gvi_boxplot_surfaceX(data)
print(fig_gvi_box)

# 4. Enregistrement haute qualité (explicite)
ggsave("Figure_GVI.pdf", plot = fig_gvi_box, width = 10, height = 8)
ggsave("Figure_GVI.png", plot = fig_gvi_box, width = 10, height = 7, dpi = 300)

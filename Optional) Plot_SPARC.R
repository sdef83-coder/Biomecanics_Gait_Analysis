setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data")

library(ggplot2)
library(dplyr)
library(tidyr)

# 1. Chargement des données
data <- read.csv2("ACP_Clustering_DATA.csv", sep = ";", dec = ".", check.names = TRUE)

# 2. Fonction SPARC COM Magnitude en version Boxplot
create_sparc_boxplot <- function(df) {
  
  # Préparation des données avec dplyr:: spécifique
  df_plot <- df %>%
    dplyr::select(AgeGroup, Surface, Mean_COM.SPARC.Magnitude..ua.) %>%
    mutate(
      Value = as.numeric(as.character(Mean_COM.SPARC.Magnitude..ua.)),
      # Traduction des étiquettes
      AgeGroup = factor(AgeGroup, 
                        levels = c("JeunesEnfants", "Enfants", "Adolescents", "Adultes"),
                        labels = c("Young Children", "Children", "Adolescents", "Adults")),
      Surface = factor(Surface, 
                       levels = c("Plat", "Medium", "High"),
                       labels = c("Even", "Medium", "High"))
    ) %>%
    filter(!is.na(Value))
  
  # Création du graphique
  ggplot(df_plot, aes(x = AgeGroup, y = Value, fill = Surface)) +
    
    # BOXPLOT : outlier.shape = NA pour éviter de doubler les points du jitter
    geom_boxplot(position = position_dodge(0.85), width = 0.7, 
                 alpha = 0.7, outlier.shape = NA, color = "black") +
    
    # Points individuels (Jitter)
    geom_jitter(aes(group = Surface),
                position = position_dodge(0.85), 
                size = 1.2, alpha = 0.4, color = "black") +
    
    # Couleurs cohérentes avec tes autres graphiques
    scale_fill_manual(values = c("Even" = "blue", "Medium" = "green", "High" = "red")) +
    
    # Labels
    labs(
      title = "Spectral Arc Length (SPARC) across age groups and surfaces",
      subtitle = "Magnitude of Centre Of Mass measurement",
      y = "SPARC Magnitude Score (ua)",
      x = "Age Group",
      fill = "Surface"
    ) +
    
    # Thème
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
}

# 3. Génération et affichage
fig_sparc_box <- create_sparc_boxplot(data)
print(fig_sparc_box)

# Sauvegarde
ggsave("SPARC_COM_Magnitude.pdf", width = 10, height = 8)
ggsave("SPARC_COM_Magnitude_Boxplot.png", plot = fig_sparc_box, width = 10, height = 7, dpi = 300)
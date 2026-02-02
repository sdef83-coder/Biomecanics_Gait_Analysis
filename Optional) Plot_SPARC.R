setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data")

library(ggplot2)
library(dplyr)

# 1. Chargement des données
data <- read.csv2("ACP_Clustering_DATA.csv", sep = ";", dec = ".", check.names = TRUE)

# 2. Préparation des données
df_plot <- data %>%
  dplyr::select(AgeGroup, Surface, Mean_COM.SPARC.Magnitude..ua.) %>%
  mutate(
    Value = as.numeric(as.character(Mean_COM.SPARC.Magnitude..ua.)),
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

# 3. Graphique : Surface en abscisse, boxplots + points
fig_sparc_surface_x <- ggplot(df_plot, aes(x = Surface, y = Value, fill = AgeGroup)) +
  
  # Boxplots
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
  
  # Couleurs des groupes d’âge
  scale_fill_manual(values = c(
    "Young Children" = "blue",
    "Children"       = "orange",
    "Adolescents"    = "green",
    "Adults"         = "purple"
  )) +
  
  # Labels
  labs(
    title = "Spectral Arc Length (SPARC) across surfaces and age groups",
    subtitle = "Magnitude of Centre Of Mass measurement",
    y = "SPARC Magnitude Score (ua)",
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

# 4. Affichage
print(fig_sparc_surface_x)

# 5. Sauvegarde
ggsave(
  "SPARC_COM_Magnitude.png",
  plot = fig_sparc_surface_x,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  "SPARC_COM_Magnitude.pdf",
  plot = fig_sparc_surface_x,
  width = 10,
  height = 8
)

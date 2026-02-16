setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/PhysicalActi_et_VAS")

# =========================================================
# 0) LIBRAIRIES
# - dplyr / tidyverse : manipulation de données (group_by, mutate, summarise, etc.)
# - readr            : fonctions d’import/export (ex: write_csv)
# - ggplot2          : visualisations (boxplots, densités, export)
# - gtsummary / gt   : tableaux descriptifs format "publication" + export (PNG/PDF)
#
# Remarque : dplyr est chargé deux fois (sans conséquence fonctionnelle, juste redondant).
# =========================================================
library(dplyr)
library(readr)
library(ggplot2)
library(gtsummary)
library(gt)
library(dplyr)
library(tidyverse)

# =========================================================
# 1) IMPORT DES DONNÉES
# - read.csv(file.choose()) : sélection manuelle d’un fichier .csv avec séparateur ";"
# - colnames() : contrôle rapide des noms de colonnes (diagnostic import)
# =========================================================
df <- read.csv(file.choose(), sep = ";")

# Vérifier le nom des colonnes
colnames(df)

# =========================================================
# 2) FACTEURS + ORDRE DES GROUPES D’ÂGE
# - Conversion de AgeGroup en facteur ordonné pour contrôler l’ordre d’affichage
#   et garantir la cohérence dans les tables/figures (Children < Adolescent < Adult)
# =========================================================
df$AgeGroup <- factor(df$AgeGroup,levels = c("Children", "Adolescent", "Adult"),ordered = TRUE)

# =========================================================
# 3) PALETTE DE COULEURS PAR GROUPE D’ÂGE
# - Définition d’une correspondance "AgeGroup -> couleur"
# - Utile pour garder une cohérence graphique entre figures
# =========================================================
cols_age <- c(
  "Children"   = "orange",  # orange
  "Adolescent" = "green",  # vert
  "Adult"      = "purple"   # violet
)

# =========================================================
# 4) STANDARDISATION : CALCUL DU Z-SCORE PAR GROUPE D’ÂGE
# Objectif : rendre les scores comparables entre groupes malgré des questionnaires différents
#            (ex: PAQ-C, PAQ-A, GPAQ) en centrant-réduisant *au sein* de chaque AgeGroup.
#
# Méthode :
# - group_by(AgeGroup) : calcul des stats séparément par groupe d’âge
# - mutate(Zscore = (RawScore - mean)/sd) : z-score intra-groupe
# - na.rm=TRUE : ignore les valeurs manquantes pour mean/sd
# - ungroup() : évite que les opérations suivantes restent "groupées"
# =========================================================
df <- df %>%
  group_by(AgeGroup) %>%
  mutate(
    Zscore = (RawScore - mean(RawScore, na.rm = TRUE)) /
      sd(RawScore, na.rm = TRUE)
  ) %>%
  ungroup()

# Inspection visuelle de la table (optionnel)
View(df)

# =========================================================
# 5) CONTRÔLE QUALITÉ DU Z-SCORE
# Objectif : vérifier que, par AgeGroup, la moyenne ~ 0 et l’écart-type ~ 1
# - summarise(mean_z, sd_z) : diagnostic attendu après centrage-réduction
# =========================================================
df %>%
  group_by(AgeGroup) %>%
  summarise(
    mean_z = mean(Zscore, na.rm = TRUE),
    sd_z   = sd(Zscore, na.rm = TRUE)
  )

# =========================================================
# 6) EXPORT DU DATASET STANDARDISÉ
# - write_csv : export en .csv (format readr) dans le répertoire courant (getwd())
# =========================================================
write_csv(df, "PhysicalActivity_Zscored.csv")

# =========================================================
# 7) TABLEAU DESCRIPTIF (RAW + Z-SCORE) PAR GROUPE D’ÂGE
# Objectif : produire un tableau "Mean ± SD" par AgeGroup, prêt à exporter.
#
# Notes techniques :
# - df_pa <- as.data.frame(df) : sécurise le type (data.frame)
# - dplyr::select : évite ambiguïté avec d’autres packages
# - tbl_summary(by=AgeGroup) : statistiques séparées par groupe
# - modify_* : mise en forme (en-têtes, gras, notes de bas de tableau)
# =========================================================

# On s'assure que df est bien considéré comme un tableau de données
# et on le renomme pour éviter le conflit avec la fonction df() de R
df_pa <- as.data.frame(df)

tab_phys_activity <- df_pa %>%
  # On utilise dplyr:: pour être sûr d'appeler la bonne fonction
  dplyr::select(AgeGroup, RawScore, Zscore) %>%
  tbl_summary(
    by = AgeGroup,
    label = list(
      RawScore ~ "Physical Activity (Raw Score)",
      Zscore   ~ "Physical Activity (Z-score)"
    ),
    statistic = list(all_continuous() ~ "{mean} ± {sd}"),
    digits = list(all_continuous() ~ 2),
    missing = "no"
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "**Variable**") %>%
  modify_spanning_header(all_stat_cols() ~ "**Age groups**") %>%
  modify_footnote(all_stat_cols() ~ "Mean ± SD; Children (PAQ-C), Adolescent (PAQ-A), Adult (GPAQ)")

# Afficher le tableau dans la console
print(tab_phys_activity)

# =========================================================
# 8) EXPORT DU TABLEAU DESCRIPTIF (GT) EN PDF + PNG
# - as_gt() : conversion gtsummary -> gt
# - gtsave() : export (PDF + PNG)
# Notes :
# - PDF : peut nécessiter Chrome/LaTeX selon config
# - PNG : expand ajoute une marge blanche autour du tableau
# =========================================================
gt_table <- tab_phys_activity %>% as_gt()

# 2. Sauvegarde en PDF
# Note : Nécessite que Chrome ou une installation LaTeX soit accessible par R
gt::gtsave(gt_table, "Table_PhysActivity.pdf")

# 3. Sauvegarde en PNG
# L'argument expand permet d'ajouter une petite marge blanche autour du tableau
gt::gtsave(gt_table, "Table_PhysActivity.png", expand = 10)

# Message de confirmation (diagnostic)
message("Tableau exporté avec succès dans : ", getwd())

# =========================================================
# 9) VISUALISATIONS
# - Redéfinition d’une palette hex (plus "publication-friendly")
# - Figure 1 : boxplot + points individuels (jitter) du z-score par AgeGroup
# - Figure 2 : densités des z-scores par AgeGroup
# - Export : PNG (600 dpi, fond blanc) + PDF (vectoriel)
# =========================================================

# Palette couleur en hex (cohérence visuelle, meilleure reproductibilité)
cols_age <- c(
  "Children"   = "#E69F00",  # orange
  "Adolescent" = "#009E73",  # vert
  "Adult"      = "#7B3294"   # violet
)

# ---------------------------------------------------------
# 9.1) BOXPLOT + POINTS INDIVIDUELS (Zscore ~ AgeGroup)
# - geom_boxplot : distribution globale (médiane, IQR), outliers masqués (outlier.shape=NA)
# - geom_jitter  : points individuels superposés avec transparence
# - geom_hline   : ligne de référence z=0
# - scale_fill_manual : couleurs fixes par groupe
# - theme_minimal + theme : mise en forme proche du style de tes autres figures
# ---------------------------------------------------------
p_box <- ggplot(df, aes(x = AgeGroup, y = Zscore, fill = AgeGroup)) +
  # Boxplot avec contour noir et légère transparence
  geom_boxplot(
    position = position_dodge(0.85), 
    width = 0.7,
    alpha = 0.7, 
    outlier.shape = NA, 
    color = "black"
  ) +
  # Points individuels : jitter avec contour noir et fond transparent
  geom_jitter(
    aes(group = AgeGroup),
    position = position_dodge(0.85),
    size = 1.2, 
    alpha = 0.3, 
    color = "black"
  ) +
  # Ligne de référence à 0
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.6, color = "grey40") +
  # Application des couleurs
  scale_fill_manual(values = cols_age) +
  labs(
    title = "Standardized Physical Activity (z-score) by age group",
    x = "Age group",
    y = "Physical Activity (z-score)",
    fill = "Age group"
  ) +
  # Thème minimaliste identique au script 11
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# Affichage dans la fenêtre Plots
p_box

# Export PNG haute résolution (publication/rapport)
ggsave(
  filename = "PA_Zscore_by_AgeGroup.png",
  plot = p_box,
  width = 18, height = 12, units = "cm",
  dpi = 600,
  bg = "white"
)

# Export PDF (vectoriel)
ggsave(
  filename = "PA_Zscore_by_AgeGroup.pdf",
  plot = p_box,
  width = 18, height = 12, units = "cm"
)

# ---------------------------------------------------------
# 9.2) DENSITÉ DES Z-SCORES PAR GROUPE
# - geom_density : forme de la distribution par AgeGroup
# - geom_vline(x=0) : référence z=0
# - Export PNG + PDF comme pour le boxplot
# ---------------------------------------------------------
p_density <- ggplot(df, aes(x = Zscore, fill = AgeGroup)) +
  geom_density(alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = cols_age) +
  labs(
    title = "Distribution of physical activity z-scores",
    x = "Physical activity (z-score)",
    y = "Density"
  ) +
  theme_minimal()

# Export PNG haute résolution
ggsave(
  filename = "PA_Zscore_Distribution.png",
  plot = p_density,
  width = 18, height = 12, units = "cm",
  dpi = 600,
  bg = "white"
)

# Export PDF (vectoriel)
ggsave(
  filename = "PA_Zscore_Distribution.pdf",
  plot = p_density,
  width = 18, height = 12, units = "cm"
)

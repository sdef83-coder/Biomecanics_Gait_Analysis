setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Fig/Clustering")

# ---------------------------
# I - Stats entre les clusters (chi², t-test, visualisation Age - Sexe - Surface)
#----------------------------

# 1. Chargement
library(tidyverse)
library(rstatix)
library(ggalluvial)
library(gridExtra)
library(grid)
library(gtable)
library(fmsb)




# 2. Importation
df_clust <- read.csv("DATA_FOR_R_GLOBAL_20260123_1409.csv", sep = ";", check.names = FALSE)
df_meta <- read.csv(file.choose(), sep = ";", check.names = FALSE)
# Aller chercher le fichier "participant.metadonnees" dans le dossier stats LMM - stats descriptives

View(df_meta)
View(df_clust)



# 3. Jointure (Lien entre l'ID et le Sexe/Age)
# On ne garde que Participant, AgeMonths et Sex du fichier meta pour éviter les doublons
df_meta_clean <- df_meta %>%
  select(Participant, AgeMonths, Sex) 

df <- left_join(df_clust, df_meta_clean, by = "Participant")


# 3.1 Préparation des variables
# On récupère les colonnes numériques de df_clust pour les boxplots
vars_marche <- colnames(df_clust)[sapply(df_clust, is.numeric)]
vars_marche <- vars_marche[!vars_marche %in% c("ClusterID", "AgeMonths")]

# Df format long pour fig descriptives
df_long <- df %>%
  pivot_longer(cols = all_of(vars_marche), 
               names_to = "Variable", 
               values_to = "Valeur")

# 3.2 Création des Boxplots Individuels avec nettoyage des noms
if(!dir.exists("Boxplots_Individuels")) dir.create("Boxplots_Individuels")

# Vos couleurs manuelles
mes_couleurs <- c("1" = "#3498db", "2" = "#e74c3c")

for (v in vars_marche) {
  
  # --- LOGIQUE DE NETTOYAGE DU NOM ---
  clean_label <- v
  
  # 1. On retire "Mean_"
  clean_label <- gsub("^Mean_", "", clean_label)
  
  # 2. Cas spécifique : Norm Gait Speed -> unité (ua)
  if (grepl("Norm Gait Speed", clean_label, ignore.case = TRUE)) {
    clean_label <- "Norm Gait Speed (ua)"
  }
  
  # 3. On remplace "CV_" par "C.V. " et on force l'unité (%)
  # (Sauf si c'est déjà traité par le cas Norm Gait Speed plus haut)
  if (startsWith(v, "CV_")) {
    clean_label <- gsub("^CV_", "C.V. ", clean_label)
    clean_label <- paste0(gsub(" \\(.*\\)", "", clean_label), " (%)")
  }
  
  # 4. On remplace "SI_" par "S.I. " et on retire TOUTE unité
  if (startsWith(v, "SI_")) {
    clean_label <- gsub("^SI_", "S.I. ", clean_label)
    clean_label <- gsub(" \\(.*\\)", "", clean_label)
  }
  
  # --- CONSTRUCTION DU GRAPHIQUE ---
  p <- ggplot(df, aes(x = factor(ClusterID), y = .data[[v]], fill = factor(ClusterID))) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1) +
    scale_fill_manual(values = mes_couleurs) + 
    labs(title = clean_label, 
         x = "Cluster ID", 
         y = clean_label, 
         fill = "Cluster") + # Note 'caption' retirée ici
    theme_classic() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
      axis.title.y = element_text(size = 10)
    )
  
  # Sauvegarde
  file_name <- gsub("[^[:alnum:]]", "_", v)
  ggsave(paste0("Boxplots_Individuels/Boxplot_", file_name, ".png"), 
         plot = p, width = 12, height = 10, units = "cm")
}
# Boxplots enregistrées dans "Clustering"




# 4 - TABLEAU RÉCAPITULATIF (Médiane [IQR]) & TESTS STATISTIQUES
# -----------------------------------------------------------

# 4.1. Sélection automatique des variables
vars_exclude <- c("Participant", "ClusterID", "Condition", "Sex", "AgeMonths", 
                  "AgeGroup", "AgeGroup.x", "AgeGroup.y", "AgeMonths.x", 
                  "AgeMonths.y", "Sex_EN", "AgeGroup_EN", "Condition_EN")

vars_gait <- colnames(df)[!colnames(df) %in% vars_exclude]

# 4.2. Calcul des Médianes et IQR (25ème et 75ème percentiles)
summary_table <- df %>%
  select(ClusterID, any_of("AgeMonths.y"), all_of(vars_gait)) %>%
  group_by(ClusterID) %>%
  summarise(across(everything(), 
                   list(med = ~median(.x, na.rm = TRUE), 
                        q25 = ~quantile(.x, 0.25, na.rm = TRUE),
                        q75 = ~quantile(.x, 0.75, na.rm = TRUE)))) %>%
  pivot_longer(cols = -ClusterID, 
               names_to = c("Variable", ".value"), 
               names_sep = "_(?=[^_]+$)")

# 4.3. Nettoyage et Formatage Médiane [IQR]
summary_formatted <- summary_table %>%
  mutate(
    CleanVar = Variable,
    CleanVar = gsub("^Mean_", "", CleanVar),
    CleanVar = if_else(grepl("Norm Gait Speed", CleanVar, ignore.case = TRUE), "Norm Gait Speed (ua)", CleanVar),
    CleanVar = if_else(grepl("^CV_", CleanVar), paste0(gsub("^CV_", "C.V. ", gsub(" \\(.*\\)", "", CleanVar)), " (%)"), CleanVar),
    CleanVar = if_else(grepl("^SI_", CleanVar), gsub("^SI_", "S.I. ", gsub(" \\(.*\\)", "", CleanVar)), CleanVar),
    # Formatage de la cellule
    Med_IQR = paste0(round(med, 2), " [", round(q25, 2), " - ", round(q75, 2), "]")
  ) %>%
  select(ClusterID, CleanVar, Med_IQR) %>%
  pivot_wider(names_from = ClusterID, values_from = Med_IQR, names_prefix = "Cluster_") %>%
  rename(Variable = CleanVar)

# 4.4 & 4.5. Tests de Wilcoxon avec Sécurité Anti-Erreur
stats_tests_clean <- df %>%
  select(ClusterID, any_of("AgeMonths.y"), all_of(vars_gait)) %>%
  pivot_longer(cols = -ClusterID, names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  wilcox_test(Value ~ ClusterID)

# SÉCURITÉ : AJOUT DES ÉTOILES SI ELLES MANQUENT ---
# Si le test n'a pas généré p.signif, on le crée manuellement
if (!"p.signif" %in% colnames(stats_tests_clean)) {
  stats_tests_clean <- stats_tests_clean %>% add_significance("p")
}

# NETTOYAGE DES NOMS (IDENTIQUE AU TABLEAU SUMMARY) ---
stats_tests_clean <- stats_tests_clean %>%
  mutate(
    Variable = gsub("^Mean_", "", Variable),
    Variable = if_else(grepl("Norm Gait Speed", Variable, ignore.case = TRUE), "Norm Gait Speed (ua)", Variable),
    Variable = if_else(grepl("^CV_", Variable), paste0(gsub("^CV_", "C.V. ", gsub(" \\(.*\\)", "", Variable)), " (%)"), Variable),
    Variable = if_else(grepl("^SI_", Variable), gsub("^SI_", "S.I. ", gsub(" \\(.*\\)", "", Variable)), Variable),
    Variable = dplyr::recode(Variable, "AgeMonths.y" = "Age (months)")
  )

# FUSION SÉCURISÉE ---
# On vérifie quelles colonnes existent réellement pour ne pas faire planter select()
cols_disponibles <- intersect(c("Variable", "p", "p.signif"), colnames(stats_tests_clean))

final_stats_table <- left_join(summary_formatted, 
                               stats_tests_clean %>% select(all_of(cols_disponibles)), 
                               by = "Variable")

# 4.6. EXPORT ET MISE EN FORME GRAPHIQUE

# On retire les variables anthropométriques et on formate les p-values
table_to_export <- final_stats_table %>%
  filter(!Variable %in% c("Height_cm", "Weight_kg", "L0_m", "IMC")) %>%
  mutate(p = format.pval(p, digits = 2, eps = 0.001))

# Thème du tableau
tt <- ttheme_default(
  core = list(fg_params=list(cex = 0.7), 
              bg_params=list(fill=c("grey95", "white"), col=NA)),
  colhead = list(fg_params=list(cex = 0.8, fontface="bold"), 
                 bg_params=list(fill="grey80"))
)

# Création du tableau graphique
g_table <- tableGrob(table_to_export, rows = NULL, theme = tt)

# 5) Ajout des annotations en bas du tableau
notes <- textGrob("Values are Median [Interquartile Range]. P-values calculated using Wilcoxon rank-sum test.\nC.V. = Coefficient of Variation (%); S.I. = Symmetry Index; (ua) = arbitrary units.", 
                  x = 0, hjust = 0, vjust = 1, 
                  gp = gpar(fontsize = 8, fontitalic = TRUE))

# Assemblage : On ajoute une ligne de 2 cm pour la note
final_plot <- gtable_add_rows(g_table, heights = unit(2, "cm"))
final_plot <- gtable_add_grob(final_plot, notes, t = nrow(final_plot), l = 1, r = ncol(final_plot))

# --- CALCUL AUTOMATIQUE DE LA HAUTEUR ---
# On compte le nombre de lignes + entête + la note (environ 0.8 cm par ligne de donnée)
hauteur_calculee <- (nrow(table_to_export) * 0.8) + 4 

# EXPORTATION AVEC GGSAVE (Remplace png() et pdf())
ggsave("Table_Stats_Descriptives.png", 
       plot = final_plot, 
       width = 22, 
       height = hauteur_calculee, 
       units = "cm", 
       dpi = 300, 
       bg = "white")

ggsave("Table_Stats_Descriptives.pdf", 
       plot = final_plot, 
       width = 22, 
       height = hauteur_calculee, 
       units = "cm", 
       bg = "white")

print(paste("Le tableau a été exporté. Hauteur utilisée :", round(hauteur_calculee, 1), "cm"))




# 5.TEST CHI² ---
tab_sexe <- table(df$Sex, df$ClusterID)
chi2_sexe <- chisq_test(tab_sexe)

tab_age <- table(df$AgeGroup, df$ClusterID)
chi2_age <- chisq_test(tab_age)

tab_surface <- table(df$Condition, df$ClusterID)
chi2_surface <- chisq_test(tab_surface)

# 5.1. Vérifier les chiffres réels pour l'Âge
table(df$AgeGroup, df$ClusterID)

# 5.2. Vérifier les chiffres réels pour le Sexe
table(df$Sex, df$ClusterID)

# 5.3. Vérifier les chiffres réels pour la surface de marche
table(df$Condition, df$ClusterID)

print("--- RÉSULTATS CHI² ---")
print(chi2_age)
print(chi2_sexe)
print(chi2_surface)

# 5.4 Exporter les données
# Regrouper les résultats des tests Chi² dans un seul dataframe
# On extrait le nom de la variable, la statistique, le ddl (n) et la p-value
chi2_results <- bind_rows(
  mutate(chi2_age, Variable = "Age Group"),
  mutate(chi2_sexe, Variable = "Sex"),
  mutate(chi2_surface, Variable = "Surface")
) %>%
  select(Variable, n, statistic, df, p) %>%
  # Ajouter les étoiles de significativité
  add_significance() %>%
  # Formater la p-value pour l'affichage
  mutate(p = format.pval(p, digits = 2, eps = 0.001))

# Création de l'objet graphique (Tableau)
# Utilisation du même thème que pour les stats descriptives pour la cohérence
tt_chi2 <- ttheme_default(
  core = list(fg_params=list(cex = 0.9)),
  colhead = list(fg_params=list(cex = 1, fontface="bold"), 
                 bg_params=list(fill="grey80"))
)

g_chi2 <- tableGrob(chi2_results, rows = NULL, theme = tt_chi2)

# EXPORTATION
# Enregistrement en PNG
png("Table_Chi2_Results.png", width = 15, height = 8, units = "cm", res = 300)
grid.draw(g_chi2)
dev.off()

# Enregistrement en PDF
pdf("Table_Chi2_Results.pdf", width = 7, height = 4)
grid.draw(g_chi2)
dev.off()

print("Le tableau des résultats Chi² a été exporté.")




# 6. Observer les migrations des participants entre les 2 clusters
# 6.1. On prépare un tableau large pour comparer Plat vs High
migration_df <- df %>%
  filter(Condition %in% c("Plat", "High")) %>%
  select(Participant, AgeGroup, Condition, ClusterID) %>%
  pivot_wider(names_from = Condition, values_from = ClusterID) %>%
  # On crée une colonne qui indique s'il y a eu migration
  mutate(Migration = if_else(Plat == High, "Stable", "Migrant")) %>%
  filter(!is.na(Migration)) # On enlève ceux qui n'ont pas fait les deux tests

# 6.2. Qui sont les migrants ?
les_migrants <- migration_df %>% filter(Migration == "Migrant")

print(paste("Nombre de migrants :", nrow(les_migrants)))
View(les_migrants)

# 6.3. Tableau récapitulatif par groupe d'âge
recap_migration <- migration_df %>%
  group_by(AgeGroup) %>%
  summarise(
    Total = n(),
    Nb_Migrants = sum(Migration == "Migrant"),
    Perc_Migrants = round(100 * Nb_Migrants / Total, 1)
  ) %>%
  arrange(desc(Perc_Migrants))

print(recap_migration)




# ___________________________________________________________________________
# II - VISUALISATION
# ___________________________________________________________________________


# --- 1. PRÉPARATION DES LABELS EN ANGLAIS ---
df <- df %>%
  mutate(
    # On renomme et on ordonne en une seule étape
    AgeGroup_EN = factor(AgeGroup, 
                         levels = c("JeunesEnfants", "Enfants", "Adolescents", "Adultes"),
                         labels = c("Young Children", "Children", "Adolescents", "Adults")),
    
    Sex_EN = factor(Sex, 
                    levels = c("F", "M"), 
                    labels = c("Female", "Male")),
    
    Condition_EN = factor(Condition, 
                          levels = c("Plat", "Medium", "High"),
                          labels = c("Even", "Medium", "High"))
  )


# --- 2. BOXPLOT DES ÂGES ---
ggplot(df, aes(x = as.factor(ClusterID), y = AgeMonths.y, fill = as.factor(ClusterID))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) + 
  scale_fill_manual(values = c("#3498db", "#e74c3c")) +
  labs(title = "Age distribution by cluster",
       x = "Cluster ID", y = "Age (months)", fill = "Cluster") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))


# --- 3. RÉPARTITIONS (PROPORTIONS) ---
# A. Age group
p1 <- ggplot(df, aes(x = AgeGroup_EN, fill = as.factor(ClusterID))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#3498db", "#e74c3c")) +
  labs(title = "Age group", x = "", y = "% of sample", fill = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1))

# B. Sex
p2 <- ggplot(df, aes(x = Sex_EN, fill = as.factor(ClusterID))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#3498db", "#e74c3c")) +
  labs(title = "Sex", x = "", y = "", fill = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

# C. Surface (Condition)
p3 <- ggplot(df, aes(x = Condition_EN, fill = as.factor(ClusterID))) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c("#3498db", "#e74c3c")) +
  labs(title = "Surface", x = "", y = "", fill = "Cluster") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))


# --- 4. AFFICHAGE COMBINÉ ---
library(patchwork)
(p1 + p2 + p3) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Distribution of clusters across Age groups, Sex and Surfaces",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )


# --- 5. MIGRATIONS DE PARTICIPANTS ENTRE CLUSTERS ---
# 1. Préparation des labels de légende avec les % calculés en Section II
# On utilise votre tableau 'recap_migration' pour extraire les valeurs exactes
labels_avec_perc <- c(
  "Young Children" = paste0("Young Children (", recap_migration$Perc_Migrants[recap_migration$AgeGroup == "JeunesEnfants"], "%)"),
  "Children"       = paste0("Children (", recap_migration$Perc_Migrants[recap_migration$AgeGroup == "Enfants"], "%)"),
  "Adolescents"    = paste0("Adolescents (", recap_migration$Perc_Migrants[recap_migration$AgeGroup == "Adolescents"], "%)"),
  "Adults"         = paste0("Adults (", recap_migration$Perc_Migrants[recap_migration$AgeGroup == "Adultes"], "%)")
)

# 2. Préparation des données pour le flux
df_flux <- migration_df %>%
  mutate(
    AgeGroup_EN = factor(AgeGroup, 
                         levels = c("JeunesEnfants", "Enfants", "Adolescents", "Adultes"),
                         labels = c("Young Children", "Children", "Adolescents", "Adults")),
    # On force l'ordre des clusters pour la clarté visuelle (1 en bas, 2 en haut)
    Plat = factor(Plat, levels = c(2, 1)),
    High = factor(High, levels = c(2, 1))
  ) %>%
  group_by(AgeGroup_EN, Plat, High) %>%
  summarise(n = n(), .groups = 'drop')

# 3. Création du graphique Alluvial
plot_migration <- ggplot(df_flux, aes(y = n, axis1 = Plat, axis2 = High)) +
  # Rubans de flux : 'discern = TRUE' aide à supprimer les messages d'avis
  geom_alluvium(aes(fill = AgeGroup_EN), width = 1/12, alpha = 0.5, discern = TRUE) +
  
  # Blocs de clusters (Strates)
  geom_stratum(width = 1/6, fill = "white", color = "black", discern = TRUE) +
  
  # Étiquettes des clusters à l'intérieur des blocs
  geom_text(stat = "stratum", aes(label = paste0("Cluster ", after_stat(stratum))), 
            size = 3.5, fontface = "bold") +
  
  # Configuration des couleurs et de la légende
  scale_fill_manual(
    values = c(
      "Young Children" = "#3498db", # Bleu
      "Children"       = "#e67e22", # Orange
      "Adolescents"    = "#27ae60", # Vert
      "Adults"         = "#8e44ad"  # Violet
    ),
    labels = labels_avec_perc
  ) +
  
  # Configuration des axes et titres
  scale_x_discrete(limits = c("Even", "High"), expand = c(.1, .1)) +
  labs(
    title = "Individual Migration from Even to High Surface",
    subtitle = "Flows represent participants shifting gait profiles (% = migration rate)",
    fill = "Age Group (Migration rate %)",
    y = "Number of participants"
  ) +
  
  # Thème et mise en page
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank(),
    legend.position = "bottom" 
  )

# Affichage et sauvegarde
print(plot_migration)


# --- 6. EXPORT DES FIGURES ---
# A. Enregistrement du Boxplot
# On le redessine rapidement pour être sûr qu'il soit le "dernier" affiché
plot_age <- ggplot(df, aes(x = as.factor(ClusterID), y = AgeMonths.y, fill = as.factor(ClusterID))) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.5) + 
  scale_fill_manual(values = c("#3498db", "#e74c3c")) +
  labs(title = "Age distribution by cluster",
       x = "Cluster ID", y = "Age (months)", fill = "Cluster") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave("Boxplot_Age_Clusters.png", plot = plot_age, width = 15, height = 12, units = "cm", dpi = 300)
ggsave("Boxplot_Age_Clusters.pdf", 
       plot = plot_age, 
       width = 15, 
       height = 12, 
       units = "cm")


# B. Enregistrement du combiné (Proportions)
plot_combined <- (p1 + p2 + p3) + 
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Distribution of clusters across Age groups, Sex and Surfaces",
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  )

ggsave("Combined_Repartition_Clusters.png", plot = plot_combined, width = 25, height = 12, units = "cm", dpi = 300)
ggsave("Combined_Repartition_Clusters.pdf", 
       plot = plot_combined, 
       width = 25, 
       height = 12, 
       units = "cm")

# B. Enregistrement du graph migration
ggsave("Migration_Alluvial_Plot.png", plot = plot_migration, width = 20, height = 15, units = "cm", dpi = 300)
ggsave("Migration_Alluvial_Final.pdf", 
       plot = plot_migration, 
       width = 22, 
       height = 16, 
       units = "cm") # Pas besoin de DPI pour le PDF car c'est du vectoriel

print("Les images ont été enregistrées dans votre dossier de travail.")




# ___________________________________________________________________________
# III - RADAR PLOT DES CLUSTERS (TOUTES SURFACES CONFONDUES)
# ___________________________________________________________________________

# --- 1. PRÉPARATION DES VARIABLES ET DOMAINES ---

# Variables présentes dans votre CSV (après standardisation des noms)
# Note : on utilise les noms EXACTS tels qu'ils apparaissent dans votre CSV
vars_radar <- c(
  "Mean_Norm Gait Speed (m.s^{-1})",
  "Mean_Norm Step length (ua)",
  "Mean_Norm WR (ua)",
  "Mean_Double support time (%)",
  "Mean_Norm Cadence (ua)",
  "Mean_COM SPARC Magnitude (ua)",
  "Mean_Norm StepWidth (ua)",
  "Mean_MoS ML HS (%L0)",
  "Mean_MoS AP HS (%L0)",
  "SI_Stride length (m)",
  "SI_Double support time (%)",
  "SI_StepWidth (cm)",
  "CV_Norm StepWidth (ua)",
  "Mean_GVI (ua)",
  "CV_Gait speed (m.s^{-1})"
)

# Vérification que toutes les variables sont présentes
vars_radar_present <- intersect(vars_radar, names(df))
if(length(vars_radar_present) < length(vars_radar)) {
  message("Variables manquantes : ", paste(setdiff(vars_radar, vars_radar_present), collapse = ", "))
}

# Labels lisibles pour l'affichage
radar_labels <- c(
  "Norm Gait Speed (ua)",
  "Norm Step length (ua)",
  "Norm WR (ua)",
  "Double support time (%)",
  "Norm Cadence (ua)",
  "COM SPARC (ua)",
  "Norm StepWidth (ua)",
  "MoS ML HS (%L0)",
  "MoS AP HS (%L0)",
  "S.I. Stride length",
  "S.I. Double support",
  "S.I. StepWidth",
  "C.V. StepWidth (%)",
  "GVI (ua)",
  "C.V. Gait speed (%)"
)

# --- 2. DÉFINITION DES DOMAINES ---
domains_vars <- list(
  PACE = c(
    "Mean_Norm Gait Speed (m.s^{-1})",
    "Mean_Norm Step length (ua)",
    "Mean_Norm WR (ua)"
  ),
  RHYTHM = c(
    "Mean_Double support time (%)",
    "Mean_Norm Cadence (ua)",
    "Mean_COM SPARC Magnitude (ua)"
  ),
  `POSTURAL CONTROL` = c(
    "Mean_Norm StepWidth (ua)",
    "Mean_MoS ML HS (%L0)",
    "Mean_MoS AP HS (%L0)"
  ),
  ASYMMETRY = c(
    "SI_Stride length (m)",
    "SI_Double support time (%)",
    "SI_StepWidth (cm)"
  ),
  VARIABILITY = c(
    "CV_Norm StepWidth (ua)",
    "Mean_GVI (ua)",
    "CV_Gait speed (m.s^{-1})"
  )
)

# Couleurs des domaines
domain_colors <- c(
  PACE = "lightblue",
  RHYTHM = "lightcoral",
  `POSTURAL CONTROL` = "palegreen",
  ASYMMETRY = "plum",
  VARIABILITY = "lightyellow"
)

# Couleurs des clusters
cluster_colors <- c("1" = "blue", "2" = "red")  # Bleu et Rouge

# --- 3. FONCTION UTILITAIRE : ASSOCIER CHAQUE VARIABLE À SON DOMAINE ---
get_domain_for_vars <- function(vars_in_radar, domains_list) {
  dom_vec <- rep(NA_character_, length(vars_in_radar))
  names(dom_vec) <- vars_in_radar
  for (d in names(domains_list)) {
    dom_vec[vars_in_radar %in% domains_list[[d]]] <- d
  }
  dom_vec
}

# --- 4. FONCTION : DESSINER LES FONDS COLORÉS PAR DOMAINE ---
draw_domain_background <- function(domains_by_var, domain_cols, alpha = 0.18, r = 1) {
  n <- length(domains_by_var)
  if (n < 3) return(invisible(NULL))
  
  # Angles des axes (radar classique: premier en haut)
  angles <- seq(0, 2*pi, length.out = n + 1)[1:n] + (pi/2)
  
  # Limites entre axes = milieux angulaires
  bounds <- angles - (pi / n)
  bounds <- c(bounds, bounds[1] + 2*pi)
  
  # Regrouper les variables contiguës d'un même domaine
  runs <- rle(domains_by_var)
  idx_end <- cumsum(runs$lengths)
  idx_start <- c(1, head(idx_end, -1) + 1)
  
  for (k in seq_along(runs$values)) {
    dom <- runs$values[k]
    if (is.na(dom)) next
    col <- domain_cols[[dom]]
    if (is.null(col) || is.na(col)) next
    
    i1 <- idx_start[k]
    i2 <- idx_end[k]
    
    # Bornes angulaires du bloc contigu
    a_start <- bounds[i1]
    a_end   <- bounds[i2 + 1]
    
    # Points du secteur
    aa <- seq(a_start, a_end, length.out = 80)
    x <- c(0, r * cos(aa), 0)
    y <- c(0, r * sin(aa), 0)
    
    polygon(
      x, y,
      col = grDevices::adjustcolor(col, alpha.f = alpha),
      border = NA
    )
  }
  
  invisible(NULL)
}

# --- 5. CALCUL DES BORNES MIN/MAX GLOBALES (NORMALISATION) ---
radar_min_max <- df %>%
  dplyr::select(dplyr::all_of(vars_radar_present)) %>%
  dplyr::summarise(dplyr::across(
    dplyr::everything(),
    list(
      min = ~min(.x, na.rm = TRUE),
      max = ~max(.x, na.rm = TRUE)
    )
  ))

mins_raw <- radar_min_max %>% dplyr::select(dplyr::ends_with("_min")) %>% unlist() %>% as.numeric()
maxs_raw <- radar_min_max %>% dplyr::select(dplyr::ends_with("_max")) %>% unlist() %>% as.numeric()

# --- 6. FONCTION PRINCIPALE : CRÉER LE RADAR PLOT ---
create_cluster_radar <- function(df_full, vars, labels) {
  
  # A) Calcul des MÉDIANES par cluster (toutes surfaces confondues)
  data_median <- df_full %>%
    dplyr::group_by(ClusterID) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(vars), ~median(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::arrange(ClusterID)
  
  # B) Extraction des données individuelles (toutes surfaces confondues)
  data_indiv <- df_full %>%
    dplyr::select(ClusterID, dplyr::all_of(vars))
  
  # C) Récupération des min/max globaux (pour normalisation)
  mins <- as.vector(radar_min_max[grep("_min$", names(radar_min_max))])
  maxs <- as.vector(radar_min_max[grep("_max$", names(radar_min_max))])
  
  # D) Normalisation des médianes (0–1)
  radar_df_median <- as.data.frame(data_median[, -1, drop = FALSE])
  
  normalized_median <- as.data.frame(lapply(seq_len(ncol(radar_df_median)), function(i) {
    denom <- (maxs[[i]] - mins[[i]])
    if (is.na(denom) || denom == 0) return(rep(0, nrow(radar_df_median)))
    (radar_df_median[, i] - mins[[i]]) / denom
  }))
  
  colnames(normalized_median) <- labels
  
  # Format attendu par fmsb::radarchart :
  # - 1ère ligne : max (=1)
  # - 2ème ligne : min (=0)
  # - lignes suivantes : données (ici = médianes par cluster)
  final_radar_median <- rbind(rep(1, length(vars)), rep(0, length(vars)), normalized_median)
  
  # E) Normalisation des individus (0–1)
  radar_df_indiv <- data_indiv[, -1, drop = FALSE]
  
  normalized_indiv <- as.data.frame(lapply(seq_len(ncol(radar_df_indiv)), function(i) {
    denom <- (maxs[[i]] - mins[[i]])
    if (is.na(denom) || denom == 0) return(rep(0, nrow(radar_df_indiv)))
    (radar_df_indiv[, i] - mins[[i]]) / denom
  }))
  
  colnames(normalized_indiv) <- labels
  
  # F) Couleurs : clusters
  colors_border_median <- cluster_colors[as.character(data_median$ClusterID)]
  colors_in_median <- grDevices::adjustcolor(colors_border_median, alpha.f = 0.20)
  
  # G) Tracé en 3 couches (du fond vers l'avant)
  # 1) Cadre (axes, grille, titre) avec polygones invisibles
  # 2) Fonds colorés par domaine
  # 3) Individus (gris, derrière)
  # 4) Médianes de cluster (couleur, devant)
  
  # =========================
  # G1) Cadre du radar
  # =========================
  
  ng <- nrow(data_median)
  transparent <- grDevices::adjustcolor("white", alpha.f = 0)
  
  fmsb::radarchart(
    final_radar_median,
    axistype = 0,
    seg = 4,
    pcol  = rep(transparent, ng),
    pfcol = rep(transparent, ng),
    plwd  = rep(0.01, ng),
    plty  = rep(1, ng),
    cglcol = "grey70", 
    cglty = 1, 
    cglwd = 0.8,
    vlcex = 0.7,
    title = "Gait Profile Comparison Between Clusters (All Surfaces)"
  )
  
  # === Préparation des variables pour les étiquettes ===
  nvar <- length(vars)
  angles <- seq(0, 2*pi, length.out = nvar + 1)[1:nvar] + (pi/2)
  
  pct <- c(0.25, 0.50, 0.75, 1.00)
  r_levels <- pct
  
  ticks_real <- sapply(seq_len(nvar), function(i) {
    mins[[i]] + pct * (maxs[[i]] - mins[[i]])
  })
  
  # =========================
  # G2) Fond coloré par domaine
  # =========================
  domains_by_var <- get_domain_for_vars(vars, domains_vars)
  
  par(new = TRUE)  # superpose sur le même repère
  draw_domain_background(
    domains_by_var = domains_by_var,
    domain_cols    = domain_colors,
    alpha          = 0.35,
    r              = 1
  )
  
  # === Affichage des étiquettes de valeurs réelles ===
  for (i in seq_len(nvar)) {
    angle <- angles[i]
    
    for (j in seq_along(pct)) {
      r <- r_levels[j]
      
      # Position du texte (légèrement décalé vers l'extérieur)
      x_pos <- r * cos(angle) * 1.05
      y_pos <- r * sin(angle) * 1.05
      
      # Valeur réelle
      val <- round(ticks_real[j, i], 2)
      
      # Afficher le texte
      text(
        x = x_pos,
        y = y_pos,
        labels = val,
        cex = 0.5,
        col = "grey20",
        font = 1
      )
    }
  }
  
  # =========================
  # G3) Individus : tracés en gris derrière
  # =========================
  indiv_col <- grDevices::adjustcolor("grey30", alpha.f = 0.10)
  
  for (i in seq_len(nrow(normalized_indiv))) {
    par(new = TRUE)
    fmsb::radarchart(
      rbind(rep(1, length(vars)), rep(0, length(vars)), normalized_indiv[i, , drop = FALSE]),
      axistype = 0,
      vlabels = rep("", length(vars)),
      pcol = indiv_col,
      pfcol = NA,
      plwd = 0.7,
      plty = 1,
      cglcol = NA,
      axislabcol = NA,
      vlcex = 0,
      seg = length(vars)
    )
  }
  
  # =========================
  # G4) Médianes : tracés colorés au premier plan
  # =========================
  par(new = TRUE)
  fmsb::radarchart(
    final_radar_median,
    axistype = 0,
    vlabels = rep("", length(vars)),
    pcol = colors_border_median,
    pfcol = colors_in_median,
    plwd = 3.5,
    plty = 1,
    cglcol = NA,
    axislabcol = NA,
    vlcex = 0,
    seg = length(vars)
  )
  
  # --- Légende
  legend(
    x = "bottom",
    legend = c("Cluster 1", "Cluster 2"),
    inset = -0.15,
    horiz = TRUE,
    bty = "n",
    pch = 20,
    col = cluster_colors,
    text.col = "black",
    cex = 0.8,
    pt.cex = 1.5,
    xpd = TRUE
  )
}

# --- 7. GÉNÉRATION DES RADAR PLOTS ---

# 7.1) Export PDF
pdf("Radar_Plot_Clusters.pdf", width = 12, height = 12)
par(mfrow = c(1, 1))

create_cluster_radar(df, vars_radar_present, radar_labels)

dev.off()

# 7.2) Export PNG haute qualité
png(filename = "Radar_Plot_Clusters.png", width = 5200, height = 5200, res = 600, type = "cairo")

op <- par(no.readonly = TRUE)

par(mfrow = c(1, 1))
par(mar = c(8, 5, 5, 5))     # bottom, left, top, right
par(oma = c(0, 0, 0, 0))
par(xaxs = "i", yaxs = "i")
par(xpd = NA)
par(cex = 0.85)

create_cluster_radar(df, vars_radar_present, radar_labels)

par(op)
dev.off()

message("Radar plot des clusters généré avec succès !")
message("Fichiers créés : Radar_Plot_Clusters.pdf et Radar_Plot_Clusters.png")
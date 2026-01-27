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
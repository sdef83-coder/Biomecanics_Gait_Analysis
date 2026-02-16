setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/PhysicalActi_et_VAS")

# =========================================================
# 0) LIBRAIRIES
# - dplyr / tidyr : manipulation de données (jointures, filtres, reshape long/large)
# - ggplot2       : visualisations (boxplots, points individuels, export figures)
# - lme4/lmerTest  : modèles linéaires mixtes + tests (Satterthwaite)
# - emmeans       : moyennes marginales estimées + comparaisons post-hoc
# - car           : fonctions utiles pour ANOVA / analyses (selon besoins)
# - performance   : indices de performance des modèles (R² de Nakagawa)
# - gtsummary/gt  : production de tableaux (mise en forme + export)
# =========================================================
library(dplyr)
library(tidyr)
library(ggplot2)
library(lme4)
library(lmerTest)
library(emmeans)
library(car)
library(performance)
library(gtsummary)
library(gt)

# =========================================================
# 1) IMPORT DES DONNÉES
# - path : chemin de travail (utile si besoin de relocaliser les exports)
# - read.csv(file.choose()) : sélection manuelle du fichier à importer (csv séparateur ;)
# - View / names : inspection rapide de la structure (colonnes, cohérence)
# =========================================================
path <- ("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/PhysicalActi_et_VAS")

df <- read.csv(file.choose(), sep = ";")

View(df)
names(df)


# =========================================================
# ------------ DECRIRE LA POPULATION ----------------------------------------------
# Objectif : décrire l’échantillon (répartition des groupes d’âge + sexes)
# =========================================================

# =========================================================
# 1) RÉPARTITION DES GROUPES D'ÂGE DANS L'ÉCHANTILLON
# - factor() : conversion en facteur ordonné pour garantir l’ordre logique (Children < Adolescent < Adult)
# - table()  : effectifs par catégorie
# - Age_prop : proportions (%) des catégories (optionnel)
# - barplot  : visualisation simple des effectifs par groupe d’âge
# =========================================================

# Transformer la colonne d'âge en facteur (pour pouvoir les manipuler) et les ordonner
df$AgeGroup <- factor(df$AgeGroup,levels = c("Children", "Adolescent", "Adult"),ordered = TRUE)

# voir la répartition de chaque âge dans mon échantillon
AGE <- table(df$AgeGroup)
AGE

# Nommer mes groupe d'âge (utile pour appeler si besoin)
age_group_name <- c("Children", "Adolescent", "Adult")

Age_prop <- table(df$AgeGroup)/length(df$AgeGroup) #Optionnel, voir la répartition en % de chaque catégorie
Age_prop

# Visualiser
barplot(AGE, col = "orange", main = "Age repartition across the included population", xlab ="Age Group", ylab = "Sample")



# =========================================================
# 2) RÉPARTITION DES SEXES DANS CHAQUE GROUPE D'ÂGE
# - table(Sex, AgeGroup) : tableau croisé effectifs Sexe x Âge
# - prop.table(..., margin = 2) : proportions par colonne (donc par groupe d’âge)
# - barplot(..., beside=TRUE) : barres côte à côte pour comparer les sexes à l’intérieur de chaque âge
# =========================================================

AGE_SEXE <- table(df$Sex,df$AgeGroup)
AGE_SEXE
PROP <- prop.table(AGE_SEXE, margin = 2)  # margin = 2 → par colonne (âge)
PROP
barplot(PROP, col = c("pink","blue"),beside = TRUE, legend.text = TRUE, main = "Sex repartition across age groups", xlab = "Age Group", ylab = "Sample")



# =========================================================
# ------------ ETUDIER LE VAS ----------------------------------------------
# Objectif : (i) décrire le VAS par âge et surface, (ii) visualiser, (iii) tester statistiquement via LMM,
#           (iv) exporter résultats (figures + CSV), (v) produire un tableau descriptif format publication.
# =========================================================

# =========================================================
# 1) MOYENNES PAR GROUPE D'ÂGE POUR CHAQUE SURFACE (VAS)
# - aggregate() : calcule la moyenne (mean) par groupe d’âge, en ignorant les NA (na.rm=TRUE)
# - Sorties : 3 data.frames (Even, Medium, High) au format "AgeGroup + moyenne"
# =========================================================
V.A.S_Even <- aggregate(V.A.S.Even ~ AgeGroup, data = df, FUN = mean, na.rm = TRUE)
V.A.S_Even

V.A.S_Medium <- aggregate(V.A.S.Medium ~ AgeGroup, data = df, FUN = mean, na.rm = TRUE)
V.A.S_Medium

V.A.S_High <- aggregate(V.A.S.High ~ AgeGroup, data = df, FUN = mean, na.rm = TRUE)
V.A.S_High



# =========================================================
# 2) CONSTRUCTION D'UNE MATRICE RÉCAPITULATIVE (SURFACES x ÂGES)
# - rbind() : empile les vecteurs de moyennes en lignes (Even/Medium/High)
# - colnames() : nomme les colonnes selon l’ordre des groupes d’âge
# - Objectif : obtenir un format compact pour inspection rapide / export éventuel
# =========================================================

# Création de la matrice avec les données VAS
VAS_matrix <- rbind(
  Even   = V.A.S_Even$V.A.S.Even,
  Medium = V.A.S_Medium$V.A.S.Medium,
  High   = V.A.S_High$V.A.S.High
)

# Nomme les colonnes après avoir coller toutes les lignes d'intérêt
colnames(VAS_matrix) <- age_group_name
VAS_matrix


# =========================================================
# 3) VISUALISATION : PASSAGE EN FORMAT LONG + BOXPLOTS + POINTS INDIVIDUELS
# - Format long : indispensable pour ggplot (et pour les modèles mixtes avec mesures répétées)
# - Participant : création d’un identifiant unique par ligne (proxy d’ID sujet si non présent)
# - pivot_longer() : transforme 3 colonnes VAS -> 1 colonne "VAS" + 1 colonne "Surface"
# - factor() : définition contrôlée de l’ordre et des labels (Even/Medium/High)
# - ggplot : boxplot (tendance centrale + dispersion) + jitter (valeurs individuelles)
# - ggsave : export de la figure en PNG et PDF dans le répertoire courant
# =========================================================

# i) S'assurer qu'on a une colonne d'identifiant unique par sujet
df <- df %>% mutate(Participant = row_number()) 

# ii) Passage au format long (1 ligne = 1 sujet x 1 condition)
df_long <- df %>%
  pivot_longer(
    cols = c(V.A.S.Even, V.A.S.Medium, V.A.S.High),
    names_to = "Surface",
    values_to = "VAS"
  ) %>%
  mutate(
    # Mise en facteurs pour le LMM et contrôle de l’ordre d’affichage des modalités
    Participant = factor(Participant), 
    Surface = factor(Surface,
                     levels = c("V.A.S.Even", "V.A.S.Medium", "V.A.S.High"),
                     labels = c("Even", "Medium", "High")),
    AgeGroup = factor(AgeGroup, levels = c("Children", "Adolescent", "Adult"))
  )

# (iii) Création des boxplots (Surface en x, VAS en y, couleur par AgeGroup)
ggplot(df_long, aes(x = Surface, y = VAS, fill = AgeGroup)) +
  
  # 1) Boxplots : médiane + IQR (outliers masqués ici pour éviter surcharge visuelle)
  geom_boxplot(
    position = position_dodge(0.85),
    width = 0.7,
    alpha = 0.7,
    outlier.shape = NA,
    color = "black"
  ) +
  
  # 2) Points individuels : jitter + dodge pour visualiser la dispersion des sujets
  geom_jitter(
    aes(group = AgeGroup),
    position = position_dodge(0.85),
    size = 1.2,
    alpha = 0.3,
    color = "black"
  ) +
  
  # 3) Palette manuelle par groupe d’âge (cohérence visuelle inter-figures)
  scale_fill_manual(
    values = c(
      "Children"   = "orange",
      "Adolescent" = "green",
      "Adult"      = "purple"
    )
  ) +
  
  # 4) Titres, axes et thème
  labs(
    title = "Visual Analogue Scale across surfaces and age groups",
    x = "Surface",
    y = "VAS (cm)",
    fill = "Age group"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  )

# Export figure PNG (qualité impression)
ggsave(
  filename = "VAS_Age_Surface.png", 
  width = 10,
  height = 7,
  dpi = 300
)

# Export figure PDF (vectoriel, utile pour publication)
ggsave(
  filename = "VAS_Figure.pdf",
  width = 10, 
  height = 7
)

message("Graphiques enregistrés en PNG et PDF dans : ", getwd())

# =========================================================
# 4) TESTS STATISTIQUES : MODÈLE LINÉAIRE MIXTE (LMM) SUR LE VAS
# - Modèle : VAS ~ Surface * AgeGroup + (1|Participant)
#   * Effets fixes : Surface, AgeGroup, et interaction Surface:AgeGroup
#   * Effet aléatoire : intercept par participant (mesures répétées par sujet)
# - ANOVA Type III : tests des effets principaux + interaction (df Satterthwaite via lmerTest)
# - Post-hocs : comparaisons Holm (contrôle du risque d’erreur de type I)
# =========================================================

# 2) Fit du Modèle Linéaire Mixte (LMM)
model_vas <- lmer(VAS ~ Surface * AgeGroup + (1|Participant), data = df_long)

# 3) ANOVA Type III (Satterthwaite)
anova_vas <- anova(model_vas, type = 3)
print("--- Résultats de l'ANOVA sur le VAS ---")
print(anova_vas)

# 4) Post-hocs : Effet GLOBAL de la Surface (tous âges confondus)
print("--- Post-hocs : Effet principal de la Surface (indépendant de l'âge) ---")
em_surface_global <- emmeans(model_vas, ~ Surface)
pairs(em_surface_global, adjust = "holm")

# 4) Post-hocs : Effet GLOBAL de l'âge (toutes surfaces confondus)
print("--- Post-hocs : Effet principal de l'AgeGroup (indépendant de la surface) ---")
em_age_global <- emmeans(model_vas, ~ AgeGroup)
pairs(em_age_global, adjust = "holm")

# 5) Post-hocs : Comparaisons par Surface au sein de chaque Groupe d'Âge
print("--- Post-hocs : Effet de la Surface par Groupe d'Âge ---")
em_surf <- emmeans(model_vas, ~ Surface | AgeGroup)
pairs(em_surf, adjust = "holm")

# 6) Post-hocs : Comparaisons par Groupe d'Âge pour chaque Surface
print("--- Post-hocs : Effet de l'Âge par Surface ---")
em_age <- emmeans(model_vas, ~ AgeGroup | Surface)
pairs(em_age, adjust = "holm")


# =========================================================
# 5) R² (NAKAGAWA) + EXPORT DES RÉSULTATS (CSV)
# - r2_nakagawa() : R² marginal (effets fixes) et conditionnel (fixes + aléatoires)
# - add_signif()  : ajoute une annotation de significativité (*, **, ***)
# - df_anova_csv  : ANOVA reformatée avec étiquettes lisibles (effets principaux + interaction)
# - df_ph_csv     : post-hocs concaténés dans un seul tableau (avec colonne "Analysis")
# - write.csv / write.table : export des résultats dans deux fichiers distincts
# =========================================================

# 1) Calcul des R-carrés (Nakagawa)
r2_vals <- r2_nakagawa(model_vas)

# 2) Fonction utilitaire : code significativité à partir des p-values
add_signif <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ".")))
}

# 2) Mise en forme ANOVA pour export : ajout d'une colonne "Effect" et des codes de significativité
df_anova_csv <- as.data.frame(anova_vas) %>%
  tibble::rownames_to_column("Effect") %>%
  mutate(
    Signif = add_signif(`Pr(>F)`),
    Effect = case_when(
      Effect == "Surface" ~ "Surface (Main Effect)",
      Effect == "AgeGroup" ~ "AgeGroup (Main Effect)",
      Effect == "Surface:AgeGroup" ~ "Interaction (Surface x AgeGroup)",
      TRUE ~ Effect
    )
  ) %>%
  relocate(Effect, Signif)

# 3) Agrégation des post-hocs dans un seul data.frame, avec provenance de l’analyse
df_ph_csv <- bind_rows(
  as.data.frame(pairs(em_surface_global, adjust = "holm")) %>% mutate(Analysis = "Global Surface"),
  as.data.frame(pairs(em_age_global, adjust = "holm"))    %>% mutate(Analysis = "Global AgeGroup"),
  as.data.frame(pairs(em_surf, adjust = "holm")) %>% mutate(Analysis = "Surface within AgeGroup"),
  as.data.frame(pairs(em_age, adjust = "holm"))           %>% mutate(Analysis = "AgeGroup within Surface")
) %>% 
  mutate(Signif = add_signif(p.value))

# 4) Export des résultats
# - Fichier A : ANOVA + R²
# - Fichier B : Post-hocs
write.csv(df_anova_csv, "VAS_ANOVA_Main_Results.csv", row.names = FALSE)

# Ajout des R² à la fin du CSV (append)
cat("\n--- Effect Size (Nakagawa) ---\n", file = "VAS_ANOVA_Main_Results.csv", append = TRUE)
write.table(data.frame(
  Metric = c("R2_Marginal", "R2_Conditional"),
  Value = c(r2_vals$R2_marginal, r2_vals$R2_conditional)
), "VAS_ANOVA_Main_Results.csv", append = TRUE, sep = ",", row.names = FALSE, col.names = TRUE)

write.csv(df_ph_csv, "VAS_PostHocs_Results.csv", row.names = FALSE)

message("Fichiers CSV enregistrés : 'VAS_ANOVA_Main_Results.csv' et 'VAS_PostHocs_Results.csv'")



# =========================================================
# ------------ TABLEAU DESCRIPTIF VAS (format avec sous-colonnes) ----------------------------------------------
# Objectif : produire un tableau descriptif "Mean ± SD" par AgeGroup et Surface
# - group_by/summarise : calcule mean, sd, n par combinaison
# - sprintf : mise en forme "moyenne ± sd" (2 décimales)
# - pivot_wider : passage au format large pour obtenir des colonnes AgeGroup_Surface
# - gt : création du tableau final avec en-têtes groupées (tab_spanner)
# - gtsave : export du tableau (PNG + PDF)
# =========================================================

# Calculer mean ± SD pour chaque combinaison AgeGroup x Surface
vas_summary <- df_long %>%
  group_by(AgeGroup, Surface) %>%
  summarise(
    mean_val = mean(VAS, na.rm = TRUE),
    sd_val = sd(VAS, na.rm = TRUE),
    n = n(),
    .groups = 'drop'
  ) %>%
  mutate(
    stat = sprintf("%.2f ± %.2f", mean_val, sd_val)
  ) %>%
  select(AgeGroup, Surface, stat)

# Pivoter pour avoir AgeGroup en colonnes et Surface en sous-colonnes
vas_wide <- vas_summary %>%
  pivot_wider(
    names_from = c(AgeGroup, Surface),
    values_from = stat,
    names_sep = "_"
  ) %>%
  mutate(Variable = "V.A.S (cm)") %>%
  select(Variable, everything())

# Créer le tableau gt avec spanning headers (groupes d’âges en en-têtes, surfaces en sous-colonnes)
gt_vas <- vas_wide %>%
  gt() %>%
  tab_header(
    title = "Descriptive values of Visual Analogue Scale according to age group and walking surface",
    subtitle = "Values are reported as Mean ± SD"
  ) %>%
  tab_spanner(
    label = "Children",
    columns = c(Children_Even, Children_Medium, Children_High)
  ) %>%
  tab_spanner(
    label = "Adolescent",
    columns = c(Adolescent_Even, Adolescent_Medium, Adolescent_High)
  ) %>%
  tab_spanner(
    label = "Adult",
    columns = c(Adult_Even, Adult_Medium, Adult_High)
  ) %>%
  cols_label(
    Variable = "",
    Children_Even = "Even",
    Children_Medium = "Medium",
    Children_High = "High",
    Adolescent_Even = "Even",
    Adolescent_Medium = "Medium",
    Adolescent_High = "High",
    Adult_Even = "Even",
    Adult_Medium = "Medium",
    Adult_High = "High"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_spanners()
  )

# Affichage
gt_vas

# Export (image)
gt::gtsave(gt_vas, "Table_VAS_Descriptive.png", vwidth = 1200, vheight = 400)

# Export (PDF)
gt::gtsave(gt_vas, "Table_VAS_Descriptive.pdf")

message("Tableau descriptif VAS créé et exporté avec succès !")

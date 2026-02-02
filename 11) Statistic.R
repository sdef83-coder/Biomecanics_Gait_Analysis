## Script en 3 parties 
# I. Stat.descriptives pop.
# II. Comparaison STP variables & tableau descriptif
# III. LMM & Post-hocs sur variables d'intérêt

## ============================================================
## I. STATISTIQUES DESCRIPTIVES POPULATION
## ============================================================
## =========================================================
## =========================================================
## Table 1 (descriptifs) à partir de participants_metadata.csv
## - Quantitatives : moyenne ± ET
## - Qualitatives  : n (%)
## - Colonnes : Global (Sexe uniquement) + par AgeGroup
## =========================================================

setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result")

library(readr)
library(dplyr)
library(stringr)
library(gtsummary)
library(gt)
library(janitor)
library(performance)

chemin <- "C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM"

# 0) Import + nettoyage des noms (snake_case)
df <- read.csv(
  file.path(chemin, "participants_metadata.csv"),
  sep = ";",
  check.names = FALSE,
  stringsAsFactors = FALSE
) %>%
  janitor::clean_names()

# 1) Typage + facteurs
df <- df %>%
  mutate(
    participant = as.character(participant),
    age_group = factor(
      age_group,
      levels = c("Jeunes Enfants","Enfants","Adolescents","Adultes"),
      labels = c("Young children", "Children", "Adolescents", "Adults")
    ),
    sex         = factor(sex, levels = c("F","M")),
    age_months  = as.numeric(age_months),
    height_cm   = as.numeric(height_cm),
    weight_kg   = as.numeric(weight_kg),
    l0_m        = as.numeric(l0_m),
    imc         = as.numeric(imc)
  )

# 2) Construction de la Table 1 en deux parties

# A. Partie Sexe (Qualitative) - SANS add_p()
tab_sex <- df %>%
  select(age_group, sex) %>%
  tbl_summary(
    by = age_group,
    label = list(sex ~ "Sex"),
    statistic = list(all_categorical() ~ "{n} ({p}%)"),
    digits = list(all_categorical() ~ c(0, 1)),
    missing = "no"
  ) %>%
  add_overall(last = FALSE, col_label = "**Global**")

# B. Partie Anthropométrique (Quantitative) + p-value
tab_others <- df %>%
  select(age_group, age_months, height_cm, weight_kg, l0_m, imc) %>%
  tbl_summary(
    by = age_group,
    label = list(
      age_months ~ "Age (months)",
      height_cm  ~ "Height (cm)",
      weight_kg  ~ "Weight (kg)",
      l0_m       ~ "L0 (m)",
      imc        ~ "BMI (kg/m²)"
    ),
    statistic = list(all_continuous() ~ "{mean} ± {sd}"),
    digits = list(all_continuous() ~ 2),
    missing = "no"
  ) %>%
  add_p(test = all_continuous() ~ "kruskal.test")

# C. Fusion et Nettoyage des indices
tab1_final <- tbl_stack(list(tab_sex, tab_others)) %>%
  bold_labels() %>%
  modify_header(
    label ~ "**Variable**",
    p.value ~ "**p-value**"
  ) %>%
  modify_spanning_header(all_stat_cols() ~ "**Age groups**") %>%
  modify_footnote(all_stat_cols() ~ "n (%), Mean ± SD")

# Affichage du résultat
tab1_final

# 3) Calculs de vérification
age_counts <- df %>% count(age_group)
sex_by_age <- df %>% count(age_group, sex)

print(age_counts)
print(sex_by_age)

# 4) Export HTML
as_gt(tab1_final) %>% gt::gtsave("Table1_participants.html")

# 5) Export pdf
# Convertir le tableau gtsummary en objet gt
gt_table <- as_gt(tab1_final)
# Sauvegarde en PDF (Qualité vectorielle parfaite)
gt::gtsave(gt_table, "Table1_participants.pdf")

print("Parti I effectuée avec succés ! On a la Table 1")




## ============================================================
## II. MOYENNE +- SD DES VARIABLES SPT
## ============================================================

# ---------------------------------------------------------
# 0) Préparation de l'environnement
# ---------------------------------------------------------
setwd('C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data')

library(tidyverse)
library(readr)
library(gt)

# ---------------------------------------------------------
# 1) Chargement du fichier de données
# ---------------------------------------------------------
file_path <- "C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data/ACP_Clustering_DATA.csv"   # adapter si besoin

# Lecture robuste (gère séparateur ; ou ,)
first_line <- readLines(file_path, n = 1, warn = FALSE)
delim <- ifelse(grepl(";", first_line), ";", ",")
df <- read_delim(file_path, delim = delim, show_col_types = FALSE)

# ---------------------------------------------------------
# 2) Standardisation des noms de colonnes (pour matcher variables_interet)
#    Objectif : uniformiser les unités et remplacer espaces/ponctuation par "_"
# ---------------------------------------------------------
standardize_names <- function(x) {
  x %>%
    # 1) trim + remplacer espaces / ponctuation par _
    str_trim() %>%
    str_replace_all("[[:space:]]+", "_") %>%
    str_replace_all("[\\-]+", "_") %>%
    str_replace_all("[,;:]+", "_") %>%
    
    # 2) unités / parenthèses -> suffixes normalisés
    str_replace_all("\\(mm\\)", "mm") %>%
    str_replace_all("\\(%L0\\)", "pL0") %>%
    str_replace_all("\\(ua\\)", "ua") %>%
    str_replace_all("\\(%\\)", "p") %>%
    str_replace_all("\\(m\\.s\\^\\{-1\\}\\)", "ms1") %>%
    str_replace_all("\\(step\\.min\\^\\{-1\\}\\)", "stepmin1") %>%
    
    # 3) enlever parenthèses restantes
    str_replace_all("[\\(\\)]", "") %>%
    
    # 4) cleanup underscores
    str_replace_all("__+", "_") %>%
    str_replace_all("_$", "")
}

names(df) <- standardize_names(names(df))

# Optionnel : vérifier les noms standardisés
# print(names(df))

# ---------------------------------------------------------
# 3) Définition des variables d’intérêt (dans l’ordre souhaité)
# ---------------------------------------------------------
variables_interet <- c(
  "NCycles_Left", "NCycles_Right",
  
  "Mean_Gait_speed_m.s^{_1}", "Mean_Norm_Gait_Speed_m.s^{_1}", "Mean_Step_length_m", "Mean_Stride_length_m",  "Mean_Norm_Step_length_ua", "Mean_WalkRatio", "Mean_Norm_WR_ua",
  
  "Mean_Double_support_time_p", "Mean_Cadence_step.min^{_1}", "Mean_Norm_Cadence_ua", "Mean_COM_SPARC_Magnitude_ua", "Mean_StepTime_s", "Mean_StanceTime_s", "Mean_SwingTime_s",
  
  "Mean_StepWidth_cm", "Mean_Norm_StepWidth_ua", "Mean_MoS_AP_HS_mm", "Mean_MoS_ML_HS_mm", "Mean_MoS_AP_Stance_mm", "Mean_MoS_ML_Stance_mm", "Mean_MoS_AP_HS_pL0", "Mean_MoS_ML_HS_pL0", "Mean_MoS_AP_Stance_pL0", "Mean_MoS_ML_Stance_pL0",
  
  "Mean_GVI_ua", "CV_Norm_StepWidth_ua", "CV_Gait_speed_m.s^{_1}",
  
  "SI_Stride_length_m", "SI_Double_support_time_p", "SI_Norm_StepWidth_ua"
)

# ---------------------------------------------------------
# 4) Identification automatique des colonnes AgeGroup et Surface
# ---------------------------------------------------------
age_candidates <- c("AgeGroup", "Groupe", "Age_Group", "AgeGroup.x", "AgeGrp")
surf_candidates <- c("Surface", "Condition", "Surf")

age_col  <- intersect(age_candidates, names(df))[1]
surf_col <- intersect(surf_candidates, names(df))[1]

if (is.na(age_col) || is.na(surf_col)) {
  stop(
    "Impossible de trouver les colonnes AgeGroup / Surface.\n",
    "Colonnes dispo: ", paste(names(df), collapse = ", "), "\n",
    "Renomme tes colonnes ou modifie age_candidates / surf_candidates."
  )
}

# ---------------------------------------------------------
# 5) Vérification des variables disponibles (présentes vs absentes)
# ---------------------------------------------------------
vars_present <- intersect(variables_interet, names(df))
vars_missing <- setdiff(variables_interet, names(df))

if (length(vars_missing) > 0) {
  message("Variables absentes dans le fichier (ignorées) :\n- ", paste(vars_missing, collapse = "\n- "))
}
if (length(vars_present) == 0) stop("Aucune variable_interet trouvée dans le fichier.")

# ---------------------------------------------------------
# 6) Passage au format long + calcul des descriptifs (moyenne, SD, n)
# ---------------------------------------------------------
df_long <- df %>%
  mutate(
    AgeGroup = as.factor(.data[[age_col]]),
    Surface  = as.factor(.data[[surf_col]])
  ) %>%
  filter(Surface %in% c("Plat", "Medium", "High")) %>%
  select(AgeGroup, Surface, all_of(vars_present)) %>%
  pivot_longer(cols = all_of(vars_present), names_to = "Variable", values_to = "Value") %>%
  mutate(Value = suppressWarnings(as.numeric(Value)))

desc <- df_long %>%
  group_by(AgeGroup, Surface, Variable) %>%
  summarise(
    n    = sum(!is.na(Value)),
    mean = mean(Value, na.rm = TRUE),
    sd   = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(cell = ifelse(
    n == 0, "",
    sprintf("%.2f \u00B1 %.2f (n=%d)", mean, sd, n)
  )) %>%
  select(AgeGroup, Surface, Variable, cell)

# ---------------------------------------------------------
# 7) Mise en forme "large" : colonnes = (AgeGroup x Surface)
# ---------------------------------------------------------
tab_wide <- desc %>%
  mutate(col = paste0(as.character(AgeGroup), "___", as.character(Surface))) %>%
  select(Variable, col, cell) %>%
  pivot_wider(names_from = col, values_from = cell) %>%
  arrange(match(Variable, variables_interet))

# ---------------------------------------------------------
# 8) Fonction d'étiquetage "publication" (affichage seulement)
#    - retire "Mean_"
#    - applique unités (mm, cm, m, s, %, %L0, ua, m/s, step/min)
#    - CV -> "C.V ... (%)"
#    - SI -> "S.I ..." (sans unité)
#    - remplace "_" par espaces
# ---------------------------------------------------------
make_pretty_label <- function(varname) {
  
  v <- varname
  
  # 1) Retirer le préfixe Mean_ dans l'affichage
  v <- str_replace(v, "^Mean_", "")
  
  # SI / CV : format publication
  is_cv <- str_detect(v, "^CV_")
  is_si <- str_detect(v, "^SI_")
  
  if (is_cv) {
    # garder uniquement "CV_<nom_variable>" et jeter tout ce qui suit (unité, même si underscores)
    v <- str_replace(v, "^(CV_[^_]+_[^_]+).*$", "\\1")
    v <- str_replace(v, "^CV_", "")
    v <- paste0("C.V ", v, " (%)")
  }
  
  if (is_si) {
    # garder uniquement "SI_<nom_variable>" (sans unité)
    v <- str_replace(v, "^(SI_[^_]+_[^_]+).*$", "\\1")
    v <- str_replace(v, "^SI_", "")
    v <- paste0("S.I ", v)
  }
  
  # 2) Unités (suffixes) -> format publication
  v <- str_replace(v, "_mm$", " (mm)")
  v <- str_replace(v, "_cm$", " (cm)")
  v <- str_replace(v, "_m$",  " (m)")
  v <- str_replace(v, "_s$",  " (s)")
  v <- str_replace(v, "_pL0$", " (%L0)")
  v <- str_replace(v, "_p$",   " (%)")
  v <- str_replace(v, "_ua$",  " (ua)")
  
  # Vitesse / cadence : gérer variantes “ms1” / “m.s^{-1}” / “m.s^{_1}”
  v <- str_replace(v, "_ms1$", " (m/s)")
  v <- str_replace(v, "_m\\.s\\^\\{-1\\}$", " (m/s)")
  v <- str_replace(v, "_m\\.s\\^\\{_1\\}$", " (m/s)")
  v <- str_replace(v, "_m\\.s\\^\\{\\-1\\}$", " (m/s)")
  v <- str_replace(v, "_m\\.s\\^\\{\\-?1\\}$", " (m/s)")
  
  v <- str_replace(v, "_stepmin1$", " (step/min)")
  v <- str_replace(v, "_step\\.min\\^\\{-1\\}$", " (step/min)")
  v <- str_replace(v, "_step\\.min\\^\\{_1\\}$", " (step/min)")
  v <- str_replace(v, "_step\\.min\\^\\{\\-1\\}$", " (step/min)")
  v <- str_replace(v, "_step\\.min\\^\\{\\-?1\\}$", " (step/min)")
  
  # 3) Exception demandée : Norm Gait Speed doit être (ua) même si m/s
  v <- str_replace(v, "^Norm Gait Speed \\(m/s\\)$", "Norm Gait Speed (ua)")
  v <- str_replace(v, "^Norm_Gait_Speed \\(m/s\\)$", "Norm Gait Speed (ua)")
  v <- str_replace(v, "^Norm_Gait_Speed$", "Norm Gait Speed (ua)")
  v <- str_replace(v, "^Norm Gait Speed$", "Norm Gait Speed (ua)")
  
  # 4) Rendre lisible : underscores -> espaces
  v <- str_replace_all(v, "_", " ")
  
  # 5) Petites mises en forme (optionnel mais utile pour article)
  v <- str_replace_all(v, "\\bAP\\b", "AP")
  v <- str_replace_all(v, "\\bML\\b", "ML")
  v <- str_replace_all(v, "\\bHS\\b", "HS")
  
  # Double support time : harmoniser
  v <- str_replace_all(v, "Double support time", "Double support time")
  v <- str_replace_all(v, "Double support", "Double support")
  
  # SI / CV : garder le préfixe explicite
  v <- str_replace(v, "^SI ", "S.I. ")
  v <- str_replace(v, "^CV ", "C.V. ")
  
  # NCycles : rendre plus propre
  v <- str_replace(v, "^NCycles Left$", "N cycles Left")
  v <- str_replace(v, "^NCycles Right$", "N cycles Right")
  
  # Trim final
  v <- str_trim(v)
  
  return(v)
}

# ---------------------------------------------------------
# 9) Forcer l'ordre des colonnes : Groupes d'âge puis Surfaces
#    (important : gt ne réordonne pas les colonnes, il faut le faire avant)
# ---------------------------------------------------------
age_order <- c("JeunesEnfants", "Enfants", "Adolescents", "Adultes")
surface_levels <- c("Plat", "Medium", "High")

wanted_cols <- as.vector(outer(age_order, surface_levels, paste, sep = "___"))
wanted_cols <- wanted_cols[wanted_cols %in% names(tab_wide)]  # garde seulement celles présentes

tab_wide <- tab_wide %>%
  select(Variable, all_of(wanted_cols))

# ---------------------------------------------------------
# 10) Construction du mapping : nom technique -> label publication
# ---------------------------------------------------------
variable_labels <- setNames(
  vapply(tab_wide$Variable, make_pretty_label, character(1)),
  tab_wide$Variable
)

# ---------------------------------------------------------
# 11) Création de la table GT + application des labels "publication"
# ---------------------------------------------------------
gt_tbl <- gt(tab_wide) %>%
  cols_label(Variable = "Variables") %>%
  tab_options(table.font.size = px(12)) %>%
  text_transform(
    locations = cells_body(columns = Variable),
    fn = function(x) unname(variable_labels[x])
  )

# ---------------------------------------------------------
# 12) Ajout des spanners : Groupes d’âge -> sous-colonnes Surface
# ---------------------------------------------------------
col_names <- setdiff(names(tab_wide), "Variable")

age_levels <- age_order[age_order %in% sub("___.*$", "", col_names)]
surface_levels <- c("Plat", "Medium", "High")

# Traduction des groupes d'âge (affichage uniquement)
age_labels_en <- c(
  JeunesEnfants = "Young Children",
  Enfants       = "Children",
  Adolescents   = "Adolescents",
  Adultes       = "Adults"
)

surface_labels_en <- c(
  Plat   = "Even",
  Medium = "Medium",
  High   = "High"
)

# (A) Renommer les sous-colonnes
for (cn in col_names) {
  surf <- sub("^.*___", "", cn)
  gt_tbl <- gt_tbl %>% cols_label(!!cn := surface_labels_en[surf])
}

# (B) Ajouter les spanners par groupe d’âge, dans l’ordre défini
for (ag in age_levels) {
  cols_ag <- paste0(ag, "___", surface_levels)
  cols_ag <- cols_ag[cols_ag %in% col_names]
  
  gt_tbl <- gt_tbl %>%
    tab_spanner(label = age_labels_en[ag], columns = all_of(cols_ag))
}

# ---------------------------------------------------------
# 13) Affichage & Export
# ---------------------------------------------------------
gt_tbl <- gt_tbl %>%
  tab_header(
    title = "Descriptive values of gait parameters according to age group and walking surface",
    subtitle = "Values are reported as Mean ± SD (n)"
  )

gt_tbl
View(gt_tbl)

gtsave(gt_tbl, "Table_Descriptive_Gait_AgeGroup_Surface.pdf")

# ---------------------------------------------------------
# 14) Construction figures de l'ensemble des variables d'intérêt
#     (Surface en abscisse)
# ---------------------------------------------------------

output_dir <- "Boxplots_Gait_Results"
if (!dir.exists(output_dir)) dir.create(output_dir)

generate_gait_boxplot <- function(var_name, data) {
  
  # Label "propre" via ta fonction existante
  pretty_title <- make_pretty_label(var_name)
  
  # Préparation des données
  df_plot <- data %>%
    select(AgeGroup, Surface, all_of(var_name)) %>%
    rename(Value = !!sym(var_name)) %>%
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
      )
    ) %>%
    filter(!is.na(Value))
  
  # Création du graphique : Surface en X, fill = AgeGroup
  p <- ggplot(df_plot, aes(x = Surface, y = Value, fill = AgeGroup)) +
    geom_boxplot(
      position = position_dodge(0.85), width = 0.7,
      alpha = 0.7, outlier.shape = NA, color = "black"
    ) +
    geom_jitter(
      aes(group = AgeGroup),
      position = position_dodge(0.85),
      size = 1.2, alpha = 0.3, color = "black"
    ) +
    # Couleurs des âges (ajuste si tu veux une autre palette)
    scale_fill_manual(values = c(
      "Young Children" = "blue",
      "Children"       = "orange",
      "Adolescents"    = "green",
      "Adults"         = "purple"
    )) +
    labs(
      title = paste(pretty_title, "across surfaces and age groups"),
      y = pretty_title,
      x = "Surface",
      fill = "Age Group"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
  
  # Sauvegarde automatique
  file_name_png <- file.path(output_dir, paste0(var_name, "_SurfaceX.png"))
  ggsave(file_name_png, plot = p, width = 8, height = 6, dpi = 300)
  
  file_name_pdf <- file.path(output_dir, paste0(var_name, "_SurfaceX.pdf"))
  ggsave(file_name_pdf, plot = p, width = 10, height = 8)
  
  return(p)
}

message("Génération des boxplots (Surface en X) en cours...")
walk(vars_present, ~generate_gait_boxplot(.x, df))
message("Terminé ! Les graphiques sont dans le dossier : ", output_dir)


# ---------------------------------------------------------
# 15) Génération des RADAR PLOTS par surface
# ---------------------------------------------------------

# 15.1) Définition des variables incluses dans le radar
vars_radar <- c(
  "Mean_Gait_speed_m.s^{_1}", "Mean_Norm_Gait_Speed_m.s^{_1}", "Mean_Step_length_m",
  "Mean_Stride_length_m", "Mean_Norm_Step_length_ua", "Mean_WalkRatio", "Mean_Norm_WR_ua",
  "Mean_Double_support_time_p", "Mean_Cadence_step.min^{_1}", "Mean_Norm_Cadence_ua",
  "Mean_COM_SPARC_Magnitude_ua", "Mean_StepWidth_cm", "Mean_Norm_StepWidth_ua",
  "Mean_MoS_AP_HS_pL0", "Mean_MoS_ML_HS_pL0", "Mean_MoS_AP_HS_mm", "Mean_MoS_ML_HS_mm", "Mean_GVI_ua",
  "CV_Norm_StepWidth_ua", "CV_Gait_speed_m.s^{_1}", "SI_Stride_length_m",
  "SI_Double_support_time_p", "SI_Norm_StepWidth_ua"
)

# 15.2) Vérification des variables présentes et préparation des labels
# On ne garde que les variables réellement disponibles dans df (après standardisation/cleaning)
vars_radar_present <- intersect(vars_radar, names(df))
# Labels lisibles (mise en forme) pour l’affichage autour du radar
radar_labels <- vapply(vars_radar_present, make_pretty_label, character(1))

# 15.3) Dossier d’export
if (!dir.exists("Radar_Plots")) dir.create("Radar_Plots")

# 15.4) Palette de couleurs par groupe d’âge (pour les moyennes) et par domaine
# Remarque : tu convertis ensuite "Adultes" -> "Adults" dans la fonction (sécurisation)
age_colors <- c(
  "Young Children" = "blue",
  "Children"       = "orange",
  "Adolescents"    = "green",
  "Adultes"        = "purple"
)

# --- (A) Définir les domaines et leurs variables (ordre non critique ici)
domains_vars <- list(
  PACE = c(
    "Mean_Gait_speed_m.s^{_1}", "Mean_Norm_Gait_Speed_m.s^{_1}",
    "Mean_Step_length_m", "Mean_Stride_length_m",
    "Mean_Norm_Step_length_ua", "Mean_WalkRatio", "Mean_Norm_WR_ua"
  ),
  RHYTHM = c(
    "Mean_Double_support_time_p", "Mean_Cadence_step.min^{_1}",
    "Mean_Norm_Cadence_ua", "Mean_COM_SPARC_Magnitude_ua"
  ),
  `POSTURAL CONTROL` = c(
    "Mean_StepWidth_cm", "Mean_Norm_StepWidth_ua",
    "Mean_MoS_AP_HS_pL0", "Mean_MoS_ML_HS_pL0",
    "Mean_MoS_AP_HS_mm", "Mean_MoS_ML_HS_mm"
  ),
  ASYMMETRY = c(
    "SI_Stride_length_m", "SI_Double_support_time_p", "SI_Norm_StepWidth_ua"
  ),
  VARIABILITY = c(
    "Mean_GVI_ua", "CV_Norm_StepWidth_ua", "CV_Gait_speed_m.s^{_1}"
  )
)

# --- (B) Couleurs des domaines (tu peux remplacer par tes hex exacts)
domain_colors <- c(
  PACE = "lightblue",
  RHYTHM = "lightcoral",
  `POSTURAL CONTROL` = "palegreen",
  ASYMMETRY = "plum",
  VARIABILITY = "lightyellow"
)

# --- (C) Fonction: associer chaque variable (dans l’ordre du radar) à son domaine
get_domain_for_vars <- function(vars_in_radar, domains_list) {
  dom_vec <- rep(NA_character_, length(vars_in_radar))
  names(dom_vec) <- vars_in_radar
  for (d in names(domains_list)) {
    dom_vec[vars_in_radar %in% domains_list[[d]]] <- d
  }
  dom_vec
}

# Dessine des secteurs (wedge) par domaine selon l'ordre des variables
draw_domain_background <- function(domains_by_var, domain_cols, alpha = 0.18, r = 1) {
  # domains_by_var: vecteur nommé ou non, longueur = nb variables, contenant le nom du domaine pour chaque variable
  n <- length(domains_by_var)
  if (n < 3) return(invisible(NULL))
  
  # Angles des axes (radar classique: premier en haut)
  angles <- seq(0, 2*pi, length.out = n + 1)[1:n] + (pi/2)
  
  # Limites entre axes = milieux angulaires
  bounds <- angles - (pi / n)
  bounds <- c(bounds, bounds[1] + 2*pi)
  
  # Pour regrouper les variables contiguës d'un même domaine
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
    
    # bornes angulaires du bloc contigu
    a_start <- bounds[i1]
    a_end   <- bounds[i2 + 1]
    
    # points du secteur
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

# 15.5) Calcul des bornes min/max globales (scaling)
# Principe :
# - on fixe les bornes min/max de chaque variable sur l’ensemble des participants
# - cela garantit une normalisation cohérente entre surfaces et groupes
radar_min_max <- df %>%
  dplyr::select(dplyr::all_of(vars_radar_present)) %>%
  dplyr::summarise(dplyr::across(
    dplyr::everything(),
    list(
      min = ~min(.x, na.rm = TRUE),
      max = ~max(.x, na.rm = TRUE)
    )
  ))
# Extraction explicite (utile pour debug/contrôle si besoin)
mins_raw <- radar_min_max %>% dplyr::select(dplyr::ends_with("_min")) %>% unlist() %>% as.numeric()
maxs_raw <- radar_min_max %>% dplyr::select(dplyr::ends_with("_max")) %>% unlist() %>% as.numeric()

# 15.6) Nettoyage graphique (pour repartir sur une base propre)
graphics.off()
par(mfrow = c(1, 1))

# 15.7) Fonction : création d’un radar plot pour UNE surface
create_surface_radar <- function(surf_name, df_full, vars, labels) {
  
  # A) Calcul des moyennes par groupe d’âge (pour la surface)
  data_avg <- df_full %>%
    dplyr::filter(Surface == surf_name) %>%
    dplyr::group_by(AgeGroup) %>%
    dplyr::summarise(dplyr::across(dplyr::all_of(vars), ~mean(.x, na.rm = TRUE)), .groups = "drop") %>%
    dplyr::mutate(
      AgeGroup = factor(
        AgeGroup,
        levels = c("JeunesEnfants", "Enfants", "Adolescents", "Adultes"),
        labels = c("Young Children", "Children", "Adolescents", "Adults")
      )
    ) %>%
    dplyr::arrange(AgeGroup)
  
  # B) Extraction des données individuelles (pour la surface)
  data_indiv <- df_full %>%
    dplyr::filter(Surface == surf_name) %>%
    dplyr::mutate(
      AgeGroup = factor(
        AgeGroup,
        levels = c("JeunesEnfants", "Enfants", "Adolescents", "Adultes"),
        labels = c("Young Children", "Children", "Adolescents", "Adults")
      )
    ) %>%
    dplyr::select(AgeGroup, dplyr::all_of(vars))
  
  # C) Récupération des min/max globaux (pour normalisation)
  mins <- as.vector(radar_min_max[grep("_min$", names(radar_min_max))])
  maxs <- as.vector(radar_min_max[grep("_max$", names(radar_min_max))])
  
  # D) Normalisation des moyennes (0–1)
  radar_df_avg <- as.data.frame(data_avg[, -1, drop = FALSE])
  
  normalized_avg <- as.data.frame(lapply(seq_len(ncol(radar_df_avg)), function(i) {
    denom <- (maxs[[i]] - mins[[i]])
    if (is.na(denom) || denom == 0) return(rep(0, nrow(radar_df_avg)))
    (radar_df_avg[, i] - mins[[i]]) / denom
  }))
  
  colnames(normalized_avg) <- labels
  
  # Format attendu par fmsb::radarchart :
  # - 1ère ligne : max (=1)
  # - 2ème ligne : min (=0)
  # - lignes suivantes : données (ici = moyennes par groupe)
  final_radar_avg <- rbind(rep(1, length(vars)), rep(0, length(vars)), normalized_avg)
  
  # E) Normalisation des individus (0–1)
  radar_df_indiv <- data_indiv[, -1, drop = FALSE]
  
  normalized_indiv <- as.data.frame(lapply(seq_len(ncol(radar_df_indiv)), function(i) {
    denom <- (maxs[[i]] - mins[[i]])
    if (is.na(denom) || denom == 0) return(rep(0, nrow(radar_df_indiv)))
    (radar_df_indiv[, i] - mins[[i]]) / denom
  }))
  
  colnames(normalized_indiv) <- labels
  
  # F) Couleurs : correction du libellé Adultes -> Adults
  age_colors_fixed <- age_colors
  if ("Adultes" %in% names(age_colors_fixed) && !("Adults" %in% names(age_colors_fixed))) {
    age_colors_fixed["Adults"] <- age_colors_fixed["Adultes"]
    age_colors_fixed <- age_colors_fixed[names(age_colors_fixed) != "Adultes"]
  }
  
  colors_border_avg <- age_colors_fixed
  colors_in_avg <- grDevices::adjustcolor(colors_border_avg, alpha.f = 0.20)
  
  # G) Tracé en 3 couches (du fond vers l’avant)
  # 1) Cadre (axes, grille, titre) avec polygones invisibles mais "valides"
  # 2) Individus (gris, derrière)
  # 3) Moyennes de groupe (couleur, devant)
  
  # =========================
  # G1) Cadre du radar
  # =========================
  
  ng <- nrow(data_avg)
  transparent <- grDevices::adjustcolor("white", alpha.f = 0)  # ou "transparent"
  
  fmsb::radarchart(
    final_radar_avg,
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
    title = paste("Gait Profile on", surf_name, "Surface across age groups")
  )
  
  # === Préparation des variables pour les étiquettes ===
  nvar <- length(vars)
  angles <- seq(0, 2*pi, length.out = nvar + 1)[1:nvar] + (pi/2)
  
  pct <- c(0.25, 0.50, 0.75, 1.00)
  r_levels <- pct
  
  ticks_real <- sapply(seq_len(nvar), function(i) {
    mins[[i]] + pct * (maxs[[i]] - mins[[i]])
  })
  
  # --- Fond coloré par domaine (AJOUT)
  domains_by_var <- get_domain_for_vars(vars, domains_vars)
  
  par(new = TRUE)  # superpose sur le même repère
  draw_domain_background(
    domains_by_var = domains_by_var,
    domain_cols    = domain_colors,
    alpha          = 0.35,
    r              = 1
  )
  
  # === AJOUT: Affichage des étiquettes de valeurs réelles ===
  # (APRÈS les fonds colorés pour qu'elles soient visibles)
  for (i in seq_len(nvar)) {
    angle <- angles[i]
    
    for (j in seq_along(pct)) {
      r <- r_levels[j]
      
      # Position du texte (légèrement décalé vers l'extérieur)
      x_pos <- r * cos(angle) * 1.05
      y_pos <- r * sin(angle) * 1.05
      
      # Valeur réelle
      val <- round(ticks_real[j, i], 2)
      
      # Afficher le texte avec fond blanc semi-transparent
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
  
  # --- G2) Individus : tracés en gris derrière
  indiv_col <- grDevices::adjustcolor("grey30", alpha.f = 0.18)
  
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
  
  # --- G3) Moyennes : tracés colorés au premier plan
  par(new = TRUE)
  fmsb::radarchart(
    final_radar_avg,
    axistype = 0,
    vlabels = rep("", length(vars)),
    pcol = colors_border_avg,
    pfcol = colors_in_avg,
    plwd = 2.2,
    plty = 1,
    cglcol = NA,
    axislabcol = NA,
    vlcex = 0,
    seg = length(vars)
  )
  
  # --- Légende
  legend(
    x = "bottom",
    legend = names(colors_border_avg),
    inset = -0.15,
    horiz = TRUE,
    bty = "n",
    pch = 20,
    col = colors_border_avg,
    text.col = "black",
    cex = 0.8,
    pt.cex = 1.5,
    xpd = TRUE
  )
}

# 15.8) Export PDF et PNGs: un radar par surface
pdf("Radar_Plots/Gait_Radar_Profiles.pdf", width = 12, height = 12)
par(mfrow = c(1, 1))

for (s in c("Plat", "Medium", "High")) {
  
  create_surface_radar(s, df, vars_radar_present, radar_labels)
}

dev.off()   # <<< fermeture DU PDF uniquement


# 15.9) Export PNG haute qualité (plus lisible) — sans modifier le PDF
# ---------------------------------------------------------
for (s in c("Plat", "Medium", "High")) {
  
  fname <- paste0("Radar_Plots/Gait_Radar_Profile_", s, ".png")
  
  png(filename = fname, width = 5200, height = 5200, res = 600, type = "cairo")
  
  op <- par(no.readonly = TRUE)
  
  par(mfrow = c(1, 1))
  
  # Marges plus petites => le radar prend plus de place
  par(mar = c(8, 5, 5, 5))     # bottom, left, top, right
  
  # Pas de marge externe (sinon ça rétrécit le plot)
  par(oma = c(0, 0, 0, 0))
  
  # Evite que R rajoute une "expansion" d'axes qui réduit visuellement le cercle
  par(xaxs = "i", yaxs = "i")
  
  # Autorise texte/legend hors zone (évite coupures)
  par(xpd = NA)
  
  # Texte légèrement plus petit si besoin (optionnel)
  par(cex = 0.85)
  
  create_surface_radar(s, df, vars_radar_present, radar_labels)
  
  par(op)
  dev.off()
}





## ============================================================
## III. LMM & POST-HOCS SUR VARIABLES D'INTÉRÊT
## ============================================================
## Structure inspirée du script 14__VAS.R :
## 1) ANOVA Type III
## 2) Post-hocs : Effet GLOBAL de la Surface (toutes surfaces confondues)
## 3) Post-hocs : Effet GLOBAL de l'Âge (tous âges confondus)
## 4) Post-hocs : Effet de la Surface par Groupe d'Âge
## 5) Post-hocs : Effet de l'Âge par Surface
## ============================================================

## 1) Chargement des packages nécessaires
library(lme4)
library(lmerTest)
library(emmeans)
library(dplyr)
library(purrr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)

## 2) Chargement et préparation des données
csv_path <- "C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data/ACP_Clustering_DATA.csv"

first_line <- readLines(csv_path, n = 1, warn = FALSE)
delim <- ifelse(grepl(";", first_line), ";", ",")
df <- read_delim(csv_path, delim = delim, show_col_types = FALSE)

# Fonction de standardisation des noms (identique à la partie II)
standardize_names <- function(x) {
  x %>%
    str_trim() %>%
    str_replace_all("[[:space:]]+", "_") %>%
    str_replace_all("[\\-]+", "_") %>%
    str_replace_all("[,;:]+", "_") %>%
    str_replace_all("\\(mm\\)", "mm") %>%
    str_replace_all("\\(%L0\\)", "pL0") %>%
    str_replace_all("\\(ua\\)", "ua") %>%
    str_replace_all("\\(%\\)", "p") %>%
    str_replace_all("\\(m\\.s\\^\\{-1\\}\\)", "ms1") %>%
    str_replace_all("\\(step\\.min\\^\\{-1\\}\\)", "stepmin1") %>%
    str_replace_all("[\\(\\)]", "") %>%
    str_replace_all("__+", "_") %>%
    str_replace_all("_$", "")
}
names(df) <- standardize_names(names(df))

## 3) Définir l'ordre des facteurs
surfaces <- c("Plat","Medium","High")
groups   <- c("JeunesEnfants","Enfants","Adolescents","Adultes")

df <- df %>%
  mutate(
    Participant = factor(Participant),
    Surface     = factor(Surface, levels = surfaces),
    AgeGroup    = factor(AgeGroup, levels = groups, ordered = TRUE)
  )

## 4) Liste des variables à tester
variables_to_test <- c(
  "Mean_Gait_speed_m.s^{_1}", "Mean_Norm_Gait_Speed_m.s^{_1}", "Mean_Step_length_m", "Mean_Stride_length_m",  "Mean_Norm_Step_length_ua", "Mean_WalkRatio", "Mean_Norm_WR_ua",
  
  "Mean_Double_support_time_p", "Mean_Cadence_step.min^{_1}", "Mean_Norm_Cadence_ua", "Mean_COM_SPARC_Magnitude_ua", "Mean_StepTime_s", "Mean_StanceTime_s", "Mean_SwingTime_s",
  
  "Mean_StepWidth_cm", "Mean_Norm_StepWidth_ua", "Mean_MoS_AP_HS_mm", "Mean_MoS_ML_HS_mm", "Mean_MoS_AP_Stance_mm", "Mean_MoS_ML_Stance_mm", "Mean_MoS_AP_HS_pL0", "Mean_MoS_ML_HS_pL0", "Mean_MoS_AP_Stance_pL0", "Mean_MoS_ML_Stance_pL0",
  
  "Mean_GVI_ua", "CV_Norm_StepWidth_ua", "CV_Gait_speed_m.s^{_1}",
  
  "SI_Stride_length_m", "SI_Double_support_time_p", "SI_Norm_StepWidth_ua"
)

## 5) Fonction pour fitter un LMM + ANOVA type III + 4 TYPES DE POST-HOCS
fit_one_variable <- function(var_name, data, p_adjust = "holm") {
  
  if(!var_name %in% names(data)) {
    return(list(
      anova = tibble::tibble(Variable = var_name, Effect = NA, df1 = NA, df2 = NA, F = NA, p = NA,
                             Note = "Variable absente du CSV"),
      ph_surface_global = NULL,
      ph_age_global = NULL,
      ph_surface_within_age = NULL,
      ph_age_within_surface = NULL
    ))
  }
  
  d <- data %>%
    dplyr::select(Participant, Surface, AgeGroup, dplyr::all_of(var_name)) %>%
    dplyr::filter(!is.na(.data[[var_name]]))
  
  fml <- stats::as.formula(paste0("`", var_name, "` ~ Surface * AgeGroup + (1|Participant)"))
  
  m <- lmerTest::lmer(fml, data = d, REML = TRUE)
  
  ## ============================================================
  ## A) ANOVA Type III
  ## ============================================================
  a <- as.data.frame(stats::anova(m, type = 3, ddf = "Satterthwaite"))
  a <- a %>% tibble::rownames_to_column("Effect")
  
  ## Normalisation robuste des noms de colonnes
  a_df <- a %>%
    dplyr::mutate(
      Variable = var_name,
      df1 = dplyr::coalesce(
        if("NumDF" %in% names(a)) .data$NumDF else NA_real_,
        if("Df" %in% names(a)) .data$Df else NA_real_
      ),
      df2 = if("DenDF" %in% names(a)) .data$DenDF else NA_real_,
      F = dplyr::coalesce(
        if("F value" %in% names(a)) .data$`F value` else NA_real_,
        if("F.value" %in% names(a)) .data$F.value else NA_real_
      ),
      p = dplyr::coalesce(
        if("Pr(>F)" %in% names(a)) .data$`Pr(>F)` else NA_real_,
        if("p.value" %in% names(a)) .data$p.value else NA_real_,
        if("P" %in% names(a)) .data$P else NA_real_
      )
    ) %>%
    dplyr::select(Variable, Effect, df1, df2, F, p)
  
  ## ============================================================
  ## B) Post-hocs : Effet GLOBAL de la Surface (toutes surfaces confondues)
  ## ============================================================
  em_surface_global <- emmeans::emmeans(m, ~ Surface)
  ph_surf_global <- as.data.frame(pairs(em_surface_global, adjust = p_adjust)) %>%
    dplyr::mutate(Variable = var_name, Analysis = "Global Surface") %>%
    dplyr::relocate(Variable, Analysis)
  
  ## ============================================================
  ## C) Post-hocs : Effet GLOBAL de l'Âge (tous âges confondus)
  ## ============================================================
  em_age_global <- emmeans::emmeans(m, ~ AgeGroup)
  ph_age_global <- as.data.frame(pairs(em_age_global, adjust = p_adjust)) %>%
    dplyr::mutate(Variable = var_name, Analysis = "Global AgeGroup") %>%
    dplyr::relocate(Variable, Analysis)
  
  ## ============================================================
  ## D) Post-hocs : Effet de la Surface par Groupe d'Âge
  ## ============================================================
  em_surf_by_age <- emmeans::emmeans(m, ~ Surface | AgeGroup)
  ph1 <- as.data.frame(pairs(em_surf_by_age, adjust = p_adjust)) %>%
    dplyr::mutate(Variable = var_name, Analysis = "Surface within AgeGroup") %>%
    dplyr::relocate(Variable, Analysis)
  
  ## ============================================================
  ## E) Post-hocs : Effet de l'Âge par Surface
  ## ============================================================
  em_age_by_surf <- emmeans::emmeans(m, ~ AgeGroup | Surface)
  ph2 <- as.data.frame(pairs(em_age_by_surf, adjust = p_adjust)) %>%
    dplyr::mutate(Variable = var_name, Analysis = "AgeGroup within Surface") %>%
    dplyr::relocate(Variable, Analysis)
  
  ## ============================================================
  ## F) Calcul des R² (Nakagawa)
  ## ============================================================
  r2_vals <- performance::r2_nakagawa(m)
  
  list(
    model = m,
    anova = a_df,
    ph_surface_global = ph_surf_global,
    ph_age_global = ph_age_global,
    ph_surface_within_age = ph1,
    ph_age_within_surface = ph2,
    r2 = tibble::tibble(
      Variable = var_name,
      R2_Marginal = r2_vals$R2_marginal,
      R2_Conditional = r2_vals$R2_conditional
    )
  )
}

## 6) Lancer sur toutes les variables
message("Analyse en cours pour ", length(variables_to_test), " variables...")
results <- map(variables_to_test, fit_one_variable, data = df, p_adjust = "holm")

## 7) Extraire les résultats
anova_all <- bind_rows(map(results, "anova"))

posthoc_surface_global <- bind_rows(compact(map(results, "ph_surface_global")))
posthoc_age_global <- bind_rows(compact(map(results, "ph_age_global")))
posthoc_surface_within_age <- bind_rows(compact(map(results, "ph_surface_within_age")))
posthoc_age_within_surface <- bind_rows(compact(map(results, "ph_age_within_surface")))

r2_all <- bind_rows(compact(map(results, "r2")))

## 8) Appliquer FDR par famille d'effets
apply_fdr <- function(tbl) {
  tbl %>%
    filter(Effect %in% c("Surface","AgeGroup","Surface:AgeGroup")) %>%
    mutate(EffectFamily = case_when(
      Effect == "Surface" ~ "Surface",
      Effect == "AgeGroup" ~ "AgeGroup",
      Effect == "Surface:AgeGroup" ~ "Interaction"
    )) %>%
    group_by(EffectFamily) %>%
    mutate(p_fdr = p.adjust(p, method = "BH")) %>%
    ungroup()
}

anova_all_fdr <- apply_fdr(anova_all)

## ============================================================
## 9) RÉCAPITULATIF POUR LA DISCUSSION
## ============================================================

# A) Variables influencées par chaque effet (p_fdr < 0.05)
listes_effets <- anova_all_fdr %>%
  filter(p_fdr < 0.05) %>%
  group_by(EffectFamily) %>%
  summarise(
    Nombre_Variables = n(),
    Variables = paste(unique(Variable), collapse = " | ")
  ) %>%
  rename(Type_Effet = EffectFamily)

# B) Tableau de maturation par surface
maturation_final <- posthoc_age_within_surface %>%
  filter(grepl("Adultes", contrast)) %>%
  mutate(Groupe = str_remove(contrast, " - Adultes") %>% str_trim()) %>%
  filter(p.value > 0.05) %>%
  group_by(Variable, Surface) %>%
  arrange(factor(Groupe, levels = c("JeunesEnfants", "Enfants", "Adolescents"))) %>%
  slice(1) %>%
  select(Variable, Surface, Groupe) %>%
  pivot_wider(names_from = Surface, values_from = Groupe)

## ============================================================
## 10) EXPORT DES RÉSULTATS
## ============================================================
out_dir <- file.path(dirname(csv_path), "R_LMM_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# A) Fichier Excel multi-onglets PRINCIPAL
xlsx_path <- file.path(out_dir, "LMM_ANOVA_and_Posthocs_COMPLET.xlsx")
wb <- createWorkbook()

# Onglet 1 : ANOVA Type III avec FDR
addWorksheet(wb, "ANOVA_TypeIII_FDR")
writeData(wb, "ANOVA_TypeIII_FDR", anova_all_fdr)

# Onglet 2 : R² (Marginal et Conditionnel)
addWorksheet(wb, "R2_Effect_Sizes")
writeData(wb, "R2_Effect_Sizes", r2_all)

# Onglet 3 : Post-hocs GLOBAL Surface
addWorksheet(wb, "Posthoc_Surface_Global")
writeData(wb, "Posthoc_Surface_Global", posthoc_surface_global)

# Onglet 4 : Post-hocs GLOBAL Âge
addWorksheet(wb, "Posthoc_Age_Global")
writeData(wb, "Posthoc_Age_Global", posthoc_age_global)

# Onglet 5 : Post-hocs Surface par Âge
addWorksheet(wb, "Posthoc_Surface_by_Age")
writeData(wb, "Posthoc_Surface_by_Age", posthoc_surface_within_age)

# Onglet 6 : Post-hocs Âge par Surface
addWorksheet(wb, "Posthoc_Age_by_Surface")
writeData(wb, "Posthoc_Age_by_Surface", posthoc_age_within_surface)

saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# B) Fichier Excel SYNTHÈSE pour la discussion
xlsx_synthese <- file.path(out_dir, "Synthese_Resultats_LMM_Discussion.xlsx")
wb_synth <- createWorkbook()

addWorksheet(wb_synth, "Liste_Effets_Significatifs")
writeData(wb_synth, "Liste_Effets_Significatifs", listes_effets)

addWorksheet(wb_synth, "Maturation_Par_Surface")
writeData(wb_synth, "Maturation_Par_Surface", maturation_final)

saveWorkbook(wb_synth, xlsx_synthese, overwrite = TRUE)

# C) Export CSV (optionnel)
write_csv(anova_all_fdr, file.path(out_dir, "ANOVA_TypeIII_with_FDR.csv"))
write_csv(r2_all, file.path(out_dir, "R2_Effect_Sizes.csv"))
write_csv(posthoc_surface_global, file.path(out_dir, "Posthoc_Surface_GLOBAL_Holm.csv"))
write_csv(posthoc_age_global, file.path(out_dir, "Posthoc_Age_GLOBAL_Holm.csv"))
write_csv(posthoc_surface_within_age, file.path(out_dir, "Posthoc_Surface_within_AgeGroup_Holm.csv"))
write_csv(posthoc_age_within_surface, file.path(out_dir, "Posthoc_AgeGroup_within_Surface_Holm.csv"))

## ============================================================
## 11) RÉCAPITULATIF CONSOLE
## ============================================================
message("\n=== ANALYSE TERMINÉE ===")
message("ANOVA Type III + FDR: ", nrow(anova_all_fdr), " lignes")
message("R² (Nakagawa) calculés pour: ", nrow(r2_all), " variables")
message("\n--- POST-HOCS (4 analyses distinctes) ---")
message("1) Effet GLOBAL Surface: ", nrow(posthoc_surface_global), " comparaisons")
message("2) Effet GLOBAL Âge: ", nrow(posthoc_age_global), " comparaisons")
message("3) Surface par Âge: ", nrow(posthoc_surface_within_age), " comparaisons")
message("4) Âge par Surface: ", nrow(posthoc_age_within_surface), " comparaisons")
message("\n--- FICHIERS EXPORTÉS ---")
message("Fichier COMPLET: ", xlsx_path)
message("Fichier SYNTHÈSE: ", xlsx_synthese)
message("\n✓ Prêt pour la rédaction de la discussion !")
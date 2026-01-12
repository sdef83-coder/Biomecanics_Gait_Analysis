## ============================================================
## LMM Surface (répété) x AgeGroup (inter) + post-hocs emmeans
## ============================================================

## 0) Packages
pkgs <- c("readr","dplyr","stringr","tibble","lme4","lmerTest","car","emmeans","purrr","openxlsx")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if(length(to_install) > 0) install.packages(to_install, dependencies = TRUE)

library(readr)
library(dplyr)
library(stringr)
library(tibble)
library(lme4)
library(lmerTest)  # Satterthwaite df + p-values pour lmer
library(car)       # Anova(type=3)
library(emmeans)   # post-hocs
library(purrr)
library(openxlsx)

## 1) Chemin vers ton CSV long
csv_path <- "C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data/DATA_all_prepared.csv"

## 2) Lire le CSV
df <- read_delim(file = csv_path, delim = ";",show_col_types = FALSE)
names(df)

## 3) Définir l'ordre des facteurs (comme ton MATLAB)
surfaces <- c("Plat","Medium","High")
groups   <- c("JeunesEnfants","Enfants","Adolescents","Adultes")

df <- df %>%
  mutate(
    Participant = factor(Participant),
    Surface     = factor(Surface, levels = surfaces),
    AgeGroup    = factor(AgeGroup, levels = groups, ordered = TRUE)
  )

## 4) Liste des variables à tester (noms EXACTS des colonnes du CSV)
variables_to_test <- c(
  "Mean_Single support time (%)",
  "Mean_Double support time (%)",
  "Mean_BaseOfSupport (cm)",
  "Mean_StepWidth (cm)",
  "Mean_Gait speed (m.s^{-1})",
  "Mean_Stride length (m)",
  "Mean_Stride time (s)",
  "Mean_WalkRatio",
  "Mean_Norm WR (ua)",
  "Mean_Cadence (step.min^{-1})",
  "Mean_Norm Step length (ua)",
  "Mean_Norm Cadence (ua)",
  "Mean_Norm StepWidth (ua)",
  "Mean_Norm Gait Speed (m.s^{-1})",
  
  "CV_Single support time (%)",
  "CV_Double support time (%)",
  "CV_BaseOfSupport (cm)",
  "CV_StepWidth (cm)",
  "CV_Gait speed (m.s^{-1})",
  "CV_Stride length (m)",
  "CV_Stride time (s)",
  "CV_WalkRatio",
  "CV_Norm WR (ua)",
  "CV_Cadence (step.min^{-1})",
  "CV_Norm Step length (ua)",
  "CV_Norm Cadence (ua)",
  "CV_Norm StepWidth (ua)",
  "CV_Norm Gait Speed (m.s^{-1})",
  
  "Mean_MoS AP HS (mm)",
  "Mean_MoS ML HS (mm)",
  "Mean_MoS AP Stance (mm)",
  "Mean_MoS ML Stance (mm)",
  "Mean_MoS AP HS (%L0)",
  "Mean_MoS ML HS (%L0)",
  "Mean_MoS AP Stance (%L0)",
  "Mean_MoS ML Stance (%L0)",
  
  "Mean_COM SPARC Magnitude (ua)",
  "Mean_COM LDLJ Magnitude (ua)",
  "Mean_STERN SPARC Magnitude (ua)",
  "Mean_STERN LDLJ Magnitude (ua)",
  
  "SI_Stride time (s)",
  "SI_Stride length (m)",
  "SI_Double support time (%)",
  "SI_Single support time (%)",
  "SI_BaseOfSupport (cm)",
  "SI_StepWidth (cm)",
  "SI_WalkRatio",
  "SI_Norm WR (ua)",
  "SI_Cadence (step.min^{-1})",
  "SI_Norm Step length (ua)",
  "SI_Norm Cadence (ua)",
  "SI_Norm StepWidth (ua)"
)

## 5) Fonction utilitaire: sécuriser les noms qui contiennent espaces et symboles
bt <- function(x) paste0("`", x, "`")  # backticks

## 6) Fonction pour fitter un LMM + ANOVA type III + post-hocs
fit_one_variable <- function(var_name, data, p_adjust = "holm") {
  
  if(!var_name %in% names(data)) {
    return(list(
      anova = tibble::tibble(Variable = var_name, Effect = NA, df1 = NA, df2 = NA, F = NA, p = NA,
                             Note = "Variable absente du CSV"),
      ph_surface_within_age = NULL,
      ph_age_within_surface = NULL
    ))
  }
  
  d <- data %>%
    dplyr::select(Participant, Surface, AgeGroup, dplyr::all_of(var_name)) %>%
    dplyr::filter(!is.na(.data[[var_name]]))
  
  fml <- stats::as.formula(paste0("`", var_name, "` ~ Surface * AgeGroup + (1|Participant)"))
  
  m <- lmerTest::lmer(fml, data = d, REML = TRUE)
  
  ## ANOVA Type III
  a <- as.data.frame(stats::anova(m, type = 3, ddf = "Satterthwaite"))
  
  # Afficher les noms de colonnes pour débugger (à retirer après)
  # print(paste("Colonnes ANOVA pour", var_name, ":", paste(names(a), collapse = ", ")))
  
  a <- a %>%
    tibble::rownames_to_column("Effect")
  
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
  
  ## Post-hocs
  ## Post-hocs
  em_surf_by_age <- emmeans::emmeans(m, ~ Surface | AgeGroup)
  ph1 <- as.data.frame(pairs(em_surf_by_age, adjust = p_adjust)) %>%
    dplyr::mutate(Variable = var_name, ContrastFamily = "Surface within AgeGroup") %>%
    dplyr::relocate(Variable, ContrastFamily)
  
  em_age_by_surf <- emmeans::emmeans(m, ~ AgeGroup | Surface)
  ph2 <- as.data.frame(pairs(em_age_by_surf, adjust = p_adjust)) %>%
    dplyr::mutate(Variable = var_name, ContrastFamily = "AgeGroup within Surface") %>%
    dplyr::relocate(Variable, ContrastFamily)
  
  list(
    model = m,
    anova = a_df,
    ph_surface_within_age = ph1,
    ph_age_within_surface = ph2
  )
}

## 7) Lancer sur toutes les variables
results <- map(variables_to_test, fit_one_variable, data = df, p_adjust = "holm")

anova_all <- bind_rows(map(results, "anova"))

posthoc_surface_within_age <- bind_rows(compact(map(results, "ph_surface_within_age")))
posthoc_age_within_surface <- bind_rows(compact(map(results, "ph_age_within_surface")))

## 8) (Option) FDR par “famille d’effets” sur les p des effets fixes (comme ton MATLAB)
## Familles: Surface / AgeGroup / Surface:AgeGroup
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

## 9) Export Excel (pratique pour papier + traçabilité)
out_dir <- file.path(dirname(csv_path), "R_LMM_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

xlsx_path <- file.path(out_dir, "LMM_ANOVA_and_Posthocs.xlsx")
wb <- createWorkbook()

addWorksheet(wb, "ANOVA_TypeIII")
writeData(wb, "ANOVA_TypeIII", anova_all_fdr)

addWorksheet(wb, "Posthoc_Surface_by_Age")
writeData(wb, "Posthoc_Surface_by_Age", posthoc_surface_within_age)

addWorksheet(wb, "Posthoc_Age_by_Surface")
writeData(wb, "Posthoc_Age_by_Surface", posthoc_age_within_surface)

saveWorkbook(wb, xlsx_path, overwrite = TRUE)

## 10) Export CSV (si tu veux aussi)
write_csv(anova_all_fdr, file.path(out_dir, "ANOVA_TypeIII_with_FDR.csv"))
write_csv(posthoc_surface_within_age, file.path(out_dir, "Posthoc_Surface_within_AgeGroup_Holm.csv"))
write_csv(posthoc_age_within_surface, file.path(out_dir, "Posthoc_AgeGroup_within_Surface_Holm.csv"))

## 11) Petit recap console (facultatif)
message("=== Terminé ===")
message("ANOVA Type III + FDR: ", nrow(anova_all_fdr), " lignes")
message("Posthoc Surface|AgeGroup: ", nrow(posthoc_surface_within_age), " comparaisons")
message("Posthoc AgeGroup|Surface: ", nrow(posthoc_age_within_surface), " comparaisons")
message("Fichier Excel: ", xlsx_path)

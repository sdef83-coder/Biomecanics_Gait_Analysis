setwd("C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data")

library(dplyr)
library(readr)

chemin <- "C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data"
outdir <- "C:/Users/silve/Desktop/DOCTORAT/UNIV MONTREAL/TRAVAUX-THESE/Surfaces_Irregulieres/Datas/Script/gaitAnalysisGUI/result/Statistical_Analysis_LMM/Prepared_Data"
  
df <- read.csv(
  file.path(chemin, "ACP_Clustering_DATA.csv"),
  sep = ";",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

names(df)
str(df)
df$Surface <- factor(df$Surface, levels = c("Plat", "Medium", "High"))



# MATRICE CORRELATION SUR PLAT
df_plat <- df %>%
  filter(Surface == "Plat") %>%
  select(where(is.numeric), -AgeMonths)

cor_plat <- cor(df_plat, use = "pairwise.complete.obs", method = "pearson")
cor_plat_abs <- abs(cor_plat)

round(cor_plat, 2)



# MATRICE CORRELATION SUR MEDIUM
df_medium <- df %>%
  filter(Surface == "Medium") %>%
  select(where(is.numeric), -AgeMonths)

cor_medium <- cor(df_plat, use = "pairwise.complete.obs", method = "pearson")
cor_medium_abs <- abs(cor_medium)

round(cor_plat, 2)



# MATRICE DE CORRELATION HIGH
df_high <- df %>%
  filter(Surface == "High") %>%
  select(where(is.numeric), -AgeMonths, -NCycles_Left, -NCycles_Right)

cor_high <- cor(df_plat, use = "pairwise.complete.obs", method = "pearson")
cor_high_abs <- abs(cor_high)

round(cor_plat, 2)


# Export en valeurs signées
write.csv(cor_plat, file = "Multicolinearity_Even.csv")
write.csv(cor_medium, file = "Multicolinearity_Medium.csv")
write.csv(cor_high, file = "Multicolinearity_High.csv")

# Export en valeurs absolues
write.csv(cor_plat_abs, file = "Multicolinearity_Even_abs.csv")
write.csv(cor_medium_abs, file = "Multicolinearity_Medium_abs.csv")
write.csv(cor_high_abs, file = "Multicolinearity_High_abs.csv")
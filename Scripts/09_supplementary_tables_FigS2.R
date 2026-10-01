###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 09 - Supplementary Tables S1 and S2, within-temperature correlations (Fig. S2)###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Collects the outputs of scripts 01-08 into two supplementary tables and computes the
## trait correlations within temperature. Run scripts 01-08 first (or 00_run_all.R).
##
## Table S1 : for every trait, cultivar and temperature the number of biological
##            replicates (n), mean, standard deviation, standard error and the Tukey
##            letter of the within-temperature comparison of cultivars.
## Table S2 : for every trait the scale of analysis, the Type II ANOVA F and P of
##            genotype, temperature and their interaction with the degrees of freedom,
##            and the Shapiro-Wilk and Levene P-values on the scale analysed.
## Fig. S2  : Spearman correlations among the 22 traits within temperature. Each trait is
##            centred within temperature (the mean of the six cultivars at that
##            temperature is subtracted) before correlation, so that the coefficients
##            describe covariation among cultivars with the temperature effect removed
##            (Fig. S1 of script 08 correlates the uncentred means). P-values use the
##            t approximation with N - k - 1 degrees of freedom, N observations and k
##            temperatures: 14 for most pairs, 9 for pairs with anthocyanin, which was
##            measured at two temperatures.
##
## Input : results/01_* to results/08_* (csv files written by scripts 01-08)
## Output: results/09_TableS1_means_SE_letters.csv, results/09_TableS2_anova.csv,
##         results/09_spearman_within_temperature_r.csv, _P.csv, _pairs.csv
##         figures/FigS2_spearman_within_temperature.tiff

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/09_supplementary_tables_FigS2.R")
main_dir <- here::here()
setwd(main_dir)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(corrplot)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)
stopifnot("Run scripts 01-08 first" = file.exists("results/08_trait_means.csv"))

## Temperature is read as text in every file (the gas exchange files store it as a number)
rd <- function(path) {
  d <- read.csv(path, check.names = FALSE, encoding = "UTF-8")
  if ("Temperature" %in% names(d)) d$Temperature <- as.character(d$Temperature)
  d
}

## Cultivar codes and names, as in the figure legends
cultivars <- c(TAR = "Mv Tarsoly", KOL = "Mv Kolompos", LUC = "Mv Lucilla",
               PIR = "Mv Pirkadat", IKV = "Mv Ikva", DAN = "Mv Dandár")
to_code <- function(g) ifelse(g %in% names(cultivars), g, names(cultivars)[match(g, cultivars)])
to_temp <- function(t) paste0(sub("[^0-9].*$", "", t), " °C")       # "20", "20C", "20°C" -> "20 °C"

## Traits in the order of the manuscript, with unit, figure and scale of analysis
trait_info <- tribble(
  ~key,            ~Trait,                              ~Unit,               ~Group,                        ~Figure,   ~Scale,
  "FW",            "Shoot fresh weight",                "g per plant",       "Growth",                      "Fig. 1",  "log",
  "Pn",            "Net photosynthetic rate (Pn)",      "µmol m-2 s-1",      "Gas exchange",                "Fig. 2A", "log",
  "gs",            "Stomatal conductance (gs)",         "mmol m-2 s-1",      "Gas exchange",                "Fig. 2B", "log",
  "E",             "Transpiration rate (E)",            "mmol m-2 s-1",      "Gas exchange",                "Fig. 2C", "log",
  "WUE",           "Water-use efficiency (WUE)",        "µmol mmol-1",       "Gas exchange",                "Fig. 2D", "log",
  "Chl_ab",        "Total chlorophyll (Chl a+b)",       "µg g-1 FW",         "Pigments",                    "Fig. 3A", "log",
  "Carotenoid",    "Carotenoids",                       "µg g-1 FW",         "Pigments",                    "Fig. 3B", "log",
  "Anthocyanin",   "Anthocyanin",                       "nmol g-1 FW",       "Pigments",                    "Fig. 3C", "log",
  "Chl_a_b_ratio", "Chlorophyll a/b ratio",             "ratio",             "Pigments",                    "Fig. 3D", "log",
  "APX",           "Ascorbate peroxidase (APX)",        "nkat g-1 FW",       "AsA-GSH enzymes and GST",     "Fig. 4A", "none",
  "DHAR",          "Dehydroascorbate reductase (DHAR)", "nkat g-1 FW",       "AsA-GSH enzymes and GST",     "Fig. 4B", "log",
  "GR",            "Glutathione reductase (GR)",        "nkat g-1 FW",       "AsA-GSH enzymes and GST",     "Fig. 4C", "log",
  "GST",           "Glutathione S-transferase (GST)",   "nkat g-1 FW",       "AsA-GSH enzymes and GST",     "Fig. 4D", "sqrt",
  "AKR",           "Aldo-keto reductase (AKR)",         "nkat g-1 FW",       "Reactive-carbonyl enzymes",   "Fig. 5A", "none",
  "AOR_AER",       "Alkenal/alkenone oxidoreductase (AOR/AER)", "nkat g-1 FW", "Reactive-carbonyl enzymes", "Fig. 5B", "sqrt",
  "GLY",           "Glyoxalase I (GLY I)",              "nkat g-1 FW",       "Reactive-carbonyl enzymes",   "Fig. 5C", "sqrt",
  "PAL",           "Phenylalanine ammonia-lyase (PAL)", "nkat g-1 FW",       "Phenylpropanoid",             "Fig. 6",  "none",
  "SPM",           "Spermine (SPM)",                    "µg g-1 FW",         "Polyamines",                  "Fig. 7",  "sqrt",
  "SPD",           "Spermidine (SPD)",                  "µg g-1 FW",         "Polyamines",                  "Fig. 7",  "sqrt",
  "CAD",           "Cadaverine (CAD)",                  "µg g-1 FW",         "Polyamines",                  "Fig. 7",  "sqrt",
  "PUT",           "Putrescine (PUT)",                  "µg g-1 FW",         "Polyamines",                  "Fig. 7",  "sqrt",
  "DAP",           "1,3-Diaminopropane (DAP)",          "µg g-1 FW",         "Polyamines",                  "Fig. 7",  "log")

######### TABLE S1 - n, MEAN, SD, SE AND TUKEY LETTER ########################################

## Means (original units) of every script, brought to the columns key, Genotype, Temperature
means_all <- bind_rows(
  rd("results/01_FW_means.csv")           %>% mutate(key = "FW"),
  rd("results/02_gas_exchange_means.csv") %>% rename(key = trait),
  rd("results/03_pigments_means.csv")     %>% rename(key = trait),
  rd("results/04_antioxidant_means.csv")  %>% rename(key = trait),
  rd("results/05_carbonyl_means.csv")     %>% rename(key = trait) %>%
    mutate(key = dplyr::recode(key, "AOR/AER" = "AOR_AER", "GLY I" = "GLY")) %>% dplyr::select(-letter),
  rd("results/06_PAL_means.csv")          %>% rename(key = trait),
  rd("results/07_polyamines_means.csv")   %>% rename(key = Polyamine) %>% dplyr::select(-letter)) %>%
  transmute(key, Code = to_code(Genotype), Temperature = to_temp(Temperature), n, mean, sd, se)

## Tukey letters of every script
letters_all <- bind_rows(
  rd("results/01_FW_tukey_letters.csv")           %>% mutate(key = "FW"),
  rd("results/02_gas_exchange_tukey_letters.csv") %>% rename(key = trait),
  rd("results/03_pigments_tukey_letters.csv")     %>% rename(key = trait),
  rd("results/04_antioxidant_tukey_letters.csv")  %>% rename(key = trait),
  rd("results/05_carbonyl_tukey_letters.csv")     %>% rename(key = trait),
  rd("results/06_PAL_tukey_letters.csv")          %>% rename(key = trait),
  rd("results/07_polyamines_tukey_letters.csv")   %>% rename(key = Polyamine)) %>%
  transmute(key, Code = to_code(Genotype), Temperature = to_temp(Temperature), letter)

table_S1 <- means_all %>%
  left_join(letters_all, by = c("key", "Code", "Temperature")) %>%
  left_join(trait_info, by = "key") %>%
  mutate(key      = factor(key, levels = trait_info$key),
         Code     = factor(Code, levels = names(cultivars)),
         Cultivar = cultivars[as.character(Code)]) %>%
  arrange(key, Temperature, Code) %>%
  transmute(`Trait group` = Group, Trait, Unit, Figure, Temperature, Code, Cultivar,
            n, Mean = signif(mean, 4), SD = signif(sd, 3), SE = signif(se, 3),
            `Tukey letter` = letter)

stopifnot(!anyNA(table_S1$`Tukey letter`), !anyNA(table_S1$Mean))
cat("Table S1:", nrow(table_S1), "rows,", length(unique(table_S1$Trait)), "traits\n")
cat("Replicates per cell, range by trait group:\n")
print(table_S1 %>% group_by(`Trait group`) %>% summarise(n_min = min(n), n_max = max(n)))
write.csv(table_S1, "results/09_TableS1_means_SE_letters.csv", row.names = FALSE, fileEncoding = "UTF-8")

######### TABLE S2 - TYPE II ANOVA AND ASSUMPTION CHECKS #####################################

anova_all <- bind_rows(
  rd("results/01_FW_anova.csv")           %>% mutate(trait = "FW"),
  rd("results/02_gas_exchange_anova.csv"),
  rd("results/03_pigments_anova.csv"),
  rd("results/04_antioxidant_anova.csv")  %>% dplyr::select(-scale),
  rd("results/05_carbonyl_anova.csv")     %>% dplyr::select(-scale),
  rd("results/06_PAL_anova.csv")          %>% dplyr::select(-scale),
  rd("results/07_polyamines_anova.csv")   %>% dplyr::select(-scale)) %>%
  rename(key = trait, F = `F value`, P = `Pr(>F)`)

anova_wide <- anova_all %>%
  mutate(term = dplyr::recode(term, "Genotype" = "G", "Temperature" = "T",
                       "Genotype:Temperature" = "GxT", "Residuals" = "res")) %>%
  dplyr::select(key, term, Df, F, P) %>%
  pivot_wider(names_from = term, values_from = c(Df, F, P)) %>%
  dplyr::select(key, Df_G, Df_T, Df_GxT, Df_res, F_G, P_G, F_T, P_T, F_GxT, P_GxT)

## Assumption tests on the scale analysed ("raw" is the label of the untransformed scale)
assump_all <- bind_rows(
  rd("results/01_FW_assumptions.csv")           %>% mutate(trait = "FW"),
  rd("results/02_gas_exchange_assumptions.csv"),
  rd("results/03_pigments_assumptions.csv"),
  rd("results/04_antioxidant_assumptions.csv"),
  rd("results/05_carbonyl_assumptions.csv"),
  rd("results/06_PAL_assumptions.csv"),
  rd("results/07_polyamines_assumptions.csv")) %>%
  rename(key = trait) %>%
  mutate(scale = dplyr::recode(scale, "raw" = "none")) %>%
  distinct(key, scale, .keep_all = TRUE)

table_S2 <- trait_info %>%
  left_join(anova_wide, by = "key") %>%
  left_join(assump_all, by = c("key", "Scale" = "scale")) %>%
  transmute(`Trait group` = Group, Trait, Figure, `Scale of analysis` = Scale,
            `df G` = Df_G, `df T` = Df_T, `df GxT` = Df_GxT, `df residual` = Df_res,
            `F genotype` = signif(F_G, 4),          `P genotype` = signif(P_G, 3),
            `F temperature` = signif(F_T, 4),       `P temperature` = signif(P_T, 3),
            `F interaction` = signif(F_GxT, 4),     `P interaction` = signif(P_GxT, 3),
            `Shapiro-Wilk P` = signif(shapiro_P, 3), `Levene P` = signif(levene_P, 3))

stopifnot(!anyNA(table_S2$`F genotype`), !anyNA(table_S2$`Shapiro-Wilk P`))
cat("\nTable S2:", nrow(table_S2), "traits\n")
print(as.data.frame(table_S2 %>% dplyr::select(Trait, `Scale of analysis`, `P genotype`,
                                               `P temperature`, `P interaction`)), digits = 3)
write.csv(table_S2, "results/09_TableS2_anova.csv", row.names = FALSE, fileEncoding = "UTF-8")

######### FIG. S2 - SPEARMAN CORRELATIONS WITHIN TEMPERATURE ##################################

display_names <- c(Anthocyanin_nmol_gFW = "Anthocyanin", Chl_ab = "Chl a+b", Chl_a_b_ratio = "Chl a/b")
relabel <- function(x) unname(ifelse(x %in% names(display_names), display_names[x], x))

all_means <- rd("results/08_trait_means.csv")
traits    <- setdiff(names(all_means), c("Genotype", "Temperature", "VPD"))

## Centre every trait within temperature: deviation of a cultivar from the mean of the six
## cultivars at that temperature
centred <- all_means %>%
  group_by(Temperature) %>%
  mutate(across(all_of(traits), ~ .x - mean(.x, na.rm = TRUE))) %>%
  ungroup()
cen_mat <- as.matrix(centred[, traits])

cor_within <- cor(cen_mat, method = "spearman", use = "pairwise.complete.obs")

## P-values: t approximation with N - k - 1 degrees of freedom (k temperature means estimated)
p_within <- matrix(NA_real_, length(traits), length(traits), dimnames = list(traits, traits))
df_mat   <- p_within
for (i in seq_along(traits)) for (j in seq_along(traits)) {
  if (i == j) next
  ok  <- complete.cases(cen_mat[, c(i, j)])
  N   <- sum(ok)
  k   <- length(unique(centred$Temperature[ok]))
  df  <- N - k - 1
  rho <- cor_within[i, j]
  p_within[i, j] <- 2 * pt(-abs(rho * sqrt(df / (1 - rho^2))), df)
  df_mat[i, j]   <- df
}
diag(p_within) <- 0

dimnames(cor_within) <- list(relabel(traits), relabel(traits))
dimnames(p_within)   <- dimnames(cor_within)
dimnames(df_mat)     <- dimnames(cor_within)

tiff("figures/FigS2_spearman_within_temperature.tiff", width = 13, height = 8.6, units = "in",
     res = 300, compression = "lzw")
corrplot.mixed(cor_within,
               lower         = "ellipse",
               upper         = "number",
               lower.col     = colorRampPalette(c("#2A9D8F", "white", "#9A031E"))(200),
               upper.col     = colorRampPalette(c("#2A9D8F", "white", "#9A031E"))(200),
               tl.col        = "black",
               tl.cex        = 0.85,
               tl.pos        = "lt",
               number.cex    = 0.65,
               number.digits = 2,
               p.mat         = p_within,
               sig.level     = 0.05,
               insig         = "blank",
               mar           = c(0, 0, 2, 0))
invisible(dev.off())
cat("\nFigure written to figures/FigS2_spearman_within_temperature.tiff\n")

## Tables: coefficients, P-values, and each pair next to its pooled coefficient of Fig. S1
pooled   <- rd("results/08_spearman_r.csv") %>% column_to_rownames("Trait") %>% as.matrix()
pooled_P <- rd("results/08_spearman_P.csv") %>% column_to_rownames("Trait") %>% as.matrix()
stopifnot(identical(rownames(pooled), rownames(cor_within)))
ut <- upper.tri(cor_within)
pairs_tab <- data.frame(
  Trait_1                = rownames(cor_within)[row(cor_within)[ut]],
  Trait_2                = colnames(cor_within)[col(cor_within)[ut]],
  rho_pooled_FigS1       = round(pooled[ut], 3),
  P_pooled               = signif(pooled_P[ut], 3),
  rho_within_temperature = round(cor_within[ut], 3),
  df                     = df_mat[ut],
  P_within_temperature   = signif(p_within[ut], 3)) %>%
  arrange(desc(abs(rho_within_temperature)))

cat("\nTrait pairs:", nrow(pairs_tab),
    "\n  significant (P < 0.05) across the pooled means, Fig. S1 :", sum(pairs_tab$P_pooled < 0.05),
    "\n  significant within temperature, Fig. S2                 :", sum(pairs_tab$P_within_temperature < 0.05),
    "\n  significant in both                                     :",
    sum(pairs_tab$P_pooled < 0.05 & pairs_tab$P_within_temperature < 0.05), "\n")

write.csv(as.data.frame(cor_within) %>% rownames_to_column("Trait"),
          "results/09_spearman_within_temperature_r.csv", row.names = FALSE)
write.csv(as.data.frame(p_within) %>% rownames_to_column("Trait"),
          "results/09_spearman_within_temperature_P.csv", row.names = FALSE)
write.csv(pairs_tab, "results/09_spearman_within_temperature_pairs.csv", row.names = FALSE)

cat("\nDone.\n")

###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 07 - Polyamines and DAP: two-way ANOVA, Tukey letters, faceted bars (Fig. 7)  ###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Free polyamine and 1,3-diaminopropane contents, ug per g fresh weight, of six Mv winter
## wheat cultivars at 20 (control), 35 and 40 degrees C: spermine (SPM), spermidine (SPD),
## cadaverine (CAD), putrescine (PUT) and 1,3-diaminopropane (DAP). One cadaverine value
## (Mv Ikva, 20 C) is missing and is carried as NA (n = 4 for that mean).
##
## Model : content ~ Genotype * Temperature (fixed effects), Type II ANOVA (car::Anova),
##         fitted separately for each polyamine.
## Scale : chosen per polyamine from Shapiro-Wilk (residuals) and Levene's test, which are
##         reported for the raw scale and for the scale used:
##           SPM  square root  (raw Shapiro P = 0.007)
##           SPD  square root  (raw passes, P = 0.14 / 0.87; the square root gives near-
##                              perfect residual diagnostics, P = 0.90 / 0.85)
##           CAD  square root  (raw Shapiro P = 0.038)
##           PUT  square root  (raw Shapiro P < 0.001; no transform normalises the residuals
##                              because three control plants have unusually high putrescine;
##                              variances are homogeneous on the square-root scale, Levene
##                              P = 0.45, and the 40 C differences are large)
##           DAP  log          (raw Shapiro P < 0.001)
## Letters: cultivars compared within each temperature on estimated marginal means
##         (emmeans), Tukey-adjusted, compact letter display with "a" = highest mean.
## Plot  : bars = means + SE in original units; facet rows = temperature, facet columns =
##         polyamine, letters beside the bars.
##
## Input : data/Polyamines/PA_polyamines.csv
##         (Genotype, Temperature, DAP_ug_gFW, PUT_ug_gFW, CAD_ug_gFW, SPD_ug_gFW, SPM_ug_gFW)
## Output: figures/Fig7_polyamines.tiff
##         results/07_polyamines_anova.csv, results/07_polyamines_assumptions.csv,
##         results/07_polyamines_tukey_letters.csv, results/07_polyamines_means.csv

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/07_polyamines.R")
main_dir <- here::here()
setwd(main_dir)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(ggh4x)
  library(car)
  library(emmeans)
  library(multcomp)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

geno_order  <- c("TAR", "KOL", "LUC", "PIR", "IKV", "DAN")
temp_order  <- c("20C", "35C", "40C")
temp_labels <- c("20°C", "35°C", "40°C")
polyamines  <- c("SPM", "SPD", "CAD", "PUT", "DAP")
poly_colors <- c(SPM = "#264653", SPD = "#5dade2", CAD = "#e76f51",
                 PUT = "#f4a261", DAP = "#2a9d8f")

## Transformation used for the model of each polyamine (see DESCRIPTION)
trans_fun <- list(SPM = sqrt, SPD = sqrt, CAD = sqrt, PUT = sqrt, DAP = log)
trans_lab <- c(SPM = "sqrt", SPD = "sqrt", CAD = "sqrt", PUT = "sqrt", DAP = "log")

######### DATA ################################################################################

df <- read.csv("data/Polyamines/PA_polyamines.csv") %>%
  mutate(Genotype    = factor(Genotype,    levels = geno_order),
         Temperature = factor(Temperature, levels = temp_order, labels = temp_labels))

cat("Replicates per genotype x temperature:\n")
print(table(df$Genotype, df$Temperature))

df_long <- df %>%
  pivot_longer(cols = ends_with("_ug_gFW"), names_to = "Polyamine", values_to = "Conc") %>%
  mutate(Polyamine = factor(sub("_ug_gFW", "", Polyamine), levels = polyamines))

######### MEANS AND SE (original units) ########################################################

summary_df <- df_long %>%
  group_by(Polyamine, Genotype, Temperature) %>%
  summarise(n = sum(!is.na(Conc)), mean = mean(Conc, na.rm = TRUE),
            sd = sd(Conc, na.rm = TRUE), se = sd / sqrt(n), .groups = "drop") %>%
  group_by(Polyamine, Genotype) %>%
  mutate(pct_change_vs_20C = 100 * (mean - mean[Temperature == "20°C"]) /
           mean[Temperature == "20°C"]) %>%
  ungroup()

######### MODEL, ASSUMPTIONS, ANOVA AND TUKEY LETTERS PER POLYAMINE ###########################

anova_out <- list(); assump_out <- list()

get_letters <- function(poly) {

  cat("\n==================== ", poly, " ====================\n", sep = "")
  d   <- df_long %>% filter(Polyamine == poly, !is.na(Conc))
  d$y_t <- trans_fun[[poly]](d$Conc)

  ## --- models on the raw scale and on the scale used ---------------------------------
  mod_raw <- lm(Conc ~ Genotype * Temperature, data = d)
  mod     <- lm(y_t  ~ Genotype * Temperature, data = d)
  lev_raw <- leveneTest(Conc ~ Genotype * Temperature, data = d)
  lev_t   <- leveneTest(y_t  ~ Genotype * Temperature, data = d)
  assumptions <- data.frame(
    trait     = poly,
    scale     = c("raw", trans_lab[[poly]]),
    shapiro_W = c(shapiro.test(residuals(mod_raw))$statistic,
                  shapiro.test(residuals(mod))$statistic),
    shapiro_P = c(shapiro.test(residuals(mod_raw))$p.value,
                  shapiro.test(residuals(mod))$p.value),
    levene_F  = c(lev_raw$`F value`[1], lev_t$`F value`[1]),
    levene_P  = c(lev_raw$`Pr(>F)`[1],  lev_t$`Pr(>F)`[1]))
  cat("\nAssumption checks (raw scale vs scale used):\n"); print(assumptions, digits = 3)
  assump_out[[poly]] <<- assumptions

  ## --- Type II ANOVA on the scale used -----------------------------------------------
  aov_tab <- Anova(mod, type = 2)
  cat("\nType II ANOVA, scale = ", trans_lab[[poly]], ":\n", sep = ""); print(aov_tab)
  anova_out[[poly]] <<- cbind(trait = poly, scale = trans_lab[[poly]],
                              term = rownames(aov_tab), as.data.frame(aov_tab))

  ## --- Tukey letters, cultivars within each temperature ------------------------------
  ## The pairwise comparisons behind the letters are Tukey-adjusted. emmeans prints a note
  ## that "tukey" was changed to "sidak": that applies only to the confidence intervals of
  ## the means shown in the table, not to the comparisons.
  cl <- as.data.frame(cld(emmeans(mod, ~ Genotype | Temperature),
                          adjust = "tukey", Letters = letters, decreasing = TRUE)) %>%
    transmute(Polyamine = factor(poly, levels = polyamines), scale = trans_lab[[poly]],
              Genotype, Temperature, emmean_model_scale = emmean, letter = trimws(.group))
  cat("\nTukey letters (a = highest mean within a temperature):\n")
  print(cl %>% arrange(Temperature, desc(emmean_model_scale)))
  cl
}

letters_df <- map_dfr(polyamines, get_letters)

plot_df <- summary_df %>%
  left_join(letters_df %>% dplyr::select(Polyamine, Genotype, Temperature, letter),
            by = c("Polyamine", "Genotype", "Temperature"))

cat("\nMeans (ug/g FW) and % change relative to own 20 C control:\n")
print(as.data.frame(plot_df %>% dplyr::select(Polyamine, Genotype, Temperature, n, mean, se,
                                              pct_change_vs_20C, letter)), digits = 3)

######### FACETED BAR PLOT ####################################################################

p <- ggplot(plot_df, aes(x = mean, y = Genotype, fill = Polyamine)) +
  geom_col() +
  geom_errorbar(aes(xmin = pmax(mean, 0), xmax = mean + se),
                width = 0.25, linewidth = 0.4, color = "black") +
  geom_text(aes(x = mean + se, label = letter),
            hjust = -0.5, size = 3.2, color = "black", fontface = "bold") +
  facet_grid2(Temperature ~ Polyamine, scales = "free_x",
              strip = strip_themed(
                background_x = elem_list_rect(fill = poly_colors),
                text_x       = elem_list_text(colour = "white", face = "bold"))) +
  scale_fill_manual(values = poly_colors, guide = "none") +
  scale_y_discrete(limits = rev) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
  labs(x = expression("Concentration ("*mu*"g g"^{-1}*" FW)"), y = NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major   = element_blank(),
        panel.grid.minor   = element_blank(),
        strip.background.y = element_rect(fill = "#ececec", colour = NA),
        strip.text         = element_text(face = "bold"),
        axis.title.x       = element_text(margin = margin(t = 8)),
        plot.margin        = margin(8, 14, 8, 8))

######### OUTPUT ##############################################################################

write.csv(bind_rows(anova_out),  "results/07_polyamines_anova.csv",         row.names = FALSE)
write.csv(bind_rows(assump_out), "results/07_polyamines_assumptions.csv",   row.names = FALSE)
write.csv(letters_df,            "results/07_polyamines_tukey_letters.csv", row.names = FALSE)
write.csv(plot_df,               "results/07_polyamines_means.csv",         row.names = FALSE)

ggsave("figures/Fig7_polyamines.tiff", p, width = 12, height = 7.1,
       dpi = 300, bg = "white", compression = "lzw")
cat("\nFigure written to figures/Fig7_polyamines.tiff\n")

###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 03 - Photosynthetic pigments: two-way ANOVA, Tukey letters, line plots (Fig. 3)###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Photosynthetic pigments of six Mv winter wheat cultivars at 20 (control), 35 and 40
## degrees C: total chlorophyll (Chl a+b), carotenoids and the chlorophyll a/b ratio at all
## three temperatures, anthocyanin at 20 and 40 degrees C only.
##
## Model : trait ~ Genotype * Temperature (fixed effects), Type II ANOVA (car::Anova).
## Scale : the four pigment traits are analysed on the log scale as a set so that the
##         panels of the figure are compared on one common scale. Shapiro-Wilk and Levene's
##         test are reported for the raw and the log scale. Chl a+b, carotenoid and the
##         a/b ratio pass on both scales; anthocyanin fails both on the raw scale (P = 0.03),
##         and on the log scale has homogeneous variances (Levene P = 0.29) with a departure
##         from normality (Shapiro P = 0.001) that does not affect the outcome (no cultivar
##         differences at either temperature under any scale of the letters shown).
## Letters: cultivars compared within each temperature on estimated marginal means
##         (emmeans), Tukey-adjusted, compact letter display with "a" = highest mean.
## Plot  : genotype x temperature means in original units joined by lines, one line per
##         cultivar, Tukey letters printed inside the markers.
##
## Input : data/Pigments/Chlorophyll_ab.csv        (Genotype, Temperature, Chl_ab_ug_gFW)
##         data/Pigments/Carotenoid.csv            (Genotype, Temperature, Carotenoid_ug_gFW)
##         data/Pigments/Anthocyanin.csv           (Genotype, Temperature, Anthocyanin_nmol_gFW)
##         data/Pigments/Chlorophyll_ab_ratio.csv  (Genotype, Temperature, Chl_a_b_ratio)
## Output: figures/Fig3_pigments.tiff
##         results/03_pigments_anova.csv, results/03_pigments_assumptions.csv,
##         results/03_pigments_tukey_letters.csv, results/03_pigments_means.csv

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/03_pigments.R")
main_dir <- here::here()
setwd(main_dir)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(car)
  library(emmeans)
  library(multcomp)
  library(patchwork)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

genotype_order <- c("LUC", "KOL", "DAN", "IKV", "PIR", "TAR")
temp_levels    <- c("20C", "35C", "40C")
temp_labels    <- c("20 °C", "35 °C", "40 °C")

## Okabe-Ito colours and one marker shape per cultivar
pal_okabe <- c(LUC = "#CC79A7", KOL = "#D55E00", DAN = "#56B4E9",
               IKV = "#009E73", PIR = "#0072B2", TAR = "#E69F00")
shp <- c(LUC = 23,   # filled diamond
         KOL = 24,   # filled triangle up
         DAN = 21,   # filled circle
         IKV = 22,   # filled square
         PIR = 25,   # filled triangle down
         TAR = 8)    # asterisk

my_theme <- theme_gray(base_size = 12) +
  theme(panel.grid.minor  = element_blank(),
        panel.grid.major  = element_blank(),
        axis.line         = element_line(color = "grey30", linewidth = 0.4),
        axis.ticks        = element_line(color = "grey30", linewidth = 0.4),
        axis.text         = element_text(color = "grey20", size = 10),
        axis.title        = element_text(color = "grey10", size = 11),
        legend.position   = "right",
        legend.title      = element_text(face = "bold", size = 10),
        legend.text       = element_text(size = 9),
        legend.key.height = unit(0.9, "lines"),
        plot.margin       = margin(8, 12, 8, 8))

######### DATA ################################################################################

chl  <- read.csv("data/Pigments/Chlorophyll_ab.csv")
car  <- read.csv("data/Pigments/Carotenoid.csv")
anth <- read.csv("data/Pigments/Anthocyanin.csv")
ab   <- read.csv("data/Pigments/Chlorophyll_ab_ratio.csv")

for (nm in c("chl", "car", "anth", "ab")) {
  d <- get(nm)
  d <- d %>% mutate(Genotype    = factor(Genotype,    levels = genotype_order),
                    Temperature = factor(Temperature, levels = temp_levels))
  assign(nm, d)
  cat("\nReplicates per genotype x temperature,", nm, ":\n"); print(table(d$Genotype, d$Temperature))
}

######### ONE PANEL: MODEL, ASSUMPTIONS, ANOVA, TUKEY LETTERS AND LINE PLOT ###################

## Results of every panel are collected in these lists and written out at the end
anova_out <- list(); assump_out <- list(); letters_out <- list(); means_out <- list()

make_line_plot <- function(data, value_col, trait, ylab, y_min, y_max, show_35 = TRUE) {

  cat("\n==================== ", trait, " ====================\n", sep = "")
  data$y <- data[[value_col]]
  data   <- data %>% filter(!is.na(y)) %>% droplevels()

  ## --- models on the raw and the log scale ---------------------------------------------
  mod_raw <- lm(y ~ Genotype * Temperature, data = data)
  mod_log <- lm(log(y) ~ Genotype * Temperature, data = data)
  lev_raw <- leveneTest(y ~ Genotype * Temperature, data = data)
  lev_log <- leveneTest(log(y) ~ Genotype * Temperature, data = data)

  assumptions <- data.frame(
    trait     = trait,
    scale     = c("raw", "log"),
    shapiro_W = c(shapiro.test(residuals(mod_raw))$statistic,
                  shapiro.test(residuals(mod_log))$statistic),
    shapiro_P = c(shapiro.test(residuals(mod_raw))$p.value,
                  shapiro.test(residuals(mod_log))$p.value),
    levene_F  = c(lev_raw$`F value`[1], lev_log$`F value`[1]),
    levene_P  = c(lev_raw$`Pr(>F)`[1],  lev_log$`Pr(>F)`[1]))
  cat("\nAssumption checks (raw vs log scale):\n"); print(assumptions, digits = 3)
  assump_out[[trait]] <<- assumptions

  ## --- Type II ANOVA on the log scale --------------------------------------------------
  mod     <- mod_log
  aov_tab <- Anova(mod, type = 2)
  cat("\nType II ANOVA, log(", trait, "):\n", sep = ""); print(aov_tab)
  anova_out[[trait]] <<- cbind(trait = trait, term = rownames(aov_tab), as.data.frame(aov_tab))

  ## --- Tukey letters, cultivars within each temperature --------------------------------
  emm <- emmeans(mod, ~ Genotype | Temperature)
  ## The pairwise comparisons behind the letters are Tukey-adjusted. emmeans prints a note
  ## that "tukey" was changed to "sidak": that applies only to the confidence intervals of
  ## the means shown in the table, not to the comparisons.
  cl <- as.data.frame(cld(emm, adjust = "tukey", Letters = letters, decreasing = TRUE)) %>%
    mutate(letter = trimws(.group)) %>%
    dplyr::select(Genotype, Temperature, emmean_log = emmean, letter)
  cat("\nTukey letters (a = highest mean within a temperature):\n")
  print(cl %>% arrange(Temperature, desc(emmean_log)))
  letters_out[[trait]] <<- cbind(trait = trait, cl)

  ## --- means, SE and % change relative to the 20 C control -----------------------------
  means_df <- data %>%
    group_by(Genotype, Temperature) %>%
    summarise(n = n(), mean = mean(y), sd = sd(y), se = sd / sqrt(n), .groups = "drop")
  ctrl   <- means_df %>% filter(Temperature == "20C") %>% dplyr::select(Genotype, ctrl = mean)
  means_df <- means_df %>% left_join(ctrl, by = "Genotype") %>%
    mutate(pct_change_vs_20C = 100 * (mean - ctrl) / ctrl) %>% dplyr::select(-ctrl)
  cat("\nMeans and % change relative to own 20 C control:\n")
  print(as.data.frame(means_df), digits = 3)
  means_out[[trait]] <<- cbind(trait = trait, means_df)

  ## --- plot ----------------------------------------------------------------------------
  plot_df <- means_df %>%
    left_join(cl %>% dplyr::select(Genotype, Temperature, letter), by = c("Genotype", "Temperature")) %>%
    mutate(x_pos = as.integer(factor(Temperature, levels = temp_levels)))

  x_breaks <- if (show_35) seq_along(temp_levels) else c(1, 3)
  x_labels <- if (show_35) temp_labels            else temp_labels[c(1, 3)]

  ggplot(plot_df, aes(x = x_pos, y = mean, color = Genotype, group = Genotype)) +
    geom_line(linewidth = 1.1, lineend = "round") +
    geom_point(aes(shape = Genotype), size = 4.5, stroke = 1.3, fill = "white") +
    geom_text(aes(label = letter), size = 2.8, color = "grey15", fontface = "bold") +
    scale_color_manual(values = pal_okabe) +
    scale_shape_manual(values = shp) +
    scale_x_continuous(breaks = x_breaks, labels = x_labels, limits = c(1, 3),
                       expand = expansion(mult = c(0.08, 0.08))) +
    coord_cartesian(ylim = c(y_min, y_max)) +
    labs(x = "Temperature", y = ylab) +
    my_theme
}

######### THE FOUR PANELS #####################################################################

p_chl  <- make_line_plot(chl,  "Chl_ab_ug_gFW",        "Chl_ab",
                         expression("Chl a+b ("*mu*"g g"^{-1}*" FW)"),      1500, 3500)
p_car  <- make_line_plot(car,  "Carotenoid_ug_gFW",    "Carotenoid",
                         expression("Carotenoid ("*mu*"g g"^{-1}*" FW)"),    430,  660)
p_anth <- make_line_plot(anth, "Anthocyanin_nmol_gFW", "Anthocyanin",
                         expression("Anthocyanin (nmol g"^{-1}*" FW)"),        0,   16, show_35 = FALSE)
p_ab   <- make_line_plot(ab,   "Chl_a_b_ratio",        "Chl_a_b_ratio",
                         "Chl a / b ratio",                                   1.8,  3.3)

combined <- (p_chl + p_car) / (p_anth + p_ab) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

######### OUTPUT ##############################################################################

write.csv(bind_rows(anova_out),   "results/03_pigments_anova.csv",         row.names = FALSE)
write.csv(bind_rows(assump_out),  "results/03_pigments_assumptions.csv",   row.names = FALSE)
write.csv(bind_rows(letters_out), "results/03_pigments_tukey_letters.csv", row.names = FALSE)
write.csv(bind_rows(means_out),   "results/03_pigments_means.csv",         row.names = FALSE)

ggsave("figures/Fig3_pigments.tiff", combined, width = 10, height = 7.7,
       dpi = 300, bg = "white", compression = "lzw")
cat("\nFigure written to figures/Fig3_pigments.tiff\n")

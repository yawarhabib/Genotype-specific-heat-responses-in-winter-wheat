###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 02 - Leaf gas exchange: two-way ANOVA, Tukey letters, 4-panel radar (Fig. 2)  ###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Leaf gas exchange of six Mv winter wheat cultivars at 20 (control), 35 and 40 degrees C,
## n = 3 plants per genotype x temperature: net photosynthetic rate (Pn), stomatal
## conductance (gs), transpiration rate (E) and water-use efficiency (WUE = Pn / E).
## Vapour-pressure deficit (VPD) is in the data file but is not analysed here (it is
## collinear with temperature).
##
## Model : trait ~ Genotype * Temperature (fixed effects), Type II ANOVA (car::Anova).
## Scale : all four traits are analysed on the log scale so that the four panels of the
##         figure are compared on one common scale. Shapiro-Wilk and Levene's test are
##         reported for the raw and the log scale; Pn, gs and E pass on both, WUE shows a
##         mild departure from normality on both scales (Shapiro P = 0.03 raw, 0.02 log)
##         with homogeneous variances (Levene P > 0.85), which the F-test tolerates.
## Letters: cultivars compared within each temperature on estimated marginal means
##         (emmeans), Tukey-adjusted, compact letter display with "a" = highest mean.
## Plot  : one radar per trait (A Pn, B gs, C E, D WUE); genotype x temperature means
##         in original units, scaled to the outer ring of each panel, one polygon per
##         temperature, Tukey letters at each vertex.
##
## Input : data/Gas_exchange/Leaf_gas_Exchange.csv
##         (Genotype = full name, Genotype_code, Temperature, Pn, E, gs, VPD, WUE)
## Output: figures/Fig2_leaf_gas_exchange.tiff
##         results/02_gas_exchange_anova.csv, results/02_gas_exchange_assumptions.csv,
##         results/02_gas_exchange_tukey_letters.csv, results/02_gas_exchange_means.csv

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/02_gas_exchange.R")
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

geno_order  <- c("TAR", "KOL", "LUC", "PIR", "IKV", "DAN")
temp_order  <- c("20", "35", "40")
temp_labels <- c("20°C", "35°C", "40°C")
pal         <- c("20" = "#3B8BC2", "35" = "#F2C14E", "40" = "#D1495B")

######### DATA ################################################################################

gas <- read.csv("data/Gas_exchange/Leaf_gas_Exchange.csv", header = TRUE, encoding = "UTF-8")

## Full cultivar names are carried in the data file (column Genotype); the code column
## (Genotype_code) fixes the plotting order
geno_order_full <- gas %>% distinct(Genotype_code, Genotype) %>%
  { .$Genotype[match(geno_order, .$Genotype_code)] }

gas <- gas %>%
  mutate(Genotype    = factor(Genotype,    levels = geno_order_full),
         Temperature = factor(Temperature, levels = temp_order))

cat("Replicates per genotype x temperature:\n")
print(table(gas$Genotype, gas$Temperature))

traits <- c("Pn", "gs", "E", "WUE")

######### RADAR GEOMETRY SHARED BY THE FOUR PANELS ############################################

n_axes  <- length(geno_order_full)
axis_df <- tibble(
  Genotype = factor(geno_order_full, levels = geno_order_full),
  axis_id  = seq_len(n_axes),
  angle    = pi / 2 + 2 * pi * (seq_len(n_axes) - 1) / n_axes)

## Ring value labels along the diagonal between the Tarsoly and Kolompos spokes
ring_angle <- pi / 2 + pi / 6

######### ONE PANEL: MODEL, ASSUMPTIONS, ANOVA, TUKEY LETTERS AND RADAR #######################

## Results of every panel are collected in these lists and written out at the end
anova_out <- list(); assump_out <- list(); letters_out <- list(); means_out <- list()

make_radar <- function(var, title_expr, show_legend = FALSE) {

  cat("\n==================== ", var, " ====================\n", sep = "")

  ## --- models on the raw and the log scale ---------------------------------------------
  f_raw   <- as.formula(paste(var, "~ Genotype * Temperature"))
  f_log   <- as.formula(paste("log(", var, ") ~ Genotype * Temperature"))
  mod_raw <- lm(f_raw, data = gas)
  mod_log <- lm(f_log, data = gas)

  lev_raw <- leveneTest(f_raw, data = gas)
  lev_log <- leveneTest(f_log, data = gas)

  assumptions <- data.frame(
    trait     = var,
    scale     = c("raw", "log"),
    shapiro_W = c(shapiro.test(residuals(mod_raw))$statistic,
                  shapiro.test(residuals(mod_log))$statistic),
    shapiro_P = c(shapiro.test(residuals(mod_raw))$p.value,
                  shapiro.test(residuals(mod_log))$p.value),
    levene_F  = c(lev_raw$`F value`[1], lev_log$`F value`[1]),
    levene_P  = c(lev_raw$`Pr(>F)`[1],  lev_log$`Pr(>F)`[1]))
  cat("\nAssumption checks (raw vs log scale):\n"); print(assumptions, digits = 3)
  assump_out[[var]] <<- assumptions

  ## --- Type II ANOVA on the log scale --------------------------------------------------
  mod     <- mod_log
  aov_tab <- Anova(mod, type = 2)
  cat("\nType II ANOVA, log(", var, "):\n", sep = ""); print(aov_tab)
  anova_out[[var]] <<- cbind(trait = var, term = rownames(aov_tab), as.data.frame(aov_tab))

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
  letters_out[[var]] <<- cbind(trait = var, cl)

  ## --- means, SE and % change relative to the 20 C control -----------------------------
  means_df <- gas %>%
    group_by(Genotype, Temperature) %>%
    summarise(n = n(), mean = mean(.data[[var]]), sd = sd(.data[[var]]),
              se = sd / sqrt(n), .groups = "drop")
  pct_df <- means_df %>%
    dplyr::select(Genotype, Temperature, mean) %>%
    pivot_wider(names_from = Temperature, values_from = mean) %>%
    mutate(pct_change_35C = 100 * (`35` - `20`) / `20`,
           pct_change_40C = 100 * (`40` - `20`) / `20`)
  cat("\nMeans and % change relative to own 20 C control:\n")
  print(as.data.frame(pct_df), digits = 3)
  means_out[[var]] <<- cbind(trait = var,
                             means_df %>% left_join(pct_df %>% dplyr::select(Genotype, pct_change_35C,
                                                                             pct_change_40C),
                                                    by = "Genotype"))

  ## --- radar geometry: means scaled to the outer ring of this panel --------------------
  max_val     <- ceiling(max(means_df$mean) * 10) / 10
  ring_levels <- pretty(c(0, max_val), n = 4)
  ring_levels <- ring_levels[ring_levels > 0 & ring_levels <= max_val]
  max_val     <- max(ring_levels)

  plot_df <- means_df %>%
    left_join(axis_df, by = "Genotype") %>%
    mutate(r = mean / max_val, x = r * cos(angle), y = r * sin(angle))

  plot_df_closed <- plot_df %>%
    group_by(Temperature) %>% arrange(axis_id) %>%
    group_modify(~ bind_rows(.x, .x[1, ])) %>% ungroup()

  ring_norm <- ring_levels / max_val
  ring_df <- expand_grid(ring = ring_norm, axis_id = seq_len(n_axes)) %>%
    left_join(axis_df %>% dplyr::select(axis_id, angle), by = "axis_id") %>%
    mutate(x = ring * cos(angle), y = ring * sin(angle)) %>%
    group_by(ring) %>% group_modify(~ bind_rows(.x, .x[1, ])) %>% ungroup()

  spoke_df <- axis_df %>% mutate(x_end = cos(angle), y_end = sin(angle))
  label_df <- axis_df %>% mutate(x = 1.22 * cos(angle), y = 1.22 * sin(angle))
  axval_df <- tibble(val = ring_levels,
                     x = ring_norm * cos(ring_angle), y = ring_norm * sin(ring_angle))

  letter_df <- plot_df %>%
    left_join(cl %>% dplyr::select(Genotype, Temperature, letter),
              by = c("Genotype", "Temperature")) %>%
    mutate(r_letter = r + 0.03,
           x_lab = r_letter * cos(angle), y_lab = r_letter * sin(angle))

  ## --- plot ----------------------------------------------------------------------------
  ggplot() +
    geom_segment(data = spoke_df, aes(x = 0, y = 0, xend = x_end, yend = y_end),
                 colour = "grey80", linewidth = 0.3) +
    geom_path(data = ring_df, aes(x = x, y = y, group = ring),
              colour = "grey80", linewidth = 0.3) +
    geom_path(data = ring_df %>% filter(ring == max(ring)), aes(x = x, y = y, group = ring),
              colour = "grey50", linewidth = 0.9) +
    geom_polygon(data = plot_df_closed,
                 aes(x = x, y = y, group = Temperature, colour = Temperature, fill = Temperature),
                 alpha = 0.30, linewidth = 0.9) +
    geom_point(data = plot_df, aes(x = x, y = y, colour = Temperature), size = 1.8) +
    ## letters justified away from the centre so that neighbouring temperatures do not overlap
    geom_text(data = letter_df, aes(x = x_lab, y = y_lab, label = letter, colour = Temperature),
              fontface = "bold", size = 3, hjust = "outward", vjust = "outward",
              show.legend = FALSE) +
    geom_text(data = label_df, aes(x = x, y = y, label = Genotype), size = 3.3) +
    geom_label(data = axval_df, aes(x = x, y = y, label = sprintf("%.2g", val)),
               size = 2.8, colour = "grey40", fill = "white", linewidth = 0,
               label.padding = unit(0.1, "lines")) +
    scale_colour_manual(values = pal, labels = temp_labels, name = NULL) +
    scale_fill_manual(values = pal, labels = temp_labels, name = NULL) +
    coord_fixed(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5)) +
    labs(title = title_expr) +
    guides(colour = "none",
           fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.5, alpha = 0.6),
                               keywidth = unit(1.2, "cm"), keyheight = unit(0.5, "cm"))) +
    theme_void(base_size = 11) +
    theme(plot.title       = element_text(hjust = 0.5, face = "bold", size = 12,
                                          margin = margin(b = 6)),
          legend.position  = if (show_legend) "bottom" else "none",
          legend.direction = "horizontal",
          legend.key       = element_blank(),
          legend.text      = element_text(size = 10),
          legend.spacing.x = unit(0.4, "cm"),
          plot.margin      = margin(5, 5, 5, 5))
}

######### THE FOUR PANELS #####################################################################

p_Pn  <- make_radar("Pn",  expression(italic(P)[n]~(mu*mol~m^{-2}~s^{-1})))
p_gs  <- make_radar("gs",  expression(italic(g)[s]~(mmol~m^{-2}~s^{-1})))
p_E   <- make_radar("E",   expression(italic(E)~(mmol~m^{-2}~s^{-1})))
p_WUE <- make_radar("WUE", expression(WUE~(mu*mol~mmol^{-1})), show_legend = TRUE)

final_figure <- (p_Pn | p_gs) / (p_E | p_WUE) +
  plot_annotation(tag_levels = "A")

######### OUTPUT ##############################################################################

write.csv(bind_rows(anova_out),   "results/02_gas_exchange_anova.csv",         row.names = FALSE)
write.csv(bind_rows(assump_out),  "results/02_gas_exchange_assumptions.csv",   row.names = FALSE)
write.csv(bind_rows(letters_out), "results/02_gas_exchange_tukey_letters.csv", row.names = FALSE)
write.csv(bind_rows(means_out),   "results/02_gas_exchange_means.csv",         row.names = FALSE)

ggsave("figures/Fig2_leaf_gas_exchange.tiff", final_figure, width = 9.5, height = 9.4,
       dpi = 300, bg = "white", compression = "lzw")
cat("\nFigure written to figures/Fig2_leaf_gas_exchange.tiff\n")

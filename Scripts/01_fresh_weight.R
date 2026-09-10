###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 01 - Shoot fresh weight: two-way ANOVA, Tukey letters and radar plot (Fig. 1) ###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Shoot fresh weight (g per plant) of six Mv winter wheat cultivars grown at 20 (control),
## 35 and 40 degrees C, n = 7 plants per genotype x temperature.
##
## Model : FW ~ Genotype * Temperature (fixed effects), Type II ANOVA (car::Anova).
## Checks: Shapiro-Wilk on residuals and Levene's test; both fail on the raw scale
##         (P = 0.002 and 0.031) and pass after log transformation (P = 0.34 and 0.85),
##         so the model is fitted to log(FW).
## Letters: cultivars compared within each temperature on estimated marginal means
##         (emmeans), Tukey-adjusted, compact letter display with "a" = highest mean.
## Plot  : radar of the raw (untransformed) genotype x temperature means in g per plant,
##         one polygon per temperature, Tukey letters at each vertex.
##
## Input : data/Fresh_weight/FW.csv (Genotype = code, Cultivar = full name, Temperature, FW_g)
## Output: figures/Fig1_shoot_fresh_weight.tiff
##         results/01_FW_anova.csv, results/01_FW_assumptions.csv,
##         results/01_FW_tukey_letters.csv, results/01_FW_means.csv

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/01_fresh_weight.R")
main_dir <- here::here()
setwd(main_dir)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(car)
  library(emmeans)
  library(multcomp)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

geno_order  <- c("TAR", "KOL", "LUC", "PIR", "IKV", "DAN")
temp_order  <- c("20C", "35C", "40C")
temp_labels <- c("20°C", "35°C", "40°C")
pal         <- c("20C" = "#3B8BC2", "35C" = "#F2C14E", "40C" = "#D1495B")

######### DATA ################################################################################

FW <- read.csv("data/Fresh_weight/FW.csv", header = TRUE, encoding = "UTF-8")

## Cultivar names are carried in the data file (column Cultivar); the code column
## (Genotype) fixes the plotting order
geno_order_full <- FW %>% distinct(Genotype, Cultivar) %>%
  { .$Cultivar[match(geno_order, .$Genotype)] }

FW <- FW %>%
  mutate(Genotype    = factor(Cultivar,    levels = geno_order_full),
         Temperature = factor(Temperature, levels = temp_order)) %>%
  dplyr::select(Genotype, Temperature, FW_g)

cat("Replicates per genotype x temperature:\n")
print(table(FW$Genotype, FW$Temperature))

######### MODEL AND ASSUMPTION CHECKS #########################################################

mod_raw <- lm(FW_g ~ Genotype * Temperature, data = FW)
mod_log <- lm(log(FW_g) ~ Genotype * Temperature, data = FW)

lev_raw <- leveneTest(FW_g ~ Genotype * Temperature, data = FW)
lev_log <- leveneTest(log(FW_g) ~ Genotype * Temperature, data = FW)

assumptions <- data.frame(
  scale     = c("raw", "log"),
  shapiro_W = c(shapiro.test(residuals(mod_raw))$statistic,
                shapiro.test(residuals(mod_log))$statistic),
  shapiro_P = c(shapiro.test(residuals(mod_raw))$p.value,
                shapiro.test(residuals(mod_log))$p.value),
  levene_F  = c(lev_raw$`F value`[1], lev_log$`F value`[1]),
  levene_P  = c(lev_raw$`Pr(>F)`[1],  lev_log$`Pr(>F)`[1]))
cat("\nAssumption checks (raw vs log scale):\n"); print(assumptions, digits = 3)
write.csv(assumptions, "results/01_FW_assumptions.csv", row.names = FALSE)

## Raw scale fails both tests, log scale passes both -> analyse log(FW)
mod     <- mod_log
aov_tab <- Anova(mod, type = 2)
cat("\nType II ANOVA, log(FW):\n"); print(aov_tab)
write.csv(cbind(term = rownames(aov_tab), as.data.frame(aov_tab)),
          "results/01_FW_anova.csv", row.names = FALSE)

######### TUKEY LETTERS (cultivars within each temperature) ###################################

emm     <- emmeans(mod, ~ Genotype | Temperature)
## The pairwise comparisons behind the letters are Tukey-adjusted. emmeans prints a note
## that "tukey" was changed to "sidak": that applies only to the confidence intervals of
## the means shown in the table, not to the comparisons.
cld_tab <- cld(emm, adjust = "tukey", Letters = letters, decreasing = TRUE)
letters_df <- as.data.frame(cld_tab) %>%
  mutate(letter = trimws(.group)) %>%
  dplyr::select(Genotype, Temperature, emmean_log = emmean, letter)
cat("\nTukey letters (a = highest mean within a temperature):\n")
print(letters_df %>% arrange(Temperature, desc(emmean_log)))
write.csv(letters_df, "results/01_FW_tukey_letters.csv", row.names = FALSE)

######### MEANS AND PERCENT CHANGE ############################################################

means_df <- FW %>%
  group_by(Genotype, Temperature) %>%
  summarise(n = n(), mean = mean(FW_g), sd = sd(FW_g), se = sd / sqrt(n), .groups = "drop")

pct_df <- means_df %>%
  dplyr::select(Genotype, Temperature, mean) %>%
  pivot_wider(names_from = Temperature, values_from = mean) %>%
  mutate(pct_change_35C = 100 * (`35C` - `20C`) / `20C`,
         pct_change_40C = 100 * (`40C` - `20C`) / `20C`)
cat("\nMeans (g) and % change relative to own 20 C control:\n")
print(as.data.frame(pct_df), digits = 3)
write.csv(means_df %>%
            left_join(pct_df %>% dplyr::select(Genotype, pct_change_35C, pct_change_40C),
                      by = "Genotype"),
          "results/01_FW_means.csv", row.names = FALSE)

######### RADAR PLOT ##########################################################################

n_axes  <- length(geno_order_full)
max_val <- 1.0                                     # outer ring, g per plant

axis_df <- tibble(
  Genotype = factor(geno_order_full, levels = geno_order_full),
  axis_id  = seq_len(n_axes),
  angle    = pi / 2 + 2 * pi * (seq_len(n_axes) - 1) / n_axes)

plot_df <- means_df %>%
  left_join(axis_df, by = "Genotype") %>%
  mutate(x = mean * cos(angle), y = mean * sin(angle))

plot_df_closed <- plot_df %>%
  group_by(Temperature) %>% arrange(axis_id) %>%
  group_modify(~ bind_rows(.x, .x[1, ])) %>% ungroup()

ring_levels <- seq(0.2, max_val, by = 0.2)
ring_df <- expand_grid(ring = ring_levels, axis_id = seq_len(n_axes)) %>%
  left_join(axis_df %>% dplyr::select(axis_id, angle), by = "axis_id") %>%
  mutate(x = ring * cos(angle), y = ring * sin(angle)) %>%
  group_by(ring) %>% group_modify(~ bind_rows(.x, .x[1, ])) %>% ungroup()

spoke_df <- axis_df %>% mutate(x_end = max_val * cos(angle), y_end = max_val * sin(angle))
label_df <- axis_df %>% mutate(x = 1.18 * max_val * cos(angle), y = 1.18 * max_val * sin(angle))
## Ring value labels along the diagonal between the Tarsoly and Kolompos spokes
ring_angle <- pi / 2 + pi / 6
axval_df <- tibble(val = ring_levels,
                   x = ring_levels * cos(ring_angle), y = ring_levels * sin(ring_angle))

letter_df <- plot_df %>%
  left_join(letters_df %>% dplyr::select(Genotype, Temperature, letter),
            by = c("Genotype", "Temperature")) %>%
  mutate(r_letter = mean + 0.03,
         x_lab = r_letter * cos(angle), y_lab = r_letter * sin(angle))

p_radar <- ggplot() +
  geom_segment(data = spoke_df, aes(x = 0, y = 0, xend = x_end, yend = y_end),
               colour = "grey80", linewidth = 0.3) +
  geom_path(data = ring_df, aes(x = x, y = y, group = ring),
            colour = "grey80", linewidth = 0.3) +
  geom_path(data = ring_df %>% filter(ring == max_val), aes(x = x, y = y, group = ring),
            colour = "grey50", linewidth = 1.0) +
  geom_polygon(data = plot_df_closed,
               aes(x = x, y = y, group = Temperature, colour = Temperature, fill = Temperature),
               alpha = 0.30, linewidth = 1.0) +
  geom_point(data = plot_df, aes(x = x, y = y, colour = Temperature), size = 1.8) +
  ## letters justified away from the centre so that neighbouring temperatures do not overlap
  geom_text(data = letter_df, aes(x = x_lab, y = y_lab, label = letter, colour = Temperature),
            fontface = "bold", size = 3.5, hjust = "outward", vjust = "outward",
            show.legend = FALSE) +
  geom_text(data = label_df, aes(x = x, y = y, label = Genotype), size = 4) +
  geom_label(data = axval_df, aes(x = x, y = y, label = sprintf("%.1f", val)),
             size = 3, colour = "grey40", fill = "white", linewidth = 0,
             label.padding = unit(0.1, "lines")) +
  scale_colour_manual(values = pal, labels = temp_labels, name = NULL) +
  scale_fill_manual(values = pal, labels = temp_labels, name = NULL) +
  coord_fixed(xlim = c(-1.4, 1.4), ylim = c(-1.4, 1.4)) +
  guides(colour = "none",
         fill = guide_legend(override.aes = list(colour = "black", linewidth = 0.5, alpha = 0.6),
                             keywidth = unit(1.2, "cm"), keyheight = unit(0.5, "cm"))) +
  theme_void(base_size = 12) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.key = element_blank(), legend.text = element_text(size = 11),
        legend.spacing.x = unit(0.4, "cm"), plot.margin = margin(10, 10, 10, 10))

ggsave("figures/Fig1_shoot_fresh_weight.tiff", p_radar, width = 9.5, height = 8.8,
       dpi = 300, bg = "white", compression = "lzw")
cat("\nFigure written to figures/Fig1_shoot_fresh_weight.tiff\n")

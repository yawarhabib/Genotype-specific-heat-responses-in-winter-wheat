###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 05 - Reactive-carbonyl enzymes: ANOVA, Tukey letters, violins, heatmap (Fig. 5)###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Activities of the reactive-carbonyl detoxification enzymes in six Mv winter wheat
## cultivars at 20 (control), 35 and 40 degrees C: aldo-keto reductase (AKR),
## alkenal/alkenone oxidoreductase (AOR/AER) and glyoxalase I (GLY I), nkat per g FW.
##
## Model : activity ~ Genotype * Temperature (fixed effects), Type II ANOVA (car::Anova).
## Scale : chosen per enzyme from Shapiro-Wilk (residuals) and Levene's test, which are
##         reported for the raw scale and for the scale used:
##           AKR      untransformed  (raw passes both tests)
##           AOR/AER  square root    (raw Levene P = 0.034)
##           GLY I    square root    (raw passes, P = 0.60 / 0.43; the activity spans a
##                                    six-fold range and the square root stabilises the
##                                    variance further, Levene P = 0.97)
## Letters: cultivars compared within each temperature on estimated marginal means
##         (emmeans), Tukey-adjusted, compact letter display with "a" = highest mean.
## Plot  : A-C violin + inner box plot + individual observations per genotype x
##         temperature with letters above each violin; D heatmap of the genotype x
##         temperature means as within-enzyme z-scores (each enzyme scaled over its 18
##         means), tiles labelled with the mean and the same Tukey letter as in A-C.
##
## Input : data/Detoxification_enzyme/AKR.csv, AOR_AER.csv, GLY.csv
##         (Genotype = code, Cultivar = full name, Temperature, <enzyme>_nkat_gFW)
## Output: figures/Fig5_carbonyl_enzymes.tiff
##         results/05_carbonyl_anova.csv, results/05_carbonyl_assumptions.csv,
##         results/05_carbonyl_tukey_letters.csv, results/05_carbonyl_means.csv

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/05_carbonyl_enzymes.R")
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
temp_order  <- c("20C", "35C", "40C")
temp_labels <- c("20C" = "20°C", "35C" = "35°C", "40C" = "40°C")
pal         <- c("20C" = "#3B8BC2", "35C" = "#F2C14E", "40C" = "#D1495B")
dodge       <- position_dodge(width = 0.8)

## Shared violin theme
violin_theme <- theme_bw(base_size = 12) +
  theme(panel.grid      = element_blank(),
        legend.position = "top",
        axis.title      = element_text(face = "bold"),
        axis.title.y    = element_text(family = "serif", face = "bold"),
        axis.text.x     = element_text(color = "black", size = 10),
        plot.tag        = element_text(size = 14))

######### PANEL A - AKR (untransformed) ########################################################

df_AKR <- read.csv("data/Detoxification_enzyme/AKR.csv", encoding = "UTF-8")
geno_order_full <- df_AKR %>% distinct(Genotype, Cultivar) %>%
  { .$Cultivar[match(geno_order, .$Genotype)] }
df_AKR <- df_AKR %>%
  transmute(Genotype    = factor(Cultivar,    levels = geno_order_full),
            Temperature = factor(Temperature, levels = temp_order),
            y = AKR_nkat_gFW, y_t = AKR_nkat_gFW)
cat("\n==================== AKR ====================\n")
cat("Replicates per genotype x temperature:\n"); print(table(df_AKR$Genotype, df_AKR$Temperature))

mod_AKR <- lm(y_t ~ Genotype * Temperature, data = df_AKR)
lev_AKR <- leveneTest(y_t ~ Genotype * Temperature, data = df_AKR)
assump_AKR <- data.frame(trait = "AKR", scale = "none",
                         shapiro_W = shapiro.test(residuals(mod_AKR))$statistic,
                         shapiro_P = shapiro.test(residuals(mod_AKR))$p.value,
                         levene_F  = lev_AKR$`F value`[1], levene_P = lev_AKR$`Pr(>F)`[1])
cat("\nAssumption checks (scale used = none):\n"); print(assump_AKR, digits = 3)
aov_AKR <- Anova(mod_AKR, type = 2)
cat("\nType II ANOVA:\n"); print(aov_AKR)

## The pairwise comparisons behind the letters are Tukey-adjusted. emmeans prints a note
## that "tukey" was changed to "sidak": that applies only to the confidence intervals of
## the means shown in the table, not to the comparisons.
cl_AKR <- as.data.frame(cld(emmeans(mod_AKR, ~ Genotype | Temperature),
                            adjust = "tukey", Letters = letters, decreasing = TRUE)) %>%
  mutate(letter = trimws(.group)) %>%
  dplyr::select(Genotype, Temperature, emmean_model_scale = emmean, letter)
cat("\nTukey letters (a = highest mean within a temperature):\n")
print(cl_AKR %>% arrange(Temperature, desc(emmean_model_scale)))

letters_AKR <- df_AKR %>%
  group_by(Genotype, Temperature) %>%
  summarise(y_top = max(density(y, adjust = 1.2)$x), .groups = "drop") %>%
  left_join(cl_AKR %>% dplyr::select(Genotype, Temperature, letter),
            by = c("Genotype", "Temperature"))

p_AKR <- ggplot(df_AKR, aes(Genotype, y, fill = Temperature)) +
  geom_vline(xintercept = seq(1.5, 5.5, 1), color = "grey70", linetype = "dashed",
             linewidth = 0.4) +
  geom_violin(position = dodge, alpha = 0.75, color = "grey25",
              trim = FALSE, scale = "width", width = 1.0, adjust = 1.2) +
  geom_boxplot(aes(group = interaction(Genotype, Temperature)),
               position = dodge, width = 0.15, fill = "gray", alpha = 0.70,
               outlier.shape = NA, color = "grey20", linewidth = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
             size = 2, alpha = 0.55, shape = 21, stroke = 0.5, color = "black",
             show.legend = FALSE) +
  geom_text(data = letters_AKR,
            aes(x = Genotype, y = y_top, label = letter, group = Temperature,
                color = Temperature),
            position = dodge, size = 4, vjust = -0.4, fontface = "bold",
            show.legend = FALSE) +
  scale_fill_manual(values = pal, labels = temp_labels) +
  scale_color_manual(values = pal, labels = temp_labels, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "", y = expression(AKR~activity~(nkat~g^{-1}~FW))) +
  violin_theme

######### PANEL B - AOR/AER (square root) #####################################################

df_AOR <- read.csv("data/Detoxification_enzyme/AOR_AER.csv", encoding = "UTF-8") %>%
  transmute(Genotype    = factor(Cultivar,    levels = geno_order_full),
            Temperature = factor(Temperature, levels = temp_order),
            y = AOR_AER_nkat_gFW, y_t = sqrt(AOR_AER_nkat_gFW))
cat("\n==================== AOR/AER ====================\n")
cat("Replicates per genotype x temperature:\n"); print(table(df_AOR$Genotype, df_AOR$Temperature))

mod_AOR_raw <- lm(y   ~ Genotype * Temperature, data = df_AOR)
mod_AOR     <- lm(y_t ~ Genotype * Temperature, data = df_AOR)
lev_AOR_raw <- leveneTest(y   ~ Genotype * Temperature, data = df_AOR)
lev_AOR     <- leveneTest(y_t ~ Genotype * Temperature, data = df_AOR)
assump_AOR <- data.frame(trait = "AOR_AER", scale = c("raw", "sqrt"),
                         shapiro_W = c(shapiro.test(residuals(mod_AOR_raw))$statistic,
                                       shapiro.test(residuals(mod_AOR))$statistic),
                         shapiro_P = c(shapiro.test(residuals(mod_AOR_raw))$p.value,
                                       shapiro.test(residuals(mod_AOR))$p.value),
                         levene_F  = c(lev_AOR_raw$`F value`[1], lev_AOR$`F value`[1]),
                         levene_P  = c(lev_AOR_raw$`Pr(>F)`[1],  lev_AOR$`Pr(>F)`[1]))
cat("\nAssumption checks (raw vs sqrt):\n"); print(assump_AOR, digits = 3)
aov_AOR <- Anova(mod_AOR, type = 2)
cat("\nType II ANOVA, sqrt scale:\n"); print(aov_AOR)

cl_AOR <- as.data.frame(cld(emmeans(mod_AOR, ~ Genotype | Temperature),
                            adjust = "tukey", Letters = letters, decreasing = TRUE)) %>%
  mutate(letter = trimws(.group)) %>%
  dplyr::select(Genotype, Temperature, emmean_model_scale = emmean, letter)
cat("\nTukey letters (a = highest mean within a temperature):\n")
print(cl_AOR %>% arrange(Temperature, desc(emmean_model_scale)))

letters_AOR <- df_AOR %>%
  group_by(Genotype, Temperature) %>%
  summarise(y_top = max(density(y, adjust = 1.2)$x), .groups = "drop") %>%
  left_join(cl_AOR %>% dplyr::select(Genotype, Temperature, letter),
            by = c("Genotype", "Temperature"))

p_AOR <- ggplot(df_AOR, aes(Genotype, y, fill = Temperature)) +
  geom_vline(xintercept = seq(1.5, 5.5, 1), color = "grey70", linetype = "dashed",
             linewidth = 0.4) +
  geom_violin(position = dodge, alpha = 0.75, color = "grey25",
              trim = FALSE, scale = "width", width = 1.0, adjust = 1.2) +
  geom_boxplot(aes(group = interaction(Genotype, Temperature)),
               position = dodge, width = 0.15, fill = "gray", alpha = 0.70,
               outlier.shape = NA, color = "grey20", linewidth = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
             size = 2, alpha = 0.55, shape = 21, stroke = 0.5, color = "black",
             show.legend = FALSE) +
  geom_text(data = letters_AOR,
            aes(x = Genotype, y = y_top, label = letter, group = Temperature,
                color = Temperature),
            position = dodge, size = 4, vjust = -0.4, fontface = "bold",
            show.legend = FALSE) +
  scale_fill_manual(values = pal, labels = temp_labels) +
  scale_color_manual(values = pal, labels = temp_labels, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "", y = expression(AOR/AER~activity~(nkat~g^{-1}~FW))) +
  violin_theme

######### PANEL C - GLY I (square root) #######################################################

df_GLY <- read.csv("data/Detoxification_enzyme/GLY.csv", encoding = "UTF-8") %>%
  transmute(Genotype    = factor(Cultivar,    levels = geno_order_full),
            Temperature = factor(Temperature, levels = temp_order),
            y = GLY_nkat_gFW, y_t = sqrt(GLY_nkat_gFW))
cat("\n==================== GLY I ====================\n")
cat("Replicates per genotype x temperature:\n"); print(table(df_GLY$Genotype, df_GLY$Temperature))

mod_GLY_raw <- lm(y   ~ Genotype * Temperature, data = df_GLY)
mod_GLY     <- lm(y_t ~ Genotype * Temperature, data = df_GLY)
lev_GLY_raw <- leveneTest(y   ~ Genotype * Temperature, data = df_GLY)
lev_GLY     <- leveneTest(y_t ~ Genotype * Temperature, data = df_GLY)
assump_GLY <- data.frame(trait = "GLY", scale = c("raw", "sqrt"),
                         shapiro_W = c(shapiro.test(residuals(mod_GLY_raw))$statistic,
                                       shapiro.test(residuals(mod_GLY))$statistic),
                         shapiro_P = c(shapiro.test(residuals(mod_GLY_raw))$p.value,
                                       shapiro.test(residuals(mod_GLY))$p.value),
                         levene_F  = c(lev_GLY_raw$`F value`[1], lev_GLY$`F value`[1]),
                         levene_P  = c(lev_GLY_raw$`Pr(>F)`[1],  lev_GLY$`Pr(>F)`[1]))
cat("\nAssumption checks (raw vs sqrt):\n"); print(assump_GLY, digits = 3)
aov_GLY <- Anova(mod_GLY, type = 2)
cat("\nType II ANOVA, sqrt scale:\n"); print(aov_GLY)

cl_GLY <- as.data.frame(cld(emmeans(mod_GLY, ~ Genotype | Temperature),
                            adjust = "tukey", Letters = letters, decreasing = TRUE)) %>%
  mutate(letter = trimws(.group)) %>%
  dplyr::select(Genotype, Temperature, emmean_model_scale = emmean, letter)
cat("\nTukey letters (a = highest mean within a temperature):\n")
print(cl_GLY %>% arrange(Temperature, desc(emmean_model_scale)))

letters_GLY <- df_GLY %>%
  group_by(Genotype, Temperature) %>%
  summarise(y_top = max(density(y, adjust = 1.2)$x), .groups = "drop") %>%
  left_join(cl_GLY %>% dplyr::select(Genotype, Temperature, letter),
            by = c("Genotype", "Temperature"))

p_GLY <- ggplot(df_GLY, aes(Genotype, y, fill = Temperature)) +
  geom_vline(xintercept = seq(1.5, 5.5, 1), color = "grey70", linetype = "dashed",
             linewidth = 0.4) +
  geom_violin(position = dodge, alpha = 0.75, color = "grey25",
              trim = FALSE, scale = "width", width = 1.0, adjust = 1.2) +
  geom_boxplot(aes(group = interaction(Genotype, Temperature)),
               position = dodge, width = 0.15, fill = "gray", alpha = 0.70,
               outlier.shape = NA, color = "grey20", linewidth = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.8),
             size = 2, alpha = 0.55, shape = 21, stroke = 0.5, color = "black",
             show.legend = FALSE) +
  geom_text(data = letters_GLY,
            aes(x = Genotype, y = y_top, label = letter, group = Temperature,
                color = Temperature),
            position = dodge, size = 4, vjust = -0.4, fontface = "bold",
            show.legend = FALSE) +
  scale_fill_manual(values = pal, labels = temp_labels) +
  scale_color_manual(values = pal, labels = temp_labels, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "", y = expression(plain("GLY I")~activity~(nkat~g^{-1}~FW))) +
  violin_theme

######### PANEL D - HEATMAP OF MEANS AS WITHIN-ENZYME Z-SCORES ################################

## Genotype x temperature means of the three enzymes with their Tukey letters
heat_df <- bind_rows(
  df_AKR %>% mutate(Trait = "AKR")     %>% left_join(cl_AKR, by = c("Genotype", "Temperature")),
  df_AOR %>% mutate(Trait = "AOR/AER") %>% left_join(cl_AOR, by = c("Genotype", "Temperature")),
  df_GLY %>% mutate(Trait = "GLY I")   %>% left_join(cl_GLY, by = c("Genotype", "Temperature"))) %>%
  group_by(Trait, Genotype, Temperature, letter) %>%
  summarise(n = n(), mean = mean(y), sd = sd(y), se = sd / sqrt(n), .groups = "drop") %>%
  group_by(Trait) %>%
  mutate(z = as.numeric(scale(mean))) %>%          # z-score over the 18 means of the enzyme
  ungroup() %>%
  mutate(Trait = factor(Trait, levels = c("AKR", "AOR/AER", "GLY I")))

p_heat <- ggplot(heat_df, aes(Trait, Genotype, fill = z)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f\n%s", mean, letter)), size = 3, colour = "gray15",
            lineheight = 0.9) +
  facet_wrap(~ Temperature, nrow = 1, labeller = labeller(Temperature = temp_labels)) +
  scale_fill_gradient2(low = "#3B8BC2", mid = "#F2C14E", high = "firebrick",
                       midpoint = 0, name = "Z-score\n(within trait)") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev(geno_order_full)) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid        = element_blank(),
        strip.background  = element_rect(fill = "gray92", colour = NA),
        strip.text        = element_text(face = "bold"),
        axis.text.x       = element_text(face = "bold", family = "serif"),
        axis.text.y       = element_text(face = "plain"),
        legend.position   = "bottom",
        legend.key.width  = unit(20, "pt"),
        legend.key.height = unit(8,  "pt"),
        plot.tag          = element_text(size = 14))

######### FOUR-PANEL FIGURE ###################################################################

panel <- (p_AKR + (p_AOR + theme(legend.position = "none"))) /
  ((p_GLY + theme(legend.position = "none")) + p_heat) +
  plot_annotation(tag_levels = "A")

######### OUTPUT ##############################################################################

## Percent change relative to the 20 C control, per enzyme, in original units
means_out <- heat_df %>%
  dplyr::select(trait = Trait, Genotype, Temperature, n, mean, sd, se, z, letter) %>%
  group_by(trait, Genotype) %>%
  mutate(pct_change_vs_20C = 100 * (mean - mean[Temperature == "20C"]) / mean[Temperature == "20C"]) %>%
  ungroup()
cat("\nMeans (nkat/g FW) and % change relative to own 20 C control:\n")
print(as.data.frame(means_out %>% dplyr::select(trait, Genotype, Temperature, mean, pct_change_vs_20C)),
      digits = 3)

write.csv(bind_rows(cbind(trait = "AKR",     scale = "none", term = rownames(aov_AKR), as.data.frame(aov_AKR)),
                    cbind(trait = "AOR_AER", scale = "sqrt", term = rownames(aov_AOR), as.data.frame(aov_AOR)),
                    cbind(trait = "GLY",     scale = "sqrt", term = rownames(aov_GLY), as.data.frame(aov_GLY))),
          "results/05_carbonyl_anova.csv", row.names = FALSE)
write.csv(bind_rows(assump_AKR, assump_AOR, assump_GLY),
          "results/05_carbonyl_assumptions.csv", row.names = FALSE)
write.csv(bind_rows(cbind(trait = "AKR",     scale = "none", cl_AKR),
                    cbind(trait = "AOR_AER", scale = "sqrt", cl_AOR),
                    cbind(trait = "GLY",     scale = "sqrt", cl_GLY)),
          "results/05_carbonyl_tukey_letters.csv", row.names = FALSE)
write.csv(means_out, "results/05_carbonyl_means.csv", row.names = FALSE)

ggsave("figures/Fig5_carbonyl_enzymes.tiff", panel, width = 11.8, height = 9.1,
       dpi = 300, bg = "white", compression = "lzw")
cat("\nFigure written to figures/Fig5_carbonyl_enzymes.tiff\n")

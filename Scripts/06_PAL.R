###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 06 - Phenylalanine ammonia-lyase: two-way ANOVA, Tukey letters, violin (Fig. 6)###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Phenylalanine ammonia-lyase (PAL) activity, nkat per g fresh weight, of six Mv winter
## wheat cultivars at 20 (control), 35 and 40 degrees C.
##
## Model : PAL ~ Genotype * Temperature (fixed effects), Type II ANOVA (car::Anova).
## Scale : untransformed. Shapiro-Wilk on the residuals and Levene's test pass on the raw
##         scale (P = 0.052 and 0.17); log and square-root transforms both fail and are
##         reported for comparison in the assumptions file.
## Letters: cultivars compared within each temperature on estimated marginal means
##         (emmeans), Tukey-adjusted, compact letter display with "a" = highest mean.
## Plot  : violin + inner box plot + individual observations per genotype x temperature,
##         letters above each violin.
##
## Input : data/PAL/PAL.csv (Genotype = code, Cultivar = full name, Temperature, PAL_nkat_gFW)
## Output: figures/Fig6_PAL.tiff
##         results/06_PAL_anova.csv, results/06_PAL_assumptions.csv,
##         results/06_PAL_tukey_letters.csv, results/06_PAL_means.csv

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/06_PAL.R")
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
temp_labels <- c("20C" = "20°C", "35C" = "35°C", "40C" = "40°C")
pal         <- c("20C" = "#3B8BC2", "35C" = "#F2C14E", "40C" = "#D1495B")
dodge       <- position_dodge(width = 0.8)

######### DATA ################################################################################

df_PAL <- read.csv("data/PAL/PAL.csv", encoding = "UTF-8")

## Cultivar names are carried in the data file (column Cultivar); the code column
## (Genotype) fixes the plotting order
geno_order_full <- df_PAL %>% distinct(Genotype, Cultivar) %>%
  { .$Cultivar[match(geno_order, .$Genotype)] }

df_PAL <- df_PAL %>%
  transmute(Genotype    = factor(Cultivar,    levels = geno_order_full),
            Temperature = factor(Temperature, levels = temp_order),
            y           = PAL_nkat_gFW)

cat("Replicates per genotype x temperature:\n")
print(table(df_PAL$Genotype, df_PAL$Temperature))

######### MODEL AND ASSUMPTION CHECKS #########################################################

mod      <- lm(y       ~ Genotype * Temperature, data = df_PAL)
mod_log  <- lm(log(y)  ~ Genotype * Temperature, data = df_PAL)
mod_sqrt <- lm(sqrt(y) ~ Genotype * Temperature, data = df_PAL)

lev      <- leveneTest(y       ~ Genotype * Temperature, data = df_PAL)
lev_log  <- leveneTest(log(y)  ~ Genotype * Temperature, data = df_PAL)
lev_sqrt <- leveneTest(sqrt(y) ~ Genotype * Temperature, data = df_PAL)

assumptions <- data.frame(
  trait     = "PAL",
  scale     = c("raw", "log", "sqrt"),
  shapiro_W = c(shapiro.test(residuals(mod))$statistic,
                shapiro.test(residuals(mod_log))$statistic,
                shapiro.test(residuals(mod_sqrt))$statistic),
  shapiro_P = c(shapiro.test(residuals(mod))$p.value,
                shapiro.test(residuals(mod_log))$p.value,
                shapiro.test(residuals(mod_sqrt))$p.value),
  levene_F  = c(lev$`F value`[1], lev_log$`F value`[1], lev_sqrt$`F value`[1]),
  levene_P  = c(lev$`Pr(>F)`[1],  lev_log$`Pr(>F)`[1],  lev_sqrt$`Pr(>F)`[1]))
cat("\nAssumption checks (raw scale used; log and sqrt for comparison):\n")
print(assumptions, digits = 3)
write.csv(assumptions, "results/06_PAL_assumptions.csv", row.names = FALSE)

## Raw scale passes both tests -> analysed untransformed
aov_tab <- Anova(mod, type = 2)
cat("\nType II ANOVA:\n"); print(aov_tab)
write.csv(cbind(trait = "PAL", scale = "none", term = rownames(aov_tab), as.data.frame(aov_tab)),
          "results/06_PAL_anova.csv", row.names = FALSE)

######### TUKEY LETTERS (cultivars within each temperature) ###################################

emm <- emmeans(mod, ~ Genotype | Temperature)
## The pairwise comparisons behind the letters are Tukey-adjusted. emmeans prints a note
## that "tukey" was changed to "sidak": that applies only to the confidence intervals of
## the means shown in the table, not to the comparisons.
cl <- as.data.frame(cld(emm, adjust = "tukey", Letters = letters, decreasing = TRUE)) %>%
  mutate(letter = trimws(.group)) %>%
  dplyr::select(Genotype, Temperature, emmean, letter)
cat("\nTukey letters (a = highest mean within a temperature):\n")
print(cl %>% arrange(Temperature, desc(emmean)))
write.csv(cbind(trait = "PAL", scale = "none", cl), "results/06_PAL_tukey_letters.csv",
          row.names = FALSE)

######### MEANS AND PERCENT CHANGE ############################################################

means_df <- df_PAL %>%
  group_by(Genotype, Temperature) %>%
  summarise(n = n(), mean = mean(y), sd = sd(y), se = sd / sqrt(n), .groups = "drop") %>%
  group_by(Genotype) %>%
  mutate(pct_change_vs_20C = 100 * (mean - mean[Temperature == "20C"]) / mean[Temperature == "20C"]) %>%
  ungroup()
cat("\nMeans (nkat/g FW) and % change relative to own 20 C control:\n")
print(as.data.frame(means_df), digits = 3)
write.csv(cbind(trait = "PAL", means_df), "results/06_PAL_means.csv", row.names = FALSE)

######### VIOLIN PLOT #########################################################################

letters_df <- df_PAL %>%
  group_by(Genotype, Temperature) %>%
  summarise(y_top = max(density(y, adjust = 1.2)$x), .groups = "drop") %>%
  left_join(cl %>% dplyr::select(Genotype, Temperature, letter),
            by = c("Genotype", "Temperature"))

p_PAL <- ggplot(df_PAL, aes(Genotype, y, fill = Temperature)) +
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
  geom_text(data = letters_df,
            aes(x = Genotype, y = y_top, label = letter, group = Temperature,
                color = Temperature),
            position = dodge, size = 4, vjust = -0.4, fontface = "bold",
            show.legend = FALSE) +
  scale_fill_manual(values = pal, labels = temp_labels) +
  scale_color_manual(values = pal, labels = temp_labels, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.10))) +
  labs(x = "", y = expression(PAL~activity~(nkat~g^{-1}~FW))) +
  theme_bw(base_size = 12) +
  theme(panel.grid      = element_blank(),
        legend.position = "top",
        axis.title      = element_text(face = "bold"),
        axis.title.y    = element_text(family = "serif", face = "bold"),
        axis.text.x     = element_text(color = "black", size = 10))

ggsave("figures/Fig6_PAL.tiff", p_PAL, width = 9.6, height = 6.1,
       dpi = 300, bg = "white", compression = "lzw")
cat("\nFigure written to figures/Fig6_PAL.tiff\n")

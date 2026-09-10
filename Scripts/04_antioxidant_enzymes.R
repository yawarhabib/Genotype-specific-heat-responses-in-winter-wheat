###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 04 - AsA-GSH enzymes and GST: two-way ANOVA, Tukey letters, violins (Fig. 4)  ###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Activities of the ascorbate-glutathione cycle enzymes and glutathione S-transferase in
## six Mv winter wheat cultivars at 20 (control), 35 and 40 degrees C: ascorbate
## peroxidase (APX), dehydroascorbate reductase (DHAR), glutathione reductase (GR) and
## glutathione S-transferase (GST), all in nkat per g fresh weight.
##
## Model : activity ~ Genotype * Temperature (fixed effects), Type II ANOVA (car::Anova).
## Scale : chosen per enzyme from Shapiro-Wilk (residuals) and Levene's test, which are
##         reported for the raw scale and for the scale used:
##           APX  untransformed  (raw passes both tests)
##           DHAR log            (raw Shapiro P = 0.053, transformed as a precaution)
##           GR   log            (raw Shapiro P = 0.011)
##           GST  square root    (raw Shapiro P = 0.031; log over-corrects and fails both)
## Letters: cultivars compared within each temperature on estimated marginal means
##         (emmeans), Tukey-adjusted, compact letter display with "a" = highest mean.
## Plot  : violin + inner box plot + individual observations per genotype x temperature,
##         letters above each violin; four panels (A APX, B DHAR, C GR, D GST).
##
## Input : data/Antioxidant_enzyme/APX.csv, DHAR.csv, GR.csv, GST.csv
##         (Genotype = code, Cultivar = full name, Temperature, <enzyme>_nkat_gFW)
## Output: figures/Fig4_antioxidant_enzymes.tiff
##         results/04_antioxidant_anova.csv, results/04_antioxidant_assumptions.csv,
##         results/04_antioxidant_tukey_letters.csv, results/04_antioxidant_means.csv

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/04_antioxidant_enzymes.R")
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

## Transformation used for the model of each enzyme (see DESCRIPTION)
trans_fun <- list(APX = identity, DHAR = log, GR = log, GST = sqrt)
trans_lab <- c(APX = "none", DHAR = "log", GR = "log", GST = "sqrt")

######### ONE PANEL: MODEL, ASSUMPTIONS, ANOVA, TUKEY LETTERS AND VIOLIN PLOT #################

## Results of every panel are collected in these lists and written out at the end
anova_out <- list(); assump_out <- list(); letters_out <- list(); means_out <- list()

make_violin <- function(trait, ylab_expr) {

  cat("\n==================== ", trait, " ====================\n", sep = "")

  ## --- data: cultivar names are carried in the file (column Cultivar) -------------------
  df <- read.csv(paste0("data/Antioxidant_enzyme/", trait, ".csv"), encoding = "UTF-8")
  geno_order_full <- df %>% distinct(Genotype, Cultivar) %>%
    { .$Cultivar[match(geno_order, .$Genotype)] }
  d <- df %>%
    transmute(Genotype    = factor(Cultivar,    levels = geno_order_full),
              Temperature = factor(Temperature, levels = temp_order),
              y           = .data[[paste0(trait, "_nkat_gFW")]])
  cat("Replicates per genotype x temperature:\n"); print(table(d$Genotype, d$Temperature))

  ## --- models on the raw scale and on the scale used -----------------------------------
  tf      <- trans_fun[[trait]]
  d$y_t   <- tf(d$y)
  mod_raw <- lm(y   ~ Genotype * Temperature, data = d)
  mod     <- lm(y_t ~ Genotype * Temperature, data = d)
  lev_raw <- leveneTest(y   ~ Genotype * Temperature, data = d)
  lev_t   <- leveneTest(y_t ~ Genotype * Temperature, data = d)

  assumptions <- data.frame(
    trait     = trait,
    scale     = c("raw", trans_lab[[trait]]),
    shapiro_W = c(shapiro.test(residuals(mod_raw))$statistic,
                  shapiro.test(residuals(mod))$statistic),
    shapiro_P = c(shapiro.test(residuals(mod_raw))$p.value,
                  shapiro.test(residuals(mod))$p.value),
    levene_F  = c(lev_raw$`F value`[1], lev_t$`F value`[1]),
    levene_P  = c(lev_raw$`Pr(>F)`[1],  lev_t$`Pr(>F)`[1]))
  cat("\nAssumption checks (raw scale vs scale used):\n"); print(assumptions, digits = 3)
  assump_out[[trait]] <<- assumptions

  ## --- Type II ANOVA on the scale used -------------------------------------------------
  aov_tab <- Anova(mod, type = 2)
  cat("\nType II ANOVA, scale = ", trans_lab[[trait]], ":\n", sep = ""); print(aov_tab)
  anova_out[[trait]] <<- cbind(trait = trait, scale = trans_lab[[trait]],
                               term = rownames(aov_tab), as.data.frame(aov_tab))

  ## --- Tukey letters, cultivars within each temperature --------------------------------
  emm <- emmeans(mod, ~ Genotype | Temperature)
  ## The pairwise comparisons behind the letters are Tukey-adjusted. emmeans prints a note
  ## that "tukey" was changed to "sidak": that applies only to the confidence intervals of
  ## the means shown in the table, not to the comparisons.
  cl <- as.data.frame(cld(emm, adjust = "tukey", Letters = letters, decreasing = TRUE)) %>%
    mutate(letter = trimws(.group)) %>%
    dplyr::select(Genotype, Temperature, emmean_model_scale = emmean, letter)
  cat("\nTukey letters (a = highest mean within a temperature):\n")
  print(cl %>% arrange(Temperature, desc(emmean_model_scale)))
  letters_out[[trait]] <<- cbind(trait = trait, scale = trans_lab[[trait]], cl)

  ## --- means, SE and % change relative to the 20 C control (original units) ------------
  means_df <- d %>%
    group_by(Genotype, Temperature) %>%
    summarise(n = n(), mean = mean(y), sd = sd(y), se = sd / sqrt(n), .groups = "drop")
  ctrl <- means_df %>% filter(Temperature == "20C") %>% dplyr::select(Genotype, ctrl = mean)
  means_df <- means_df %>% left_join(ctrl, by = "Genotype") %>%
    mutate(pct_change_vs_20C = 100 * (mean - ctrl) / ctrl) %>% dplyr::select(-ctrl)
  cat("\nMeans (nkat/g FW) and % change relative to own 20 C control:\n")
  print(as.data.frame(means_df), digits = 3)
  means_out[[trait]] <<- cbind(trait = trait, means_df)

  ## --- letter position: top of each violin -----------------------------------------------
  letters_df <- d %>%
    group_by(Genotype, Temperature) %>%
    summarise(y_top = max(density(y, adjust = 1.2)$x), .groups = "drop") %>%
    left_join(cl %>% dplyr::select(Genotype, Temperature, letter),
              by = c("Genotype", "Temperature"))

  ## --- plot ----------------------------------------------------------------------------
  ggplot(d, aes(Genotype, y, fill = Temperature)) +
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
    labs(x = "", y = ylab_expr) +
    theme_bw(base_size = 12) +
    theme(panel.grid   = element_blank(),
          legend.position = "top",
          axis.title   = element_text(face = "bold"),
          axis.title.y = element_text(family = "serif", face = "bold"),
          axis.text.x  = element_text(color = "black", size = 10),
          plot.tag     = element_text(size = 14))
}

######### THE FOUR PANELS #####################################################################

p_APX  <- make_violin("APX",  expression(APX~activity~(nkat~g^{-1}~FW)))
p_DHAR <- make_violin("DHAR", expression(DHAR~activity~(nkat~g^{-1}~FW)))
p_GR   <- make_violin("GR",   expression(GR~activity~(nkat~g^{-1}~FW)))
p_GST  <- make_violin("GST",  expression(GST~activity~(nkat~g^{-1}~FW)))

## 2 x 2 panels with one shared legend on top. patchwork places a collected legend
## according to the session theme, so that is set for the assembly and restored afterwards.
old_theme <- theme_set(theme_get() + theme(legend.position = "top",
                                           legend.justification = "center"))
panel <- (p_APX + p_DHAR) / (p_GR + p_GST) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A")

######### OUTPUT ##############################################################################

write.csv(bind_rows(anova_out),   "results/04_antioxidant_anova.csv",         row.names = FALSE)
write.csv(bind_rows(assump_out),  "results/04_antioxidant_assumptions.csv",   row.names = FALSE)
write.csv(bind_rows(letters_out), "results/04_antioxidant_tukey_letters.csv", row.names = FALSE)
write.csv(bind_rows(means_out),   "results/04_antioxidant_means.csv",         row.names = FALSE)

ggsave("figures/Fig4_antioxidant_enzymes.tiff", panel, width = 11.8, height = 9.1,
       dpi = 300, bg = "white", compression = "lzw")
theme_set(old_theme)
cat("\nFigure written to figures/Fig4_antioxidant_enzymes.tiff\n")

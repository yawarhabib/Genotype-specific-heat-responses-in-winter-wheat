###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 08 - PCA biplot (Fig. 8), STI heatmap (Fig. 9), Spearman matrix (Fig. S1)     ###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Integrative analyses on the genotype x temperature means (6 cultivars x 3 temperatures
## = 18 rows) of every trait measured in scripts 01-07: shoot fresh weight, gas exchange,
## pigments, AsA-GSH enzymes and GST, reactive-carbonyl enzymes, PAL and polyamines.
##
## PCA (Fig. 8) : prcomp on the standardised means of 22 traits. Vapour-pressure deficit
##                (collinear with temperature) and anthocyanin (measured at 20 and 40 C
##                only) are excluded so that all 18 samples enter the ordination. Biplot of
##                sample scores (68 % normal ellipses per temperature) and trait loadings,
##                loading arrows coloured by functional group and rescaled for display
##                only; the inset shows the variance explained by the first ten components.
## STI (Fig. 9) : trait-wise stress-tolerance index after Fernandez (1992),
##                  STI = (Yp x Ys) / mean(Yp)^2
##                with Yp the 20 C mean and Ys the 35 or 40 C mean of a trait in a cultivar
##                and mean(Yp) the control mean of that trait over the six cultivars.
##                Computed for every trait except VPD and anthocyanin, standardised to
##                within-trait z-scores, and shown as a heatmap with rows split by stress
##                level (35 vs 40 C) and clustered within each block, columns clustered
##                over all rows; Ward's method (ward.D2) on Euclidean distances.
## Spearman (Fig. S1): rank correlations among all traits except VPD over the 18 means,
##                pairwise-complete observations (anthocyanin has 12), significance from
##                corrplot::cor.mtest; only P < 0.05 values are printed.
##
## Input : all 15 data files under data/ (see scripts 01-07 for their structure)
## Output: figures/Fig8_PCA_biplot.tiff, figures/Fig9_STI_heatmap.tiff,
##         figures/FigS1_spearman_correlation.tiff
##         results/08_trait_means.csv          (the 18 genotype x temperature means)
##         results/08_PCA_scores.csv, 08_PCA_loadings.csv, 08_PCA_variance.csv
##         results/08_spearman_r.csv, 08_spearman_P.csv, 08_spearman_pairs.csv
##         results/08_STI_long.csv, 08_STI_matrix.csv, 08_STI_zscores.csv,
##         results/08_STI_row_order.csv, 08_STI_column_order.csv

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/08_multivariate_STI.R")
main_dir <- here::here()
setwd(main_dir)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(tibble)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(ggnewscale)      # two colour scales in one ggplot (arrows vs points)
  library(ComplexHeatmap)
  library(circlize)
  library(corrplot)
})

dir.create("figures", showWarnings = FALSE)
dir.create("results", showWarnings = FALSE)

######### SHARED DEFINITIONS ##################################################################

## Functional group of every trait (loading arrow colours in the biplot)
trait_group <- c(
  APX = "Antioxidant", DHAR = "Antioxidant", GR = "Antioxidant", GST = "Antioxidant",
  AKR = "RCF-detoxification", AOR_AER = "RCF-detoxification", GLY = "RCF-detoxification",
  PUT = "Polyamine", SPD = "Polyamine", SPM = "Polyamine", DAP = "Polyamine", CAD = "Polyamine",
  Pn = "Gas exchange", E = "Gas exchange", gs = "Gas exchange", WUE = "Gas exchange",
  Chl_ab = "Pigment", Chl_a_b_ratio = "Pigment", Carotenoid = "Pigment",
  Anthocyanin_nmol_gFW = "Pigment",
  FW = "Other", PAL = "Other")

group_colors <- c("Antioxidant"        = "#1b9e77",
                  "RCF-detoxification" = "#d95f02",
                  "Polyamine"          = "#7570b3",
                  "Gas exchange"       = "#e7298a",
                  "Pigment"            = "#66a61e",
                  "Other"              = "#8c8c8c")

## Display labels: data column name -> label shown on the plots
display_names <- c(Anthocyanin_nmol_gFW = "Anthocyanin",
                   Chl_ab               = "Chl a+b",
                   Chl_a_b_ratio        = "Chl a/b")
relabel <- function(x) unname(ifelse(x %in% names(display_names), display_names[x], x))

temp_colors <- c("20C" = "#2E8B57", "35C" = "#3f6d94", "40C" = "#e23b3b")
geno_shapes <- c(DAN = 16, IKV = 17, KOL = 15, LUC = 18, PIR = 8, TAR = 4)

######### 1. LOAD ALL DATA BLOCKS AND AVERAGE PER GENOTYPE x TEMPERATURE ######################

rd <- function(path) read.csv(path, encoding = "UTF-8")

APX     <- rd("data/Antioxidant_enzyme/APX.csv")
DHAR    <- rd("data/Antioxidant_enzyme/DHAR.csv")
GR      <- rd("data/Antioxidant_enzyme/GR.csv")
GST     <- rd("data/Antioxidant_enzyme/GST.csv")
AKR     <- rd("data/Detoxification_enzyme/AKR.csv")
AOR_AER <- rd("data/Detoxification_enzyme/AOR_AER.csv")
GLY     <- rd("data/Detoxification_enzyme/GLY.csv")
PA      <- rd("data/Polyamines/PA_polyamines.csv")
FW      <- rd("data/Fresh_weight/FW.csv")
PAL     <- rd("data/PAL/PAL.csv")
Antho   <- rd("data/Pigments/Anthocyanin.csv")
Carot   <- rd("data/Pigments/Carotenoid.csv")
Chl     <- rd("data/Pigments/Chlorophyll_ab.csv")
Chl_ab  <- rd("data/Pigments/Chlorophyll_ab_ratio.csv")

## Gas exchange: use the genotype code and write the temperature as "20C" etc.
Gas <- rd("data/Gas_exchange/Leaf_gas_Exchange.csv") %>%
  transmute(Genotype = Genotype_code, Temperature = paste0(Temperature, "C"),
            Pn, E, gs, VPD, WUE)

## Mean of every numeric column per genotype x temperature (text columns are ignored)
agg_block <- function(d) {
  d %>%
    group_by(Genotype, Temperature) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE)), .groups = "drop")
}

all_means <- list(APX, DHAR, GR, GST, AKR, AOR_AER, GLY, Gas, PA,
                  Antho, Carot, Chl, Chl_ab, FW, PAL) %>%
  map(agg_block) %>%
  reduce(full_join, by = c("Genotype", "Temperature")) %>%
  rename_with(~ str_remove(., "_(nkat_gFW|ug_gFW|g)$"), -c(Genotype, Temperature)) %>%
  arrange(Temperature, Genotype)

cat("Genotype x temperature means:", nrow(all_means), "rows,",
    ncol(all_means) - 2, "traits\n")
print(names(all_means))
write.csv(all_means, "results/08_trait_means.csv", row.names = FALSE)

######### 2. PCA ##############################################################################

## VPD is confounded with temperature; anthocyanin has no 35 C value (NA would drop
## six samples). Both stay in all_means for the STI and correlation steps that can use them.
mat <- all_means %>% dplyr::select(-Genotype, -Temperature, -VPD, -Anthocyanin_nmol_gFW)
stopifnot(all(is.finite(as.matrix(mat))))
stopifnot(all(apply(mat, 2, sd) > 0))

pca     <- prcomp(mat, scale. = TRUE)
var_exp <- round(summary(pca)$importance["Proportion of Variance", ] * 100, 1)
cat("\nVariance explained (%):\n"); print(var_exp)

## The sign of a principal component is arbitrary. Orient the axes so that the control
## (20 C) samples score positive on PC1 and the 35 C samples positive on PC2, the
## orientation described in the Results; scores and loadings are flipped together.
flip <- c(PC1 = sign(mean(pca$x[all_means$Temperature == "20C", "PC1"])),
          PC2 = sign(mean(pca$x[all_means$Temperature == "35C", "PC2"])))
for (pc in names(flip)) {
  pca$x[, pc]        <- pca$x[, pc]        * flip[[pc]]
  pca$rotation[, pc] <- pca$rotation[, pc] * flip[[pc]]
}

scores <- as.data.frame(pca$x) %>%
  mutate(Genotype = all_means$Genotype, Temperature = all_means$Temperature)

## Loading arrows: eigenvectors rescaled for display only. The longest arrow reaches
## 30 % of the smaller span of the score cloud. The true loadings are written to results.
arrow_frac <- 0.30
loadings   <- as.data.frame(pca$rotation) %>%
  rownames_to_column("Trait") %>%
  mutate(Group = trait_group[Trait], Label = relabel(Trait))
score_span <- min(diff(range(scores$PC1)), diff(range(scores$PC2)))
cur_max    <- max(sqrt(loadings$PC1^2 + loadings$PC2^2))
loadings   <- loadings %>%
  mutate(PC1 = PC1 / cur_max * score_span * arrow_frac,
         PC2 = PC2 / cur_max * score_span * arrow_frac)
stopifnot(!any(is.na(loadings$Group)))

var_df10 <- tibble(PC  = factor(paste0("PC", seq_along(var_exp)),
                                levels = paste0("PC", seq_along(var_exp))),
                   Var = var_exp) %>% slice(1:10)

## Results tables
write.csv(scores %>% relocate(Genotype, Temperature), "results/08_PCA_scores.csv",
          row.names = FALSE)
write.csv(as.data.frame(pca$rotation) %>% rownames_to_column("Trait") %>%
            mutate(Label = relabel(Trait), Group = trait_group[Trait]) %>%
            relocate(Trait, Label, Group),
          "results/08_PCA_loadings.csv", row.names = FALSE)
write.csv(data.frame(PC           = paste0("PC", seq_along(pca$sdev)),
                     SD           = pca$sdev,
                     Eigenvalue   = pca$sdev^2,
                     Prop_Var_pct = round(pca$sdev^2 / sum(pca$sdev^2) * 100, 2),
                     Cum_Var_pct  = round(cumsum(pca$sdev^2) / sum(pca$sdev^2) * 100, 2)),
          "results/08_PCA_variance.csv", row.names = FALSE)

######### 3. FIG. 8 - PCA BIPLOT WITH SCREE INSET #############################################

p_main <- ggplot() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey75") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey75") +
  stat_ellipse(data = scores, aes(PC1, PC2, fill = Temperature),
               geom = "polygon", alpha = 0.15, color = NA, type = "norm", level = 0.68) +
  scale_fill_manual(values = temp_colors, guide = "none") +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2, color = Group),
               arrow = arrow(length = unit(0.22, "cm")), linewidth = 0.6) +
  geom_text_repel(data = loadings, aes(PC1, PC2, label = Label, color = Group),
                  fontface = "bold", size = 3.6, seed = 1, show.legend = FALSE) +
  scale_color_manual(values = group_colors, name = "Trait group",
                     guide = guide_legend(order = 3)) +
  new_scale_color() +
  geom_point(data = scores, aes(PC1, PC2, color = Temperature, shape = Genotype),
             size = 4.5, stroke = 1.2) +
  scale_color_manual(values = temp_colors, name = "Temperature",
                     guide = guide_legend(order = 1)) +
  scale_shape_manual(values = geno_shapes, name = "Genotype",
                     guide = guide_legend(order = 2)) +
  geom_text_repel(data = scores,
                  aes(PC1, PC2, label = paste(Genotype, Temperature, sep = "|")),
                  size = 3, color = "grey25", box.padding = 0.4, seed = 2) +
  labs(x = paste0("PC1 (", var_exp[1], "%)"), y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

p_scree <- ggplot(var_df10, aes(PC, Var, group = 1)) +
  geom_col(aes(fill = Var), width = 0.7, show.legend = FALSE) +
  scale_fill_gradient(low = "#bcd4c6", high = "#3f6d94") +
  geom_line(color = "grey40", linewidth = 0.5) +
  geom_point(color = "grey25", size = 1.4) +
  geom_text(aes(label = ifelse(Var >= 2, paste0(Var, "%"), "")),
            vjust = -0.6, size = 2.3, color = "grey20") +
  scale_y_continuous(limits = c(0, max(var_df10$Var) * 1.22),
                     expand = expansion(mult = c(0, 0.05))) +
  labs(x = NULL, y = "Variance (%)", title = "Scree") +
  theme_classic(base_size = 12) +
  theme(plot.title      = element_text(face = "bold", hjust = 0.5),
        plot.background = element_rect(fill = "white", color = NA))

p_pca <- p_main +
  inset_element(p_scree, left = 0.60, bottom = 0.72, right = 0.999, top = 0.99,
                align_to = "panel")

ggsave("figures/Fig8_PCA_biplot.tiff", p_pca, width = 12, height = 9.25,
       dpi = 300, bg = "white", compression = "lzw")
cat("\nFigure written to figures/Fig8_PCA_biplot.tiff\n")

######### 4. STRESS-TOLERANCE INDEX ###########################################################

drop_traits <- c("VPD", "Anthocyanin_nmol_gFW")

ctrl_long <- all_means %>%
  filter(Temperature == "20C") %>%
  pivot_longer(-c(Genotype, Temperature), names_to = "Trait", values_to = "Yp") %>%
  filter(!Trait %in% drop_traits) %>%
  dplyr::select(Genotype, Trait, Yp)

trait_meanYp <- ctrl_long %>% group_by(Trait) %>% summarise(meanYp = mean(Yp), .groups = "drop")

sti_df <- all_means %>%
  filter(Temperature %in% c("35C", "40C")) %>%
  pivot_longer(-c(Genotype, Temperature), names_to = "Trait", values_to = "Ys") %>%
  filter(!Trait %in% drop_traits) %>%
  left_join(ctrl_long,    by = c("Genotype", "Trait")) %>%
  left_join(trait_meanYp, by = "Trait") %>%
  mutate(STI    = (Yp * Ys) / (meanYp^2),
         row_id = paste(Genotype, Temperature, sep = "|"))

sti_mat <- sti_df %>%
  dplyr::select(row_id, Trait, STI) %>%
  pivot_wider(names_from = Trait, values_from = STI) %>%
  column_to_rownames("row_id") %>%
  as.matrix()
stopifnot(!anyNA(sti_mat))

mat_z   <- scale(sti_mat)                                  # z-score per trait (column)
ann_row <- data.frame(Stress = factor(ifelse(grepl("35C$", rownames(sti_mat)), "35C", "40C")),
                      row.names = rownames(sti_mat))

######### 5. FIG. 9 - STI HEATMAP #############################################################

ht <- Heatmap(mat_z,
              row_split                 = ann_row$Stress,
              cluster_row_slices        = FALSE,
              clustering_method_rows    = "ward.D2",
              clustering_method_columns = "ward.D2",
              column_labels             = relabel(colnames(mat_z)),
              col = colorRamp2(c(-2, 0, 2), c("#2A9D8F", "white", "#9A031E")),
              left_annotation = rowAnnotation(
                Stress = ann_row$Stress,
                col    = list(Stress = c("35C" = "#577590", "40C" = "#f94144"))),
              name = "z(STI)")

tiff("figures/Fig9_STI_heatmap.tiff", width = 13, height = 8.6, units = "in", res = 300,
     compression = "lzw")
ht_drawn <- draw(ht)
invisible(dev.off())
cat("Figure written to figures/Fig9_STI_heatmap.tiff\n")

## Cluster order as drawn, rows flattened block by block (top to bottom)
ro <- row_order(ht_drawn); co <- column_order(ht_drawn)
row_order_tbl <- imap_dfr(ro, function(idx, slice)
  tibble(Stress = slice, Rank_in_block = seq_along(idx), row_id = rownames(mat_z)[idx])) %>%
  mutate(Overall_rank = row_number()) %>% relocate(Overall_rank)
col_order_tbl <- tibble(Rank = seq_along(co), Trait = colnames(mat_z)[co],
                        Label = relabel(colnames(mat_z)[co]))

write.csv(sti_df %>% mutate(Trait_label = relabel(Trait)) %>%
            dplyr::select(Genotype, Temperature, Trait, Trait_label, Yp, Ys, meanYp, STI),
          "results/08_STI_long.csv", row.names = FALSE)
write.csv(as.data.frame(sti_mat) %>% rownames_to_column("row_id"),
          "results/08_STI_matrix.csv", row.names = FALSE)
write.csv(as.data.frame(mat_z) %>% rownames_to_column("row_id"),
          "results/08_STI_zscores.csv", row.names = FALSE)
write.csv(row_order_tbl, "results/08_STI_row_order.csv",    row.names = FALSE)
write.csv(col_order_tbl, "results/08_STI_column_order.csv", row.names = FALSE)

######### 6. FIG. S1 - SPEARMAN CORRELATION MATRIX ############################################

trait_mat <- all_means %>% dplyr::select(-Genotype, -Temperature, -VPD) %>% as.matrix()
cor_res   <- cor(trait_mat, method = "spearman", use = "pairwise.complete.obs")
## cor.test gives exact P-values for pairs without tied ranks and the asymptotic
## approximation where ties occur; its warning about the latter is silenced here
p_mat     <- suppressWarnings(cor.mtest(trait_mat, conf.level = 0.95, method = "spearman"))$p
dimnames(cor_res) <- list(relabel(rownames(cor_res)), relabel(colnames(cor_res)))
dimnames(p_mat)   <- dimnames(cor_res)

tiff("figures/FigS1_spearman_correlation.tiff", width = 13, height = 8.6, units = "in",
     res = 300, compression = "lzw")
corrplot.mixed(cor_res,
               lower         = "ellipse",
               upper         = "number",
               lower.col     = colorRampPalette(c("#2A9D8F", "white", "#9A031E"))(200),
               upper.col     = colorRampPalette(c("#2A9D8F", "white", "#9A031E"))(200),
               tl.col        = "black",
               tl.cex        = 0.85,
               tl.pos        = "lt",
               number.cex    = 0.65,
               number.digits = 2,
               p.mat         = p_mat,
               sig.level     = 0.05,
               insig         = "blank",
               mar           = c(0, 0, 2, 0))
invisible(dev.off())
cat("Figure written to figures/FigS1_spearman_correlation.tiff\n")

ut <- upper.tri(cor_res)
write.csv(as.data.frame(cor_res) %>% rownames_to_column("Trait"),
          "results/08_spearman_r.csv", row.names = FALSE)
write.csv(as.data.frame(p_mat) %>% rownames_to_column("Trait"),
          "results/08_spearman_P.csv", row.names = FALSE)
write.csv(data.frame(Trait_1 = rownames(cor_res)[row(cor_res)[ut]],
                     Trait_2 = colnames(cor_res)[col(cor_res)[ut]],
                     rho     = round(cor_res[ut], 3),
                     P       = signif(p_mat[ut], 3)) %>% arrange(desc(abs(rho))),
          "results/08_spearman_pairs.csv", row.names = FALSE)

cat("\nDone.\n")

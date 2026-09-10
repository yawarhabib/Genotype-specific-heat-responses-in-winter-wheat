###############################################################################################
###   Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying   ###
###                    heat tolerance in six winter wheat cultivars                         ###
###                                       ----------                                        ###
###   Script 00 - Run the complete analysis (scripts 01 to 08) and record the session      ###
###                                       ----------                                        ###
###   Author: Yawar Habib                                                                   ###
###############################################################################################

######### DESCRIPTION #########################################################################

## Sources the eight analysis scripts in order. Each script is self-contained and can also
## be run on its own; this runner only fixes the order and records the R session used.
##
##   01_fresh_weight.R          Fig. 1   shoot fresh weight
##   02_gas_exchange.R          Fig. 2   Pn, gs, E, WUE
##   03_pigments.R              Fig. 3   Chl a+b, carotenoid, anthocyanin, Chl a/b
##   04_antioxidant_enzymes.R   Fig. 4   APX, DHAR, GR, GST
##   05_carbonyl_enzymes.R      Fig. 5   AKR, AOR/AER, GLY I and their z-score heatmap
##   06_PAL.R                   Fig. 6   PAL
##   07_polyamines.R            Fig. 7   SPM, SPD, CAD, PUT, DAP
##   08_multivariate_STI.R      Fig. 8, Fig. 9, Fig. S1   PCA, STI heatmap, Spearman matrix
##
## Input : data/ (15 csv files)
## Output: figures/ (10 tiff files), results/ (csv tables), results/sessionInfo.txt

######### INITIALISATION OF THE WORKING SPACE ##################################################

graphics.off()
rm(list = ls())

here::i_am("Scripts/00_run_all.R")
main_dir <- here::here()
setwd(main_dir)

######### RUN #################################################################################

scripts <- c("01_fresh_weight.R", "02_gas_exchange.R", "03_pigments.R",
             "04_antioxidant_enzymes.R", "05_carbonyl_enzymes.R", "06_PAL.R",
             "07_polyamines.R", "08_multivariate_STI.R")

for (s in scripts) {
  cat("\n\n############################  ", s, "  ############################\n\n", sep = "")
  source(file.path("Scripts", s), echo = FALSE)
  setwd(main_dir)                                   # every script sets the same root
}

######### SESSION INFORMATION #################################################################

## The scripts above clear the workspace; reload only what the session record needs
suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(purrr); library(stringr); library(tibble)
  library(ggplot2); library(ggrepel); library(patchwork); library(ggnewscale); library(ggh4x)
  library(car); library(emmeans); library(multcomp)
  library(ComplexHeatmap); library(circlize); library(corrplot)
})
writeLines(c(paste("Run on", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
             capture.output(sessionInfo())),
           "results/sessionInfo.txt")
cat("\nAll scripts finished. Session recorded in results/sessionInfo.txt\n")

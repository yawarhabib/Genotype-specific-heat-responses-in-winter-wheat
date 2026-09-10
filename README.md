# Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying heat tolerance in six winter wheat cultivars

R code and data for the statistical analysis and figures of the manuscript. Repository: https://github.com/yawarhabib/Genotype-specific-heat-responses-in-winter-wheat

> Rahman A., Habib Y., Jara Quispe A., Elguera J.P., Krecsák H., Khan I., Pál M., Balla K., Gechev T., Majláth I.
> *Genotype-specific antioxidant, reactive-carbonyl and polyamine responses underlying heat tolerance in six winter wheat cultivars* (manuscript under revision).

Six Martonvásár winter wheat cultivars (Mv Tarsoly, Mv Kolompos, Mv Lucilla, Mv Pirkadat, Mv Ikva, Mv Dandár) were grown at 20 °C (control), 35 °C and 40 °C and profiled for shoot biomass, leaf gas exchange, photosynthetic pigments, ascorbate–glutathione enzymes, reactive-carbonyl detoxification enzymes, phenylalanine ammonia-lyase and free polyamines. Every figure and every number quoted in the Results is produced by the scripts in this repository from the data files as provided.

## Repository layout

```
data/          raw measurements, one folder per trait group, one csv per file
Scripts/       numbered R scripts, one per figure (00 runs them all)
```

Running the scripts creates two further folders: `figures/` with the ten figures (300 dpi TIFF) and `results/` with the ANOVA tables, assumption checks, Tukey letters, means and the PCA, STI and correlation tables (csv), plus `sessionInfo.txt`.

### Data

| Folder | File(s) | Traits | Replicates per genotype × temperature |
|---|---|---|---|
| `Fresh_weight` | `FW.csv` | shoot fresh weight (g per plant) | 7 |
| `Gas_exchange` | `Leaf_gas_Exchange.csv` | Pn, E, gs, VPD, WUE | 3 |
| `Pigments` | `Chlorophyll_ab.csv`, `Carotenoid.csv`, `Chlorophyll_ab_ratio.csv`, `Anthocyanin.csv` | Chl a+b, carotenoid, Chl a/b, anthocyanin (20 and 40 °C only) | 4 (anthocyanin 3–4) |
| `Antioxidant_enzyme` | `APX.csv`, `DHAR.csv`, `GR.csv`, `GST.csv` | enzyme activities, nkat g⁻¹ FW | 4 (APX 4–5) |
| `Detoxification_enzyme` | `AKR.csv`, `AOR_AER.csv`, `GLY.csv` | enzyme activities, nkat g⁻¹ FW | 4 (35 °C: 5) |
| `PAL` | `PAL.csv` | PAL activity, nkat g⁻¹ FW | 4–5 |
| `Polyamines` | `PA_polyamines.csv` | DAP, PUT, CAD, SPD, SPM, µg g⁻¹ FW | 5 |

Every file has a `Genotype` column with the cultivar code (TAR, KOL, LUC, PIR, IKV, DAN) and a `Temperature` column (`20C`, `35C`, `40C`). Files whose figures show full cultivar names also carry them in a `Cultivar` column; the gas exchange file carries the full name in `Genotype` and the code in `Genotype_code`. One cadaverine value (Mv Ikva, 20 °C) is missing and is recorded as `NA`.

### Scripts

| Script | Figure | Content |
|---|---|---|
| `01_fresh_weight.R` | Fig. 1 | shoot fresh weight, radar plot |
| `02_gas_exchange.R` | Fig. 2 | Pn, gs, E, WUE, four radar panels |
| `03_pigments.R` | Fig. 3 | pigments, line plots |
| `04_antioxidant_enzymes.R` | Fig. 4 | APX, DHAR, GR, GST, violin panels |
| `05_carbonyl_enzymes.R` | Fig. 5 | AKR, AOR/AER, GLY I, violin panels and z-score heatmap |
| `06_PAL.R` | Fig. 6 | PAL, violin plot |
| `07_polyamines.R` | Fig. 7 | polyamines, faceted bar plot |
| `08_multivariate_STI.R` | Fig. 8, Fig. 9, Fig. S1 | PCA biplot, stress-tolerance-index heatmap, Spearman correlation matrix |
| `00_run_all.R` | — | runs 01 to 08 in order and writes `results/sessionInfo.txt` |

Each script starts with a header describing the model, the scale on which the trait was analysed and why, and its inputs and outputs. The statistical procedure is the same throughout: a two-way linear model (genotype × temperature), Type II ANOVA (`car::Anova`), Shapiro–Wilk and Levene checks on the raw scale and on the scale used, and Tukey-adjusted comparisons of the cultivars within each temperature on estimated marginal means (`emmeans`), summarised as compact letter displays in which `a` marks the highest mean.

## How to run

1. Clone or download the repository and open `Genotype-specific-heat-responses-in-winter-wheat.Rproj` in RStudio, or start R in the repository root. The scripts locate the root with the `here` package, so they can be run from any working directory inside the project.
2. Install the packages listed below if needed.
3. Run `Scripts/00_run_all.R`, or any single script. Figures are written to `figures/` and tables to `results/`.

The analysis was run with R 4.4.0 on Windows; `00_run_all.R` records the exact package versions of a run in `results/sessionInfo.txt`.

### Packages

`here`, `dplyr`, `tidyr`, `purrr`, `stringr`, `tibble`, `ggplot2`, `patchwork`, `ggrepel`, `ggnewscale`, `ggh4x`, `car`, `emmeans`, `multcomp`, `corrplot`, `circlize`, `ComplexHeatmap` (Bioconductor).

```r
install.packages(c("here", "dplyr", "tidyr", "purrr", "stringr", "tibble", "ggplot2",
                   "patchwork", "ggrepel", "ggnewscale", "ggh4x", "car", "emmeans",
                   "multcomp", "corrplot", "circlize", "BiocManager"))
BiocManager::install("ComplexHeatmap")
```

## Licence

Code: [MIT](LICENSE). Data (`data/`): [CC BY 4.0](data/LICENSE).

## Contact

Yawar Habib (analysis and code) — Center for Bio and Medical Technologies, Skolkovo Institute of Science and Technology, Moscow, Russia. y.habib@skoltech.ru · GitHub: [@yawarhabib](https://github.com/yawarhabib)

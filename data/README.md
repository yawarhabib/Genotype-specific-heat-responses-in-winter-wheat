# Data

Raw measurements of six Martonvásár winter wheat cultivars grown at 20 °C (control), 35 °C and 40 °C, as analysed in

Rahman A., Habib Y., Jara Quispe A., Elguera J.P., Krecsák H., Khan I., Pál M., Balla K., Gechev T., Majláth I. (2026). Genotype-specific antioxidant, reactive-carbonyl, and polyamine responses underlying heat tolerance in six winter wheat cultivars. *Frontiers in Plant Science* 17:1973386. https://doi.org/10.3389/fpls.2026.1973386

One folder per trait group, one csv per file, one row per plant. Licence: CC BY 4.0 (see `LICENSE` in this folder).

## Files

| File | Rows | Columns | Unit |
|---|---|---|---|
| `Fresh_weight/FW.csv` | 126 | Genotype, Cultivar, Temperature, FW_g | g per plant |
| `Gas_exchange/Leaf_gas_Exchange.csv` | 54 | Genotype, Genotype_code, Temperature, Pn, E, gs, VPD, WUE | Pn µmol m⁻² s⁻¹; E and gs mmol m⁻² s⁻¹; VPD kPa; WUE µmol mmol⁻¹ |
| `Pigments/Chlorophyll_ab.csv` | 72 | Genotype, Temperature, Chl_ab_ug_gFW | µg g⁻¹ FW |
| `Pigments/Carotenoid.csv` | 72 | Genotype, Temperature, Carotenoid_ug_gFW | µg g⁻¹ FW |
| `Pigments/Chlorophyll_ab_ratio.csv` | 72 | Genotype, Temperature, Chl_a_b_ratio | ratio |
| `Pigments/Anthocyanin.csv` | 45 | Genotype, Temperature, Anthocyanin_nmol_gFW | nmol g⁻¹ FW |
| `Antioxidant_enzyme/APX.csv` | 78 | Genotype, Cultivar, Temperature, APX_nkat_gFW | nkat g⁻¹ FW |
| `Antioxidant_enzyme/DHAR.csv` | 72 | Genotype, Cultivar, Temperature, DHAR_nkat_gFW | nkat g⁻¹ FW |
| `Antioxidant_enzyme/GR.csv` | 72 | Genotype, Cultivar, Temperature, GR_nkat_gFW | nkat g⁻¹ FW |
| `Antioxidant_enzyme/GST.csv` | 72 | Genotype, Cultivar, Temperature, GST_nkat_gFW | nkat g⁻¹ FW |
| `Detoxification_enzyme/AKR.csv` | 80 | Genotype, Cultivar, Temperature, AKR_nkat_gFW | nkat g⁻¹ FW |
| `Detoxification_enzyme/AOR_AER.csv` | 72 | Genotype, Cultivar, Temperature, AOR_AER_nkat_gFW | nkat g⁻¹ FW |
| `Detoxification_enzyme/GLY.csv` | 80 | Genotype, Cultivar, Temperature, GLY_nkat_gFW | nkat g⁻¹ FW |
| `PAL/PAL.csv` | 78 | Genotype, Cultivar, Temperature, PAL_nkat_gFW | nkat g⁻¹ FW |
| `Polyamines/PA_polyamines.csv` | 90 | Genotype, Temperature, DAP_ug_gFW, PUT_ug_gFW, CAD_ug_gFW, SPD_ug_gFW, SPM_ug_gFW | µg g⁻¹ FW |

## Columns

| Column | Content |
|---|---|
| `Genotype` | cultivar code: TAR (Mv Tarsoly), KOL (Mv Kolompos), LUC (Mv Lucilla), PIR (Mv Pirkadat), IKV (Mv Ikva), DAN (Mv Dandár). In the gas exchange file this column holds the full cultivar name and `Genotype_code` holds the code |
| `Cultivar` | full cultivar name, present in the files whose figures show full names |
| `Temperature` | growth temperature, `20C`, `35C`, `40C`; the gas exchange file codes it as `20`, `35`, `40` |
| trait columns | one measurement per plant, named as trait and unit, for example `APX_nkat_gFW` |

Abbreviations: Pn net photosynthetic rate, E transpiration rate, gs stomatal conductance, VPD vapour-pressure deficit, WUE water-use efficiency (Pn / E), APX ascorbate peroxidase, DHAR dehydroascorbate reductase, GR glutathione reductase, GST glutathione S-transferase, AKR aldo-keto reductase, AOR/AER alkenal/one reductase, GLY glyoxalase I, PAL phenylalanine ammonia-lyase, DAP 1,3-diaminopropane, PUT putrescine, CAD cadaverine, SPD spermidine, SPM spermine, FW fresh weight.

## Notes

- Anthocyanin was measured at 20 °C and 40 °C only.
- One cadaverine value (Mv Ikva, 20 °C) is missing and is recorded as `NA`.
- Replicates per genotype × temperature are listed in the repository README.

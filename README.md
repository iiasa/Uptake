# Climate Mitigation Pathways Reduce Wildfire Carbon Losses in European Forests

**UPTAKE Project, IIASA**

This repository provides a reproducible R workflow to quantify how climate‑mitigation pathways influence wildfire‑driven carbon losses in European forests (EU27+UK) and how these losses interact with harvested and energy biomass. The pipeline integrates **G4M** biomass stocks with **FLAM** burned area, builds annual forest‑type (afforested vs. managed) proportions, computes annual and spatial biomass loss (2015–2100), and generates publication‑ready figures and exportable CSVs.

---

## Table of Contents

- [Scope & Key Questions](#scope--key-questions)
- [Repository Structure](#repository-structure)
- [Data Inputs](#data-inputs)
- [Units & Conventions](#units--conventions)
- [Software Requirements](#software-requirements)
- [Quick Start](#quick-start)
- [Script Reference](#script-reference)
- [Outputs](#outputs)
- [Quality Control](#quality-control)
- [Reproducibility Notes](#reproducibility-notes)
- [Citing](#citing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## Scope & Key Questions

1. **Carbon impacts of fire under mitigation**  
   How do different mitigation pathways alter wildfire‑driven biomass (carbon) losses across EU27+UK?

2. **Allocation across forest types**  
   How are losses distributed between **afforested (AFF)** and **managed (MGD)** forests over time?

3. **Downstream implications**  
   What are the implications for **harvest volumes**, **residues**, and **energy biomass** when accounting for fire?

---

## Repository Structure

```
.
├─ 1_Uptake_data_sorting.R
├─ 2_Uptake_harvest_res_plot.R
├─ 3_Uptake_proportion_estimates.R
├─ 4_Uptake_burn_area_biomass_loss.R
├─ 5‑Final_plots.R
├─ 5_1‑Final_plots_spatial_maps_total_loss_remaining.R
├─ 5_2‑Final_plots_comparisons.R
├─ 5_3‑Final_plots_BIOMASS_BURNT_LOSS_REMAINING.R
├─ 5_4‑Final_plots_burned_area.R
└─ README.md
```

> **Run order:** `1 → 3 → 4 → (2) → 5`  
> Script **2** is independent (harvest/residues analytics) and can be run once inputs are available. The `5_*` helpers offer alternative figure sets.

---

## Data Inputs

- **G4M biomass tables (UPTAKE export)**  
  Files like:  
  `bioclimaDetails_t1_UPTAKE_31072025_Ref_Mit{2p6,4p5,7p0}_RCPref_NONE_Pco2_0.txt`

- **FLAM burned area (monthly GeoTIFFs, later annualized)**  
  Paths like:  
  `FLAM_Output/SSP126_A_burn_monthly/*.tif`, `.../SSP245...`, `.../SSP370...`

- **Scenario harvest/residues (UPTAKE CSV)**  
  `tabs_gui_ACCREU_v01_18062025_scen_csv_alln1.csv`

- **Country energy shares (Excel)**  
  `Energy_share_country.xlsx`

> Configure the root paths at the top of each script (`root_dir`, `g4m_dir`, `flam_dir`, `prop_dir`, `output_dir`). Replace system‑ or drive‑specific paths (e.g. `H:/…`) with your local setup.

---

## Units & Conventions

- **CRS:** EPSG:4326 (longitude/latitude)
- **Biomass inputs:** tC per hectare (tC ha⁻¹) from G4M; converted to **tC** when aggregated by forest area
- **Burned area:** Monthly FLAM rasters; aggregated to **annual km²** for time‑series analyses
- **Outputs:** Biomass loss in **tC**; large‑scale reports may use **MtC** or **GtC** with unit labels
- **Forest type proportions:** **AFF** and **MGD**, with `AFF + MGD = 1` within forested areas
- **Scenario mapping:** `Mit2p6 → SSP126`, `Mit4p5 → SSP245`, `Mit7p0 → SSP370`

---

## Software Requirements

- **R ≥ 4.2**
- **Required R packages:**  
  `dplyr, readr, stringr, purrr, tidyr, glue, tibble, data.table, readxl, raster, terra, sp, sf, rnaturalearth, exactextractr, ggplot2, ggtext, ggrepel, ggridges, patchwork, viridis, scales, gganimate, transformr`

```r
pkgs <- c(
  "dplyr", "readr", "stringr", "purrr", "tidyr", "glue", "tibble",
  "data.table", "readxl", "raster", "terra", "sp", "sf", "rnaturalearth",
  "exactextractr", "ggplot2", "ggtext", "ggrepel", "ggridges",
  "patchwork", "viridis", "scales", "gganimate", "transformr"
)
install.packages(setdiff(pkgs, installed.packages()[,1]))
```

- **GIS / raster handling notes**  
  - Proportion stacks should be resampled with **nearest neighbour** to maintain forest/non‑forest mask  
  - Biomass density rasters resampled **bilinearly** where needed to match loss grids

---

## Quick Start

1. Set paths in scripts (`root_dir`, `g4m_dir`, `flam_dir`, `prop_dir`, `output_dir`) to your local directories  
2. Run **`1_Uptake_data_sorting.R`** – process & aggregate G4M biomass data  
3. Run **`3_Uptake_proportion_estimates.R`** – generate annual AFF/MGD proportion stacks; perform sanity checks  
4. Run **`4_Uptake_burn_area_biomass_loss.R`** – compute yearly biomass loss via FLAM burned area & proportions  
5. *(Optional)* Run **`2_Uptake_harvest_res_plot.R`** for harvest & residues analyses  
6. Run **`5‑Final_plots.R`** (and/or `5_1`–`5_4`) – generate publication‑ready figures & CSVs

---

## Script Reference

| Script | Purpose | Key Inputs | Key Outputs |
|---|---|---|---|
| **1_Uptake_data_sorting.R** | Extract & aggregate G4M biomass (tC ha⁻¹), compute totals; optionally create per‑scenario rasters | G4M biomass tables | `G4M_biomass_totals_tc.csv`, scenario biomass rasters |
| **2_Uptake_harvest_res_plot.R** | Analyze harvest & residues for EU27+UK; produce relevant figures | Harvest/residues CSV | CSVs + PNG figures |
| **3_Uptake_proportion_estimates.R** | Build annual forest type proportion stacks (AFF vs MGD); check sums | G4M data + forest masks | AFF / MGD proportion stacks; summary stats & diagnostics |
| **4_Uptake_burn_area_biomass_loss.R** | Compute annual biomass loss by forest type using FLAM and proportion stacks | FLAM rasters, G4M biomass, proportion stacks | Loss rasters; annual results CSV; scenario comparisons |
| **5‑Final_plots.R** | Create main publication plots; export visualization data | Annual results, loss & biomass rasters | Final figures; CSVs for plotting |
| **5_1‑Final_plots_spatial_maps_total_loss_remaining.R** | Spatial maps: total losses and remaining biomass | Loss & biomass grids | High resolution spatial figure outputs |
| **5_2‑Final_plots_comparisons.R** | Compare scenarios over time | Scenario results | Comparative figures (time series etc.) |
| **5_3‑Final_plots_BIOMASS_BURNT_LOSS_REMAINING.R** | Combined views: burned area, loss, and residual biomass | Loss & biomass data | Composite visual summaries |
| **5_4‑Final_plots_burned_area.R** | Burned area time series & area summary figures | FLAM annualized burned area | Burned area plots & tables |

---

## Outputs

Typical output files and folder structure (customizable via `output_dir`):

```
Uptake/
├─ G4M_biomass_totals_tc.csv
├─ G4M_biomass_totals_GtC.csv
├─ {scenario}/Biomass_rasters/
│   ├─ Afforested_Biomass_<year>.tif
│   └─ Managed_Biomass_<year>.tif
├─ Area_Proportions/{scenario}/
│   ├─ {scenario}_afforested_proportion_stack.tif
│   ├─ {scenario}_managed_proportion_stack.tif
│   └─ Summary_stats & diagnostics
├─ Biomass_Loss_Results/{scenario}/
│   ├─ Afforested_Loss_<year>.tif
│   ├─ Managed_Loss_<year>.tif
│   ├─ Total_Loss_<year>.tif
│   ├─ annual_results.csv
│   └─ scenario comparison figures and CSVs
├─ Visualizations/
│   ├─ PNG figures
│   └─ CSVs for plotting
├─ Fire_Impact_Analysis/
│   └─ detailed result tables & supplementary figures
```

---

## Quality Control

- **Forest‑only normalization**: check that `AFF + MGD ≈ 1` in all forest grid cells; flag deviations beyond tolerance  
- **Mass balance checks**: ensure biomass loss never exceeds available biomass  
- **Unit consistency** in all outputs; units clearly labelled in figures & CSVs

---

## Reproducibility Notes

- All steps produce CSV exports of major intermediate / final results so figures can be reproduced without rerunning computationally heavy raster operations  
- Resampling and interpolation methods are fixed in scripts  
- Session and package versions should be recorded via `sessionInfo()` at end of runs  

---

## Citing

The correct author order is:

**Mahdi (Andre) Nakhavali, Johanna San‑Pedro, Andrey Krasovskiy, Fulvio di Fulvio, Mykola Gusti, Florian Kraxner, Petr Havlík**

Please cite this work as:

> Nakhavali, M. (Andre), San‑Pedro, J., Krasovskiy, A., di Fulvio, F., Gusti, M., Kraxner, F., & Havlík, P. (2025). *Climate Mitigation Pathways Reduce Wildfire Carbon Losses in European Forests*. UPTAKE Project, IIASA. GitHub repository.

BibTeX:

```bibtex
@misc{Nakhavali2025UptakeWildfireCarbonLoss,
  author       = {Nakhavali, Mahdi (Andre) and San‑Pedro, Johanna and Krasovskiy, Andrey and di Fulvio, Fulvio and Gusti, Mykola and Kraxner, Florian and Havlík, Petr},
  title        = {Climate Mitigation Pathways Reduce Wildfire Carbon Losses in European Forests},
  howpublished = {GitHub repository, UPTAKE Project, IIASA},
  year         = {2025},
  url          = {https://github.com/iiasa/Uptake}
}
```

---

## License

- **Code**: MIT  
- **Documentation & Figures**: CC BY 4.0  

Ensure a `LICENSE` file is present in the root.

---

## Acknowledgements

- Developed as part of **UPTAKE Project**, at **International Institute for Applied Systems Analysis (IIASA)**  
- Thanks to the G4M and FLAM model teams for data frameworks  
- Data inputs: EU27+UK forest shapefiles, country energy shares, harvest/residues scenario datasets  
- Visualization and GIS support via R packages and community tools  

# Climate Mitigation Pathways Reduce Wildfire Carbon Losses in European Forests
**UPTAKE Project · FLAM × G4M integration**

Quantifies how climate-mitigation pathways affect wildfire-driven carbon losses in European forests (EU27+UK) and downstream implications for harvested and energy biomass. The workflow integrates **G4M** biomass stocks with **FLAM** burned area, builds forest-type proportions, computes annual and spatial biomass loss (2015–2100), and generates publication-quality figures and CSVs for reproducibility.

---

## Table of Contents
- [Overview](#overview)
- [Repository Structure](#repository-structure)
- [Inputs & Data Sources](#inputs--data-sources)
- [Units & Conventions](#units--conventions)
- [Software Requirements](#software-requirements)
- [Quick Start](#quick-start)
- [Scripts](#scripts)
- [Outputs](#outputs)
- [Quality Control & Sanity Checks](#quality-control--sanity-checks)
- [Reproducibility Notes](#reproducibility-notes)
- [Citing](#citing)
- [License](#license)
- [Acknowledgements](#acknowledgements)

---

## Overview
This repository provides a reproducible R pipeline to estimate wildfire-induced biomass (carbon) losses under three mitigation pathways and to assess potential impacts on harvested and energy biomass supply. It covers:

- Reading and aggregating G4M biomass tables
- Constructing annual afforested/managed area proportions
- Splitting FLAM burned area by forest type
- Computing annual and spatial biomass loss (tC)
- Aggregating EU27+UK harvest and residue volumes
- Estimating before/after-fire energy biomass using country-specific energy shares
- Producing publication-ready figures and exportable CSVs

---

## Repository Structure
```
.
├─ 1_Uptake_data_sorting.R
├─ 2_Uptake_harvest_rest_plot.R
├─ 3_Uptake_proportion_estimates.R
├─ 4_Uptake_burn_area_biomass_loss.R
├─ 5-Final_plots.R
├─ 6_harvest_plots.R
└─ README.md
```

> **Run order:** 1 → 3 → 4 → (2) → 5 → 6  
> (Script 2 is independent for harvest/residues analytics and can be run any time after inputs are available.)

---

## Inputs & Data Sources
- **G4M biomass tables** (UPTAKE export):  
  `bioclimaDetails_t1_UPTAKE_31072025_Ref_Mit{2p6,4p5,7p0}_RCPref_NONE_Pco2_0.txt`
- **FLAM burned area** (monthly GeoTIFFs):  
  `FLAM_Output/SSP126_A_burn_monthly/*.tif`, `.../SSP245...`, `.../SSP370...`
- **Scenario harvest/residues** (UPTAKE CSV):  
  `tabs_gui_ACCREU_v01_18062025_scen_csv_alln1.csv`
- **Country energy shares** (Excel):  
  `Energy_share_country.xlsx`

> Paths are defined at the top of each script (e.g., `root_dir`, `g4m_dir`, `flam_dir`, `output_dir`).  
> Replace Windows-style `H:/...` with your environment paths as needed.

---

## Units & Conventions
- Coordinates/CRS: **EPSG:4326** (lon/lat).
- G4M biomass inputs: **tC ha⁻¹** (per-ha); converted to **tC** using forest area where appropriate.
- Burned area (FLAM): monthly rasters; aggregated to **annual km²** in analysis steps (scripts handle conversions).
- Biomass loss outputs: **tC**; figures often reported in **GtC** or **MtC** (see plot labels/CSV units).
- Forest type proportions: afforested (AFF) and managed (MGD); **AFF + MGD = 1** **inside forest**; **0** outside.
- Scenario mapping used in loss calculations:  
  `Mit2p6 → SSP126`, `Mit4p5 → SSP245`, `Mit7p0 → SSP370`.

---

## Software Requirements
- **R ≥ 4.2**
- Packages (superset across scripts):
  ```
  dplyr, readr, stringr, purrr, tidyr, glue, tibble, data.table, readxl,
  raster, terra, sp, sf, rnaturalearth, exactextractr,
  ggplot2, ggtext, ggrepel, ggridges, patchwork, viridis, scales,
  gganimate, transformr
  ```
Install:
```r
pkgs <- c(
  "dplyr","readr","stringr","purrr","tidyr","glue","tibble","data.table","readxl",
  "raster","terra","sp","sf","rnaturalearth","exactextractr",
  "ggplot2","ggtext","ggrepel","ggridges","patchwork","viridis","scales",
  "gganimate","transformr"
)
install.packages(setdiff(pkgs, installed.packages()[,1]))
```

**GIS Notes**
- Proportion stacks are resampled via **nearest** (to avoid mixing forest/non-forest).
- Biomass densities are resampled **bilinearly** where appropriate to match loss grids.

---

## Quick Start
1. **Configure paths** at the top of each script (`root_dir`, `g4m_dir`, `flam_dir`, `prop_dir`, `output_dir`).
2. Run **`1_Uptake_data_sorting.R`** to produce yearly totals and (optionally) per-year biomass rasters.
3. Run **`3_Uptake_proportion_estimates.R`** to create annual AFF/MGD proportion stacks (1990–2100) and sanity checks.
4. Run **`4_Uptake_burn_area_biomass_loss.R`** to compute annual biomass loss, save rasters, and create core figures.
5. (Optional) Run **`2_Uptake_harvest_rest_plot.R`** for harvest/residue analytics.
6. Run **`5-Final_plots.R`** to generate publication-ready visualizations and CSVs.
7. Run **`6_harvest_plots.R`** to estimate before/after-fire energy biomass with country-specific energy shares and produce summary plots.

---

## Scripts
| Script | Purpose | Key Inputs | Key Outputs |
|---|---|---|---|
| `1_Uptake_data_sorting.R` | Extract and aggregate G4M biomass; write totals; optional per-year rasters | G4M tables (`bioclimaDetails_...`) | `G4M_biomass_totals_tc.csv`, `G4M_biomass_totals_GtC.csv`, `<scenario>/*Biom_*.tif` |
| `2_Uptake_harvest_rest_plot.R` | EU27+UK harvest & residue analysis and figures | `tabs_gui_ACCREU_...csv` | `EU_harvest_residues_*` CSVs + PNGs |
| `3_Uptake_proportion_estimates.R` | Build annual AFF/MGD proportion stacks; enforce `AFF+MGD=1` inside forest | G4M tables | `<scenario>_afforested_proportion_stack.tif`, `<scenario>_managed_proportion_stack.tif`, sanity PNGs/CSV |
| `4_Uptake_burn_area_biomass_loss.R` | Compute annual biomass loss (tC) by splitting FLAM burn with AFF/MGD proportions; spatial & temporal plots | Proportion stacks, FLAM annualized burn, G4M rasters | `Aff/Managed/Total_{Burn, Loss}_<year>.tif`, `annual_results.csv`, multi-scenario PNGs |
| `5-Final_plots.R` | Publication-ready composite figures; exports plot CSVs | Scenario `annual_results.csv`, loss rasters | `Visualizations/*.png`, `Visualization_Data_CSV/*.csv` |
| `6_harvest_plots.R` | Impact on harvest & energy biomass (country energy shares) + QC | `EU_harvest_residues_selected_scenarios.csv`, loss rasters, G4M rasters, `Energy_share_country.xlsx` | `Fire_Impact_Analysis/*.png`, `*Results.csv` (detailed & summary) |

---

## Outputs
Default output roots (customizable in scripts):

```
H:/Uptake/
  G4M_biomass_totals_tc.csv
  G4M_biomass_totals_GtC.csv
  Mit{2p6,4p5,7p0}/*Biom_*.tif

H:/Uptake/Area_Proportions/<scenario>/
  <scenario>_afforested_proportion_stack.tif
  <scenario>_managed_proportion_stack.tif
  Summary_Stats_forestOnly.csv
  Sanity_*.png

H:/UPTAKE/Biomass_Loss_Results/
  <scenario>/Afforested_{Burn,Loss}_<year>.tif
  <scenario>/Managed_{Burn,Loss}_<year>.tif
  <scenario>/Total_{Burn,Loss}_<year>.tif
  <scenario>/annual_results.csv
  scenario_comparison.png
  cumulative_loss.png
  Visualizations/*.png
  Visualization_Data_CSV/*.csv

H:/UPTAKE/
  EU_harvest_residues_*.csv
  EU_biomass_comparison_creative.png

H:/UPTAKE/Biomass_Loss_Results/Fire_Impact_Analysis/
  Fire_Impact_Results.csv
  Detailed_Fire_Impact_Results.csv
  Summary_Fire_Impact_Results.csv
  Harvest_Volume_Impact.png
  Energy_Biomass_Difference.png
  Combined_Fire_Impact_Analysis.png
  Reduction_Percentage_Scatter.png
  Cumulative_Loss_Radial.png
```

---

## Quality Control & Sanity Checks
- **Forest-only normalization:** `3_Uptake_proportion_estimates.R` reports mean(AFF+MGD)≈1 inside forest; fraction of invalid cells (`|AFF+MGD−1| > tol`) by year.
- **Mass balance:** `6_harvest_plots.R` ensures fire loss **≤** total biomass and warns if outside expected ranges.
- **Unit consistency:** `5-Final_plots.R` applies explicit unit corrections in derived plots/CSVs; y-axes and captions state units.

---

## Reproducibility Notes
- All critical visualizations have matching **CSV exports** (`Visualization_Data_CSV` and `Fire_Impact_Analysis`).
- Randomness: none by default.
- Resampling methods are fixed and documented in code comments.
- Please track exact versions of input datasets and package session info for archival purposes:
  ```r
  sessionInfo()
  ```

---

## Citing
If you use this code or derived data/products, please cite:

> **Nakhavali, M.A.** et al. *Climate Mitigation Pathways Reduce Wildfire Carbon Losses in European Forests* (UPTAKE project). GitHub repository, 2025.  
> FLAM and G4M model references as per their respective documentation.

**BibTeX (template):**
```bibtex
@misc{Nakhavali2025UptakeWildfireCarbonLoss,
  author    = {Nakhavali, Mahdi (Andre) and Gusti, Mykola and di Fulvio, Fulvio and Krasovskiy, Andrey and San-Pedro, Johanna and Kraxner, Florian and Havl{'\i}k, Petr},
  title     = {Climate Mitigation Pathways Reduce Wildfire Carbon Losses in European Forests},
  howpublished = {GitHub},
  year      = {2025},
  note      = {UPTAKE project; FLAM × G4M integration},
  url       = {https://github.com/<org>/<repo>}
}
```

---

## License
Specify a license to enable reuse (recommended: MIT for code, CC-BY 4.0 for figures/derived CSVs).  
Example:
- Code: `MIT`
- Documentation & figures: `CC BY 4.0`

Add `LICENSE` and/or `LICENSE-Docs` files accordingly.

---

## Acknowledgements
- **UPTAKE Project** and contributing institutions.
- **G4M** (Global Forest Model) and **FLAM** teams for model frameworks and data.
- EU27+UK country shapefiles: *Natural Earth* via `{rnaturalearth}`.

For questions or collaboration, please open a GitHub Issue or contact the maintainers.

---

# Chinook Salmon Marine Consumption in Central-South Chile

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Code and data supporting the manuscript:

> **[Title]** — Hans [Apellido], [Co-autores]. *Journal of Applied Ecology*, [año].  
> DOI: [pending]

---

## Overview

This repository contains the R analysis pipeline for estimating the marine consumption of invasive Chinook salmon (*Oncorhynchus tshawytscha*) in central-south Chile. The pipeline integrates:

- Biological data (age, weight, energy density) from commercial fishing monitoring
- Stomach content and stable isotope analysis (?¹³C, ?¹?N, ?³?S) for diet reconstruction
- Bayesian diet mixing models (COSIMMR)
- Individual-level bioenergetic modeling (Fish Bioenergetics 4.0)
- Population scaling to estimate total prey consumption by river and spatial grid

Study rivers: **Toltén, Imperial, Valdivia, Río Bueno**  
Study period: **[años]**

---

## Repository structure

```
repo-chinook-consumption/
??? scripts/                        # R analysis scripts (run in numerical order)
?   ??? 00-extra-analysis.R         # Data cleaning and sample size summary
?   ??? 01-energy-density-weight-by-age.R
?   ??? 02-population-estimation-from-nb.R
?   ??? 03-age-structure-whole-population.R
?   ??? 04-total-population-estimation-ocean.R
?   ??? 05-diet-composition.R
?   ??? 06.0-descriptive-analysis-isotopes.R
?   ??? 06.1-descriptive-analysis-isotopes-2.R
?   ??? 06.2-mixing_models_cosimmr.R
?   ??? 06.3-reading-mixing-models.R
?   ??? 06.4-saving-diet-parameters-muscle.R
?   ??? 07-temperature-explore.R
?   ??? 08-age-consumption-model-individual-musculo.R
?   ??? 09-consumption-by-river-final-muscle.R
?   ??? 10-consumption-by-grid-final.R
?   ??? Graficos.R                  # All publication figures and tables
??? data/
?   ??? data_raw/
?   ?   ??? biological-data/        # Diet proportions, prey parameters, age composition
?   ?   ??? diet-isotopes/          # Stable isotope data and mixing model inputs
?   ?   ??? energy-density/         # Energy density measurements (salmon and prey)
?   ?   ??? river-population/       # Population estimates by river
?   ?   ??? bioenergetic-model/     # Individual consumption results (CSV)
?   ?   ??? population-consumption/ # Population-scaled consumption results
?   ?   ??? temperature/            # Daily SST time series (processed from Copernicus)
?   ??? spatial/                    # Coastal polygon (60 nm) and river points
??? README.md
```

> **Note:** Raw biological sampling data (`data_chinook_cleaned.csv`) and population data from
> Espinoza (2023) are not included in this repository pending data sharing agreements.
> See **Data availability** below.

---

## Analysis workflow

Run scripts in the following order from the repository root:

| Step | Script | Input | Output |
|------|--------|-------|--------|
| 1 | `00-extra-analysis.R` | Raw biological DB | `data_chinook_cleaned.csv` |
| 2 | `01-energy-density-weight-by-age.R` | `DE_salmones.csv` | `energy-density-weight-by-age.csv` |
| 3 | `02-population-estimation-from-nb.R` | `nb_rivers.csv`, `poblacion_espinoza.csv` | `total_population_from_nb.csv` |
| 4 | `04-total-population-estimation-ocean.R` | Step 3 output | `poblacion_mar_global.csv` |
| 5 | `05-diet-composition.R` | `data_chinook_cleaned.csv` | `diet_proportion_by_age.csv` |
| 6 | `06.0–06.4` (isotopes) | `stable-isotopes.csv` | `diet_proportion_by_age_ISOTOPES.csv` |
| 7 | `07-temperature-explore.R` | Copernicus NetCDF (see below) | `temperature_time_serie_daily.csv` |
| 8 | `08-age-consumption-model-individual-musculo.R` | Steps 1–7 outputs | `modelos_total_by_age_ISOTOPES.rds` |
| 9 | `09-consumption-by-river-final-muscle.R` | Step 8 + population | `consumo_poblacional_por_rio_detallado.csv` |
| 10 | `10-consumption-by-grid-final.R` | Steps 8–9 + temperature grids | `resultado_consumo_grillas_ISOTOPES.rds` |
| 11 | `Graficos.R` | All previous outputs | Figures and tables |

---

## Requirements

### R packages

Install all dependencies using `renv`:

```r
install.packages("renv")
renv::restore()
```

Key packages: `tidyverse`, `cosimmr`, `fb4package`, `terra`, `ncdf4`, `sf`, `ggplot2`, `TMB`

### Temperature data (Copernicus — not included)

Script `07-temperature-explore.R` requires GLORYS12 daily SST NetCDF files (2010–2024).  
Download free from: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030

Place files as: `data/data_raw/temperature/nc_files/mercator_YYYY.nc`

---

## Data availability

| Dataset | Status | Source |
|---------|--------|--------|
| Stable isotopes | Included (`data/data_raw/diet-isotopes/`) | This study |
| Energy density | Included (`data/data_raw/energy-density/`) | This study |
| Diet proportions | Included (`data/data_raw/biological-data/`) | Derived — this study |
| Population estimates (Nb) | Included (`data/data_raw/river-population/`) | This study |
| Population estimates (Espinoza) | Not included | Espinoza et al. (2023) — contact authors |
| Individual biological data | Not included | INVASAL monitoring program — contact [institucion] |
| Prey catch data (artisanal/industrial) | Not included | SERNAPESCA — available at [url] |
| SST temperature (GLORYS12) | Not included | Copernicus Marine Service — see above |
| Large model outputs (.rds) | Not included (reproducible) | Run pipeline from included inputs |

---

## Citation

If you use this code or data, please cite:

```
[Apellido], H. et al. ([año]). [Título]. Journal of Applied Ecology. DOI: [pending]
```

---

## License

Code: [MIT License](LICENSE)  
Data (included CSVs): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

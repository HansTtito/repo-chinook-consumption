# Bioenergetic assessment of non-native Chinook salmon predation on anchovy and sardine in south-central Chile

[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/)

Code and data supporting the manuscript:

> **Bioenergetic assessment of non-native Chinook salmon (*Oncorhynchus tshawytscha*) predation on anchovy and sardine in south-central Chile**
>
> Hans Ttito, Billy Ernst, David Beauchamp, Rachelle Johnson, Pablo Fierro, Felipe Docmac, Carolina Parada, Konrad Gorski, Jaime Tapia, Daniel Gomez-Uchida.
>
> *Manuscript in preparation*, 2026.

---

## Overview

We developed an age-structured bioenergetic model to estimate marine prey consumption by naturalized Chinook salmon (*Oncorhynchus tshawytscha*) populations associated with four major river systems in south-central Chile across nine years of marine residence (2014�2022). Over the study period, Chinook salmon consumed an estimated **94,515 tonnes** of prey, with sardine (*Strangomera bentincki*, 34.8%) and anchovy (*Engraulis ringens*, 22.8%) comprising the majority of prey biomass. Annual consumption varied 17-fold (1,297�22,114 t), driven primarily by interannual variation in reconstructed at-sea abundance.

The analysis integrates:

- Individual biological data (age, weight) from commercial fishing monitoring
- Stable isotope analysis (delta-13C, delta-15N, delta-34S) for diet reconstruction via Bayesian mixing models (COSIMMR)
- Energy density measurements fitted to nonlinear growth models
- Population estimates from mark-recapture (Nb) methods
- Individual-level bioenergetic modeling (Fish Bioenergetics 4.0, `fb4package`)
- Population scaling by river and spatial grid

Study rivers: **Tolten, Imperial, Valdivia, Rio Bueno**  
Study period: **marine residence 2014�2022; spawning seasons 2014/15 to 2022/23**

---

## Repository structure

```
repo-chinook-consumption/
??? scripts/                              # R analysis scripts (run in order)
?   ??? 01-energy-density-weight-by-age.R
?   ??? 02-population-estimation-from-nb.R
?   ??? 04-total-population-estimation-ocean.R
?   ??? 06.0-descriptive-analysis-isotopes.R
?   ??? 06.1-descriptive-analysis-isotopes-2.R
?   ??? 06.2-mixing_models_cosimmr.R
?   ??? 06.3-reading-mixing-models.R
?   ??? 06.4-saving-diet-parameters-muscle.R
?   ??? 07-temperature-explore.R
?   ??? 08-age-consumption-model-individual-musculo.R
?   ??? 09-consumption-by-river-final-muscle.R
?   ??? 10-consumption-by-grid-final.R
??? data/
?   ??? data_raw/
?   ?   ??? biological-data/    # chinook_model_inputs.csv, diet proportions, prey parameters
?   ?   ??? diet-isotopes/      # Stable isotope data (muscle, prey sources)
?   ?   ??? energy-density/     # Energy density measurements and NLS model outputs
?   ?   ??? river-population/   # Population estimates by river
?   ?   ??? bioenergetic-model/ # Individual consumption results (CSV)
?   ?   ??? population-consumption/ # Population-scaled consumption results
?   ?   ??? temperature/        # Daily SST time series (CSV; NetCDF not included)
?   ??? spatial/                # Coastal polygon (60 nm) and river points
??? README.md
??? LICENSE
??? .gitignore
```

---

## Analysis workflow

Most intermediate data files are **already provided as CSVs** in `data/`. This means you can run each script independently using the provided inputs, or run the full pipeline end to end.

> **Set your working directory to the repository root before running any script.**
> In R: `setwd("path/to/repo-chinook-consumption")`

| Step | Script | Key inputs | Key outputs | Pre-computed? |
|------|--------|-----------|-------------|---------------|
| 1 | `01-energy-density-weight-by-age.R` | `DE_salmones.csv`, `chinook_model_inputs.csv` | `energy-density-weight-by-age.csv` | Yes |
| 2 | `02-population-estimation-from-nb.R` | `nb_rivers.csv`, `poblacion_espinoza.csv` | `total_population_from_nb.csv` | Yes |
| 3 | `04-total-population-estimation-ocean.R` | `age-structure-population-by-river.csv` | `poblacion_mar_global.csv`, `poblacion_total_mar_por_rio.csv` | Yes |
| 4 | `06.0`, `06.1` | `stable-isotopes.csv` | Descriptive isotope statistics | -- |
| 5 | `06.2-mixing_models_cosimmr.R` | `stable-isotopes.csv`, `datos_rios.csv` | `cosimmr_model_8_muscle_5sources_cov.rds` (in `output/`) | **No -- must run** |
| 6 | `06.3-reading-mixing-models.R` | Output from step 5 | Model comparison diagnostics | -- |
| 7 | `06.4-saving-diet-parameters-muscle.R` | Output from step 5, `diet_proportion_by_age.csv` | `diet_proportion_by_age_ISOTOPES.csv` | Yes |
| 8 | `07-temperature-explore.R` | Copernicus GLORYS12 NetCDF (see below) | `temperature_time_serie_daily.csv` | Yes (CSV provided) |
| 9 | `08-age-consumption-model-individual-musculo.R` | Steps 1, 3, 7, 8 outputs + `chinook_model_inputs.csv` | `resultados_total_consumption_by_age_ISOTOPES.csv`, `modelos_total_by_age_ISOTOPES.rds` | Yes (CSV provided) |
| 10 | `09-consumption-by-river-final-muscle.R` | Steps 3 + 9 outputs | `consumo_poblacional_*_ISOTOPES.csv` | Yes |
| 11 | `10-consumption-by-grid-final.R` | Step 9 outputs + full temperature RDS | `resultado_consumo_grillas_ISOTOPES.rds` | Requires external data |

### Notes on the workflow

- **Step 5 (COSIMMR)** is the only step that must be run to regenerate its output. The `.rds` model file is not included in the repository (large binary file). After running `06.2`, the remaining scripts can use the saved model.
- **Step 11 (spatial grid)** requires the full processed temperature dataset (`total_data_temperature.rds`) generated from Copernicus NetCDF files (see below). The script is included for full transparency but this output is not strictly needed to reproduce the main consumption estimates.
- All other intermediate CSVs are provided; you can start the pipeline at any step.

---

## Requirements

### R packages

```r
# Core analysis
install.packages(c("tidyverse", "ggplot2", "janitor", "nlstools",
                   "cosimmr", "sf", "terra", "ncdf4", "patchwork",
                   "flextable", "officer", "box"))

# Fish Bioenergetics 4.0 (from GitHub)
remotes::install_github("Gaudeamus013/fb4package")
```

Key packages: `cosimmr` (Bayesian mixing models), `fb4package` (Fish Bioenergetics 4.0), `terra`/`ncdf4` (spatial/temperature processing), `sf` (spatial data).

### Temperature data (Copernicus � not included)

Script `07-temperature-explore.R` requires GLORYS12 daily SST NetCDF files (2010�2024).  
Download free (registration required) from:  
<https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030>

Place files as: `data/data_raw/temperature/nc_files/mercator_YYYY.nc`

The processed daily SST time series (`temperature_time_serie_daily.csv`) is already included, so **script 07 only needs to be re-run if you want to regenerate it from the raw NetCDF files**.

---

## Data availability

| Dataset | Status | Location | Notes |
|---------|--------|----------|-------|
| Biological model inputs (age, weight, season) | **Included** | `data_raw/biological-data/chinook_model_inputs.csv` | Reduced from full database (7 columns) |
| Diet proportions (stomach content) | **Included** | `data_raw/biological-data/diet_proportion_by_age.csv` | Derived - this study |
| Diet proportions (stable isotopes) | **Included** | `data_raw/biological-data/diet_proportion_by_age_ISOTOPES.csv` | Output of COSIMMR (script 06.4) |
| Stable isotope data | **Included** | `data_raw/diet-isotopes/stable-isotopes.csv` | Muscle tissue only |
| Energy density | **Included** | `data_raw/energy-density/` | Salmon and prey |
| Population estimates (Nb mark-recapture) | **Included** | `data_raw/river-population/` | This study |
| Population estimates (Espinoza) | **Included** | `data_raw/river-population/poblacion_espinoza.csv` | Espinoza Henriquez (2023) |
| SST time series (daily, CSV) | **Included** | `data_raw/temperature/temperature_time_serie_daily.csv` | Processed from Copernicus |
| Individual consumption results | **Included** | `data_raw/bioenergetic-model/` | Output of script 08 |
| Population consumption results | **Included** | `data_raw/population-consumption/` | Output of scripts 09, 10 |
| Full biological database | **Not included** | � | INVASAL monitoring program; contact corresponding author |
| SST NetCDF (GLORYS12) | **Not included** | � | Copernicus Marine Service � see above |
| Large model outputs (.rds) | **Not included** | � | Reproducible by running the pipeline |

---

## Citation

If you use this code or data, please cite:

```
Ttito, H., Ernst, B., Beauchamp, D., Johnson, R., Fierro, P., Docmac, F., Parada, C.,
Gorski, K., Tapia, J., Gomez-Uchida, D. (2026). Bioenergetic assessment of non-native
Chinook salmon (Oncorhynchus tshawytscha) predation on anchovy and sardine in
south-central Chile. Manuscript in preparation.
```

> This repository will be updated with the full citation once the manuscript is submitted and accepted.

---

## License

Code: [MIT License](LICENSE)  
Data (included CSVs): [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

---

## Contact

Hans Ttito � kvttitos@gmail.com  
Universidad de Concepci�n, Concepci�n, Chile

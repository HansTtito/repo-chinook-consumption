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

We developed an age-structured bioenergetic model to estimate marine prey consumption by naturalized Chinook salmon (*Oncorhynchus tshawytscha*) populations associated with four major river systems in south-central Chile across nine years of marine residence (2014-2022). Over the study period, Chinook salmon consumed an estimated **94,515 tonnes** of prey, with sardine (*Strangomera bentincki*, 34.8%) and anchovy (*Engraulis ringens*, 22.8%) comprising the majority of prey biomass. Annual consumption varied 17-fold (1,297-22,114 t), driven primarily by interannual variation in reconstructed at-sea abundance.

The pipeline integrates:

- Biological data (age, weight, energy density) from commercial fishing monitoring
- Stomach content and stable isotope analysis (delta-13C, delta-15N, delta-34S) for diet reconstruction
- Bayesian diet mixing models (COSIMMR)
- Individual-level bioenergetic modeling (Fish Bioenergetics 4.0)
- Population scaling to estimate total prey consumption by river and spatial grid

Study rivers: **Tolten, Imperial, Valdivia, Rio Bueno**
Study period: **2014-2022 (marine residence); spawning seasons 2014/15 to 2022/23**

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

> **Note:** Raw biological sampling data and population data from Espinoza (2023) are not
> included pending data sharing agreements. See **Data availability** below.
> The file `chinook_model_inputs.csv` (year, month, age, weight, species) provides the
> minimum biological inputs required to run the bioenergetic model (script 08).

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
| 6 | `06.0-06.4` (isotopes) | `stable-isotopes.csv` | `diet_proportion_by_age_ISOTOPES.csv` |
| 7 | `07-temperature-explore.R` | Copernicus NetCDF (see below) | `temperature_time_serie_daily.csv` |
| 8 | `08-age-consumption-model-individual-musculo.R` | Steps 1-7 outputs | `modelos_total_by_age_ISOTOPES.rds` |
| 9 | `09-consumption-by-river-final-muscle.R` | Step 8 + population | `consumo_poblacional_por_rio_detallado.csv` |
| 10 | `10-consumption-by-grid-final.R` | Steps 8-9 + temperature grids | `resultado_consumo_grillas_ISOTOPES.rds` |
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

### Temperature data (Copernicus - not included)

Script `07-temperature-explore.R` requires GLORYS12 daily SST NetCDF files (2010-2024).
Download free from: https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030

Place files as: `data/data_raw/temperature/nc_files/mercator_YYYY.nc`

---

## Data availability

| Dataset | Status | Source |
|---------|--------|--------|
| Stable isotopes | Included (`data/data_raw/diet-isotopes/`) | This study |
| Energy density | Included (`data/data_raw/energy-density/`) | This study |
| Diet proportions | Included (`data/data_raw/biological-data/`) | Derived - this study |
| Population estimates (Nb) | Included (`data/data_raw/river-population/`) | This study |
| Biological model inputs (age, weight) | Included (`chinook_model_inputs.csv`) | This study |
| Population estimates (Espinoza) | Not included | Espinoza et al. (2023) - contact authors |
| Full biological database | Not included | INVASAL monitoring program - contact corresponding author |
| Prey catch data (artisanal/industrial) | Not included | SERNAPESCA - available at sernapesca.cl |
| SST temperature (GLORYS12) | Not included | Copernicus Marine Service - see above |
| Large model outputs (.rds) | Not included (reproducible) | Run pipeline from included inputs |

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

Hans Ttito - kvttitos@gmail.com
Universidad de Concepcion, Concepcion, Chile

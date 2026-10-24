# EcopathIMTA

DOI: [10.5281/zenodo.23615487](https://doi.org/10.5281/zenodo.23615487)

Data and R code for:

Yong, Z.-J., Lin, W.-J., Wu, Y.-C., Lee, Y.-C., Lin, H.-J. Trophic structure and plant-based feed substitution in integrated multi-trophic aquaculture: An Ecopath modeling approach.

The code builds Rpath models of a seawater (mud crab, red algae) and a freshwater (sea bass, tilapia, duckweed) IMTA system in Taiwan (2024), runs the red algae and duckweed feed-substitution scenarios, and produces Tables 1-5, S1-S3 and S5 and Figures 1-3 and S2.

## Running the code

Tested with R 4.6.1. Package versions are in `renv.lock` (Rpath 1.1.0 and enaR 3.1.1 are installed from GitHub).

```r
renv::restore()
source("R/EcopathIMTA.R")
```

Results are written to `output/`. A full run takes about 10 minutes. `Sys.setenv(N_MC = 20)` gives a quicker run with fewer Monte Carlo draws for Table S3.

Two Rpath 1.1.0 functions are replaced in the code: `MTI()` reads biomass from `y$BB`, which does not exist, and `rpath()` adds every group's detrital input to every detritus group when computing detritus EE.

## Files

- `R/EcopathIMTA.R`: main script
- `R/00_load_inputs.R`: reads `data/`
- `R/EcopathIMTA_checks.R`: checks of the indices against a small hand-calculated model
- `R/EcopathIMTA_uncertainty.R`: parameter uncertainty (Table S3)
- `R/provenance/30_Ecopath_raw_data_processing.R`: how `data/` was derived from the raw records (not runnable without them)
- `data/ecopath_parameters.csv`: B, P/B, Q/B, catch and detrital input by group (g WW m-2, per culture period)
- `data/red_algae_duckweed_production.csv`: net production of red algae and duckweed from light-dark bottle incubations (mg O2 g-1 WW d-1)
- `data/ponds_2024.csv`: pond areas (m2)
- `data/derived_constants.csv`: culture period lengths and sea bass pond feed input

## Data availability

The raw field records belong to the Fisheries Research Institute, Ministry of Agriculture, Taiwan, and are available from the corresponding author on reasonable request with FRI's permission.

## License

Code: MIT. Data in `data/`: CC BY 4.0.

# Setup -----

source("R/00_R_preparation.R")
source("R/05a_load_raw_data.R")
source("R/05b_process_environmental.R")
source("R/05c_process_do.R")
source("R/05e_process_carbon.R")

## Constants -----

SW_AREA_2023 <- 1550
SW_AREA_2024 <- 2000
FW_T1_AREA <- 186
FW_T2_AREA <- 197
DAYLIGHT_HR <- 12
INVERT_PQ_RATIO <- 0.2

BREY_A  <- -0.473
BREY_B1 <- 1.007
BREY_B2 <- -0.274

AFDW_WW_POLYCHAETE <- 0.160
AFDW_WW_AMPHIPOD   <- 0.160
AFDW_WW_BIVALVE    <- 0.058
AFDW_WW_PROSOBRANCH <- 0.075
AFDW_WW_OPISTHOBRANCH <- 0.138
AFDW_WW_INSECT     <- 0.160
afdw_ww_family <- c(
  Capitellidae = AFDW_WW_POLYCHAETE, Nereididae = AFDW_WW_POLYCHAETE, Goniadidae = AFDW_WW_POLYCHAETE,
  Syllidae = AFDW_WW_POLYCHAETE, Sabellidae = AFDW_WW_POLYCHAETE,
  Corophiidae = AFDW_WW_AMPHIPOD,
  Potamididae = AFDW_WW_PROSOBRANCH,
  Scaphandridae = AFDW_WW_OPISTHOBRANCH, Cylichnidae = AFDW_WW_OPISTHOBRANCH,
  Tellinidae = AFDW_WW_BIVALVE, Laternulidae = AFDW_WW_BIVALVE,
  Chironomidae = AFDW_WW_INSECT, Hydrophilidae = AFDW_WW_INSECT
)
QB_MUD_CRAB <- 11.5
QB_SEA_BASS <- 3.35
QB_TILAPIA <- 24.53
ZOO_PB <- 67
ZOO_QB <- 280

## Helpers -----

calc_new_site <- function(water, group) {

  case_when(
    water == "sw" ~ "sw_T",
    water == "fw" & group == "T1" ~ "fw_T1",
    water == "fw" & group == "T2" ~ "fw_T2"
  )
}

separate_site_columns <- function(df) {
  df |>
    separate_wider_delim(cols = site, delim = "_",
                         names = c("water", "type", "species", "group"),
                         cols_remove = FALSE)
}

split_light_dark <- function(df, rate_col, group_cols) {
  df_light <- df |>
    filter(condition == "light") |>
    mutate(nep_hr = .data[[rate_col]]) |>
    select(-condition)

  df_dark <- df |>
    filter(condition == "dark") |>
    mutate(r_hr = if_else(.data[[rate_col]] > 0, 0, .data[[rate_col]])) |>
    select(all_of(c(group_cols, "rep", "r_hr")))

  df_respiration <- df_dark |>
    summarise(mean_r_hr = mean(r_hr, na.rm = TRUE), .by = all_of(group_cols))

  df_light |>
    left_join(df_respiration, join_by(!!!syms(group_cols))) |>
    left_join(df_dark, join_by(!!!syms(c(group_cols, "rep"))))
}

calc_daily_gpp_nep <- function(df, daylight_hr = DAYLIGHT_HR) {
  df |> mutate(
    daylight_hr = daylight_hr,
    nep_daytime = nep_hr * daylight_hr,
    r_daytime = mean_r_hr * daylight_hr,
    r_day = mean_r_hr * 24,
    gpp = nep_daytime - r_daytime,
    gpp = if_else(gpp < 0, 0, gpp),
    r_day = if_else(r_day > 0, 0, r_day),
    nep = gpp + r_day
  )
}

summarise_producer <- function(df) {
  df |>
    separate_site_columns() |>
    summarise(
      biomass      = mean(biomass, na.rm = TRUE),
      production_c = mean(production_c, na.rm = TRUE),
      biomass_c    = mean(biomass_c, na.rm = TRUE),
      area_m2      = mean(area_m2, na.rm = TRUE),
      .by = c(year, sampling, water, species, group)
    ) |>
    mutate(pb_ratio = production_c / biomass_c) |>
    summarise(
      pb_ratio = mean(pb_ratio, na.rm = TRUE),
      biomass  = mean(biomass, na.rm = TRUE),
      area_m2  = mean(area_m2, na.rm = TRUE),
      .by = c(year, water, species, group)
    )
}

calc_weighted_summary <- function(df) {
  df |>
    mutate(new_site = calc_new_site(water, group)) |>
    mutate(
      total_area_m2 = sum(area_m2),
      weighted_biomass = (area_m2 / total_area_m2) * biomass,
      .by = c(year, water, new_site)
    ) |>
    summarise(
      avg_weighted_biomass = sum(weighted_biomass, na.rm = TRUE),
      avg_weighted_pb = mean(pb_ratio, na.rm = TRUE),
      .by = c(year, water, new_site)
    ) |>
    as.data.frame()
}

## Data loading -----

info_area <- info |>
  filter(type == "imta") |>
  select(year, site, area_m2) |>
  mutate(year = as.character(year))

df_bod <- read.xlsx("data/raw/raw_bod.xlsx")

df_bod <- df_bod |>
  mutate(
    sampling_date = as.Date(sampling_date, origin = "1899-12-30"),
    year = as.character(lubridate::year(sampling_date)),
    site = if_else(site == "fw_imta_plant", "fw_imta_plant_T1", site)
  )

df_info <- info |>
  filter(type == "imta") |>
  select(year, site, depth_m, area_m2, stock_animals_ind, day_culture,
         stock_animals_total_kg, harvest_animals_total_kg) |>
  mutate(
    year = as.character(year),
    day_culture = as.numeric(day_culture),
    stock_animals_total_kg = as.numeric(stock_animals_total_kg),
    harvest_animals_total_kg = as.numeric(harvest_animals_total_kg)
  )

df_bod <- df_bod |>
  left_join(df_info, join_by(year, site))

# Culture period conversion -----

.pond_days <- info |>
  filter(type == "imta", year == 2024, group %in% c("T", "T1", "T2")) |>
  transmute(
    water, species,
    stock_date   = as.Date(as.numeric(stock_date),   origin = "1899-12-30"),
    harvest_date = as.Date(as.numeric(harvest_date), origin = "1899-12-30"),
    day_culture  = as.numeric(day_culture)
  )

envelope_days <- .pond_days |>
  summarise(envelope = as.integer(max(harvest_date) - min(stock_date)), .by = water)

group_days <- .pond_days |>
  mutate(func_group = case_when(
    water == "sw" & species == "crab"  ~ "Mud_crab",
    water == "fw" & species == "fish"  ~ "Sea_bass",
    water == "fw" & species == "plant" ~ "Tilapia"
  )) |>
  filter(!is.na(func_group)) |>
  summarise(group_yr = mean(day_culture) / 365, .by = c(water, func_group))

stopifnot(nrow(group_days) == 3)

culture_conv <- envelope_days |>
  mutate(day_culture = envelope, culture_yr = envelope / 365)

animal_conv <- group_days

# Consumers -----

## Cultured animals -----

### mud crab -----

sw_imta_crab_T_2024 <- tibble::tribble(
  ~site,           ~sampling_date, ~sampling_days, ~mean_weight_g,
  "sw_imta_crab_T", "2024-04-22",   0,               0.1322,
  "sw_imta_crab_T", "2024-05-23",   31,              9.133333,
  "sw_imta_crab_T", "2024-06-20",   60,              29.4,
  "sw_imta_crab_T", "2024-07-19",   89,              68.325523,
  "sw_imta_crab_T", "2024-08-15",   116,             116.425,
  "sw_imta_crab_T", "2024-09-11",   143,             211.4,
  "sw_imta_crab_T", "2024-10-17",   176,             269.755
) |>
  dplyr::mutate(sampling_date = as.Date(sampling_date),
                weight_gain_per_day = (mean_weight_g - dplyr::lag(mean_weight_g)) /
                  (sampling_days - dplyr::lag(sampling_days))) |>
  summarise(ind_weight_g = mean(mean_weight_g),
            ind_weight_gain_per_day = mean(weight_gain_per_day, na.rm = T),
            .by = "site") |>
  left_join(df_info |>
              filter(year == "2024", site == "sw_imta_crab_T") |>
              select(year, site, day_culture, stock_animals_ind, area_m2,
                     stock_animals_total_kg, harvest_animals_total_kg), join_by(site)) |>
  mutate(total_biomass_g = ind_weight_g * stock_animals_ind,
         biomass = total_biomass_g / area_m2,
         production = (harvest_animals_total_kg - stock_animals_total_kg) * kg_to_g / area_m2,
         pb_ratio = production / biomass
         )

### tilapia/seabass -----

fw_imta_fish_T1 <- tibble::tribble(
  ~site,                  ~sampling_date, ~sampling_days, ~mean_weight_g,
  "fw_imta_fish_T1",   "2024-03-27",   0,              10.06,
  "fw_imta_fish_T1",   "2024-04-22",   26,             21.133,
  "fw_imta_fish_T1",   "2024-05-30",   64,             74.85,
  "fw_imta_fish_T1",   "2024-06-25",   90,             158.1,
  "fw_imta_fish_T1",   "2024-07-31",   126,            268.467,
  "fw_imta_fish_T1",   "2024-08-28",   154,            390.767,
  "fw_imta_fish_T1",   "2024-10-22",   209,            629.117
) |>
  dplyr::mutate(
    sampling_date       = as.Date(sampling_date),
    weight_gain_per_day = (mean_weight_g - dplyr::lag(mean_weight_g)) / (sampling_days - dplyr::lag(sampling_days))
  )

fw_imta_fish_T2 <- tibble::tribble(
  ~site,                  ~sampling_date, ~sampling_days, ~mean_weight_g,
  "fw_imta_fish_T2",   "2024-03-27",   0,              10.06,
  "fw_imta_fish_T2",   "2024-04-22",   26,             21.95,
  "fw_imta_fish_T2",   "2024-05-30",   64,             80.767,
  "fw_imta_fish_T2",   "2024-06-25",   90,             156.667,
  "fw_imta_fish_T2",   "2024-07-31",   126,            298.467,
  "fw_imta_fish_T2",   "2024-08-28",   154,            338.617,
  "fw_imta_fish_T2",   "2024-10-22",   209,            621.917
) |>
  dplyr::mutate(
    sampling_date       = as.Date(sampling_date),
    weight_gain_per_day = (mean_weight_g - dplyr::lag(mean_weight_g)) / (sampling_days - dplyr::lag(sampling_days))
  )

fw_imta_plant_T1 <- tibble::tribble(
  ~site,                  ~sampling_date, ~sampling_days, ~mean_weight_g,
  "fw_imta_plant_T1",   "2024-03-29",   0,              0.82,
  "fw_imta_plant_T1",   "2024-04-22",   24,             8.261,
  "fw_imta_plant_T1",   "2024-05-30",   62,             51.05,
  "fw_imta_plant_T1",   "2024-06-25",   88,             105.85,
  "fw_imta_plant_T1",   "2024-07-31",   124,            200.183,
  "fw_imta_plant_T1",   "2024-08-28",   152,            280.617,
  "fw_imta_plant_T1",   "2024-10-22",   207,            453.4
) |>
  dplyr::mutate(
    sampling_date       = as.Date(sampling_date),
    weight_gain_per_day = (mean_weight_g - dplyr::lag(mean_weight_g)) / (sampling_days - dplyr::lag(sampling_days))
  )

fw_imta_plant_T2 <- tibble::tribble(
  ~site,                  ~sampling_date, ~sampling_days, ~mean_weight_g,
  "fw_imta_plant_T2",   "2024-03-29",   0,              0.82,
  "fw_imta_plant_T2",   "2024-04-22",   24,             9.361,
  "fw_imta_plant_T2",   "2024-05-30",   62,             56.683,
  "fw_imta_plant_T2",   "2024-06-25",   88,             105.696,
  "fw_imta_plant_T2",   "2024-07-31",   124,            199.667,
  "fw_imta_plant_T2",   "2024-08-28",   152,            298.033,
  "fw_imta_plant_T2",   "2024-10-22",   207,            524.317
) |>
  dplyr::mutate(
    sampling_date       = as.Date(sampling_date),
    weight_gain_per_day = (mean_weight_g - dplyr::lag(mean_weight_g)) / (sampling_days - dplyr::lag(sampling_days))
  )

fw_imta_2024 <- rbind(fw_imta_fish_T1, fw_imta_fish_T2, fw_imta_plant_T1, fw_imta_plant_T2)

fw_imta_2024 <- fw_imta_2024 |>
  summarise(ind_weight_g = mean(mean_weight_g),
                            ind_weight_gain_per_day = mean(weight_gain_per_day, na.rm = T),
                            .by = "site") |>
  left_join(df_info |>
              filter(year == "2024", str_detect(site, "fw_imta_")) |>
              select(year, site, day_culture, stock_animals_ind, area_m2,
                     stock_animals_total_kg, harvest_animals_total_kg), join_by(site)) |>
  mutate(total_biomass_g = ind_weight_g * stock_animals_ind,
         biomass = total_biomass_g / area_m2,
         production = (harvest_animals_total_kg - stock_animals_total_kg) * kg_to_g / area_m2,
         pb_ratio = production / biomass
  )

### all -----

df_animals <- info |>
  filter(year == "2024", type == "imta", group %in% c("T", "T1", "T2")) |>
  select(year, water, species, group, site,
         culture_species, area_m2, contains("stock"), contains("harvest"))

animals_bpb <- df_animals |>
  mutate(
    new_site = calc_new_site(water, group),
    biomass_kg = (as.numeric(stock_animals_total_kg) + as.numeric(harvest_animals_total_kg)) / 2,
    biomass_g = biomass_kg * kg_to_g,
    biomass = biomass_g / area_m2,

    biomass = case_when(
      water == "sw" & species == "crab" ~ sw_imta_crab_T_2024$biomass,
      water == "fw" & species == "fish" & group == "T1" ~ fw_imta_2024$biomass[1],
      water == "fw" & species == "fish" & group == "T2" ~ fw_imta_2024$biomass[2],
      water == "fw" & species == "plant" & group == "T1" ~ fw_imta_2024$biomass[3],
      water == "fw" & species == "plant" & group == "T2" ~ fw_imta_2024$biomass[4],
      TRUE ~ biomass),

    production_kg = as.numeric(harvest_animals_total_kg) - as.numeric(stock_animals_total_kg),
    production_g = production_kg * kg_to_g,
    pb_ratio = production_g / biomass_g,

    pb_ratio = case_when(
      water == "sw" & species == "crab" ~ sw_imta_crab_T_2024$pb_ratio,
      water == "fw" & species == "fish" & group == "T1" ~ fw_imta_2024$pb_ratio[1],
      water == "fw" & species == "fish" & group == "T2" ~ fw_imta_2024$pb_ratio[2],
      water == "fw" & species == "plant" & group == "T1" ~ fw_imta_2024$pb_ratio[3],
      water == "fw" & species == "plant" & group == "T2" ~ fw_imta_2024$pb_ratio[4],
      TRUE ~ pb_ratio),

    biomass = ifelse(water == "sw" & species == "plant", 0, biomass),
    pb_ratio = ifelse(water == "sw" & species == "plant", 0, pb_ratio),

    catch = (as.numeric(harvest_animals_total_kg) - as.numeric(stock_animals_total_kg)) * kg_to_g / area_m2
  ) |>
  select(year, water, species, new_site, biomass, pb_ratio, catch, area_m2) |>
  as.data.frame()

tab_animals <- animals_bpb |>
  mutate(
    total_area_m2 = sum(area_m2),
    weighted_biomass = (area_m2 / total_area_m2) * biomass,
    weighted_catch = (area_m2 / total_area_m2) * catch,
    .by = c(year, water, new_site)
  ) |>
  filter(year != "2023", biomass != 0)

## Invertebrates -----

df_benthic <- read.xlsx("data/raw/raw_invertebrates.xlsx", detectDates = TRUE)

benthic_corer_radius_m <- 5 / 100

df_benthic <- df_benthic |>
  separate_site_columns() |>
  mutate(
    year = as.character(year),
    w_biomass = w_dry - w_plate,
    sampling_area_m2 = pi * (benthic_corer_radius_m^2),
    biomass = (w_biomass * dw_to_ww_organicTissue) / sampling_area_m2,
    biomass_unit = "g WW m-2",
    func_groups = case_when(
      family == "Capitellidae" ~ "Detritivores",
      family == "Nereididae" ~ "Omnivores",
      family == "Goniadidae" ~ "Carnivores",
      family == "Chironomidae" ~ "Omnivores",
      family == "Syllidae" ~ "Omnivores",
      family == "Corophiidae" ~ "Detritivores",
      family == "Scaphandridae" ~ "Carnivores",
      family == "Potamididae" ~ "Detritivores",
      family == "Hydrophilidae" ~ "Carnivores",
      family == "Sabellidae" ~ "Filter feeders",
      family == "Laternulidae" ~ "Filter feeders",
      family == "Cylichnidae" ~ "Carnivores",
      family == "Tellinidae" ~ "Detritivores",
      .default = NA_character_
    ),
    abundance = n / sampling_area_m2,
    afdw_ww = unname(afdw_ww_family[family]),
    biomass_afdw = biomass * afdw_ww,
    w_ind_afdw = biomass_afdw / abundance,
    pb_brey = 10^BREY_A * biomass_afdw^(BREY_B1 - 1) * w_ind_afdw^BREY_B2,
    production = pb_brey * biomass
  )
stopifnot(!anyNA(df_benthic$afdw_ww))

df_benthic <- df_benthic |>
  left_join(df_info |> select(year, site, day_culture), join_by(year, site)) |>
  left_join(culture_conv |> select(water, envelope_yr = culture_yr), join_by(water)) |>
  mutate(culture_yr = if_else(water == "sw", envelope_yr, day_culture / 365),
         production = production * culture_yr) |>
  select(-envelope_yr)

df_benthic_weighted <- df_benthic |>
  summarise(
    total_weighted_biomass = sum(biomass),
    total_weighted_production = sum(production),
    .by = c(year, sampling_date, site, rep, func_groups)
  )

stopifnot(all(df_benthic_weighted$rep %in% 1:3))
corer_grid <- merge(
  distinct(df_benthic_weighted, year, sampling_date, site),
  distinct(df_benthic_weighted, year, site, func_groups),
  by = c("year", "site")
) |>
  merge(data.frame(rep = 1:3))
df_benthic_weighted <- corer_grid |>
  left_join(df_benthic_weighted, by = c("year", "sampling_date", "site", "func_groups", "rep")) |>
  mutate(
    total_weighted_biomass = coalesce(total_weighted_biomass, 0),
    total_weighted_production = coalesce(total_weighted_production, 0)
  )

df_benthic_avg <- df_benthic_weighted |>
  mutate(
    pb_ratio = ifelse(total_weighted_biomass > 0, total_weighted_production / total_weighted_biomass, NA_real_),
    pq_ratio = INVERT_PQ_RATIO,
    qb_ratio = pb_ratio * (1 / pq_ratio)
  ) |>
  rename(biomass = total_weighted_biomass) |>
  summarise(
    biomass = mean(biomass),
    pb_ratio = mean(pb_ratio, na.rm = TRUE),
    qb_ratio = mean(qb_ratio, na.rm = TRUE),
    .by = c(year, sampling_date, site, func_groups)
  ) |>
  summarise(
    biomass = mean(biomass),
    pb_ratio = mean(pb_ratio, na.rm = TRUE),
    qb_ratio = mean(qb_ratio, na.rm = TRUE),
    .by = c(year, site, func_groups)
  ) |>
  separate_site_columns() |>
  left_join(info_area, join_by(year, site))

tab_invert <- df_benthic_avg |>
  mutate(new_site = calc_new_site(water, group)) |>
  mutate(
    total_area_m2 = sum(area_m2),
    weighted_biomass = (area_m2 / total_area_m2) * biomass,
    .by = c(year, water, new_site, func_groups)
  ) |>
  summarise(
    biomass = sum(weighted_biomass),
    pb_ratio = mean(pb_ratio),
    qb_ratio = mean(qb_ratio),
    .by = c(year, water, new_site, func_groups)
  ) |>
  as.data.frame()

## Zooplankton -----

df_zoo <- read.xlsx("data/raw/raw_zooplankton.xlsx") |>
  mutate(year = as.character(year),
         sampling_date = as.Date(sampling_date, origin = "1899-12-30"))

df_zoo <- df_zoo |>
  left_join(df_info, join_by(year, site))

zoo_bm <- df_zoo |>
  mutate(
    zooplankton_ww_g = wet_g - mesh_weight_g,
    zooplankton_dw_g = dry_g - plate_g,
    filter_volume_m3 = filter_volume_L / 1000,
    pond_volume_m3 = area_m2 * depth_m,
    zoo_total_ww_g = (zooplankton_ww_g / filter_volume_m3) * pond_volume_m3,
    zoo_total_dw_g = (zooplankton_dw_g / filter_volume_m3) * pond_volume_m3,
    biomass_ww = zoo_total_ww_g / area_m2,
    biomass_dw = zoo_total_dw_g / area_m2
  ) |>
  summarise(
    biomass_ww = mean(biomass_ww, na.rm = TRUE),
    biomass_dw = mean(biomass_dw, na.rm = TRUE),
    .by = c(year, sampling_date, site)
  ) |>
  mutate(biomass_ww_from_dw = biomass_dw * dw_to_ww_organicTissue)

zoo_bm <- zoo_bm |>
  separate_site_columns() |>
  left_join(info_area, join_by(year, site))

tab_zooplankton <- zoo_bm |>
  mutate(new_site = calc_new_site(water, group)) |>
  mutate(
    total_area_m2 = sum(area_m2),
    weighted_biomass = (area_m2 / total_area_m2) * biomass_ww,
    .by = c(year, water, new_site)
  ) |>
  summarise(avg_weighted_biomass = sum(weighted_biomass), .by = c(year, water, new_site)) |>
  as.data.frame()

# Primary producers -----

## Phytoplankton -----

filter_chla <- function(input_water, input_year, input_sampling, input_sampling_new) {
  envr |>
    filter(year == input_year, water == input_water, type == "imta", sampling == input_sampling) |>
    mutate(sampling = input_sampling_new, func_group = "phytoplankton") |>
    summarise(chla_conc = mean(envr_chla, na.rm = TRUE), .by = c(year, site, sampling, func_group))
}

phyto_chla_from_cb <- do.call(rbind, lapply(list(
  c("sw", "2023", "3rd", "pre"),
  c("sw", "2023", "6th", "1st"),
  c("sw", "2023", "8th", "2nd"),
  c("sw", "2024", "7th", "1st"),
  c("sw", "2024", "8th", "2nd"),
  c("fw", "2023", "3rd", "pre"),
  c("fw", "2023", "4th", "1st"),
  c("fw", "2023", "5th", "2nd"),
  c("fw", "2024", "6th", "1st")
), function(params) {
  filter_chla(params[1], params[2], params[3], params[4])
}))

phyto_chla_rest <- process_chla |>
  filter(str_detect(sampling, "food web"), func_group == "phytoplankton") |>
  unite(col = "site", water, type, species, group, sep = "_", remove = FALSE) |>
  mutate(
    chla_content_µgml = calc_chla_trichromatic(wl_630, wl_647, wl_664, wl_750),
    chla_conc = (chla_content_µgml * acetone_mL) / filter_volume,
    sampling = case_when(
      sampling == "food web_1" ~ "1st",
      sampling == "food web_2" ~ "2nd",
      sampling == "food web_3" ~ "3rd"
    ),
    year = as.character(year)
  ) |>
  summarise(chla_conc = mean(chla_conc, na.rm = TRUE), .by = c(year, site, sampling, func_group))

phyto_chla <- bind_rows(phyto_chla_from_cb, phyto_chla_rest) |> filter(chla_conc >= 0)
rm(phyto_chla_rest, phyto_chla_from_cb)

phyto <- df_bod |>
  filter(func_group == "phytoplankton") |>
  mutate(
    volume_L = 0.3,
    do_change_mgL = end_do_mgL - start_do_mgL,
    incubation_hr = as.numeric(end_time - start_time) * 24,
    total_do_change_mg = do_change_mgL * volume_L,
    do_change_rate = do_change_mgL / incubation_hr,
    rate_unit = "mg O2 L-1 hr-1"
  ) |>
  left_join(phyto_chla, join_by(year, site, sampling, func_group))

phyto_light <- split_light_dark(phyto, "do_change_rate",
                                c("year", "sampling", "sampling_date", "site", "func_group"))

phyto_light |>
  calc_daily_gpp_nep() |>
  filter(nep < 0)

phyto_do <- phyto_light |>
  calc_daily_gpp_nep() |>
  mutate(
    nep = pmax(nep, 0),
    production_co2 = nep * mg_to_g * (1/o2_molar_mass) * (1/pq),
    season_days = culture_conv$day_culture[match(substr(site, 1, 2), culture_conv$water)],
    production_c = (production_co2 * co2_mass * co2e_to_c * depth_m * season_days) / L_to_m3,
    biomass_c = (chla_conc * phyto_chla_to_c * µg_to_g * depth_m) / L_to_m3,
    biomass = biomass_c * c_to_dw_producer * dw_to_ww_producer,
    production_unit = "g C m-2 culture-period-1",
    biomass_unit = "g WW m-2"
  )

phyto_do_change <- phyto |>
  summarise(phyto_do_rate_mgL_hr = mean(do_change_mgL / incubation_hr), .by = c(year, sampling, sampling_date, site, condition))

df_phyto <- summarise_producer(phyto_do)
tab_phyto <- calc_weighted_summary(df_phyto)

## Benthic algae -----

benthic_chla_2023 <- df_bod |>
  filter(year == "2023", func_group == "benthic") |>
  mutate(
    acetone_mL = 10,
    sample_area_cm2 = pi * (0.45)^2 * 2,
    chla_content_µgml = calc_chla_trichromatic(wl_630, wl_647, wl_664, wl_750),
    benthic_chla_µg = chla_content_µgml * acetone_mL,
    chla_conc = benthic_chla_µg / sample_area_cm2
  ) |>
  select(year, site, chla_conc, sampling, func_group, rep, condition)

benthic_chla_2024 <- process_chla |>
  filter(str_detect(sampling, "food web"), func_group == "benthic") |>
  unite(col = "site", water, type, species, group, sep = "_", remove = FALSE) |>
  mutate(
    year = as.character(year),
    chla_content_µgml = calc_chla_trichromatic(wl_630, wl_647, wl_664, wl_750),
    chla_conc = (chla_content_µgml * acetone_mL) / sample_area,
    sampling = case_when(
      sampling == "food web_1" ~ "1st",
      sampling == "food web_2" ~ "2nd"
    )
  ) |>
  rename(condition = con) |>
  select(year, site, chla_conc, sampling, func_group, rep, condition)

benthic_chla <- bind_rows(benthic_chla_2023, benthic_chla_2024)
rm(benthic_chla_2023, benthic_chla_2024)

benthic_bod <- df_bod |>
  filter(func_group == "benthic") |>
  left_join(phyto_do_change, join_by(year, sampling, sampling_date, site, condition)) |>
  left_join(benthic_chla, join_by(year, sampling, func_group, site, rep, condition)) |>
  mutate(
    incubation_area_cm2 = 38,
    volume_L = (benthic_depth_cm * incubation_area_cm2) * cm3_to_m3 * m3_to_L,
    do_change_mgL = end_do_mgL - start_do_mgL,
    incubation_hr = as.numeric(end_time - start_time) * 24,
    all_do_change_mg = do_change_mgL * volume_L,
    phyto_do_change_mg = phyto_do_rate_mgL_hr * incubation_hr * volume_L,
    benthic_do_change_mg = all_do_change_mg - phyto_do_change_mg,
    do_change_rate = benthic_do_change_mg / ((incubation_area_cm2 * cm2_to_m2) * incubation_hr),
    rate_unit = "mg O2 m-2 hr-1"
  )

benthic_bod_light <- split_light_dark(benthic_bod, "do_change_rate",
                                      c("year", "sampling", "sampling_date", "site", "func_group"))

benthic_bod_do <- benthic_bod_light |>
  calc_daily_gpp_nep() |>
  mutate(
    nep = pmax(nep, 0),
    season_days = culture_conv$day_culture[match(substr(site, 1, 2), culture_conv$water)],
    production_c = nep * mg_to_g * (1/o2_molar_mass) * (1/pq) * co2_mass * co2e_to_c * season_days,
    biomass_c = (chla_conc * benthic_chla_to_c * µg_to_g) / cm2_to_m2,
    biomass = biomass_c * c_to_dw_producer * dw_to_ww_producer,
    production_unit = "g C m-2 culture-period-1",
    biomass_unit = "g WW m-2"
  )

ba_summary <- summarise_producer(benthic_bod_do)
tab_benthic <- calc_weighted_summary(ba_summary)

## Periphyton -----

peri_chla <- process_chla |>
  filter(str_detect(sampling, "food web"), func_group == "periphyton") |>
  unite(col = "site", water, type, species, group, sep = "_", remove = FALSE) |>
  mutate(
    year = as.character(year),
    sample_area = (5 * 11) + (6 * 11),
    chla_content_µgml = calc_chla_trichromatic(wl_630, wl_647, wl_664, wl_750),
    chla_conc = (chla_content_µgml * acetone_mL) / sample_area * (periphyton_total_volume / periphyton_filter_volume),
    sampling = case_when(
      sampling == "food web_2" ~ "2nd",
      sampling == "food web_3" ~ "3rd"
    )
  ) |>
  rename(condition = con) |>
  summarise(
    chla_conc = mean(chla_conc),
    .by = c(site, sampling, func_group)
  )

peri <- df_bod |>
  filter(str_detect(site, "fw"), func_group == "periphyton", sampling_date == as.Date('2024-11-07')) |>
  left_join(phyto_do_change, join_by(year, sampling, sampling_date, site, condition)) |>
  left_join(peri_chla, join_by(site, sampling, func_group)) |>
  mutate(
    volume_L = case_when(year == "2023" ~ 0.3, year == "2024" ~ 1.3),
    incubation_area_cm2 = case_when(year == "2023" ~ 16, year == "2024" ~ (5 * 11) + (6 * 11)),
    do_change_mgL = end_do_mgL - start_do_mgL,
    incubation_hr = as.numeric(end_time - start_time) * 24,
    total_do_change_mg = do_change_mgL * volume_L,
    phyto_do_change_mg = phyto_do_rate_mgL_hr * incubation_hr * volume_L,
    peri_do_change_mg = total_do_change_mg - phyto_do_change_mg,
    do_change_rate = peri_do_change_mg / ((incubation_area_cm2 * cm2_to_m2) * incubation_hr),
    rate_unit = "mg O2 m-2 hr-1"
  )

peri_light <- split_light_dark(peri, "do_change_rate",
                               c("year", "sampling", "sampling_date", "site", "func_group"))

peri_do <- peri_light |>
  calc_daily_gpp_nep() |>
  mutate(
    nep = pmax(nep, 0),
    season_days = culture_conv$day_culture[match(substr(site, 1, 2), culture_conv$water)],
    production_c = nep * mg_to_g * (1/o2_molar_mass) * (1/pq) * co2_mass * co2e_to_c * season_days,
    biomass_c = (chla_conc * peri_chla_to_c * µg_to_g) / cm2_to_m2,
    biomass = biomass_c * c_to_dw_producer * dw_to_ww_producer,
    production_unit = "g C m-2 culture-period-1",
    biomass_unit = "g WW m-2"
  )

df_peri <- summarise_producer(peri_do)
tab_peri <- calc_weighted_summary(df_peri)

## Macroalgae / plants -----

sw_imta_plant_T_2024 <- tibble::tribble(
  ~site,            ~stocking_date, ~harvest_date, ~stocking_kg, ~harvest_kg,
  "sw_imta_plant_T", "2024-03-27",   "2024-04-24",  50,           64.580,
  "sw_imta_plant_T", "2024-04-24",   "2024-05-22",  50,           43.525,
  "sw_imta_plant_T", "2024-05-22",   "2024-06-19",  50,           62.505,
  "sw_imta_plant_T", "2024-06-19",   "2024-07-16",  50,           61.355,
  "sw_imta_plant_T", "2024-07-16",   "2024-08-14",  50,           66.011,
  "sw_imta_plant_T", "2024-08-14",   "2024-09-12",  50,           61.711,
  "sw_imta_plant_T", "2024-09-12",   "2024-10-15",  50,           56.656
) |>
  dplyr::mutate(
    stocking_date      = as.Date(stocking_date),
    harvest_date       = as.Date(harvest_date),
    culture_days       = as.integer(harvest_date - stocking_date),
    weight_gain_per_day = (harvest_kg - stocking_kg) / culture_days) |>
  dplyr::summarise(
    weight_gain_per_day = mean(weight_gain_per_day),
    mean_weight_kg      = mean(c(rbind(stocking_kg, harvest_kg))),
    mean_harvest_kg     = sum(harvest_kg) - sum(stocking_kg),
    .by = "site"
  ) |>
  left_join(df_info |>
              filter(year == "2024", site == "sw_imta_plant_T") |>
              select(year, site, day_culture, area_m2), join_by(site)) |>
  mutate(total_biomass_g = mean_weight_kg * kg_to_g,
         biomass = total_biomass_g / area_m2,
         total_production_g = weight_gain_per_day * day_culture * kg_to_g,
         production = total_production_g / area_m2,
         pb_ratio = production / biomass
  ) |>
  select(year, site, biomass, pb_ratio, mean_harvest_kg)

df_bod |>
  filter(func_group == "macroalgae") |>
  left_join(phyto_do_change, join_by(year, sampling, sampling_date, site, condition)) |>
  mutate(
    volume_L = 1.3,
    biomass_wet_g = case_when(
      year == "2023" ~ biomass_dry_g * dw_to_ww_organicTissue,
      year == "2024" ~ biomass_wet_g
    )) |>
  select(sampling_date, site, biomass_wet_g)

macro <- df_bod |>
  filter(func_group == "macroalgae") |>
  left_join(phyto_do_change, join_by(year, sampling, sampling_date, site, condition)) |>
  mutate(
    volume_L = 1.3,
    biomass_wet_g = case_when(
      year == "2023" ~ biomass_dry_g * dw_to_ww_organicTissue,
      year == "2024" ~ biomass_wet_g
    ),
    do_change_mgL = end_do_mgL - start_do_mgL,
    incubation_hr = as.numeric(end_time - start_time) * 24,
    all_do_change_mg = do_change_mgL * volume_L,
    phyto_do_change_mg = phyto_do_rate_mgL_hr * incubation_hr * volume_L,
    macro_do_change_mg = all_do_change_mg - phyto_do_change_mg,
    do_change_rate = macro_do_change_mg / (biomass_wet_g * incubation_hr),
    rate_unit = "mg O2 g-1 hr-1"
  )

macro_light <- split_light_dark(macro, "do_change_rate",
                                c("year", "sampling", "sampling_date", "site", "func_group"))

info_plant <- info |>
  filter(str_detect(site, "imta") & str_detect(site, "plant") & str_detect(site, "T")) |>
  select(year, site, stock_plant_total_kg, harvest_plant_total_kg) |>
  mutate(
    year = as.character(year),
    harvest_plant_total_kg = ifelse(
      harvest_plant_total_kg == "ND", stock_plant_total_kg, as.numeric(harvest_plant_total_kg)
    )
  )

macro_do <- macro_light |>
  left_join(info_plant, join_by(year, site)) |>
  calc_daily_gpp_nep() |>
  mutate(
    nep = pmax(nep, 0),
    mean_biomass = (stock_plant_total_kg + harvest_plant_total_kg) / 2,
    biomass = (mean_biomass * kg_to_g) / area_m2,

    biomass = ifelse(
      year == "2024" & site == "sw_imta_plant_T", sw_imta_plant_T_2024$biomass, biomass
    ),
    days_2024 = df_info$day_culture[match(paste("2024", site), paste(df_info$year, df_info$site))],
    production_c = nep * mg_to_g * (1/o2_molar_mass) * (1/pq) * co2_mass * co2e_to_c * biomass * days_2024,
    biomass_c = biomass / c_to_ww_plant,
    biomass_unit = "g WW m-2",
    production_unit = "g m-2 culture period-1"
  )

df_macro <- summarise_producer(macro_do) |>
  left_join(
    info_plant |>
      mutate(
        group = str_extract(site, "T\\d*$"),
        harvest_plant_total_kg = as.numeric(harvest_plant_total_kg)
      ) |>
      select(year, group, harvest_plant_total_kg),
    join_by(year, group)
  ) |>
  mutate(catch = ifelse(
    year == "2024" & water == "sw",
    sw_imta_plant_T_2024$mean_harvest_kg * kg_to_g / area_m2,
    as.numeric(harvest_plant_total_kg) * kg_to_g / area_m2
  )) |>
  select(-harvest_plant_total_kg)

zero_rows_macro_plant <- data.frame(
  year = as.character(c(2023, 2023, 2023, 2024, 2024, 2024)),
  water = c("sw", "fw", "fw", "sw", "fw", "fw"),
  species = NA,
  group = c("T", "T1", "T2", "T", "T1", "T2"),
  biomass = 0,
  pb_ratio = 0,
  area_m2 = c(1200, 135, 146, 1000, 135, 146),
  catch = 0
)
df_macro <- rbind(df_macro, zero_rows_macro_plant)

tab_plant_pb <- calc_weighted_summary(filter(df_macro, !is.na(species))) |> select(-avg_weighted_biomass)
tab_plant_b <- calc_weighted_summary(df_macro) |> select(-avg_weighted_pb)

tab_plant_catch <- df_macro |>
  mutate(new_site = calc_new_site(water, group)) |>
  mutate(
    total_area_m2  = sum(area_m2),
    weighted_catch = (area_m2 / total_area_m2) * catch,
    .by = c(year, water, new_site)
  ) |>
  summarise(avg_weighted_catch = sum(weighted_catch, na.rm = TRUE),
            .by = c(year, water, new_site))

# Detritus -----

## Water column -----

df_doc <- doc |>
  filter(type == "imta", str_detect(group, "T"), source == "water") |>
  mutate(
    pond_volume_L = depth_m * area_m2 * m3_to_L,
    total_doc_mg = doc * pond_volume_L,
    doc_mg_m2 = total_doc_mg / area_m2,
    .by = c(year, site)
  ) |>
  summarise(doc_mg_m2 = mean(doc_mg_m2, na.rm = TRUE), .by = c(year, site))

df_poc <- poc |>
  filter(type == "imta", str_detect(group, "T"), source == "water") |>
  mutate(
    pond_volume_L = depth_m * area_m2 * m3_to_L,
    total_poc_mg = poc * pond_volume_L,
    poc_mg_m2 = total_poc_mg / area_m2,
    .by = c(year, site)
  ) |>
  summarise(poc_mg_m2 = mean(poc_mg_m2, na.rm = TRUE), .by = c(year, site))

df_water <- df_doc |>
  left_join(df_poc, join_by(year, site)) |>
  mutate(
    detritus_doc = doc_mg_m2 * mg_to_g * oc_to_ww_organicTissue,
    detritus_poc = poc_mg_m2 * mg_to_g * oc_to_ww_organicTissue,
    detritus_unit = "g WW m-2"
  )

water_det <- df_water |>
  separate_site_columns() |>
  left_join(info_area, join_by(year, site))

tab_water <- water_det |>
  mutate(new_site = calc_new_site(water, group)) |>
  mutate(
    total_area_m2 = sum(area_m2),
    weighted_detritus_doc = (area_m2 / total_area_m2) * detritus_doc,
    weighted_detritus_poc = (area_m2 / total_area_m2) * detritus_poc,
    .by = c(year, water, new_site)
  ) |>
  summarise(
    avg_weighted_detritus_doc = sum(weighted_detritus_doc),
    avg_weighted_detritus_poc = sum(weighted_detritus_poc),
    .by = c(year, water, new_site)
  ) |>
  as.data.frame()

## Sediment -----

crab_dbd <- soc |>
  filter(year == "2024", species == "crab", group == "T") |>
  summarise(envr_dbd = mean(envr_dbd, na.rm = TRUE)) |>
  pull(envr_dbd)

df_soc <- soc |>
  filter(type == "imta", str_detect(group, "T")) |>
  mutate(
    envr_dbd = if_else(year == "2023" & str_detect(sampling, "foodweb"), crab_dbd, envr_dbd),
    soc_g_cm2 = 5 * envr_dbd * (C_percent / 100),
    detritus_soil = soc_g_cm2 * (1 / cm2_to_m2) * oc_to_ww_organicTissue,
    detritus_unit = "g WW m-2",
    new_site = calc_new_site(water, group)
  ) |>
  summarise(
    area_m2 = mean(area_m2),
    detritus_soil = mean(detritus_soil),
    .by = c(year, water, site, new_site)
  )

zero_rows_soil <- data.frame(
  year = as.character(c(2024, 2024)),
  site = c("fw_imta_plant_T1", "fw_imta_plant_T2"),
  water = c("fw", "fw"),
  new_site = c("fw_T1", "fw_T2"),
  area_m2 = c(51, 51),
  detritus_soil = c(0, 0)
)

df_soc <- rbind(df_soc, zero_rows_soil)

tab_soil <- df_soc |>
  mutate(
    total_area_m2 = sum(area_m2),
    weighted_detritus_soil = (area_m2 / total_area_m2) * detritus_soil,
    .by = c(year, water, new_site)
  ) |>
  summarise(avg_weighted_detritus_soil = sum(weighted_detritus_soil),
            .by = c(year, water, new_site)) |>
  as.data.frame()

## Commercial feed -----

tab_feed <- info |>
  filter(type == "imta", str_detect(group, "T")) |>
  mutate(
    detritus_feed = (as.numeric(feed_input_rate) * kg_to_g) / area_m2,
    detritus_unit = "g WW m-2",
    detritus_feed_import = (as.numeric(feed_input_kg) * kg_to_g) / area_m2,
    new_site = calc_new_site(water, group)
  ) |>
  mutate(
    total_area_m2 = sum(area_m2),
    weighted_detritus_feed = (area_m2 / total_area_m2) * detritus_feed,
    weighted_detritus_feed_import = (area_m2 / total_area_m2) * detritus_feed_import,
    .by = c(year, water, new_site)
  ) |>
  summarise(avg_weighted_detritus_feed = sum(weighted_detritus_feed),
            avg_weighted_detritus_feed_import = sum(weighted_detritus_feed_import),
            .by = c(year, water, new_site)) |>
  as.data.frame()

fw_fish_feed_detinput <- info |>
  filter(type == "imta", water == "fw", str_detect(group, "T"), year == 2024) |>
  mutate(
    detritus_feed_import = (as.numeric(feed_input_kg) * kg_to_g) / area_m2,
    new_site = calc_new_site(water, group),
    detritus_feed_import = if_else(species == "fish", detritus_feed_import, 0)
  ) |>
  mutate(
    total_area_m2 = sum(area_m2),
    weighted_feed_import = (area_m2 / total_area_m2) * detritus_feed_import,
    .by = new_site
  ) |>
  summarise(avg_weighted_feed_import = sum(weighted_feed_import), .by = new_site) |>
  pull(avg_weighted_feed_import) |>
  mean()

# Ecopath assembly -----

ecopath_animals <- tab_animals |>
  filter(year == "2024") |>
  mutate(func_group = case_when(
    species == "crab" ~ "Mud_crab",
    species == "fish" ~ "Sea_bass",
    species == "plant" ~ "Tilapia"
  )) |>
  select(water, new_site, weighted_biomass, pb_ratio, weighted_catch, func_group) |>
  rename(B = weighted_biomass,
         PB = pb_ratio,
         Catch = weighted_catch) |>
  summarise(PB = sum(B * PB) / sum(B),
            B = mean(B),
            Catch = mean(Catch),
            .by = c(water, func_group)) |>
  left_join(animal_conv, join_by(water, func_group)) |>
  mutate(QB = case_when(
    func_group == "Mud_crab" ~ QB_MUD_CRAB,
    func_group == "Sea_bass" ~ QB_SEA_BASS,
    func_group == "Tilapia" ~ QB_TILAPIA
  ) * group_yr,
  type = 0) |>
  select(-group_yr)

ecopath_inverts <- tab_invert |>
  filter(year == "2024") |>
  rename(B = biomass,
         PB = pb_ratio,
         QB = qb_ratio,
         func_group = func_groups) |>
  summarise(B = mean(B),
            PB = mean(PB),
            QB = mean(QB),
            .by = c(water, func_group)) |>
  mutate(type = 0,
         func_group = case_when(
           func_group == "Omnivores" ~ "Omnivores_inverts",
           func_group == "Carnivores" ~ "Carnivores_inverts",
           func_group == "Detritivores" ~ "Detritivores_inverts",
           func_group == "Filter feeders" ~ "Filter_feeders",
         ))

ecopath_zoo <- tab_zooplankton |>
  filter(year == "2024") |>
  mutate(func_group = 'Zooplankton') |>
  rename(B = avg_weighted_biomass) |>
  summarise(B = mean(B),
            .by = c(water, func_group)) |>
  left_join(culture_conv |> select(water, culture_yr), join_by(water)) |>
  mutate(PB = ZOO_PB * culture_yr,
         QB = ZOO_QB * culture_yr) |>
  mutate(type = 0) |>
  select(-culture_yr)

ecopath_phyto <- tab_phyto |>
  mutate(func_group = 'Phytoplankton') |>
  rename(B = avg_weighted_biomass,
         PB = avg_weighted_pb) |>
  summarise(B = mean(B[year == "2024"]),
            PB = mean(PB),
            .by = c(water, func_group)) |>
  mutate(QB = NA,
         type = 1)

ecopath_benthic <- tab_benthic |>
  mutate(func_group = 'Benthic_algae') |>
  rename(B = avg_weighted_biomass,
         PB = avg_weighted_pb) |>
  summarise(B = mean(B[year == "2024"]),
            PB = mean(PB),
            .by = c(water, func_group)) |>
  mutate(QB = NA,
         type = 1)

ecopath_peri <- tab_peri |>
  filter(year == "2024") |>
  mutate(func_group = 'Periphyton') |>
  rename(B = avg_weighted_biomass,
         PB = avg_weighted_pb) |>
  summarise(B = mean(B),
            PB = mean(PB),
            .by = c(water, func_group)) |>
  mutate(QB = NA,
         type = 1)

ecopath_plant <- tab_plant_b |>
  filter(year == "2024") |>
  mutate(func_group = ifelse(water == "sw", "Red_algae", "Duckweed")) |>
  rename(B = avg_weighted_biomass) |>
  left_join(tab_plant_catch |> filter(year == "2024") |> select(-year),
            join_by(water, new_site)) |>
  left_join(
    tab_plant_pb |>
      mutate(func_group = ifelse(water == "sw", "Red_algae", "Duckweed")) |>
      rename(PB = avg_weighted_pb) |>
      summarise(PB = mean(PB), .by = c(water, func_group)),
    join_by(water, func_group)
  ) |>
  summarise(B     = mean(B),
            PB    = mean(PB),
            Catch = mean(avg_weighted_catch),
            .by = c(water, func_group)) |>
  mutate(QB = NA,
         type = 1)

fw_duck_bc <- info |>
  filter(year == 2024, type == "imta", water == "fw", species == "plant", group %in% c("T1", "T2")) |>
  transmute(water, group, area_m2,
            biomass = (as.numeric(stock_plant_total_kg) + as.numeric(harvest_plant_total_kg)) / 2 * kg_to_g / area_m2,
            catch = (as.numeric(harvest_plant_total_kg) - as.numeric(stock_plant_total_kg)) * kg_to_g / area_m2) |>
  bind_rows(filter(zero_rows_macro_plant, year == "2024", water == "fw") |> select(water, group, area_m2, biomass, catch)) |>
  mutate(new_site = calc_new_site(water, group)) |>
  summarise(B = sum(area_m2 / sum(area_m2) * biomass), Catch = sum(area_m2 / sum(area_m2) * catch), .by = new_site) |>
  summarise(B = mean(B), Catch = mean(Catch))
ecopath_plant <- ecopath_plant |>
  mutate(B = if_else(water == "fw", fw_duck_bc$B, B),
         Catch = if_else(water == "fw", fw_duck_bc$Catch, Catch))

ecopath_water_det <- tab_water |>
  filter(year == "2024") |>
  mutate(func_group = 'Water_detritus',
         B = avg_weighted_detritus_doc + avg_weighted_detritus_poc) |>
  summarise(B = mean(B),
            .by = c(water, func_group)) |>
  mutate(PB = NA,
         QB = NA,
         type = 2)

ecopath_soil_det <- tab_soil |>
  filter(year == "2024") |>
  mutate(func_group = 'Sediment_detritus',
         B = avg_weighted_detritus_soil) |>
  summarise(B = mean(B),
            .by = c(water, func_group)) |>
  mutate(PB = NA,
         QB = NA,
         type = 2)

ecopath_feed <- tab_feed |>
  filter(year == "2024") |>
  mutate(func_group = 'Commercial_feed',
         B = avg_weighted_detritus_feed,
         DetInput = avg_weighted_detritus_feed_import) |>
  summarise(B = mean(B),
            DetInput = mean(DetInput),
            .by = c(water, func_group)) |>
  mutate(PB = NA,
         QB = NA,
         type = 2)

ecopath_df <- bind_rows(
  ecopath_animals,
  ecopath_inverts,
  ecopath_zoo,
  ecopath_phyto,
  ecopath_benthic,
  ecopath_peri,
  ecopath_plant,
  ecopath_water_det,
  ecopath_soil_det,
  ecopath_feed
) |>
  mutate(Catch = replace_na(Catch, 0))

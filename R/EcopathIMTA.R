# Setup -----

library(tidyverse)
library(Rpath)
library(patchwork)
library(scales)
library(enaR)
source(here::here("R", "00_load_inputs.R"))

# Culture period lengths -----
SW_CULTURE_DAYS <- const[["sw_culture_days"]]
FW_CULTURE_DAYS <- const[["fw_culture_days"]]
SW_CULTURE_YR <- SW_CULTURE_DAYS / 365
FW_CULTURE_YR <- FW_CULTURE_DAYS / 365

macro_days <- function(water) {
  d <- unique(macro_do$day_culture[stringr::str_detect(macro_do$site, paste0("^", water)) &
    macro_do$year == "2024"])
  stopifnot(length(d) == 1, !is.na(d))
  d
}
SW_MACRO_DAYS <- macro_days("sw")
FW_MACRO_DAYS <- macro_days("fw")

# Feed emission factors (kg CO2e per kg feed) -----

EF_SW_LOW <- 0.20
EF_SW_HIGH <- 1.90
EF_FW_LOW <- 0.82
EF_FW_HIGH <- 2.40
EF_FW_CANDIDATES <- c(0.82, 1.02, 1.07, 1.33, 1.70, 2.40)
EF_SW_MID <- (EF_SW_LOW + EF_SW_HIGH) / 2
EF_FW_MID <- median(EF_FW_CANDIDATES)

sys_area_2024 <- info |>
  filter(year == 2024, type == "imta", str_detect(group, "^T")) |>
  mutate(area_m2 = as.numeric(area_m2)) |>
  summarise(pair_area_m2 = sum(area_m2, na.rm = TRUE), .by = c(water, group)) |>
  summarise(system_area_m2 = mean(pair_area_m2), .by = water)
SW_SYSTEM_AREA_M2 <- sys_area_2024$system_area_m2[sys_area_2024$water == "sw"]
FW_SYSTEM_AREA_M2 <- sys_area_2024$system_area_m2[sys_area_2024$water == "fw"]

fleet <- data.frame(
  func_group = "Fleets",
  B = NA,
  PB = NA,
  Catch = NA,
  QB = NA,
  type = 3,
  DetInput = NA
)

# Field diet fractions -----

DC_MUD_FEED <- 0.387
DC_MUD_SED <- 0.4022
DC_TIL_FEED <- 0.15
DC_TIL_ZOO <- 0.01
DC_TIL_PHYTO <- 0.35
DC_TIL_PERI <- 0.26
DC_TIL_SED_FIELD <- round(1 - DC_TIL_FEED - DC_TIL_ZOO - DC_TIL_PHYTO - DC_TIL_PERI, 10)

SW_FATE_WATER <- c(
  Mud_crab = 0, Detritivores_inverts = 0, Carnivores_inverts = 0, Omnivores_inverts = 0,
  Filter_feeders = 0, Zooplankton = 1, Phytoplankton = 0.6, Benthic_algae = 0,
  Red_algae = 0.5, Water_detritus = 1, Sediment_detritus = 0, Commercial_feed = 0, Fleets = 0
)
SW_FATE_SED <- c(
  Mud_crab = 1, Detritivores_inverts = 1, Carnivores_inverts = 1, Omnivores_inverts = 1,
  Filter_feeders = 1, Zooplankton = 0, Phytoplankton = 0.4, Benthic_algae = 1,
  Red_algae = 0.5, Water_detritus = 0, Sediment_detritus = 1, Commercial_feed = 0, Fleets = 0
)
FW_FATE_WATER <- c(
  Sea_bass = 0.5, Tilapia = 0.5, Zooplankton = 1, Benthic_invert = 0, Phytoplankton = 1,
  Periphyton = 0, Duckweed = 0.5, Water_detritus = 1, Sediment_detritus = 0,
  Commercial_feed = 0, Fleets = 0
)
FW_FATE_SED <- c(
  Sea_bass = 0.5, Tilapia = 0.5, Zooplankton = 0, Benthic_invert = 1, Phytoplankton = 0,
  Periphyton = 1, Duckweed = 0.5, Water_detritus = 0, Sediment_detritus = 1,
  Commercial_feed = 0, Fleets = 0
)

# Group settings -----
UNASSIM_LOW <- c("Sea_bass", "Mud_crab")
apply_static_cols <- function(p, fate_water, fate_sed, detinput) {
  p$model[, BioAcc := ifelse(Type == 3, NA_real_, 0)]
  p$model[, Unassim := ifelse(Type == 0, ifelse(Group %in% UNASSIM_LOW, 0.2, 0.4), ifelse(Type == 3, NA_real_, 0))]
  p$model[, Water_detritus := unname(fate_water[Group])]
  p$model[, Sediment_detritus := unname(fate_sed[Group])]
  p$model[, Commercial_feed := as.numeric(Group == "Commercial_feed")]
  p$model[Group == "Commercial_feed", DetInput := detinput]
  invisible(p)
}

# Detritus EE -----
detritus_ee <- function(model, params) {
  liv <- which(model$type %in% c(0, 1))
  det <- which(model$type == 2)
  fate <- model$DetFate
  fate[is.na(fate)] <- 0
  off_diag <- fate[det, , drop = FALSE]
  diag(off_diag) <- 0
  stopifnot(all(off_diag == 0))
  zero_na <- function(x) replace(x, is.na(x), 0)
  B <- model$Biomass
  loss <- (1 - model$EE[liv]) * B[liv] * model$PB[liv] + B[liv] * zero_na(model$QB[liv]) * model$Unassim[liv]
  det_input <- zero_na(params$model$DetInput[match(model$Group[det], params$model$Group)])
  inflow <- colSums(loss * fate[liv, , drop = FALSE]) + det_input
  cons <- which(model$type == 0)
  eaten <- as.vector(zero_na(model$DC[det, match(cons, liv), drop = FALSE]) %*% (B[cons] * model$QB[cons]))
  eaten / inflow
}

# Mixed trophic impact -----
mti_src <- deparse(body(Rpath::MTI))
stopifnot(sum(grepl("bio <- y$BB", mti_src, fixed = TRUE)) == 1)
mti_rpath_patched <- Rpath::MTI
body(mti_rpath_patched) <- str2lang(paste(sub("bio <- y$BB", "bio <- y$Biomass", mti_src, fixed = TRUE), collapse = "\n"))
mti_fixed <- function(Rpath, Rpath.params) {
  mti <- mti_rpath_patched(Rpath, Rpath.params)
  det <- which(Rpath.params$model$Type == 2)
  mti[cbind(det, det)] <- 0
  mti
}

# Run model -----

run_rpath_model <- function(params, name) {
  warnings_caught <- character(0)
  withCallingHandlers(
    check.rpath.params(params),
    warning = function(w) {
      warnings_caught <<- c(warnings_caught, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  if (length(warnings_caught) > 0) {
    message(name, " check failed:\n", paste("-", warnings_caught, collapse = "\n"))
    return(NULL)
  }
  model <- rpath(params, eco.name = name)
  model$EE[model$type == 2] <- detritus_ee(model, params)
  catch_grp <- params$model$Group[!is.na(params$model$Fleets) & params$model$Fleets > 0]
  unbalanced <- names(which(model$EE > 1 + ifelse(names(model$EE) %in% catch_grp, 1e-9, 0)))
  if (length(unbalanced) > 0) {
    message(name, " UNBALANCED: ", paste(unbalanced, collapse = ", "))
  } else {
    message(name, ": balanced")
  }
  attr(model, "unbalanced") <- unbalanced
  model
}

is_balanced <- function(m) !is.null(m) && length(attr(m, "unbalanced")) == 0

catch_vec <- function(g) {
  v <- replace(g$Catch, is.na(g$Catch), 0)
  v[length(v)] <- NA
  v
}

# Diet blocks -----

apply_mudcrab_base_diet <- function(diet_dt) {
  diet_dt[Group == "Benthic_algae", Mud_crab := 0.100]
  diet_dt[Group == "Zooplankton", Mud_crab := 0.100]
  diet_dt[Group == "Omnivores_inverts", Mud_crab := 0.00534]
  diet_dt[Group == "Carnivores_inverts", Mud_crab := 0.00027]
  diet_dt[Group == "Filter_feeders", Mud_crab := 0.00007]
  diet_dt[Group == "Detritivores_inverts", Mud_crab := 0.00512]
  invisible(diet_dt)
}

apply_sw_invert_diet <- function(diet_dt) {
  diet_dt[Group == "Phytoplankton", Zooplankton := 0.50]
  diet_dt[Group == "Water_detritus", Zooplankton := 0.50]
  diet_dt[Group == "Sediment_detritus", Detritivores_inverts := 1.00]
  diet_dt[Group == "Phytoplankton", Omnivores_inverts := 0.05]
  diet_dt[Group == "Benthic_algae", Omnivores_inverts := 0.45]
  diet_dt[Group == "Sediment_detritus", Omnivores_inverts := 0.50]
  diet_dt[Group == "Omnivores_inverts", Carnivores_inverts := 0.250]
  diet_dt[Group == "Carnivores_inverts", Carnivores_inverts := 0.020]
  diet_dt[Group == "Zooplankton", Carnivores_inverts := 0.250]
  diet_dt[Group == "Filter_feeders", Carnivores_inverts := 0.040]
  diet_dt[Group == "Detritivores_inverts", Carnivores_inverts := 0.010]
  diet_dt[Group == "Sediment_detritus", Carnivores_inverts := 0.430]
  diet_dt[Group == "Phytoplankton", Filter_feeders := 0.02]
  diet_dt[Group == "Water_detritus", Filter_feeders := 0.98]
  invisible(diet_dt)
}

apply_fw_diet <- function(diet_dt, dc_feed_til, dc_duck, dc_sed_til) {

  diet_dt[Group == "Commercial_feed", Sea_bass := 0.686]
  diet_dt[Group == "Zooplankton", Sea_bass := 0.020]
  diet_dt[Group == "Periphyton", Sea_bass := 0.100]
  diet_dt[Group == "Benthic_invert", Sea_bass := 0.194]
  diet_dt[Group == "Import", Sea_bass := 0]

  diet_dt[Group == "Sediment_detritus", Benthic_invert := 0.70]
  diet_dt[Group == "Water_detritus", Benthic_invert := 0.20]
  diet_dt[Group == "Phytoplankton", Benthic_invert := 0.10]

  diet_dt[Group == "Commercial_feed", Tilapia := dc_feed_til]
  diet_dt[Group == "Duckweed", Tilapia := dc_duck]
  diet_dt[Group == "Zooplankton", Tilapia := DC_TIL_ZOO]
  diet_dt[Group == "Phytoplankton", Tilapia := DC_TIL_PHYTO]
  diet_dt[Group == "Periphyton", Tilapia := DC_TIL_PERI]
  diet_dt[Group == "Sediment_detritus", Tilapia := dc_sed_til]
  diet_dt[Group == "Phytoplankton", Zooplankton := 0.40]
  diet_dt[Group == "Water_detritus", Zooplankton := 0.60]
  invisible(diet_dt)
}

# Seawater field model -----

sw_field_groups <- ecopath_df |>
  filter(water == "sw") |>
  select(-water)
sw_field_groups <- rbind(sw_field_groups, fleet)
sw_feed_detinput <- ecopath_df |>
  filter(water == "sw", func_group == "Commercial_feed") |>
  pull(DetInput)

print(sw_field_groups)

sw_field_rpath <- create.rpath.params(
  group = sw_field_groups$func_group,
  type  = sw_field_groups$type
)

n_sw <- nrow(sw_field_rpath$model)
sw_field_rpath$model[, Biomass := sw_field_groups$B]
sw_field_rpath$model[, PB := sw_field_groups$PB]
sw_field_rpath$model[, QB := sw_field_groups$QB]

sw_catch <- replace(sw_field_groups$Catch, is.na(sw_field_groups$Catch), 0)
sw_catch[n_sw] <- NA
sw_field_rpath$model[, Fleets := sw_catch]

apply_static_cols(sw_field_rpath, SW_FATE_WATER, SW_FATE_SED, sw_feed_detinput)

sw_field_rpath$diet[Group == "Commercial_feed", Mud_crab := DC_MUD_FEED]
apply_mudcrab_base_diet(sw_field_rpath$diet)
sw_field_rpath$diet[Group == "Sediment_detritus", Mud_crab := DC_MUD_SED]

apply_sw_invert_diet(sw_field_rpath$diet)

sw_field_model <- run_rpath_model(sw_field_rpath, "Seawater_Field")
sw_field_model

# Freshwater field model -----

invert_dw <- c(Oligochaeta = 48.51, Chironomidae = 3.86, Mollusca = 0.87)
invert_ww_ratio_lo <- c(5.00, 7.04, 7.04)
invert_ww_ratio_hi <- c(6.25, 7.04, 7.04)
invert_B_pond <- mean(c(sum(invert_dw * invert_ww_ratio_lo), sum(invert_dw * invert_ww_ratio_hi)))
invert_seabass_area <- info |>
  filter(year == 2024, water == "fw", type == "imta", species == "fish", str_detect(group, "^T")) |>
  summarise(a = mean(as.numeric(area_m2))) |>
  pull(a)
invert_B <- invert_B_pond * invert_seabass_area / FW_SYSTEM_AREA_M2
invert_PB <- 3.6 * FW_CULTURE_DAYS / 365
invert_QB <- invert_PB / 0.3
benthic_invert_row <- tibble(
  func_group = "Benthic_invert",
  B = invert_B,
  PB = invert_PB,
  Catch = 0,
  QB = invert_QB,
  type = 0,
  DetInput = NA_real_
)

fw_field_groups <- ecopath_df |>
  filter(water == "fw") |>
  select(-water)
fw_field_groups <- add_row(fw_field_groups, benthic_invert_row, .after = 3)
fw_field_groups <- rbind(fw_field_groups, fleet)
fw_feed_detinput <- ecopath_df |>
  filter(water == "fw", func_group == "Commercial_feed") |>
  pull(DetInput)

print(fw_field_groups)

fw_field_rpath <- create.rpath.params(
  group = fw_field_groups$func_group,
  type  = fw_field_groups$type
)

n_fw <- nrow(fw_field_rpath$model)
fw_field_rpath$model[, Biomass := fw_field_groups$B]
fw_field_rpath$model[, PB := fw_field_groups$PB]
fw_field_rpath$model[, QB := fw_field_groups$QB]

fw_catch <- replace(fw_field_groups$Catch, is.na(fw_field_groups$Catch), 0)
fw_catch[n_fw] <- NA
fw_field_rpath$model[, Fleets := fw_catch]

apply_static_cols(fw_field_rpath, FW_FATE_WATER, FW_FATE_SED, fw_feed_detinput)

apply_fw_diet(fw_field_rpath$diet,
  dc_feed_til = DC_TIL_FEED, dc_duck = 0, dc_sed_til = DC_TIL_SED_FIELD
)

fw_field_model <- run_rpath_model(fw_field_rpath, "Freshwater_Field")
fw_field_model

# Gradient scenarios -----

# FW gradient builder -----
build_fw_grad <- function(b_duck, dc_duck, dc_feed_til, name,
                          groups = fw_field_groups, dc_pool = DC_TIL_FEED, diet_fun = NULL) {
  detinput <- fw_fish_feed_detinput +
    (fw_feed_detinput - fw_fish_feed_detinput) * (dc_feed_til / dc_pool)
  dc_sed_til <- round(1 - dc_feed_til - dc_duck - DC_TIL_ZOO - DC_TIL_PHYTO - DC_TIL_PERI, 10)

  p <- create.rpath.params(group = groups$func_group, type = groups$type)
  p$model[, Biomass := groups$B]
  p$model[, PB := groups$PB]
  p$model[, QB := groups$QB]
  p$model[, Fleets := catch_vec(groups)]
  p$model[Group == "Duckweed", Fleets := 0]
  p$model[Group == "Duckweed", Biomass := b_duck]
  apply_static_cols(p, FW_FATE_WATER, FW_FATE_SED, detinput)

  apply_fw_diet(p$diet, dc_feed_til, dc_duck, dc_sed_til)
  if (!is.null(diet_fun)) diet_fun(p$diet)

  list(params = p, model = run_rpath_model(p, name))
}

GRAD_MULT <- c(1, 1.5, 2, 2.5)
fw_b_steps <- fw_field_groups$B[fw_field_groups$func_group == "Duckweed"] * GRAD_MULT

# SW gradient builder -----
build_sw_grad <- function(b_algae, dc_algae, dc_feed_mud, name,
                          groups = sw_field_groups, dc_pool = DC_MUD_FEED, diet_fun = NULL) {
  detinput <- sw_feed_detinput * (dc_feed_mud / dc_pool)

  g <- groups
  g[g$func_group == "Red_algae", "B"] <- b_algae
  g[g$func_group == "Red_algae", "Catch"] <- 0

  sw_catch_g <- replace(g$Catch, is.na(g$Catch), 0)
  sw_catch_g[n_sw] <- NA

  p <- create.rpath.params(group = g$func_group, type = g$type)
  p$model[, Biomass := g$B]
  p$model[, PB := g$PB]
  p$model[, QB := g$QB]
  p$model[, Fleets := sw_catch_g]
  apply_static_cols(p, SW_FATE_WATER, SW_FATE_SED, detinput)

  p$diet[Group == "Commercial_feed", Mud_crab := dc_feed_mud]
  p$diet[Group == "Red_algae", Mud_crab := dc_algae]
  apply_mudcrab_base_diet(p$diet)
  p$diet[Group == "Sediment_detritus", Mud_crab := round(DC_MUD_SED + (DC_MUD_FEED - dc_pool), 10)]

  apply_sw_invert_diet(p$diet)
  if (!is.null(diet_fun)) diet_fun(p$diet)

  list(params = p, model = run_rpath_model(p, name))
}

sw_b_steps <- sw_field_groups$B[sw_field_groups$func_group == "Red_algae"] * GRAD_MULT

# Realistic gradient scenarios -----
TARGET_EE_FW <- fw_field_model$EE[fw_field_model$Group == "Duckweed"]
TARGET_EE_SW <- sw_field_model$EE[sw_field_model$Group == "Red_algae"]

RATIO_DUCK_PELLET <- 9.0
RATIO_ALGAE_FISH <- 22.3

Q_tilapia <- fw_field_groups$B[fw_field_groups$func_group == "Tilapia"] *
  fw_field_groups$QB[fw_field_groups$func_group == "Tilapia"]
PB_duckweed <- fw_field_groups$PB[fw_field_groups$func_group == "Duckweed"]
Q_mudcrab <- sw_field_groups$B[sw_field_groups$func_group == "Mud_crab"] *
  sw_field_groups$QB[sw_field_groups$func_group == "Mud_crab"]
PB_RED_ALGAE <- sw_field_groups$PB[sw_field_groups$func_group == "Red_algae"]

realistic_dc <- function(b_actual, target_ee, ratio, pb, q, cap, pool) {

  b_input <- b_actual / ratio
  dc_macro <- min(max(target_ee * pb * b_input / q, 0), cap)
  list(dc_macro = dc_macro, dc_feed = pool - dc_macro, displaced_g = dc_macro * q, b_input = b_input)
}

fw_dc_args <- list(target_ee = TARGET_EE_FW, pb = PB_duckweed, q = Q_tilapia, cap = DC_TIL_FEED, pool = DC_TIL_FEED)
sw_dc_args <- list(target_ee = TARGET_EE_SW, pb = PB_RED_ALGAE, q = Q_mudcrab, cap = DC_MUD_FEED, pool = DC_MUD_FEED)
dc_steps <- function(b_steps, dc_args, ratio) {
  lapply(b_steps, \(b) do.call(realistic_dc, c(list(b_actual = b, ratio = ratio), dc_args)))
}
fw_real <- dc_steps(fw_b_steps, fw_dc_args, RATIO_DUCK_PELLET)
sw_real <- dc_steps(sw_b_steps, sw_dc_args, RATIO_ALGAE_FISH)

stopifnot(
  all(vapply(fw_real, \(r) r$dc_macro * Q_tilapia * RATIO_DUCK_PELLET, numeric(1)) <=
    TARGET_EE_FW * PB_duckweed * fw_b_steps * (1 + 1e-9)),
  all(vapply(sw_real, \(r) r$dc_macro * Q_mudcrab * RATIO_ALGAE_FISH, numeric(1)) <=
    TARGET_EE_SW * PB_RED_ALGAE * sw_b_steps * (1 + 1e-9))
)

fw_grad_1x_realistic <- build_fw_grad(fw_real[[1]]$b_input, fw_real[[1]]$dc_macro, fw_real[[1]]$dc_feed, "FW_grad_1x_realistic")
fw_grad_1.5x_realistic <- build_fw_grad(fw_real[[2]]$b_input, fw_real[[2]]$dc_macro, fw_real[[2]]$dc_feed, "FW_grad_1.5x_realistic")
fw_grad_2x_realistic <- build_fw_grad(fw_real[[3]]$b_input, fw_real[[3]]$dc_macro, fw_real[[3]]$dc_feed, "FW_grad_2x_realistic")
fw_grad_2.5x_realistic <- build_fw_grad(fw_real[[4]]$b_input, fw_real[[4]]$dc_macro, fw_real[[4]]$dc_feed, "FW_grad_2.5x_realistic")

sw_grad_1x_realistic <- build_sw_grad(sw_real[[1]]$b_input, sw_real[[1]]$dc_macro, sw_real[[1]]$dc_feed, "SW_grad_1x_realistic")
sw_grad_1.5x_realistic <- build_sw_grad(sw_real[[2]]$b_input, sw_real[[2]]$dc_macro, sw_real[[2]]$dc_feed, "SW_grad_1.5x_realistic")
sw_grad_2x_realistic <- build_sw_grad(sw_real[[3]]$b_input, sw_real[[3]]$dc_macro, sw_real[[3]]$dc_feed, "SW_grad_2x_realistic")
sw_grad_2.5x_realistic <- build_sw_grad(sw_real[[4]]$b_input, sw_real[[4]]$dc_macro, sw_real[[4]]$dc_feed, "SW_grad_2.5x_realistic")

# Save and summarize -----

dir.create(here::here("output"), showWarnings = FALSE)

models_main <- list(sw_field = sw_field_model, fw_field = fw_field_model)
params_main <- list(sw_field = sw_field_rpath, fw_field = fw_field_rpath)

grads <- mget(paste0(rep(c("sw", "fw"), each = 4), "_grad_", c("1x", "1.5x", "2x", "2.5x"), "_realistic"))
models_list <- c(models_main, lapply(grads, `[[`, "model"))
params_list <- c(params_main, lapply(grads, `[[`, "params"))

saveRDS(models_list, here::here("output/rpath_models_final.rds"))
stopifnot(length(models_list) == 10, all(vapply(models_list, is_balanced, logical(1))))

# Helper functions -----

calc_ecosystem_indices <- function(model, params, name, macro_real = NULL) {
  m <- model
  B <- m$Biomass
  PB <- m$PB
  QB <- m$QB
  EE <- m$EE
  GE <- m$GE
  U <- m$Unassim
  TL <- m$TL
  type <- m$type
  DC <- m$DC
  living <- type == 0
  producer <- type == 1
  detritus <- type == 2
  alive <- living | producer
  consumer_idx <- which(living)
  producer_rows <- which(producer)
  detritus_rows <- which(detritus)
  P <- B * PB

  if (!is.null(macro_real)) {
    idx <- which(m$Group == macro_real$group)
    B[idx] <- macro_real$b_real
    P[idx] <- macro_real$b_real * macro_real$pb_real
  }
  Q <- ifelse(!is.na(QB), B * QB, 0)
  R <- Q * (1 - GE - U)
  FD <- Q * U + PB * (1 - EE) * B
  exports <- rowSums(m$Landings, na.rm = TRUE)

  n_cons_cols <- length(consumer_idx)
  herb_flow <- 0
  detri_flow <- 0
  for (j in seq_len(n_cons_cols)) {
    qi <- Q[consumer_idx[j]]
    dc_col <- DC[, m$Group[consumer_idx[j]]]
    herb_flow <- herb_flow + qi * sum(dc_col[producer_rows], na.rm = TRUE)
    detri_flow <- detri_flow + qi * sum(dc_col[detritus_rows], na.rm = TRUE)
  }

  det_import <- sum(params$model$DetInput[match(m$Group[detritus], params$model$Group)], na.rm = TRUE)
  TST <- sum(Q[alive] + R[alive] + exports[alive] + FD[alive], na.rm = TRUE) + det_import
  TP <- sum(P[alive], na.rm = TRUE)
  TPP <- sum(P[producer], na.rm = TRUE)
  TC <- sum(Q[alive], na.rm = TRUE)
  TR <- sum(R[alive], na.rm = TRUE)
  TFD <- sum(FD[alive], na.rm = TRUE)
  TB <- sum(B[alive], na.rm = TRUE)
  NSP <- TPP - TR
  Total_Catch <- sum(exports, na.rm = TRUE)
  tl_hi <- alive & TL > 2
  mTLco <- weighted.mean(TL[tl_hi], B[tl_hi], na.rm = TRUE)
  n_living <- m$NUM_LIVING
  CI <- sum(DC > 0, na.rm = TRUE) / n_living^2
  soi_per <- vapply(seq_len(n_cons_cols), function(j) {
    frac <- DC[, m$Group[consumer_idx[j]]]
    nonzero <- !is.na(frac) & frac > 0
    if (sum(nonzero) < 2) {
      return(0)
    }
    mu <- sum(frac[nonzero] * TL[nonzero])
    sum(frac[nonzero] * (TL[nonzero] - mu)^2)
  }, numeric(1))
  log_q <- log(Q[consumer_idx] / min(Q[consumer_idx]))
  SOI <- sum(soi_per * log_q) / sum(log_q)
  feed_in <- params$model[Group == "Commercial_feed", DetInput]

  if (length(feed_in) == 0 || is.na(feed_in)) feed_in <- 0
  data.frame(
    Model = name,
    Total_Biomass = round(TB, 2),
    Total_Production = round(TP, 2),
    Total_Consumption = round(TC, 2),
    Total_Respiration = round(TR, 2),
    Total_Flow_Detritus = round(TFD, 2),
    Total_System_Throughput = round(TST, 2),
    Net_System_Production = round(NSP, 2),
    Total_Catch = round(Total_Catch, 2),
    Total_Primary_Production = round(TPP, 2),
    PR_ratio = round(TPP / TR, 3),
    mTLco = round(mTLco, 3),
    Connectance_Index = round(CI, 3),
    System_Omnivory_Index = round(SOI, 4),
    Herb_flow = round(herb_flow, 2),
    Detri_flow = round(detri_flow, 2),
    DH_ratio = round(detri_flow / herb_flow, 3),
    Feed_DetInput = round(feed_in, 1),
    Feed_dependency = round(feed_in / TST, 4),
    Q_TST = round(TC / TST, 4),
    Ex_TST = round(Total_Catch / TST, 4),
    PP_TST = round(TPP / TST, 4),
    stringsAsFactors = FALSE
  )
}

extract_diet_matrix <- function(model, params) {
  m <- model
  DC <- m$DC
  groups <- m$Group
  living_names <- groups[m$type %in% c(0, 1)]
  consumer_names <- groups[m$type == 0]

  dc_df <- as.data.frame(DC)
  colnames(dc_df) <- living_names

  import_row <- params$diet[Group == "Import", .SD, .SDcols = consumer_names]
  import_cons <- as.numeric(import_row[1, ])
  import_full <- setNames(
    c(import_cons, rep(0, length(living_names) - length(consumer_names))),
    living_names
  )

  dc_full <- rbind(dc_df, as.data.frame(t(import_full)))
  rownames(dc_full) <- c(groups, "Import")
  round(dc_full, 4)
}

extract_mortality <- function(model, name) {
  m <- model
  PB <- m$PB
  EE <- m$EE
  type <- m$type
  B <- m$Biomass
  M0 <- PB * (1 - EE)
  F <- rowSums(m$Landings, na.rm = TRUE) / ifelse(B > 0, B, NA)
  M2 <- PB * EE - F
  alive <- type %in% c(0, 1)
  data.frame(
    Model = name, Group = m$Group[alive], TL = round(m$TL[alive], 3),
    M0 = round(M0[alive], 4), M2 = round(pmax(M2[alive], 0), 4),
    F = round(F[alive], 4), PB = round(PB[alive], 4),
    stringsAsFactors = FALSE
  )
}

calc_keystoneness <- function(model, params, name) {
  mti_mat <- mti_fixed(model, params)
  all_groups <- model$Group
  dimnames(mti_mat) <- list(all_groups, all_groups)
  living_grp <- all_groups[model$type %in% c(0, 1)]
  mti_sub <- mti_mat[living_grp, living_grp]
  B_living <- model$Biomass[model$type %in% c(0, 1)]
  B_total <- sum(B_living, na.rm = TRUE)
  diag(mti_sub) <- 0
  eps <- sqrt(rowSums(mti_sub^2, na.rm = TRUE))
  p <- B_living / B_total
  KS <- log10(eps * (1 - p))
  data.frame(
    Model = name,
    Group = living_grp,
    TL = model$TL[model$type %in% c(0, 1)],
    epsilon = round(eps, 4),
    rel_total_impact = round(eps / max(eps), 4),
    p_B = round(p, 4),
    KS = round(KS, 4),
    stringsAsFactors = FALSE
  )
}

calc_network_indices <- function(model, params, name, return_spine = FALSE) {
  incl_idx <- which(model$type %in% c(0, 1, 2))
  n_incl <- length(incl_idx)
  incl_grp <- as.character(model$Group[incl_idx])
  incl_type <- model$type[incl_idx]
  liv <- which(incl_type %in% c(0, 1))
  det <- which(incl_type == 2)

  zero_na <- function(x) replace(x, is.na(x), 0)
  B <- zero_na(model$Biomass[incl_idx])
  PB <- zero_na(model$PB[incl_idx])
  QB <- zero_na(model$QB[incl_idx])
  U <- zero_na(model$Unassim[incl_idx])
  EE <- zero_na(model$EE[incl_idx])
  P <- B * PB
  Q <- B * QB
  catch <- rowSums(model$Landings[incl_idx, , drop = FALSE], na.rm = TRUE)
  fate <- model$DetFate[incl_idx, , drop = FALSE]
  fate[is.na(fate)] <- 0
  det_input <- zero_na(params$model$DetInput[match(incl_grp, params$model$Group)])

  flow_mat <- matrix(0, n_incl, n_incl, dimnames = list(incl_grp, incl_grp))
  cons <- which(incl_type == 0)
  DC <- model$DC
  dc_col <- match(cons, liv)
  flow_mat[, cons] <- sweep(zero_na(DC[seq_len(n_incl), dc_col, drop = FALSE]), 2, Q[cons], "*")
  import_in <- rep(0, n_incl)
  import_in[cons] <- Q[cons] * zero_na(DC[nrow(DC), dc_col])
  loss <- (1 - EE[liv]) * P[liv] + Q[liv] * U[liv]
  flow_mat[liv, det] <- flow_mat[liv, det] + loss * fate[liv, ]
  det_in <- colSums(det_input * fate) + colSums(flow_mat[liv, det, drop = FALSE])
  det_out <- rowSums(flow_mat[det, , drop = FALSE])
  det_unused <- pmax(det_in - det_out, 0)
  pass <- det_unused * fate[det, , drop = FALSE]
  diag(pass) <- 0
  flow_mat[det, det] <- flow_mat[det, det] + pass

  input_vec <- import_in + ifelse(incl_type == 1, P, 0)
  input_vec[det] <- colSums(det_input * fate)
  resp_vec <- ifelse(incl_type == 0, Q * (1 - U) - P, 0)
  export_vec <- catch
  export_vec[det] <- export_vec[det] + input_vec[det] + colSums(flow_mat[, det, drop = FALSE]) - rowSums(flow_mat[det, , drop = FALSE])

  storage_vec <- pmax(B, 1e-6)
  living_vec <- incl_type %in% c(0, 1)

  net <- enaR::pack(
    flow = flow_mat, input = input_vec, export = export_vec,
    respiration = resp_vec, storage = storage_vec, living = living_vec
  )
  flow_res <- enaR::enaFlow(net, balance.override = TRUE)
  asc_res <- enaR::enaAscendency(net)
  tro_res <- enaR::enaTroAgg(net, balance.override = TRUE)

  harvest_by_level <- as.vector(tro_res$A[seq_along(liv), seq_along(liv), drop = FALSE] %*% catch[liv])
  ls <- tro_res$LS
  thr1_ewe <- sum(P[incl_type == 1]) + sum(flow_mat[det, liv]) + sum(det_input[det])
  te <- (harvest_by_level[1:2] + ls[2:3]) / c(thr1_ewe, ls[2])
  te_enar <- tro_res$TE[1:2]
  fci <- unname(flow_res$ns[, "FCI"])
  a_c <- unname(asc_res[, "ASC.CAP"])
  o_c <- unname(asc_res[, "OH.CAP"])

  if (return_spine) {
    return(list(tro = tro_res, harvest_by_level = harvest_by_level,
                biomass_by_level = as.vector(tro_res$A[seq_along(liv), seq_along(liv), drop = FALSE] %*% B[liv]),
                thr1_ewe = thr1_ewe))
  }

  ewe <- ewe_fci_ac(flow_mat, import_in + ifelse(incl_type == 2, det_input, 0), catch, resp_vec)

  data.frame(
    Model = name, FCI = round(ewe[["FCI"]], 4), A_C = round(ewe[["A_C"]], 4), O_C = round(ewe[["O_C"]], 4),
    TE_I = round(te[1], 4), TE_II = round(te[2], 4),
    TE_I_enaR = round(te_enar[1], 4), TE_II_enaR = round(te_enar[2], 4),
    FCI_enaR = round(fci, 4), A_C_enaR = round(a_c, 4), O_C_enaR = round(o_c, 4),
    stringsAsFactors = FALSE
  )
}

ewe_fci_ac <- function(flow, imp, ex, resp) {
  Q <- rowSums(flow) + ex + resp
  Tin <- colSums(flow) + imp
  tru <- sum(Q) + sum(imp)
  term <- function(f, den) ifelse(f > 0 & den > 0, f * log(f * tru / den), 0)
  ovh <- function(f, den) ifelse(f > 0 & den > 0, -f * log(f^2 / den), 0)
  dens <- list(list(imp, sum(imp) * Q), list(flow, outer(Q, Tin)), list(ex, sum(ex) * Tin), list(resp, sum(resp) * Tin))
  asc <- sum(sapply(dens, function(d) sum(term(d[[1]], d[[2]]))))
  oh <- sum(sapply(dens, function(d) sum(ovh(d[[1]], d[[2]]))))
  f_all <- c(imp, flow, ex, resp)
  f_all <- f_all[f_all > 0]
  cap <- -sum(f_all * log(f_all / tru))
  n <- solve(diag(length(Q)) - flow / ifelse(Q > 0, Q, 1))
  cycled <- sum(ifelse(diag(n) > 0, Q * (diag(n) - 1) / diag(n), 0))
  c(FCI = cycled / tru, A_C = asc / cap, O_C = oh / cap)
}

# Table 1 -----

table1_raw <- do.call(rbind, lapply(names(models_main), function(nm) {
  m <- models_main[[nm]]
  tmp <- tempfile(fileext = ".csv")
  write.Rpath(m, file = tmp)
  df <- read.csv(tmp)
  df$Model <- nm
  df[, c("Model", "Group", "type", "TL", "Biomass", "PB", "QB", "EE", "GE", "Removals")]
}))

table1 <- table1_raw |>
  filter(type %in% c(0, 1)) |>
  select(Model, Group, TL, Biomass, PB, QB, EE, GE, Removals) |>
  mutate(
    Period_days = ifelse(startsWith(Model, "sw"), SW_CULTURE_DAYS, FW_CULTURE_DAYS),
    Unit        = "g WW m-2 period-1",
    across(where(is.numeric), \(x) round(x, 3))
  ) |>
  select(Model, Group, Period_days, Unit, TL, Biomass, PB, QB, EE, GE, Removals)

writexl::write_xlsx(table1, here::here("output/Table1_basic_estimates.xlsx"))

# Tables 2 and 3 -----

table_labels <- c(sw_field = "Table2", fw_field = "Table3")
system_labels <- c(sw_field = "seawater", fw_field = "freshwater")

for (nm in names(models_main)) {
  if (!nm %in% names(params_main)) next
  m <- models_main[[nm]]
  dc_df <- extract_diet_matrix(m, params_main[[nm]])
  cons_cols <- as.character(m$Group[m$type == 0])
  dc_ms <- cbind(Prey = rownames(dc_df), dc_df[, cons_cols, drop = FALSE])
  nonzero <- rowSums(dc_ms[, -1, drop = FALSE], na.rm = TRUE) > 0
  dc_ms <- dc_ms[nonzero, ]
  fname <- paste0(here::here("output/"), table_labels[[nm]], "_diet_", system_labels[[nm]], ".xlsx")
  writexl::write_xlsx(as.data.frame(dc_ms), fname)
}

# Tables 4 and 5 -----

macro_real_lookup <- list(
  sw_field = NULL,
  fw_field = NULL,
  sw_grad_1x_realistic = list(group = "Red_algae", b_real = sw_b_steps[1], pb_real = PB_RED_ALGAE),
  sw_grad_1.5x_realistic = list(group = "Red_algae", b_real = sw_b_steps[2], pb_real = PB_RED_ALGAE),
  sw_grad_2x_realistic = list(group = "Red_algae", b_real = sw_b_steps[3], pb_real = PB_RED_ALGAE),
  sw_grad_2.5x_realistic = list(group = "Red_algae", b_real = sw_b_steps[4], pb_real = PB_RED_ALGAE),
  fw_grad_1x_realistic = list(group = "Duckweed", b_real = fw_b_steps[1], pb_real = PB_duckweed),
  fw_grad_1.5x_realistic = list(group = "Duckweed", b_real = fw_b_steps[2], pb_real = PB_duckweed),
  fw_grad_2x_realistic = list(group = "Duckweed", b_real = fw_b_steps[3], pb_real = PB_duckweed),
  fw_grad_2.5x_realistic = list(group = "Duckweed", b_real = fw_b_steps[4], pb_real = PB_duckweed)
)

table5_raw <- do.call(rbind, lapply(names(models_list), function(nm) {
  p <- if (nm %in% names(params_list)) params_list[[nm]] else NULL
  calc_ecosystem_indices(models_list[[nm]], p, nm, macro_real_lookup[[nm]])
}))

redalgae.duckweed_b <- c(
  sw_field = sw_b_steps[1],
  fw_field = fw_b_steps[1],
  sw_grad_1x_realistic = sw_b_steps[1], sw_grad_1.5x_realistic = sw_b_steps[2],
  sw_grad_2x_realistic = sw_b_steps[3], sw_grad_2.5x_realistic = sw_b_steps[4],
  fw_grad_1x_realistic = fw_b_steps[1], fw_grad_1.5x_realistic = fw_b_steps[2],
  fw_grad_2x_realistic = fw_b_steps[3], fw_grad_2.5x_realistic = fw_b_steps[4]
)

feed_supply_cut_g <- function(dc_macro, is_sw) ifelse(is_sw,
  sw_feed_detinput * dc_macro / DC_MUD_FEED,
  (fw_feed_detinput - fw_fish_feed_detinput) * dc_macro / DC_TIL_FEED)
redalgae.duckweed_supply_cut_g <- c(
  sw_field = 0, fw_field = 0,
  setNames(vapply(sw_real, \(r) feed_supply_cut_g(r$dc_macro, TRUE), numeric(1)), paste0("sw_grad_", GRAD_MULT, "x_realistic")),
  setNames(vapply(fw_real, \(r) feed_supply_cut_g(r$dc_macro, FALSE), numeric(1)), paste0("fw_grad_", GRAD_MULT, "x_realistic"))
)

redalgae.duckweed_displaced_g <- c(
  sw_field = 0, fw_field = 0,
  sw_grad_1x_realistic = sw_real[[1]]$displaced_g, sw_grad_1.5x_realistic = sw_real[[2]]$displaced_g,
  sw_grad_2x_realistic = sw_real[[3]]$displaced_g, sw_grad_2.5x_realistic = sw_real[[4]]$displaced_g,
  fw_grad_1x_realistic = fw_real[[1]]$displaced_g, fw_grad_1.5x_realistic = fw_real[[2]]$displaced_g,
  fw_grad_2x_realistic = fw_real[[3]]$displaced_g, fw_grad_2.5x_realistic = fw_real[[4]]$displaced_g
)

table5 <- table5_raw |>
  mutate(
    RedAlgae_Duckweed_B = round(redalgae.duckweed_b[Model], 2),
    System = ifelse(startsWith(Model, "sw"), "Seawater", "Freshwater"),
    Scenario = sub("_realistic$", "", sub("^(sw|fw)_", "", Model)),
    area_m2 = ifelse(startsWith(Model, "sw"), SW_SYSTEM_AREA_M2, FW_SYSTEM_AREA_M2),
    ef_low = ifelse(startsWith(Model, "sw"), EF_SW_LOW, EF_FW_LOW),
    ef_high = ifelse(startsWith(Model, "sw"), EF_SW_HIGH, EF_FW_HIGH),
    ef_mid = ifelse(startsWith(Model, "sw"), EF_SW_MID, EF_FW_MID),
    feed_avoided_kg = redalgae.duckweed_supply_cut_g[Model] * area_m2 / 1000,
    GHG_saved_kg_CO2e = round(feed_avoided_kg * ef_mid, 1),
    GHG_low_kg_CO2e = round(feed_avoided_kg * ef_low, 1),
    GHG_high_kg_CO2e = round(feed_avoided_kg * ef_high, 1)
  )

network_idx <- do.call(rbind, lapply(names(models_list), function(nm) {
  calc_network_indices(models_list[[nm]], params_list[[nm]], nm)
}))

table5 <- table5 |>
  left_join(network_idx, by = "Model") |>
  select(
    Model, System, Scenario, RedAlgae_Duckweed_B,
    Feed_DetInput, Total_Biomass, Total_Catch,
    Total_System_Throughput, PR_ratio, DH_ratio,
    Q_TST, Ex_TST, PP_TST,
    TE_I, TE_II, TE_I_enaR, TE_II_enaR,
    FCI, A_C, O_C, System_Omnivory_Index, Connectance_Index, mTLco,
    GHG_saved_kg_CO2e, GHG_low_kg_CO2e, GHG_high_kg_CO2e
  )

writexl::write_xlsx(table5, here::here("output/Table5_ecosystem_indices.xlsx"))
print(table5)

tableS1 <- network_idx |>
  filter(Model %in% names(models_main)) |>
  transmute(Model,
    TE_I.EwE = TE_I, TE_I.enaR = TE_I_enaR, TE_II.EwE = TE_II, TE_II.enaR = TE_II_enaR,
    FCI.EwE = FCI, FCI.enaR = FCI_enaR, A_C.EwE = A_C, A_C.enaR = A_C_enaR, O_C.EwE = O_C, O_C.enaR = O_C_enaR
  ) |>
  pivot_longer(-Model, names_to = c("Index", ".value"), names_sep = "\\.") |>
  mutate(System = ifelse(startsWith(Model, "sw"), "Seawater", "Freshwater")) |>
  select(System, Index, EwE, enaR)
writexl::write_xlsx(tableS1, here::here("output/TableS1_EwE_vs_enaR.xlsx"))
print(tableS1)

table4_gradient <- table5 |>
  mutate(
    field_feed_ref = ifelse(System == "Seawater", sw_feed_detinput, fw_feed_detinput),
    Feed_pct_reduction = round(redalgae.duckweed_supply_cut_g[Model] / field_feed_ref * 100, 2)
  ) |>
  select(
    Model, System, Scenario, RedAlgae_Duckweed_B,
    Feed_DetInput,
    Feed_pct_reduction,
    GHG_saved_kg_CO2e, GHG_low_kg_CO2e, GHG_high_kg_CO2e
  )
writexl::write_xlsx(table4_gradient, here::here("output/Table4_gradient.xlsx"))

# Figure 1 -----

drop_import <- function(model) {
  if ("Import" %in% rownames(model$DC)) {
    model$DC <- model$DC[rownames(model$DC) != "Import", , drop = FALSE]
  }
  model
}

get_lane_map <- function(model) {
  pm <- Rpath:::summarize.for.webplot(drop_import(model), fleets = FALSE)$pointmap
  setNames(rank(pm$x.pos, ties.method = "first"), pm$Group)
}

foodweb_labels <- c(
  Mud_crab = "Mud crab", Carnivores_inverts = "Carniv. inv.", Omnivores_inverts = "Omniv. inv.",
  Detritivores_inverts = "Detriv. inv.", Filter_feeders = "Filter fdr.", Zooplankton = "Zooplankton",
  Red_algae = "Red algae", Phytoplankton = "Phytoplankton", Benthic_algae = "Benthic algae",
  Commercial_feed = "Comm. feed", Sediment_detritus = "Sed. detritus", Water_detritus = "Water det.",
  Sea_bass = "Sea bass", Tilapia = "Tilapia", Duckweed = "Duckweed", Periphyton = "Periphyton",
  Benthic_invert = "Benthic inv."
)

plot_foodweb_biomass <- function(model, scenario_label = "", lane_map = NULL) {
  web.info <- Rpath:::summarize.for.webplot(drop_import(model), fleets = FALSE)
  pointmap <- web.info$pointmap
  if (!is.null(lane_map)) pointmap$x.pos <- lane_map[pointmap$Group]

  groups <- model$Group
  type_vec <- model$type
  consumer_idx <- which(type_vec == 0)
  Q <- ifelse(!is.na(model$QB), model$Biomass * model$QB, 0)
  DC <- drop_import(model)$DC
  dc_type <- type_vec[match(rownames(DC), groups)]

  conn_rows <- lapply(seq_along(consumer_idx), function(j) {
    pred_i <- consumer_idx[j]
    dc_col <- DC[, groups[pred_i]]
    keep <- !is.na(dc_col) & dc_col > 0 & dc_type != 3
    prey_i <- match(rownames(DC)[keep], groups)
    if (length(prey_i) == 0) {
      return(NULL)
    }
    data.frame(
      predator = groups[pred_i], prey = groups[prey_i],
      flow = Q[pred_i] * dc_col[keep],
      stringsAsFactors = FALSE
    )
  })
  connections <- do.call(rbind, Filter(Negate(is.null), conn_rows))

  pm <- pointmap[, c("Group", "x.pos", "TL")]
  connections <- merge(connections, pm, by.x = "predator", by.y = "Group")
  names(connections)[names(connections) == "x.pos"] <- "pred.x"
  names(connections)[names(connections) == "TL"] <- "pred.y"
  connections <- merge(connections, pm,
    by.x = "prey", by.y = "Group",
    suffixes = c("", ".prey")
  )
  names(connections)[names(connections) == "x.pos"] <- "prey.x"
  names(connections)[names(connections) == "TL"] <- "prey.y"
  connections$log_flow <- log10(pmax(connections$flow, 1e-3))
  connections$flow_type <- factor(
    ifelse(
      connections$prey == "Commercial_feed", "feed",
      ifelse(connections$prey %in% c("Red_algae", "Duckweed"), "redalgae.duckweed", "other")
    ),
    levels = c("feed", "redalgae.duckweed", "other")
  )
  missing_flow <- setdiff(levels(connections$flow_type), as.character(connections$flow_type))
  if (length(missing_flow) > 0) {
    connections <- rbind(connections, data.frame(
      prey = NA_character_, predator = NA_character_, flow = NA_real_,
      pred.x = 0, pred.y = -100, prey.x = 0, prey.y = -100,
      log_flow = 0, flow_type = factor(missing_flow, levels = levels(connections$flow_type))
    ))
  }

  type_cols <- c("0" = "#0072B2", "1" = "#009E73", "2" = "#9AA0A6", "feed" = "#E69F00")
  pointmap$logB <- log10(pmax(pointmap$Biomass, 0.001))
  pointmap$type <- ifelse(pointmap$Group == "Commercial_feed", "feed", as.character(pointmap$type))
  pointmap$label <- foodweb_labels[pointmap$Group]

  ggplot() +
    geom_hline(
      yintercept = 1:3, linetype = "dotted",
      colour = "grey88", linewidth = 0.3
    ) +
    geom_segment(
      data = connections,
      aes(
        x = prey.x, xend = pred.x, y = prey.y, yend = pred.y,
        linewidth = log_flow, colour = flow_type
      ),
      alpha = 0.75,
      arrow = arrow(length = unit(0.09, "cm"), type = "closed")
    ) +
    scale_linewidth_continuous(range = c(0.15, 1.8), guide = "none") +
    scale_colour_manual(
      name = "Flow",
      values = c("feed" = "#E69F00", "redalgae.duckweed" = "#009E73", "other" = "#C4C9CF"),
      labels = c("feed" = "Commercial feed", "redalgae.duckweed" = "Red algae/Duckweed", "other" = "Other"),
      drop = FALSE,
      guide = guide_legend(override.aes = list(linewidth = 1.3))
    ) +
    geom_point(
      data = pointmap,
      aes(x = x.pos, y = TL, size = logB, fill = type),
      pch = 21, colour = "white", stroke = 0.6
    ) +
    ggrepel::geom_text_repel(
      data = subset(pointmap, type == "0"),
      aes(x = x.pos, y = TL, label = label),
      size = 2.4,
      colour = "#333333",
      family = "Century Gothic",
      fontface = "bold",
      seed = 42,
      force = 20,
      force_pull = 0.4,
      nudge_y = 0.35,
      max.time = 3,
      max.iter = 200000,
      xlim = grDevices::extendrange(pointmap$x.pos, f = 0.22),
      ylim = c(1.55, 4.3),
      min.segment.length = 0,
      segment.size = 0.25,
      segment.color = "grey50",
      box.padding = 0.7,
      point.padding = 0.3,
      max.overlaps = Inf,
      bg.color = "white",
      bg.r = 0.12
    ) +
    ggrepel::geom_text_repel(
      data = subset(pointmap, type != "0"),
      aes(x = x.pos, y = TL, label = label),
      size = 2.4,
      colour = "#333333",
      family = "Century Gothic",
      fontface = "bold",
      seed = 42,
      force = 20,
      force_pull = 0.4,
      nudge_y = -0.35,
      max.time = 3,
      max.iter = 200000,
      xlim = grDevices::extendrange(pointmap$x.pos, f = 0.22),
      ylim = c(-1.1, 1.45),
      min.segment.length = 0,
      segment.size = 0.25,
      segment.color = "grey50",
      box.padding = 0.7,
      point.padding = 0.3,
      max.overlaps = Inf,
      bg.color = "white",
      bg.r = 0.12
    ) +
    scale_fill_manual(
      name = "Type",
      values = type_cols,
      labels = c("0" = "Consumer", "1" = "Producer", "2" = "Detritus", "feed" = "Commercial feed"),
      breaks = c("0", "1", "2", "feed")
    ) +
    scale_size_continuous(range = c(3, 12), guide = "none") +
    scale_x_continuous(expand = expansion(mult = 0.22)) +
    coord_cartesian(ylim = c(-1.1, 4.1)) +
    labs(x = NULL, y = "Trophic level", subtitle = scenario_label) +
    theme_classic(base_size = 13, base_family = "Century Gothic") +
    theme(
      text = element_text(face = "bold"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      plot.subtitle = element_text(size = 10, colour = "grey40", face = "bold")
    )
}

sw_lanes <- get_lane_map(models_main[["sw_field"]])
fw_lanes <- get_lane_map(models_main[["fw_field"]])

field_plots <- list(
  plot_foodweb_biomass(models_main[["sw_field"]], lane_map = sw_lanes),
  plot_foodweb_biomass(models_main[["fw_field"]], lane_map = fw_lanes)
)
grad_plots <- list(
  plot_foodweb_biomass(models_list[["sw_grad_2.5x_realistic"]],
    scenario_label = "2.5× red algae", lane_map = sw_lanes
  ),
  plot_foodweb_biomass(models_list[["fw_grad_2.5x_realistic"]],
    scenario_label = "2.5× duckweed", lane_map = fw_lanes
  )
)

fig1 <- wrap_plots(c(field_plots, grad_plots), ncol = 2, nrow = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 12, family = "Century Gothic")
  )

ggsave(here::here("output/Figure1_foodweb_diagrams.png"), fig1, width = 250, height = 150, units = "mm", dpi = 300)
ggsave(here::here("output/Figure1_foodweb_diagrams.pdf"), fig1, width = 250, height = 150, units = "mm", device = cairo_pdf)

# Figure 2 -----
spine_col <- list(harv = "#2a78d6", resp = "#eb6834", ink = "#0b0b0b", ink2 = "#52514e",
                  box = "#f4f4f2", det = "#f1ebe1", line = "#52514e")
spine_lw <- function(v) pmin(0.3 + 0.24 * log10(pmax(v, 1)), 1.4)
spine_fmt <- function(v) vapply(v, function(z) {
  if (z >= 100) format(round(z), big.mark = ",")
  else if (z >= 10) formatC(z, format = "f", digits = 1)
  else if (z >= 0.01) formatC(z, format = "f", digits = 2)
  else "<0.01"
}, "")

spine_data <- function(nm) {
  sp <- calc_network_indices(models_main[[nm]], params_main[[nm]], nm, return_spine = TRUE)
  ls <- sp$tro$LS; h <- sp$harvest_by_level
  sp$te <- c((h[1] + ls[2]) / sp$thr1_ewe, (h[2:3] + ls[3:4]) / ls[2:3])
  sp
}

gg_spine <- function(sp, tag, nlev_axis) {
  tro <- sp$tro; ns <- as.data.frame(tro$ns); C <- spine_col
  nlev <- sum(tro$LS > 0)
  bw <- 2.0; bh <- 1.5; gap <- 1.6; yT <- 4.2; yD <- 1.6; yBus <- 0.25; tsz <- 2.8
  xk <- (bw + gap) * (seq_len(nlev) - 1) + 1.5
  xD <- xk[1]; top <- yT + bh / 2; bot <- yT - bh / 2
  ar <- arrow(length = unit(1.5, "mm"), type = "closed", angle = 25)
  S <- list(); Tx <- list()
  sg <- function(x, y, xe, ye, col = C$line, w = 0.4, arrow = FALSE)
    S[[length(S) + 1]] <<- data.frame(x, y, xe, ye, col, w, arrow)
  tx <- function(x, y, l, h = 0.5, v = 0.5, col = C$ink2, face = "plain")
    Tx[[length(Tx) + 1]] <<- data.frame(x, y, l, h, v, col, face)
  for (k in seq_len(nlev)) {
    if (k <= length(sp$te)) tx(xk[k] + bw / 2 - 0.1, top - 0.12, sprintf("%.2f%%", 100 * sp$te[k]), 1, 1, C$ink, "bold")
    tx(xk[k] - bw / 2 + 0.1, bot + 0.1, spine_fmt(sp$biomass_by_level[k]), 0, 0)
    hk <- sp$harvest_by_level[k]
    if (hk > 0) { sg(xk[k] + 0.35, top, xk[k] + 0.35, top + 0.75, C$harv, spine_lw(hk), TRUE)
                  tx(xk[k] + 0.35, top + 0.95, spine_fmt(hk), 0.5, 0) }
    if (k > 1 && tro$CR[k] > 0) { sg(xk[k] - 0.35, bot, xk[k] - 0.35, bot - 0.75, C$resp, spine_lw(tro$CR[k]), TRUE)
                                  tx(xk[k] - 0.35, bot - 0.95, spine_fmt(tro$CR[k]), 0.5, 1) }
    if (k < nlev) { sg(xk[k] + bw / 2, yT, xk[k + 1] - bw / 2, yT, C$ink, spine_lw(tro$GC[k + 1]), TRUE)
                    tx((xk[k] + xk[k + 1]) / 2, yT + 0.14, spine_fmt(tro$GC[k + 1]), 0.5, 0, C$ink) }
    if (k > 1) { xr <- xk[k] + bw / 2 + 0.3; w <- spine_lw(tro$RDP[k])
                 sg(xk[k] + bw / 2, yT - 0.35, xr, yT - 0.35, w = w); sg(xr, yT - 0.35, xr, yBus, w = w)
                 tx(xr - 0.1, yBus + 0.12, spine_fmt(tro$RDP[k]), 1, 0) }
  }
  sg(xD, bot, xD, yD + bh / 2, C$line, spine_lw(tro$RDP[1]), TRUE)
  tx(xD + 0.12, (bot + yD + bh / 2) / 2, spine_fmt(tro$RDP[1]), 0, 0.5)
  sg(xD - bw / 2 - 0.9, yD + bh / 2 + 0.5, xD - bw / 2, yD + bh / 2 - 0.1, C$ink, spine_lw(ns$DetritalInput), TRUE)
  tx(xD - bw / 2 - 0.95, yD + bh / 2 + 0.55, spine_fmt(ns$DetritalInput), 1, 0)
  xj <- xk[2] - bw / 2 - 0.35; wd <- spine_lw(ns$Detritivory)
  sg(xD + bw / 2, yD, xj, yD, C$ink, wd); sg(xj, yD, xj, bot + 0.3, C$ink, wd)
  sg(xj, bot + 0.3, xk[2] - bw / 2, bot + 0.3, C$ink, wd, TRUE)
  tx(xj + 0.1, yD + 0.12, spine_fmt(ns$Detritivory), 0, 0, C$ink)
  xr_last <- xk[nlev] + bw / 2 + 0.3
  sg(xr_last, yBus, xD, yBus); sg(xD, yBus, xD, yD - bh / 2, arrow = TRUE)
  S <- do.call(rbind, S); Tx <- do.call(rbind, Tx)
  boxes <- data.frame(x = c(xk, xD), y = c(rep(yT, nlev), yD),
                      lab = c("P", as.character(as.roman(2:nlev)), "D"), fill = c(rep(C$box, nlev), C$det))
  x_right <- (bw + gap) * (max(nlev, nlev_axis) - 1) + 1.5 + bw / 2 + 0.9
  ggplot() +
    geom_segment(data = subset(S, !arrow), aes(x, y, xend = xe, yend = ye, colour = I(col), linewidth = I(w)),
                 lineend = "butt", linejoin = "mitre") +
    geom_segment(data = subset(S, arrow), aes(x, y, xend = xe, yend = ye, colour = I(col), linewidth = I(w)),
                 arrow = ar, lineend = "butt", linejoin = "mitre") +
    geom_rect(data = boxes, aes(xmin = x - bw / 2, xmax = x + bw / 2, ymin = y - bh / 2, ymax = y + bh / 2,
                                fill = I(fill)), colour = C$ink, linewidth = 0.35) +
    geom_text(data = boxes, aes(x, y + 0.05, label = lab), family = "Century Gothic", size = 4.4, colour = C$ink) +
    geom_text(data = Tx, aes(x, y, label = l, hjust = h, vjust = v, colour = I(col), fontface = face),
              family = "Century Gothic", size = tsz) +
    coord_fixed(xlim = c(-1.6, x_right), ylim = c(-0.1, top + 1.4), expand = FALSE, clip = "off") +
    labs(tag = tag) + theme_void(base_family = "Century Gothic") +
    theme(plot.tag = element_text(face = "bold", size = 10, colour = C$ink), plot.tag.position = c(0, 1),
          plot.margin = margin(4, 6, 4, 6))
}

gg_spine_key <- function() {
  C <- spine_col; ar <- arrow(length = unit(1.5, "mm"), type = "closed", angle = 25)
  t <- function(x, y, l, ...) annotate("text", x, y, label = l, family = "Century Gothic", size = 2.8, colour = C$ink, ...)
  ggplot() +
    geom_rect(aes(xmin = 0, xmax = 2, ymin = 0, ymax = 1.5), fill = C$box, colour = C$ink, linewidth = 0.35) +
    annotate("text", 1, 0.75, label = "TL", family = "Century Gothic", size = 4.4, colour = C$ink) +
    t(1.92, 1.38, "TE %", hjust = 1, vjust = 1, fontface = "bold") + t(0.08, 0.1, "Biomass", hjust = 0, vjust = 0) +
    annotate("segment", x = 1.35, y = 1.5, xend = 1.35, yend = 2.2, colour = C$harv, linewidth = 0.9, arrow = ar) +
    t(1.35, 2.3, "Harvest", vjust = 0) +
    annotate("segment", x = 0.65, y = 0, xend = 0.65, yend = -0.7, colour = C$resp, linewidth = 0.9, arrow = ar) +
    t(0.65, -0.8, "Respiration", vjust = 1) +
    annotate("segment", x = 2, y = 0.75, xend = 2.8, yend = 0.75, colour = C$line, arrow = ar) +
    t(2.9, 0.75, "Predation", hjust = 0) +
    annotate("segment", x = 2, y = 0.4, xend = 2.3, yend = 0.4, colour = C$line) +
    annotate("segment", x = 2.3, y = 0.4, xend = 2.3, yend = -0.3, colour = C$line) +
    t(2.4, -0.3, "Return to detritus", hjust = 0) +
    annotate("segment", x = -0.8, y = 1.9, xend = 0, yend = 1.45, colour = C$line, arrow = ar) +
    t(-0.85, 1.95, "Import", hjust = 1, vjust = 0) +
    coord_fixed(xlim = c(-2, 5), ylim = c(-1.3, 2.8), expand = FALSE, clip = "off") + theme_void()
}

spine <- list(sw = spine_data("sw_field"), fw = spine_data("fw_field"))

fig2 <- (gg_spine(spine$sw, "(a)", 4) / gg_spine(spine$fw, "(b)", 4)) | gg_spine_key()
fig2 <- fig2 + plot_layout(widths = c(1, 0.36))
ggsave(here::here("output/Figure2_lindeman_spine.png"), fig2, width = 10.8, height = 6.6, dpi = 300, bg = "white")
ggsave(here::here("output/Figure2_lindeman_spine.pdf"), fig2, width = 10.8, height = 6.6, device = cairo_pdf)

# Figure 3 -----
plot_mti_filtered <- function(model, params, scenario_label = "") {
  mti_mat <- mti_fixed(model, params)
  all_groups <- model$Group
  dimnames(mti_mat) <- list(all_groups, all_groups)
  tl_vals <- setNames(model$TL, all_groups)

  living_grp <- all_groups[model$type %in% c(0, 1)]
  mti_sub <- mti_mat[living_grp, living_grp]

  order_asc <- names(sort(tl_vals[living_grp]))
  order_desc <- names(sort(tl_vals[living_grp], decreasing = TRUE))

  mti_long <- as.data.frame(as.table(mti_sub)) |>
    setNames(c("Impacting", "Impacted", "MTI")) |>
    mutate(
      MTI = as.numeric(MTI),
      Impacting = factor(unname(foodweb_labels[as.character(Impacting)]),
        levels = unname(foodweb_labels[order_asc])
      ),
      Impacted = factor(unname(foodweb_labels[as.character(Impacted)]),
        levels = unname(foodweb_labels[order_desc])
      ),
    )

  ggplot(mti_long, aes(x = Impacting, y = Impacted)) +
    geom_tile(aes(fill = MTI), colour = "#ececea", linewidth = 0.3) +
    geom_text(
      data = \(d) dplyr::filter(d, abs(MTI) >= 0.20),
      aes(label = sub("^(-?)0\\.", "\\1.", sprintf("%.2f", MTI))),
      size = 2.2, colour = "grey15",
      family = "Century Gothic", fontface = "bold"
    ) +
    scale_fill_gradient2(
      low = "#D55E00", mid = "white", high = "#0072B2",
      midpoint = 0, name = "MTI",
      limits = c(-1, 1), oob = scales::squish,
      na.value = "grey96"
    ) +
    coord_equal() +
    labs(subtitle = scenario_label, x = "Impacting group", y = "Impacted group") +
    theme_classic(base_size = 10, base_family = "Century Gothic") +
    theme(
      text = element_text(face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "bold"),
      axis.text.y = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 10, colour = "grey40", face = "bold")
    )
}

mti_grad_models <- list(
  sw_grad_2.5x_realistic = models_list[["sw_grad_2.5x_realistic"]],
  fw_grad_2.5x_realistic = models_list[["fw_grad_2.5x_realistic"]]
)
mti_grad_params <- list(
  sw_grad_2.5x_realistic = params_list[["sw_grad_2.5x_realistic"]],
  fw_grad_2.5x_realistic = params_list[["fw_grad_2.5x_realistic"]]
)
mti_grad_labels <- c(
  sw_grad_2.5x_realistic = "Seawater 2.5× red algae",
  fw_grad_2.5x_realistic = "Freshwater 2.5× duckweed"
)

mti_plots <- lapply(names(mti_grad_models), function(nm) {
  plot_mti_filtered(mti_grad_models[[nm]], mti_grad_params[[nm]], scenario_label = mti_grad_labels[[nm]])
})

mti_field_plots <- list(
  plot_mti_filtered(models_main$sw_field, params_main$sw_field, scenario_label = "Seawater field"),
  plot_mti_filtered(models_main$fw_field, params_main$fw_field, scenario_label = "Freshwater field")
)
fig3 <- wrap_plots(c(mti_field_plots, mti_plots), ncol = 2, nrow = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 12, family = "Century Gothic"))

ggsave(here::here("output/Figure3_MTI_heatmaps.png"), fig3, width = 260, height = 270, units = "mm", dpi = 300)
ggsave(here::here("output/Figure3_MTI_heatmaps.pdf"), fig3, width = 260, height = 270, units = "mm", device = cairo_pdf)

# MTI matrices -----
for (nm in names(models_main)) {
  if (!nm %in% names(params_main)) next
  mti_mat <- mti_fixed(models_main[[nm]], params_main[[nm]])
  groups <- models_main[[nm]]$Group
  mti_df <- as.data.frame(round(mti_mat, 4))
  colnames(mti_df) <- groups
  rownames(mti_df) <- groups
  mti_df <- cbind(Group = rownames(mti_df), mti_df)
  writexl::write_xlsx(mti_df, paste0(here::here("output/MTI_matrix_"), nm, ".xlsx"))
}

# Table S2 -----

pedigree_sw <- data.frame(
  Model = "sw_field",
  Group = c(
    "Mud_crab", "Carnivores_inverts", "Omnivores_inverts",
    "Detritivores_inverts", "Filter_feeders", "Zooplankton",
    "Red_algae", "Phytoplankton", "Benthic_algae"
  ),
  B = c(1.0, 0.7, 0.7, 0.7, 0.7, 0.4, 1.0, 0.4, 0.4),
  PB = c(1.0, 0.5, 0.5, 0.5, 0.5, 0.2, 1.0, 1.0, 1.0),
  QB = c(0.6, 0.2, 0.2, 0.2, 0.2, 0.2, NA, NA, NA),
  Diet = c(0.5, 0.2, 0.2, 0.2, 0.2, 0.5, NA, NA, NA),
  Catch = c(1.0, NA, NA, NA, NA, NA, 1.0, NA, NA),
  stringsAsFactors = FALSE
)

pedigree_fw <- data.frame(
  Model = "fw_field",
  Group = c(
    "Sea_bass", "Tilapia", "Zooplankton",
    "Duckweed", "Phytoplankton", "Periphyton", "Benthic_invert"
  ),
  B = c(1.0, 1.0, 0.4, 1.0, 0.4, 0.4, 0.0),
  PB = c(1.0, 1.0, 0.2, 0.8, 1.0, 1.0, 0.6),
  QB = c(0.5, 0.2, 0.2, NA, NA, NA, 0.2),
  Diet = c(0.5, 0.7, 0.5, NA, NA, NA, 0.0),
  Catch = c(1.0, 1.0, NA, 1.0, NA, NA, NA),
  stringsAsFactors = FALSE
)

pedigree_all <- rbind(pedigree_sw, pedigree_fw)

pedigree_all$Mean_group <- apply(
  pedigree_all[, c("B", "PB", "QB", "Diet", "Catch")],
  1, function(x) round(mean(x, na.rm = TRUE), 3)
)

mean_pedigree <- data.frame(
  Model = c("sw_field", "fw_field"),
  Mean_pedigree = c(
    round(mean(
      unlist(pedigree_all[
        pedigree_all$Model == "sw_field",
        c("B", "PB", "QB", "Diet", "Catch")
      ]),
      na.rm = TRUE
    ), 3),
    round(mean(
      unlist(pedigree_all[
        pedigree_all$Model == "fw_field",
        c("B", "PB", "QB", "Diet", "Catch")
      ]),
      na.rm = TRUE
    ), 3)
  )
)

writexl::write_xlsx(
  list(
    SW_field = pedigree_all[
      pedigree_all$Model == "sw_field",
      c(
        "Group", "B", "PB", "QB", "Diet", "Catch", "Mean_group"
      )
    ],
    FW_field = pedigree_all[
      pedigree_all$Model == "fw_field",
      c(
        "Group", "B", "PB", "QB", "Diet", "Catch", "Mean_group"
      )
    ]
  ),
  here::here("output/TableS2_pedigree.xlsx")
)
print(mean_pedigree)

# Sampling ranges -----

aydin_range_pb <- c("1" = 0.1, "0.8" = 0.2, "0.7" = 0.3, "0.6" = 0.4, "0.5" = 0.5, "0.2" = 0.6, "0.1" = 0.7, "0" = 0.8)
aydin_range_qb <- aydin_range_pb
aydin_range_b <- c(0.1, 0.1, 0.5, 0.5, 0.5, 0.8, 0.8, 0.8)
aydin_range_diet <- c(0.1, 0.3, 0.5, 0.6, 0.7, 0.8, 0.8, 0.8)

aydin_rank <- data.frame(
  Group = c(
    "Mud_crab", "Carnivores_inverts", "Omnivores_inverts", "Detritivores_inverts", "Filter_feeders",
    "Zooplankton", "Red_algae", "Phytoplankton", "Benthic_algae",
    "Sea_bass", "Tilapia", "Duckweed", "Periphyton", "Benthic_invert"
  ),
  Model = c(rep("sw_field", 9), rep("fw_field", 5)),
  B_rank = c(2, 3, 3, 3, 3, 5, 2, 5, 5, 2, 2, 2, 5, 6),
  Diet_rank = c(5, 6, 6, 6, 6, 5, NA, NA, NA, 5, 3, NA, NA, 7),
  stringsAsFactors = FALSE
)
aydin_rank <- rbind(
  aydin_rank,
  data.frame(Group = c("Zooplankton", "Phytoplankton"), Model = "fw_field", B_rank = c(5, 5), Diet_rank = c(5, NA))
)

pedigree_ranges <- merge(
  pedigree_all[, c("Model", "Group", "B", "PB", "QB", "Diet")],
  aydin_rank,
  by = c("Model", "Group"), all.x = TRUE, suffixes = c("_cw", "")
)
pedigree_ranges$PB_range <- unname(aydin_range_pb[as.character(pedigree_ranges$PB)])
pedigree_ranges$QB_range <- unname(aydin_range_qb[as.character(pedigree_ranges$QB)])
pedigree_ranges$B_range <- aydin_range_b[pedigree_ranges$B_rank]
pedigree_ranges$Diet_range <- aydin_range_diet[pedigree_ranges$Diet_rank]
pedigree_ranges <- pedigree_ranges[, c(
  "Model", "Group", "B_rank", "B_range", "PB", "PB_range", "QB", "QB_range", "Diet_rank", "Diet_range"
)]
names(pedigree_ranges)[names(pedigree_ranges) %in% c("PB", "QB")] <- c("PB_cw_score", "QB_cw_score")

writexl::write_xlsx(pedigree_ranges, here::here("output/TableS3_pedigree_ranges.xlsx"))
print(pedigree_ranges)

# Table S5 -----

ideal_real <- c(
  list(NULL, NULL),
  dc_steps(sw_b_steps, sw_dc_args, 1),
  dc_steps(fw_b_steps, fw_dc_args, 1)
)
ideal_pick <- function(field, f) setNames(
  vapply(seq_along(ideal_real), \(i) if (is.null(ideal_real[[i]])) field[i] else ideal_real[[i]][[f]], numeric(1)),
  names(redalgae.duckweed_displaced_g)
)
ideal_displaced_g <- ideal_pick(c(0, 0), "displaced_g")
ideal_dc_macro <- ideal_pick(c(0, 0), "dc_macro")
ideal_dc_feed <- ideal_pick(c(DC_MUD_FEED, DC_TIL_FEED), "dc_feed")

diet_frac <- function(nm) {
  is_sw <- startsWith(nm, "sw")
  pred <- if (is_sw) "Mud_crab" else "Tilapia"
  prey_macro <- if (is_sw) "Red_algae" else "Duckweed"
  dc <- extract_diet_matrix(models_list[[nm]], params_list[[nm]])
  data.frame(
    Model = nm,
    Commercial_feed = dc["Commercial_feed", pred],
    RedAlgae_Duckweed = dc[prey_macro, pred]
  )
}
diet_raw <- do.call(rbind, lapply(names(models_list), diet_frac))

tableS5 <- table4_gradient |>
  left_join(diet_raw, by = "Model") |>
  mutate(
    is_sw = System == "Seawater",
    ideal_cut_g = feed_supply_cut_g(ideal_dc_macro[Model], is_sw),
    Feed_pct_reduction_idealized = round(ideal_cut_g / ifelse(is_sw, sw_feed_detinput, fw_feed_detinput) * 100, 2),
    ideal_feed_kg = ideal_cut_g * ifelse(is_sw, SW_SYSTEM_AREA_M2, FW_SYSTEM_AREA_M2) / 1000,
    GHG_saved_idealized = round(ideal_feed_kg * ifelse(is_sw, EF_SW_MID, EF_FW_MID), 1),
    GHG_low_idealized = round(ideal_feed_kg * ifelse(is_sw, EF_SW_LOW, EF_FW_LOW), 1),
    GHG_high_idealized = round(ideal_feed_kg * ifelse(is_sw, EF_SW_HIGH, EF_FW_HIGH), 1),
    Commercial_feed_idealized = round(ideal_dc_feed[Model], 4),
    RedAlgae_Duckweed_idealized = round(ideal_dc_macro[Model], 4),
    Scenario = ifelse(Scenario == "field", "Field", sub("x$", "×", sub("grad_", "", Scenario)))
  ) |>
  select(
    System, Scenario, RedAlgae_Duckweed_B,
    Feed_pct_reduction_idealized, Feed_pct_reduction,
    GHG_saved_idealized, GHG_low_idealized, GHG_high_idealized,
    GHG_saved_kg_CO2e, GHG_low_kg_CO2e, GHG_high_kg_CO2e,
    Commercial_feed_idealized, RedAlgae_Duckweed_idealized,
    Commercial_feed, RedAlgae_Duckweed
  ) |>
  rename(
    Commercial_feed_realistic = Commercial_feed,
    RedAlgae_Duckweed_realistic = RedAlgae_Duckweed,
    Feed_pct_reduction_realistic = Feed_pct_reduction,
    GHG_saved_realistic_kg_CO2e = GHG_saved_kg_CO2e
  )
writexl::write_xlsx(tableS5, here::here("output/TableS5_idealized_vs_realistic.xlsx"))

# Figure S2 -----
ks_df <- do.call(rbind, lapply(names(models_main), \(nm) calc_keystoneness(models_main[[nm]], params_main[[nm]], nm)))
ks_df$System <- ifelse(startsWith(ks_df$Model, "sw"), "Seawater IMTA", "Freshwater IMTA")
ks_df$Group_label <- unname(foodweb_labels[ks_df$Group])

sys_cols <- c("Seawater IMTA" = "#0072B2", "Freshwater IMTA" = "#D55E00")

plot_keystoneness <- function(df, colour) {
  ggplot(df, aes(x = p_B, y = epsilon, label = Group_label)) +
    geom_point(size = 3.5, pch = 19, colour = colour) +
    ggrepel::geom_text_repel(
      colour = colour, size = 3, family = "Century Gothic", fontface = "bold", seed = 42,
      force = 3, box.padding = 0.4, min.segment.length = 0,
      segment.size = 0.25, max.overlaps = Inf
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.15))) +
    labs(x = "Relative biomass (p)", y = "Trophic impact (ε)") +
    theme_classic(base_size = 11, base_family = "Century Gothic") +
    theme(text = element_text(face = "bold"))
}

ks_plots <- lapply(names(sys_cols), function(sys) {
  plot_keystoneness(dplyr::filter(ks_df, System == sys), sys_cols[[sys]])
})

figS2 <- wrap_plots(ks_plots, ncol = 2) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")") &
  theme(plot.tag = element_text(face = "bold", size = 12, family = "Century Gothic"))

ggsave(here::here("output/FigureS2_keystoneness.png"), figS2, width = 190, height = 100, units = "mm", dpi = 300)

# GHG per kg -----

avoided_per_kg <- table4_gradient |>
  filter(Scenario != "field") |>
  mutate(
    water = ifelse(System == "Seawater", "sw", "fw"),
    plant_kg = RedAlgae_Duckweed_B * ifelse(water == "sw", SW_SYSTEM_AREA_M2, FW_SYSTEM_AREA_M2) / 1000,
    mid = GHG_saved_kg_CO2e / plant_kg, high = GHG_high_kg_CO2e / plant_kg
  )
avoided_per_kg <- summarise(avoided_per_kg, avoided_mid = mean(mid), avoided_high = mean(high), .by = water)

ghg_per_kg <- macro_do |>
  mutate(water = str_extract(site, "^sw|^fw")) |>
  summarise(nep = mean(nep, na.rm = TRUE), .by = c(water, year, site)) |>
  summarise(nep = mean(nep), .by = water) |>
  mutate(
    System = ifelse(water == "sw", "Seawater", "Freshwater"),
    day_culture = ifelse(water == "sw", SW_MACRO_DAYS, FW_MACRO_DAYS),
    CO2e_removed_per_kg = nep * mg_to_g * (1 / o2_molar_mass) * (1 / pq) *
      co2_mass * day_culture
  ) |>
  left_join(avoided_per_kg, by = "water") |>
  mutate(Removed_to_avoided = CO2e_removed_per_kg / avoided_mid) |>
  select(System, CO2e_removed_per_kg, CO2e_avoided_per_kg = avoided_mid, CO2e_avoided_high_per_kg = avoided_high,
         Removed_to_avoided)

print(ghg_per_kg)

# Mortality breakdown -----

mortality_tab <- do.call(rbind, lapply(names(models_main), function(nm) {
  extract_mortality(models_main[[nm]], nm)
}))

writexl::write_xlsx(mortality_tab, here::here("output/Mortality_breakdown.xlsx"))

# Checks -----
source(here::here("R", "EcopathIMTA_checks.R"))

# Uncertainty -----
source(here::here("R", "EcopathIMTA_uncertainty.R"))

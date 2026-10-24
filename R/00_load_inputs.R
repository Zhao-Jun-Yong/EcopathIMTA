read_input <- function(file) read.csv(here::here("data", file), stringsAsFactors = FALSE)

ecopath_df <- read_input("ecopath_parameters.csv")

macro_do <- read_input("red_algae_duckweed_production.csv")

info <- read_input("ponds_2024.csv")

constants <- read_input("derived_constants.csv")
const <- setNames(constants$value, constants$constant)
fw_fish_feed_detinput <- const[["fw_seabass_feed"]]

pq <- 1.25
co2_mass <- 44.01
o2_molar_mass <- 32
mg_to_g <- 1 / 1000

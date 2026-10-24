# Toy food chain -----
toy <- create.rpath.params(
  group = c("Producer", "Herbivore", "Carnivore", "Detritus", "Fleet"),
  type = c(1, 0, 0, 2, 3)
)
toy$model[, `:=`(
  Biomass = c(10, 5, 1, 1, NA), PB = c(100, 10, 5, NA, NA), QB = c(NA, 40, 20, NA, NA),
  EE = NA_real_, ProdCons = NA_real_, BioAcc = c(0, 0, 0, 0, NA), Unassim = c(0, 0.2, 0.2, 0, NA),
  DetInput = c(NA, NA, NA, 0, NA), Detritus = c(1, 1, 1, 1, 0), Fleet = c(0, 5, 2, 0, NA), Fleet.disc = c(0, 0, 0, 0, NA)
)]
toy$diet[, Herbivore := c(0.7, 0, 0, 0.3, 0)]
toy$diet[, Carnivore := c(0, 1, 0, 0, 0)]
toy$stanzas$NStanzaGroups <- 0
toy_model <- run_rpath_model(toy, "toy")
stopifnot(is_balanced(toy_model))

toy_net <- calc_network_indices(toy_model, toy, "toy")
toy_eco <- calc_ecosystem_indices(toy_model, toy, "toy")

toy_flow <- matrix(0, 4, 4, dimnames = list(toy$model$Group[1:4], toy$model$Group[1:4]))
toy_flow["Producer", c("Herbivore", "Detritus")] <- c(140, 860)
toy_flow["Herbivore", c("Carnivore", "Detritus")] <- c(20, 65)
toy_flow["Carnivore", "Detritus"] <- 7
toy_flow["Detritus", "Herbivore"] <- 60
toy_T <- c(1000, 200, 20, 932)
toy_N <- solve(diag(4) - toy_flow / toy_T)
toy_fci <- sum((diag(toy_N) - 1) / diag(toy_N) * toy_T) / sum(toy_T)

toy_DC <- matrix(0, 5, 5)
toy_DC[1, 2] <- 0.7; toy_DC[4, 2] <- 0.3; toy_DC[2, 3] <- 1; toy_DC[2, 5] <- 5 / 7; toy_DC[3, 5] <- 2 / 7
toy_FC <- matrix(0, 5, 5)
toy_FC[1, 2] <- 1; toy_FC[4, 2] <- 1; toy_FC[2, 3] <- 20 / 25; toy_FC[2, 5] <- 5 / 25; toy_FC[3, 5] <- 1
toy_M <- solve(diag(5) - (toy_DC - t(toy_FC))) - diag(5)
toy_M[4, 4] <- 0
toy_ML <- toy_M[1:3, 1:3]
diag(toy_ML) <- 0
toy_ks <- log10(sqrt(rowSums(toy_ML^2)) * (1 - c(10, 5, 1) / 16))

ewe_A <- matrix(0, 7, 7)
ewe_A[1, 2] <- 140; ewe_A[1, 4] <- 860; ewe_A[2, 3] <- 20; ewe_A[2, 4] <- 65; ewe_A[3, 4] <- 7; ewe_A[4, 2] <- 60
ewe_A[2, 6] <- 5; ewe_A[3, 6] <- 2; ewe_A[2, 7] <- 110; ewe_A[3, 7] <- 11
ewe_out <- rowSums(ewe_A); ewe_in <- colSums(ewe_A)
ewe_N <- solve(diag(4) - ewe_A[1:4, 1:4] / ewe_out[1:4])
toy_fci_ewe <- sum(ewe_out[1:4] * (diag(ewe_N) - 1) / diag(ewe_N)) / 1280
ewe_asc <- 0
for (k in which(ewe_A > 0)) {
  i <- row(ewe_A)[k]; j <- col(ewe_A)[k]; f <- ewe_A[k]
  den <- if (j <= 4) ewe_out[i] * ewe_in[j] else ewe_in[j] * ewe_in[i]
  ewe_asc <- ewe_asc + f * log(f * 1280 / den)
}
toy_ac_ewe <- ewe_asc / -sum(ewe_A[ewe_A > 0] * log(ewe_A[ewe_A > 0] / 1280))

toy_expected <- c(
  TE_I = 200 / 1060, TE_II = 0.125, TE_I_enaR = 200 / 1072, TE_II_enaR = 0.1, FCI_enaR = toy_fci,
  FCI = toy_fci_ewe, A_C = toy_ac_ewe, CI = 3 / 9, TST = 1280, DH = 60 / 140,
  det_EE = 60 / 932, EE_P = 0.14, EE_H = 0.5, EE_C = 0.4, TL_C = 3, mTLco = 3, MTI = toy_M, KS = toy_ks
)
toy_got <- c(
  TE_I = toy_net$TE_I, TE_II = toy_net$TE_II, TE_I_enaR = toy_net$TE_I_enaR, TE_II_enaR = toy_net$TE_II_enaR,
  FCI_enaR = toy_net$FCI_enaR, FCI = toy_net$FCI, A_C = toy_net$A_C, CI = toy_eco$Connectance_Index, TST = toy_eco$Total_System_Throughput,
  DH = toy_eco$DH_ratio,
  det_EE = unname(toy_model$EE[toy_model$Group == "Detritus"]),
  EE_P = unname(toy_model$EE[1]), EE_H = unname(toy_model$EE[2]), EE_C = unname(toy_model$EE[3]), TL_C = unname(toy_model$TL[toy_model$Group == "Carnivore"]),
  mTLco = toy_eco$mTLco, MTI = unname(mti_fixed(toy_model, toy)), KS = calc_keystoneness(toy_model, toy, "toy")$KS
)
toy_tol <- c(
  TE_I = 5e-5, TE_II = 5e-5, TE_I_enaR = 5e-5, TE_II_enaR = 5e-5, FCI_enaR = 5e-5, FCI = 5e-5, A_C = 5e-5,
  CI = 5e-4, TST = 0.05, DH = 5e-4, det_EE = 1e-9, EE_P = 1e-9, EE_H = 1e-9, EE_C = 1e-9, TL_C = 1e-9,
  mTLco = 5e-4, setNames(rep(1e-9, 25), paste0("MTI", 1:25)), setNames(rep(5e-5, 3), paste0("KS", 1:3))
)
stopifnot(identical(names(toy_got), names(toy_expected)), identical(names(toy_tol), names(toy_expected)))
toy_fail <- abs(toy_got - toy_expected) > toy_tol
if (any(toy_fail)) {
  stop("toy food chain check failed:\n", paste0(names(toy_got)[toy_fail], ": got ", signif(toy_got[toy_fail], 6),
    ", expected ", signif(toy_expected[toy_fail], 6), collapse = "\n"))
}
message("toy food chain check passed")
for (nm in names(models_main)) {
  m <- models_main[[nm]]
  cons <- which(m$type == 0)
  q <- m$Biomass[cons] * m$QB[cons]
  num <- den <- 0
  for (k in seq_along(cons)) {
    dc <- m$DC[, m$Group[cons[k]]]
    dc[is.na(dc)] <- 0
    tl <- m$TL[match(names(dc), m$Group)]
    tl[is.na(tl)] <- 0
    oi <- sum(dc * (tl - sum(dc * tl))^2)
    num <- num + log(q[k] / min(q)) * oi
    den <- den + log(q[k] / min(q))
  }
  got <- calc_ecosystem_indices(m, params_main[[nm]], nm)$System_Omnivory_Index
  if (abs(got - num / den) > 5e-5) stop(nm, " SOI ", got, " != EwE form ", num / den)
}
message("SOI check passed")
for (nm in names(models_main)) {
  m <- models_main[[nm]]
  p <- params_main[[nm]]
  z <- function(x) ifelse(is.na(x), 0, x)
  liv <- which(m$type %in% 0:1)
  cons <- which(m$type == 0)
  q <- sum(m$Biomass[cons] * m$QB[cons])
  resp <- sum(m$Biomass[cons] * m$QB[cons] * (1 - m$Unassim[cons]) - m$Biomass[cons] * m$PB[cons])
  to_det <- sum((1 - m$EE[liv]) * m$Biomass[liv] * m$PB[liv] + m$Biomass[liv] * z(m$QB[liv]) * m$Unassim[liv])
  det_imp <- sum(z(p$model$DetInput[match(m$Group[m$type == 2], p$model$Group)]))
  hand <- q + sum(m$Landings, na.rm = TRUE) + to_det + det_imp + resp
  got <- calc_ecosystem_indices(m, p, nm)$Total_System_Throughput
  if (abs(got - hand) > 0.05) stop(nm, " TST ", got, " != EwE form ", hand)
}
message("TST check passed")

# Detritus EE -----
for (nm in names(models_main)) {
  m <- models_main[[nm]]
  p <- params_main[[nm]]
  z <- function(x) ifelse(is.na(x), 0, x)
  liv <- which(m$type %in% c(0, 1))
  rpath_raw <- rpath(p, eco.name = nm)$EE
  for (d in which(m$type == 2)) {
    dcol <- match(d, which(m$type == 2))
    inflow <- z(p$model$DetInput[p$model$Group == m$Group[d]])
    for (i in liv) {
      lost <- (1 - m$EE[i]) * m$Biomass[i] * m$PB[i] + m$Biomass[i] * z(m$QB[i]) * m$Unassim[i]
      inflow <- inflow + lost * z(m$DetFate[i, dcol])
    }
    eaten <- 0
    for (j in which(m$type == 0)) eaten <- eaten + m$Biomass[j] * m$QB[j] * z(m$DC[d, match(j, liv)])
    hand <- eaten / inflow
    if (abs(hand - m$EE[d]) > 1e-9) stop(nm, " ", m$Group[d], ": detritus EE used ", m$EE[d], " != hand budget ", hand)
    msg <- if (abs(rpath_raw[d] - hand) <= 1e-6) {
      "Rpath agrees"
    } else {
      "Rpath differs"
    }
    message(sprintf("detritus EE %s %s: %.4f; Rpath %.4f (%s)", nm, m$Group[d], hand, rpath_raw[d], msg))
  }
}

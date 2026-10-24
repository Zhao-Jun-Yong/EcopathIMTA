# Parameter uncertainty -----

N_MC <- as.integer(Sys.getenv("N_MC", 500))
EE_TOLERANCE <- 1e-9
MAX_ATTEMPTS <- 400000

quiet_rpath <- function(x) {
  utils::capture.output(result <- suppressMessages(x))
  result
}
with_params <- function(model, params) {
  if (!is.null(model)) attr(model, "params") <- params
  model
}

# Model builders -----
build_field_sw <- function(groups, diet_fun = NULL, name = "SW_unc_field") {
  p <- create.rpath.params(group = groups$func_group, type = groups$type)
  p$model[, Biomass := groups$B]
  p$model[, PB := groups$PB]
  p$model[, QB := groups$QB]
  p$model[, Fleets := catch_vec(groups)]
  apply_static_cols(p, SW_FATE_WATER, SW_FATE_SED, sw_feed_detinput)
  p$diet[Group == "Commercial_feed", Mud_crab := DC_MUD_FEED]
  apply_mudcrab_base_diet(p$diet)
  p$diet[Group == "Sediment_detritus", Mud_crab := DC_MUD_SED]
  apply_sw_invert_diet(p$diet)
  if (!is.null(diet_fun)) diet_fun(p$diet)
  with_params(run_rpath_model(p, name), p)
}
build_field_fw <- function(groups, diet_fun = NULL, name = "FW_unc_field") {
  p <- create.rpath.params(group = groups$func_group, type = groups$type)
  p$model[, Biomass := groups$B]
  p$model[, PB := groups$PB]
  p$model[, QB := groups$QB]
  p$model[, Fleets := catch_vec(groups)]
  apply_static_cols(p, FW_FATE_WATER, FW_FATE_SED, fw_feed_detinput)
  apply_fw_diet(p$diet, dc_feed_til = DC_TIL_FEED, dc_duck = 0, dc_sed_til = DC_TIL_SED_FIELD)
  if (!is.null(diet_fun)) diet_fun(p$diet)
  with_params(run_rpath_model(p, name), p)
}

# Indices -----
UNC_INDICES <- c("Total_System_Throughput", "PR_ratio", "DH_ratio", "Total_Biomass", "TE_I", "TE_II", "FCI", "A_C", "System_Omnivory_Index")
unc_indices <- function(model) {
  eco <- calc_ecosystem_indices(model, attr(model, "params"), "x", NULL)
  net <- calc_network_indices(model, attr(model, "params"), "x")
  vals <- c(
    unlist(eco[c("Total_System_Throughput", "PR_ratio", "DH_ratio", "Total_Biomass", "System_Omnivory_Index")]),
    unlist(net[c("TE_I", "TE_II", "FCI", "A_C")])
  )
  vals[UNC_INDICES]
}

# Setup -----
make_setup <- function(model, groups, field) {
  stopifnot(all(groups$type %in% 0:3))
  sampled <- groups$func_group[groups$type %in% 0:1]
  pedigree <- pedigree_ranges[pedigree_ranges$Model == model & pedigree_ranges$Group %in% sampled, ]
  stopifnot(all(sampled %in% pedigree$Group))
  base <- quiet_rpath(field(groups))
  base_diet <- attr(base, "params")$diet
  consumer_names <- setdiff(names(base_diet), "Group")
  diet_range <- setNames(pedigree$Diet_range, pedigree$Group)
  diet_template <- lapply(consumer_names, function(cn) {
    d <- base_diet[[cn]]
    rows <- which(!is.na(d) & d > 0)
    list(rows = rows, base = d[rows], width = unname(diet_range[cn]))
  })
  names(diet_template) <- consumer_names
  diet_template <- diet_template[!vapply(diet_template, function(t) is.na(t$width) || length(t$base) < 1, TRUE)]
  diet_matrix <- base$DC
  prey_names <- if (!is.null(rownames(diet_matrix))) rownames(diet_matrix) else base$Group[seq_len(nrow(diet_matrix))]
  is_eaten <- setNames(rowSums(diet_matrix > 0, na.rm = TRUE) > 0, prey_names)
  prescreen_groups <- groups$func_group[groups$type %in% 0:1 & !is.na(groups$Catch) & groups$Catch > 0 & !(is_eaten[groups$func_group] %in% TRUE)]
  consumers <- groups$func_group[groups$type %in% 0:1 & !is.na(groups$QB)]
  catch_groups <- groups$func_group[groups$type %in% 0:1 & !is.na(groups$Catch) & groups$Catch > 0]
  bounds <- do.call(rbind, lapply(seq_len(nrow(pedigree)), function(i) {
    grp <- pedigree$Group[i]
    out <- data.frame(group = grp, param = "B", lo = 1 - pedigree$B_range[i], hi = 1 + pedigree$B_range[i])
    out <- rbind(out, data.frame(group = grp, param = if (grp %in% catch_groups) "Y" else "PB", lo = 1 - pedigree$PB_range[i], hi = 1 + pedigree$PB_range[i]))
    if (grp %in% consumers) {
      out <- rbind(out, data.frame(group = grp, param = "GE", lo = 1 - pedigree$QB_range[i], hi = 1 + pedigree$QB_range[i]))
    }
    out
  }))
  bounds$name <- paste(bounds$group, bounds$param, sep = "|")
  stopifnot(!anyNA(bounds$lo) && !anyNA(bounds$hi))
  model_table <- attr(base, "params")$model
  ge_limit <- 1 - setNames(model_table$Unassim, model_table$Group)[consumers]
  list(
    model = model, groups = groups, field = field, pedigree = pedigree, base = base, diet_template = diet_template,
    prescreen_groups = prescreen_groups, consumers = consumers, catch_groups = catch_groups, bounds = bounds, ge_limit = ge_limit
  )
}
unc_setups <- list(
  sw = make_setup("sw_field", sw_field_groups, build_field_sw),
  fw = make_setup("fw_field", fw_field_groups, build_field_fw)
)

# Draws -----
draw_multipliers <- function(setup) setNames(runif(nrow(setup$bounds), setup$bounds$lo, setup$bounds$hi), setup$bounds$name)
draw_diet <- function(setup) {
  lapply(setup$diet_template, function(t) {
    x <- runif(length(t$base), t$base * (1 - t$width), t$base * (1 + t$width))
    x / sum(x) * sum(t$base)
  })
}
make_diet_fun <- function(setup, diet_draw) {
  function(diet_table) {
    for (cn in names(diet_draw)) data.table::set(diet_table, i = setup$diet_template[[cn]]$rows, j = cn, value = diet_draw[[cn]])
    invisible(diet_table)
  }
}
apply_multipliers <- function(setup, mult) {
  groups <- setup$groups
  for (i in which(groups$type %in% 0:1)) {
    grp <- groups$func_group[i]
    base_pb <- groups$PB[i]
    base_qb <- groups$QB[i]
    groups$B[i] <- groups$B[i] * mult[paste(grp, "B", sep = "|")]
    if (grp %in% setup$catch_groups) {
      groups$Catch[i] <- groups$Catch[i] * mult[paste(grp, "Y", sep = "|")]
      groups$PB[i] <- base_pb * mult[paste(grp, "Y", sep = "|")] / mult[paste(grp, "B", sep = "|")]
    } else {
      groups$PB[i] <- base_pb * mult[paste(grp, "PB", sep = "|")]
    }
    if (!is.na(base_qb)) groups$QB[i] <- groups$PB[i] / ((base_pb / base_qb) * mult[paste(grp, "GE", sep = "|")])
  }
  groups
}
prescreen_failures <- function(setup, groups) {
  i <- match(setup$prescreen_groups, groups$func_group)
  tol <- ifelse(setup$prescreen_groups %in% setup$catch_groups, EE_TOLERANCE, 0)
  setup$prescreen_groups[groups$Catch[i] / (groups$B[i] * groups$PB[i]) > 1 + tol]
}

energy_failures <- function(setup, groups) {
  i <- match(setup$consumers, groups$func_group)
  setup$consumers[groups$PB[i] / groups$QB[i] > setup$ge_limit]
}

with_worst <- function(res, group, ratio, rule) {
  k <- which.max(ratio)
  res$bind_group <- group[k]
  res$bind_rule <- rule
  res$worst_ratio <- unname(ratio[k])
  res
}

run_attempt <- function(setup, mult, diet_draw) {
  groups <- apply_multipliers(setup, mult)
  res <- list(status = NA_integer_, indices = rep(NA_real_, length(UNC_INDICES)),
    bind_group = NA_character_, bind_rule = NA_character_, worst_ratio = NA_real_)
  failing <- prescreen_failures(setup, groups)
  if (length(failing)) {
    res$status <- 0L
    i <- match(failing, groups$func_group)
    return(with_worst(res, failing, groups$Catch[i] / (groups$B[i] * groups$PB[i]), "EE>1 (pre-screen)"))
  }
  negative_r <- energy_failures(setup, groups)
  if (length(negative_r)) {
    res$status <- 3L
    i <- match(negative_r, groups$func_group)
    return(with_worst(res, negative_r, groups$PB[i] / groups$QB[i] / setup$ge_limit[negative_r], "R<0"))
  }
  field_model <- quiet_rpath(setup$field(groups, make_diet_fun(setup, diet_draw)))
  if (is.null(field_model)) {
    res$status <- 1L
    res$bind_group <- res$bind_rule <- "check_failed"
    return(res)
  }
  ee_tol <- ifelse(field_model$Group %in% setup$catch_groups, EE_TOLERANCE, 0)
  over_one <- field_model$Group[which(field_model$EE > 1 + ee_tol)]
  if (length(over_one)) {
    res$status <- 1L
    return(with_worst(res, over_one, field_model$EE[match(over_one, field_model$Group)], "EE>1"))
  }
  res$status <- 2L
  res$indices <- unname(unc_indices(field_model))
  res
}

# Monte Carlo -----
mc_run <- function(system_name, target = N_MC) {
  setup <- unc_setups[[system_name]]
  attempts <- list()
  n_balanced <- 0
  while (n_balanced < target && length(attempts) < MAX_ATTEMPTS) {
    mult <- draw_multipliers(setup)
    diet_draw <- draw_diet(setup)
    attempt <- run_attempt(setup, mult, diet_draw)
    attempts[[length(attempts) + 1]] <- attempt
    if (attempt$status == 2L) n_balanced <- n_balanced + 1
  }
  if (n_balanced < target) warning(system_name, ": ", n_balanced, " of ", target, " balanced draws", call. = FALSE)
  list(system = system_name,
    status = vapply(attempts, `[[`, integer(1), "status"),
    bind_group = vapply(attempts, `[[`, character(1), "bind_group"),
    bind_rule = vapply(attempts, `[[`, character(1), "bind_rule"),
    worst_ratio = vapply(attempts, `[[`, numeric(1), "worst_ratio"),
    indices = do.call(rbind, lapply(attempts, `[[`, "indices")))
}
set.seed(20261002)
mc_results <- lapply(setNames(names(unc_setups), names(unc_setups)), mc_run)

quant3 <- function(x) quantile(x, c(0.025, 0.5, 0.975), na.rm = TRUE)
base_indices <- lapply(unc_setups, function(setup) setNames(unc_indices(setup$base), UNC_INDICES))

mc_accept <- do.call(rbind, lapply(mc_results, function(run) {
  n_attempts <- length(run$status)
  n_balanced <- sum(run$status == 2L)
  ci <- binom.test(n_balanced, n_attempts)$conf.int
  data.frame(system = run$system, attempts = n_attempts, prescreen_rejected = sum(run$status == 0L), full_model_unbalanced = sum(run$status == 1L),
    energy_rejected = sum(run$status == 3L), balanced = n_balanced, accept = n_balanced / n_attempts, accept_lo = ci[1], accept_hi = ci[2])
}))
print(mc_accept, row.names = FALSE, digits = 3)

mc_index <- do.call(rbind, lapply(mc_results, function(run) {
  kept <- run$status == 2L
  do.call(rbind, lapply(seq_along(UNC_INDICES), function(j) {
    values <- run$indices[kept, j]
    data.frame(system = run$system, index = UNC_INDICES[j], n_balanced = length(values),
      baseline = base_indices[[run$system]][j],
      lo = unname(quant3(values)[1]), med = unname(quant3(values)[2]), hi = unname(quant3(values)[3]),
      baseline_percentile = round(100 * mean(values <= base_indices[[run$system]][j]), 1))
  }))
}))
print(mc_index, row.names = FALSE, digits = 3)

mc_binding <- do.call(rbind, lapply(mc_results, function(run) {
  rejected <- run$status != 2L
  d <- data.frame(group = run$bind_group[rejected], rule = run$bind_rule[rejected], ratio = run$worst_ratio[rejected])
  out <- do.call(rbind, lapply(split(d, list(d$group, d$rule), drop = TRUE), function(g) {
    data.frame(group = g$group[1], rule = g$rule[1], attempts_binding = nrow(g), median_ratio = median(g$ratio))
  }))
  data.frame(system = run$system, out, share_of_attempts = out$attempts_binding / length(run$status), row.names = NULL)
}))
print(mc_binding[order(mc_binding$system, -mc_binding$attempts_binding), ], row.names = FALSE, digits = 3)

writexl::write_xlsx(
  list(
    MC_accept = mc_accept, MC_field_indices = mc_index, MC_binding_constraint = mc_binding,
    Ranges_SW = unc_setups$sw$bounds, Ranges_FW = unc_setups$fw$bounds,
    Run_info = data.frame(
      item = c("software", "target balanced draws per system", "max attempts", "seed"),
      value = c(paste0("R ", getRversion(), "; Rpath ", as.character(utils::packageVersion("Rpath")), "; enaR ", as.character(utils::packageVersion("enaR"))),
        N_MC, MAX_ATTEMPTS, 20261002)
    )
  ),
  here::here("output/TableS3_uncertainty.xlsx")
)

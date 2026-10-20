############################################################
### Export Table 1 values: five fixed-rate measles waves   ###
############################################################

### Run from the project root containing fitData/.
### Output is a wide CSV with one row per table quantity and
### separate estimate/support-min/support-max cells for LaTeX formatting.

if (!requireNamespace("deSolve", quietly = TRUE)) {
  stop("Package 'deSolve' is required to calculate FWHM.")
}

base_path <- getwd()
fitData_dir <- file.path(base_path, "fitData")
out_dir <- file.path(fitData_dir, "K1_table_fixedRates_fiveWaves")
table_dir <- file.path(out_dir, "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

nu_fixed <- 0.100
gamma_fixed <- 0.125
rate_tol <- 1e-8
delta_BIC_threshold <- 2
fwhm_horizon_days <- 730
fwhm_dt_days <- 0.10

fit_series <- list(
  Ontario_2025 = list(label = "Ontario 2025", files = c("K1_incidence_Ontario_fitR0S0csv.csv", "K1_incidence_Ontario_fitR0S0.csv")),
  Alberta_2025 = list(label = "Alberta 2025", files = "K1_incidence_Alberta_fitR0S0_2025.csv"),
  Alberta_2026 = list(label = "Alberta 2026", files = "K1_incidence_Alberta_fitR0S0_2026.csv"),
  Manitoba_2025 = list(label = "Manitoba 2025", files = "K1_incidence_Manitoba_fitR0S0_2025.csv"),
  Manitoba_2026 = list(label = "Manitoba 2026", files = "K1_incidence_Manitoba_fitR0S0_2026.csv")
)
series_order <- names(fit_series)

numify <- function(x) suppressWarnings(as.numeric(as.character(x)))

first_existing <- function(files, label) {
  candidates <- file.path(fitData_dir, files)
  found <- candidates[file.exists(candidates)]
  if (!length(found)) stop("Missing ", label, ". Checked: ", paste(candidates, collapse = "; "))
  found[1]
}

find_col <- function(df, candidates) {
  hit <- candidates[candidates %in% names(df)]
  if (!length(hit)) stop("Missing required column: ", paste(candidates, collapse = ", "))
  hit[1]
}

get_se <- function(row, parameter_name, fallback_index) {
  parameter_columns <- grep("^se_.*parameter[0-9]+$", names(row), value = TRUE)
  for (parameter_column in parameter_columns) {
    if (identical(as.character(row[[parameter_column]][1]), parameter_name)) {
      index <- sub("^.*parameter([0-9]+)$", "\\1", parameter_column)
      se_column <- grep(paste0("^se_.*se", index, "$"), names(row), value = TRUE)
      if (length(se_column)) return(numify(row[[se_column[1]]][1]))
    }
  }
  se_column <- grep(paste0("^se_.*se", fallback_index, "$"), names(row), value = TRUE)
  if (length(se_column)) numify(row[[se_column[1]]][1]) else NA_real_
}

initial_growth_rate <- function(Reff0) {
  0.5 * (-(nu_fixed + gamma_fixed) +
           sqrt((nu_fixed - gamma_fixed)^2 + 4 * nu_fixed * gamma_fixed * Reff0))
}

simulate_incidence <- function(Reff0, S0, I0, times) {
  beta <- Reff0 * gamma_fixed / S0
  rhs <- function(t, y, p) {
    with(as.list(c(y, p)), {
      incidence <- beta * S * I
      list(c(-incidence, incidence - nu * E1, nu * E1 - gamma * I, gamma * I))
    })
  }
  sol <- as.data.frame(deSolve::ode(
    y = c(S = S0, E1 = 0, I = I0, R = 0), times = times, func = rhs,
    parms = c(beta = beta, nu = nu_fixed, gamma = gamma_fixed), method = "lsoda"
  ))
  data.frame(time = sol$time, incidence = pmax(0, beta * sol$S * sol$I))
}

interpolate_crossing <- function(t1, y1, t2, y2, target) {
  t1 + (target - y1) * (t2 - t1) / (y2 - y1)
}

full_width_half_maximum <- function(Reff0, S0, I0) {
  times <- seq(0, fwhm_horizon_days, by = fwhm_dt_days)
  fit <- simulate_incidence(Reff0, S0, I0, times)
  peak_index <- which.max(fit$incidence)
  half_max <- fit$incidence[peak_index] / 2
  above <- which(fit$incidence >= half_max)
  first_above <- min(above)
  last_above <- max(above)
  if (first_above == 1 || last_above == nrow(fit)) {
    return(NA_real_)
  }
  left <- interpolate_crossing(
    fit$time[first_above - 1], fit$incidence[first_above - 1],
    fit$time[first_above], fit$incidence[first_above], half_max
  )
  right <- interpolate_crossing(
    fit$time[last_above], fit$incidence[last_above],
    fit$time[last_above + 1], fit$incidence[last_above + 1], half_max
  )
  right - left
}

summarise_series <- function(id, cfg) {
  file <- first_existing(cfg$files, paste0(id, " fit results"))
  dat <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  if ("success" %in% names(dat)) dat <- dat[tolower(as.character(dat$success)) %in% "true", , drop = FALSE]

  bic_column <- find_col(dat, c("ll_importanceSampling.BIC", "ll_linearization.BIC", "BIC"))
  required <- c("nu_fixed", "gamma_fixed", "I0_fixed", "est_ReproNum_pop", "est_S0_pop")
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Missing columns in ", file, ": ", paste(missing, collapse = ", "))

  dat$nu <- numify(dat$nu_fixed)
  dat$gamma <- numify(dat$gamma_fixed)
  dat$I0 <- numify(dat$I0_fixed)
  dat$BIC <- numify(dat[[bic_column]])
  dat$Reff0 <- numify(dat$est_ReproNum_pop)
  dat$S0 <- numify(dat$est_S0_pop)
  dat <- dat[
    abs(dat$nu - nu_fixed) <= rate_tol & abs(dat$gamma - gamma_fixed) <= rate_tol &
      is.finite(dat$I0) & is.finite(dat$BIC) & is.finite(dat$Reff0) & is.finite(dat$S0),
    , drop = FALSE
  ]
  if (!nrow(dat)) stop("No valid rows at nu=", nu_fixed, ", gamma=", gamma_fixed, " in ", file)

  dat <- dat[order(dat$I0, dat$BIC), , drop = FALSE]
  profile <- do.call(rbind, lapply(split(dat, dat$I0), function(x) x[which.min(x$BIC), , drop = FALSE]))
  profile <- profile[order(profile$I0), , drop = FALSE]
  profile$delta_BIC <- profile$BIC - min(profile$BIC)
  profile$delta_BIC_lt_2 <- profile$delta_BIC < delta_BIC_threshold
  profile$Reff0_SE <- vapply(seq_len(nrow(profile)), function(i) get_se(profile[i, , drop = FALSE], "ReproNum_pop", 1), numeric(1))
  profile$S0_SE <- vapply(seq_len(nrow(profile)), function(i) get_se(profile[i, , drop = FALSE], "S0_pop", 2), numeric(1))
  profile$Reff0_RSE_percent <- 100 * profile$Reff0_SE / profile$Reff0
  profile$S0_RSE_percent <- 100 * profile$S0_SE / profile$S0

  best <- profile[which.min(profile$BIC), , drop = FALSE]
  supported <- profile[profile$delta_BIC_lt_2, , drop = FALSE]
  if (!nrow(supported)) stop("No Delta BIC < 2 support rows in ", file)

  growth_rate <- initial_growth_rate(best$Reff0[1])
  doubling_time <- if (growth_rate > 0) log(2) / growth_rate else NA_real_
  fwhm_days <- full_width_half_maximum(best$Reff0[1], best$S0[1], best$I0[1])

  list(
    profile = profile,
    values = list(
      label = cfg$label, source_file = basename(file),
      Reff0 = best$Reff0[1], Reff0_support_min = min(supported$Reff0), Reff0_support_max = max(supported$Reff0),
      Reff0_RSE_percent = best$Reff0_RSE_percent[1],
      Reff0_RSE_percent_support_min = min(supported$Reff0_RSE_percent), Reff0_RSE_percent_support_max = max(supported$Reff0_RSE_percent),
      S0 = best$S0[1], S0_support_min = min(supported$S0), S0_support_max = max(supported$S0),
      S0_RSE_percent = best$S0_RSE_percent[1],
      S0_RSE_percent_support_min = min(supported$S0_RSE_percent), S0_RSE_percent_support_max = max(supported$S0_RSE_percent),
      I0 = best$I0[1], I0_support_min = min(supported$I0), I0_support_max = max(supported$I0),
      minimum_BIC = best$BIC[1],
      initial_growth_rate_per_day = growth_rate,
      doubling_time_days = doubling_time,
      FWHM_days = fwhm_days
    )
  )
}

make_wide_table <- function(results) {
  quantities <- data.frame(
    quantity_key = c("Reff0", "Reff0_RSE_percent", "S0", "S0_RSE_percent", "I0", "initial_growth_rate_per_day", "doubling_time_days", "FWHM_days", "minimum_BIC"),
    quantity = c("R_eff,0", "RSE of R_eff,0", "S0", "RSE of S0", "I0", "Initial exponential growth rate, r", "Doubling time", "Full width at half maximum", "Minimum BIC"),
    unit = c("", "percent", "individuals", "percent", "individuals", "day^-1", "days", "days", ""),
    interval_type = c("Delta BIC < 2 profile support range", "Delta BIC < 2 profile support range",
                      "Delta BIC < 2 profile support range", "Delta BIC < 2 profile support range",
                      "Delta BIC < 2 profile support range", "best-fit trajectory only", "best-fit trajectory only", "best-fit trajectory only", "not applicable"),
    stringsAsFactors = FALSE
  )

  output <- quantities
  for (id in series_order) {
    values <- results[[id]]$values
    prefix <- id
    output[[paste0(prefix, "_estimate")]] <- vapply(output$quantity_key, function(key) values[[key]], numeric(1))
    output[[paste0(prefix, "_support_min")]] <- vapply(output$quantity_key, function(key) {
      min_key <- paste0(key, "_support_min")
      if (min_key %in% names(values)) values[[min_key]] else NA_real_
    }, numeric(1))
    output[[paste0(prefix, "_support_max")]] <- vapply(output$quantity_key, function(key) {
      max_key <- paste0(key, "_support_max")
      if (max_key %in% names(values)) values[[max_key]] else NA_real_
    }, numeric(1))
  }
  output$quantity_key <- NULL
  output
}

results <- lapply(names(fit_series), function(id) summarise_series(id, fit_series[[id]]))
names(results) <- names(fit_series)
wide_table <- make_wide_table(results)

output_file <- file.path(table_dir, "K1_fixedRates_fiveWaves_latex_table_values.csv")
write.csv(wide_table, output_file, row.names = FALSE, na = "")
for (id in names(results)) {
  write.csv(results[[id]]$profile,
            file.path(table_dir, paste0("K1_", id, "_fixedRates_I0_profile_for_table.csv")),
            row.names = FALSE)
}

print(wide_table)
message("Table values written to: ", output_file)

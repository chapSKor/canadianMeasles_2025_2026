# Weekly supplementary table: 45 rows, five waves x nine (L,D) pairs.
# Period labels remain in DAYS; matching rates are nu=7/L, gamma=7/D per WEEK.
# ONE SCRIPT: exports supplementary AND main-table CSVs, including r, doubling time and FWHM.
# Supplementary parameter table for five waves x nine period combinations.
# Run from the project root containing fitData/. Requires deSolve for trajectories; no Monolix runs.
# Reads saved CSV scan RESULTS, not the R scripts used to run the scans.
# IMPORTANT: pool only runs of the SAME wave/data/model/error specification.
# Keep superseded or incompatible runs outside fitData/, or set files_by_wave.
# All automatically matched inputs are listed in the exported manifest.

if (!requireNamespace("deSolve", quietly = TRUE)) {
  stop("Install deSolve first: install.packages('deSolve')")
}
# Source from the project root containing fitData/.
base_path <- getwd()
fit_dir <- file.path(base_path, "fitData")
out_dir <- file.path(fit_dir, "K1_supplement_period_sensitivity_wkly")
bic_column <- "ll_importanceSampling.BIC" # One consistent likelihood method.
rate_tol <- 1e-8
support_cutoff <- 2
periods <- expand.grid(latent_days = 9:11, infectious_days = 7:9)
periods <- periods[order(periods$latent_days, periods$infectious_days), , drop = FALSE]
waves <- c("Ontario_2025", "Alberta_2025", "Alberta_2026",
           "Manitoba_2025", "Manitoba_2026")

# NULL: discover coarse/fine CSVs in fitData (nonrecursive).
# Or give explicit filenames relative to fitData to control which runs to pool.
files_by_wave <- setNames(vector("list", length(waves)), waves)
# Example:
# files_by_wave$Alberta_2026 <- c("K1_incidence_Alberta_fitR0S0_2026_wkly.csv",
#                                "K1_incidence_Alberta_fitR0S0_2026_wkly_FineTuneI0.csv")

# Exact weekly filename patterns: never mix daily results or grid definitions.
# Both _wkly_FineTuneI0 and _FineTuneI0_wkly are recognized, case-insensitively.
weekly_suffix <- "(_wkly|_wkly_FineTuneI0|_FineTuneI0_wkly)\\.csv$"
patterns <- setNames(vapply(waves, function(wave) {
  province <- sub("_[0-9]+$", "", wave)
  year <- sub(".*_", "", wave)
  stem <- if (province == "Ontario") "K1_incidence_Ontario_fitR0S0(csv|_2025)?"
          else paste0("K1_incidence_", province, "_fitR0S0_", year)
  paste0("^", stem, weekly_suffix)
}, ""), waves)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))

# Resolve SE/RSE by the parameter NAME, never by an assumed column index.
# Monolix exported RSE is already in percent. Prefer stochastic approximation.
rse_percent <- function(row, parameter, estimate) {
  pc <- grep("^se_.*parameter[0-9]+$", names(row), value = TRUE)
  pc <- pc[order(!grepl("stochasticApproximation", pc), pc)]
  for (name in pc) {
    if (!identical(as.character(row[[name]][1]), parameter)) next
    rse_name <- sub("parameter([0-9]+)$", "rse\\1", name)
    se_name <- sub("parameter([0-9]+)$", "se\\1", name)
    if (rse_name %in% names(row)) {
      value <- num(row[[rse_name]][1])
      if (is.finite(value) && value >= 0) return(value)
    }
    if (se_name %in% names(row)) {
      value <- num(row[[se_name]][1])
      if (is.finite(value) && value >= 0) return(100 * value / estimate)
    }
  }
  NA_real_
}

# All ODE calculations use WEEKS. Reported durations are converted to days.
# Incidence remains beta*S*I, matching the existing fitted observation model.
# FWHM is measured around the global peak of each independent fitted wave;
# display-only joins used in Figure 1A are not used here.
fwhm_dt_days <- 0.1
fwhm_initial_horizon_days <- 730
fwhm_max_horizon_days <- 5840
trajectory_metrics <- function(Reff0, S0, I0, nu, gamma) {
  r_week <- (sqrt((nu - gamma)^2 + 4 * nu * gamma * Reff0) - nu - gamma) / 2
  beta <- Reff0 * gamma / S0
  rhs <- function(t, y, p) {
    inc <- beta * y[1] * y[3]
    list(c(-inc, inc - nu * y[2], nu * y[2] - gamma * y[3], gamma * y[3]))
  }
  horizon <- fwhm_initial_horizon_days
  repeat {
    times <- seq(0, horizon / 7, by = fwhm_dt_days / 7)
    fit <- as.data.frame(deSolve::ode(
      y = c(S = S0, E1 = 0, I = I0, R = 0), times = times,
      func = rhs, parms = NULL, method = "lsoda", rtol = 1e-9, atol = 1e-9))
    inc <- pmax(0, beta * fit$S * fit$I)
    if (length(inc) != length(times) || any(!is.finite(inc))) stop("Invalid ODE trajectory")
    peak <- which.max(inc)
    half <- inc[peak] / 2
    # Find the nearest crossings bracketing this peak, not separate early transients.
    left_candidates <- which(seq_along(inc) < peak & inc < half)
    right_candidates <- which(seq_along(inc) > peak & inc < half)
    if (length(right_candidates) || horizon >= fwhm_max_horizon_days) break
    horizon <- min(2 * horizon, fwhm_max_horizon_days)
  }
  cross <- function(i, j) times[i] + (half - inc[i]) *
    (times[j] - times[i]) / (inc[j] - inc[i])
  left <- if (length(left_candidates)) {
    i <- max(left_candidates); cross(i, i + 1L)
  } else NA_real_
  right <- if (length(right_candidates)) {
    j <- min(right_candidates); cross(j - 1L, j)
  } else NA_real_
  status <- if (half <= 0) "no_positive_peak" else if (peak == 1L) {
    "peak_at_wave_origin"
  } else if (!is.finite(right)) "right_crossing_not_reached" else if (!is.finite(left)) {
    "left_crossing_before_wave_origin"
  } else "ok"
  width <- if (status == "ok") right - left else NA_real_
  data.frame(initial_growth_rate_per_week = r_week,
    initial_growth_rate_per_day = r_week / 7,
    doubling_time_weeks = if (r_week > 0) log(2) / r_week else NA_real_,
    doubling_time_days = if (r_week > 0) 7 * log(2) / r_week else NA_real_,
    FWHM_weeks = width, FWHM_days = 7 * width, FWHM_status = status,
    peak_time_days = 7 * times[peak], peak_incidence_per_week = inc[peak],
    half_peak_left_days = 7 * left, half_peak_right_days = 7 * right,
    simulated_horizon_days = horizon)
}

manifest <- list()
read_results <- function(file, wave) {
  dat <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  if (!bic_column %in% names(dat)) {
    # Grid definitions have no estimates/likelihoods and are not fit results.
    likelihood_columns <- grep("BIC|est_ReproNum", names(dat), value = TRUE)
    if (length(likelihood_columns)) stop("Incompatible likelihood columns in ", file,
                                         "; expected ", bic_column)
    manifest[[length(manifest) + 1L]] <<- data.frame(wave = wave, file = file,
      status = "skipped: no likelihood results", input_rows = nrow(dat), valid_rows = 0)
    return(NULL)
  }
  required <- c("success", "nu_fixed", "gamma_fixed", "I0_fixed",
                "est_ReproNum_pop", "est_S0_pop")
  if (!all(required %in% names(dat))) stop("Missing result columns in ", file)
  expected_province <- sub("_[0-9]+$", "", wave)
  if ("province" %in% names(dat) && any(!is.na(dat$province) &
      as.character(dat$province) != expected_province)) stop("Province mismatch in ", file)
  if ("model" %in% names(dat) && any(!is.na(dat$model) & dat$model != "K1")) {
    stop("Non-K1 rows in ", file)
  }
  out <- data.frame(nu = num(dat$nu_fixed), gamma = num(dat$gamma_fixed),
    I0 = num(dat$I0_fixed), S0 = num(dat$est_S0_pop),
    Reff0 = num(dat$est_ReproNum_pop), BIC = num(dat[[bic_column]]))
  keep <- tolower(trimws(as.character(dat$success))) %in% c("true", "1") &
    apply(out, 1, function(z) all(is.finite(z))) &
    out$nu > 0 & out$gamma > 0 & out$I0 > 0 & out$S0 > 0 & out$Reff0 > 0
  out$source_row <- seq_len(nrow(out))
  out <- out[keep, , drop = FALSE]
  out$source_file <- rep(file, nrow(out))
  out$S0_RSE_percent <- vapply(out$source_row, function(i)
    rse_percent(dat[i, , drop = FALSE], "S0_pop", num(dat$est_S0_pop[i])), numeric(1))
  out$Reff0_RSE_percent <- vapply(out$source_row, function(i)
    rse_percent(dat[i, , drop = FALSE], "ReproNum_pop", num(dat$est_ReproNum_pop[i])), numeric(1))
  manifest[[length(manifest) + 1L]] <<- data.frame(wave = wave, file = file,
    status = "included", input_rows = nrow(dat), valid_rows = nrow(out))
  out
}

tables <- profiles <- list()
for (wave in waves) {
  files <- files_by_wave[[wave]]
  if (is.null(files)) {
    files <- list.files(fit_dir, pattern = patterns[[wave]], full.names = TRUE, ignore.case = TRUE)
    # Do not accidentally assign explicitly dated later Ontario runs to 2025.
    if (wave == "Ontario_2025") files <- files[!grepl("_20(2[6-9]|[3-9][0-9])", basename(files))]
  } else files <- file.path(fit_dir, files)
  if (!length(files) || any(!file.exists(files))) stop("Missing scan files for ", wave)
  if (any(!grepl(patterns[[wave]], basename(files), ignore.case = TRUE))) {
    stop("Input filenames must identify weekly results for ", wave,
         "; use the recognized _wkly.csv or weekly FineTuneI0 patterns, not _grid.csv.")
  }
  files <- sort(unique(normalizePath(files, mustWork = TRUE)))
  message(wave, ": pooling ", paste(basename(files), collapse = ", "))
  scans <- do.call(rbind, lapply(files, read_results, wave = wave))
  if (is.null(scans) || !nrow(scans)) stop("No valid scan results for ", wave)
  for (j in seq_len(nrow(periods))) {
    L <- periods$latent_days[j]; D <- periods$infectious_days[j]
    slice <- scans[abs(scans$nu - 7/L) < rate_tol & abs(scans$gamma - 7/D) < rate_tol, , drop = FALSE]
    if (!nrow(slice)) stop("Missing weekly-rate slice: ", wave, ", L=", L, ", D=", D,
                          " days; expected nu=", 7/L, " and gamma=", 7/D, " week^-1")
    # Retain the minimum BIC per I0 across repeated coarse/fine fits.
    # Exact ties prefer complete RSEs, then file path and source row.
    missing_rse <- !is.finite(slice$S0_RSE_percent) | !is.finite(slice$Reff0_RSE_percent)
    slice <- slice[order(slice$I0, slice$BIC, missing_rse, slice$source_file, slice$source_row), , drop = FALSE]
    profile <- slice[!duplicated(slice$I0), , drop = FALSE]
    profile$delta_BIC <- profile$BIC - min(profile$BIC)
    profile$accepted <- profile$delta_BIC <= support_cutoff
    best <- profile[which.min(profile$BIC), , drop = FALSE]
    accepted <- profile$I0[profile$accepted]
    lower_edge <- min(accepted) == min(profile$I0)
    upper_edge <- max(accepted) == max(profile$I0)
    # A range is an envelope over evaluated values, not an interpolated interval.
    runs <- rle(profile$accepted)
    disjoint <- sum(runs$values) > 1
    tied <- profile$I0[profile$BIC == min(profile$BIC)]
    tables[[length(tables) + 1L]] <- data.frame(
      wave = wave, province = sub("_[0-9]+$", "", wave), year = sub(".*_", "", wave),
      latent_days = L, infectious_days = D, nu = 7/L, gamma = 7/D, rate_units = "week^-1",
      reference = L == 10 && D == 8,
      S0 = best$S0, S0_RSE_percent = best$S0_RSE_percent,
      Reff0 = best$Reff0, Reff0_RSE_percent = best$Reff0_RSE_percent,
      I0_minBIC = best$I0, I0_support_min = min(accepted), I0_support_max = max(accepted),
      minimum_BIC = best$BIC, I0_accepted_values = paste(accepted, collapse = ";"),
      I0_tied_minima = paste(tied, collapse = ";"),
      support_hits_lower_scan_boundary = lower_edge,
      support_hits_upper_scan_boundary = upper_edge, support_disjoint_on_evaluated_grid = disjoint,
      examined_I0_min = min(profile$I0), examined_I0_max = max(profile$I0),
      n_distinct_I0 = nrow(profile), max_I0_grid_gap = if (nrow(profile) > 1) max(diff(profile$I0)) else NA_real_,
      n_result_files_in_slice = length(unique(slice$source_file)),
      best_source_file = best$source_file, best_source_row = best$source_row)
    profiles[[length(profiles) + 1L]] <- cbind(data.frame(wave = wave, latent_days = L, infectious_days = D), profile)
    if (lower_edge || upper_edge) warning(wave, " L=", L, " D=", D, ": support reaches a scanned boundary.", call. = FALSE)
    if (disjoint) warning(wave, " L=", L, " D=", D, ": support has separated groups; inspect accepted values.", call. = FALSE)
    if (!is.finite(best$S0_RSE_percent) || !is.finite(best$Reff0_RSE_percent)) {
      warning(wave, " L=", L, " D=", D, ": missing RSE on the minimum-BIC row; exported as NA.", call. = FALSE)
    }
  }
}

result <- do.call(rbind, tables)
message("Computing growth, doubling time and FWHM for each selected fit...")
metrics <- lapply(seq_len(nrow(result)), function(i) {
  x <- result[i, ]
  m <- trajectory_metrics(x$Reff0, x$S0, x$I0_minBIC, x$nu, x$gamma)
  if (m$FWHM_status != "ok") warning(x$wave, " L=", x$latent_days,
    " D=", x$infectious_days, ": FWHM exported as NA (", m$FWHM_status, ")", call. = FALSE)
  m
})
result <- cbind(result, do.call(rbind, metrics))
# Compare only with the weekly supplementary-figure selections in this directory.
# Pooling fine runs can change best fits: the table may require a figure update.
figure_csv <- file.path(out_dir, "FigureS_K1_period_sensitivity_fiveWaves_selected_fits.csv")
if (file.exists(figure_csv)) {
  old <- read.csv(figure_csv, check.names = FALSE)
  needed <- c("wave", "latent_days", "infectious_days", "I0", "S0", "Reff0", "BIC")
  if (all(needed %in% names(old))) {
    key <- function(x) paste(x$wave, x$latent_days, x$infectious_days, sep = ":")
    idx <- match(key(result), key(old))
    result$differs_from_saved_figure <- NA
    hit <- which(!is.na(idx))
    result$differs_from_saved_figure[hit] <- vapply(hit, function(i) {
      a <- unlist(result[i, c("I0_minBIC", "S0", "Reff0", "minimum_BIC")], use.names = FALSE)
      b <- unlist(old[idx[i], c("I0", "S0", "Reff0", "BIC")], use.names = FALSE)
      any(!is.finite(b)) || any(abs(a - b) > 1e-7 * pmax(1, abs(b)))
    }, logical(1))
    if (any(result$differs_from_saved_figure %in% TRUE)) warning(
      "Pooled selections differ from the saved supplementary figure. Regenerate its curves from the pooled best fits.", call. = FALSE)
  }
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
# Publication-facing table: exactly one row per wave and period combination.
# The range is the envelope of evaluated I0 values with BIC - min(BIC) <= 2
# within that combination; audit files retain the exact accepted values.
publication_table <- data.frame(
  province = result$province,
  year = result$year,
  latent_period_days = result$latent_days,
  infectious_period_days = result$infectious_days,
  minimum_BIC = result$minimum_BIC,
  I0_at_minimum_BIC = result$I0_minBIC,
  accepted_I0_range = paste(result$I0_support_min, result$I0_support_max, sep = "--"),
  S0_at_minimum_BIC = result$S0,
  S0_RSE_percent = result$S0_RSE_percent,
  Reff0_at_minimum_BIC = result$Reff0,
  Reff0_RSE_percent = result$Reff0_RSE_percent,
  initial_growth_rate_per_day = result$initial_growth_rate_per_day,
  doubling_time_days = result$doubling_time_days,
  FWHM_days = result$FWHM_days,
  FWHM_status = result$FWHM_status
)
stopifnot(nrow(publication_table) == 45L)
write.csv(publication_table, file.path(out_dir, "TableS_K1_period_sensitivity.csv"), row.names = FALSE, na = "NA")
write.csv(result, file.path(out_dir, "TableS_K1_period_sensitivity_audit.csv"), row.names = FALSE, na = "NA")
write.csv(do.call(rbind, profiles), file.path(out_dir, "TableS_K1_period_sensitivity_profiles.csv"), row.names = FALSE, na = "NA")
write.csv(do.call(rbind, manifest), file.path(out_dir, "TableS_K1_period_sensitivity_inputs.csv"), row.names = FALSE)
message("Saved ", nrow(result), " table rows plus profile and input audits to ", out_dir)

# MAIN TABLE: reference periods L=10 and D=8 days, selected from the same pooled fits.
main_dir <- file.path(fit_dir, "K1_table_fixedRates_fiveWaves_wkly", "tables")
dir.create(main_dir, recursive = TRUE, showWarnings = FALSE)
main <- result[result$reference, , drop = FALSE]
main <- main[match(waves, main$wave), , drop = FALSE]
stopifnot(nrow(main) == 5L, !anyNA(main$wave))
main_publication <- publication_table[result$reference, , drop = FALSE]
write.csv(main_publication, file.path(main_dir, "K1_fixedRates_fiveWaves_summary.csv"),
          row.names = FALSE, na = "NA")
write.csv(main, file.path(main_dir, "K1_fixedRates_fiveWaves_audit.csv"), row.names = FALSE, na = "NA")

keys <- c("Reff0", "Reff0_RSE_percent", "S0", "S0_RSE_percent", "I0_minBIC",
          "initial_growth_rate_per_day", "doubling_time_days", "FWHM_days", "minimum_BIC")
wide <- data.frame(
  quantity = c("R_eff,0", "RSE of R_eff,0", "S0", "RSE of S0", "I0",
    "Initial exponential growth rate, r", "Doubling time", "Full width at half maximum", "Minimum BIC"),
  unit = c("dimensionless", "percent", "individuals", "percent", "individuals",
           "day^-1", "days", "days", "dimensionless"),
  interval_type = c(rep("Delta BIC <= 2 profile support range", 5),
    "selected-fit early exponential rate", "selected-fit early doubling time",
    "selected-fit trajectory only", "not applicable"))
all_profiles <- do.call(rbind, profiles)
finite_range <- function(x) {
  x <- x[is.finite(x)]
  if (length(x)) range(x) else c(NA_real_, NA_real_)
}
for (id in waves) {
  best <- main[main$wave == id, , drop = FALSE]
  prof <- all_profiles[all_profiles$wave == id & all_profiles$latent_days == 10 &
                        all_profiles$infectious_days == 8, , drop = FALSE]
  support <- prof[prof$accepted, , drop = FALSE]
  wide[[paste0(id, "_estimate")]] <- vapply(keys, function(k) best[[k]], numeric(1))
  ranges <- vapply(seq_along(keys), function(j) {
    if (j > 5) return(c(NA_real_, NA_real_))
    k <- if (keys[j] == "I0_minBIC") "I0" else keys[j]
    finite_range(support[[k]])
  }, numeric(2))
  wide[[paste0(id, "_support_min")]] <- ranges[1, ]
  wide[[paste0(id, "_support_max")]] <- ranges[2, ]
  write.csv(prof, file.path(main_dir, paste0("K1_", id, "_fixedRates_I0_profile_for_table.csv")),
            row.names = FALSE, na = "NA")
}
write.csv(wide, file.path(main_dir, "K1_fixedRates_fiveWaves_latex_table_values.csv"),
          row.names = FALSE, na = "NA")
message("Main table (5 waves, including r, doubling time and FWHM): ", main_dir)
message("Supplementary table (45 rows, including all three quantities): ", out_dir)
message("FWHM_status explains any missing width; no interval is invented when a crossing is absent.")

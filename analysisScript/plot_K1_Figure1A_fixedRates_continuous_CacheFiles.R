# All model times/rates use weeks. Output beta*S*I is the fitted instantaneous
# rate (cases/week); this update does not change it to integrated weekly counts.
############################################################
### Figure 1A: fixed-rate, province-continuous K1 fits    ###
############################################################

if (!requireNamespace("deSolve", quietly = TRUE)) {
  stop("Package 'deSolve' is required.")
}

base_path <- getwd()
fitData_dir <- file.path(base_path, "fitData")
out_root <- file.path(fitData_dir, "K1_figure1A_continuous_wkly")
plot_dir <- file.path(out_root, "plots")
table_dir <- file.path(out_root, "tables")
cache_dir <- file.path(out_root, "cache")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

nu_fixed <- 7 / 10  # week^-1
gamma_fixed <- 7 / 8  # week^-1
rate_tol <- 1e-8
delta_BIC_threshold <- 2
n_uncertainty_sims <- 500
uncertainty_seed <- 20260917
dense_dt_weeks <- 0.10 / 7
uncertainty_dt_weeks <- 1.00 / 7
recompute_uncertainty <- FALSE

province_colours <- list(
  Ontario = list(line = "#2F6FAE", fill = "#BFD7F1", pch = 16, lwd = 1.3),
  Alberta = list(line = "#D98C2B", fill = "#F2C38B", pch = 0, lwd = 2.1),
  Manitoba = list(line = "#3F8F5A", fill = "#B7DDB8", pch = 17, lwd = 1.5)
)

### The grid-definition CSVs are intentionally not used because they do not
### contain likelihoods or parameter estimates.
fits <- list(
  Ontario_2025 = list(
    province = "Ontario", wave = "2025",
    fit_file = c("K1_incidence_Ontario_fitR0S0csv_wkly.csv",
                 "K1_incidence_Ontario_fitR0S0_wkly.csv"),
    raw_file = "Ontario/MeaslesOntData_Province_2025.csv"
  ),
  Alberta_2025 = list(
    province = "Alberta", wave = "2025",
    fit_file = "K1_incidence_Alberta_fitR0S0_2025_wkly.csv",
    raw_file = "Alberta/MeaslesAlbData_Province_2025.csv"
  ),
  Alberta_2026 = list(
    province = "Alberta", wave = "2026",
    fit_file = "K1_incidence_Alberta_fitR0S0_2026_wkly.csv",
    raw_file = "Alberta/MeaslesAlbData_Province_2026_revised.csv"
  ),
  Manitoba_2025 = list(
    province = "Manitoba", wave = "2025",
    fit_file = "K1_incidence_Manitoba_fitR0S0_2025_wkly.csv",
    raw_file = "Manitoba/MeaslesManData_Province_2025.csv"
  ),
  Manitoba_2026 = list(
    province = "Manitoba", wave = "2026",
    fit_file = "K1_incidence_Manitoba_fitR0S0_2026_wkly.csv",
    raw_file = "Manitoba/MeaslesManData_Province_2026.csv"
  )
)


# Weekly input utilities. Pool coarse + optional fine RESULTS for the same wave.
# Grid definitions and old daily fits are never read.
weekly_scan <- function(root, province, year) {
  stems <- paste0("K1_incidence_", province, "_fitR0S0_", year)
  if (province == "Ontario" && as.character(year) == "2025")
    stems <- c(stems, "K1_incidence_Ontario_fitR0S0", "K1_incidence_Ontario_fitR0S0csv")
  wanted <- unlist(lapply(stems, function(s) c(paste0(s, "_wkly.csv"),
    paste0(s, "_wkly_FineTuneI0.csv"), paste0(s, "_FineTuneI0_wkly.csv"))))
  available <- list.files(root, full.names = TRUE)
  paths <- available[tolower(basename(available)) %in% tolower(wanted)]
  if (!length(paths)) stop("No weekly results for ", province, " ", year, " in ", root)
  frames <- lapply(paths, function(p) {
    d <- read.csv(p, check.names = FALSE, stringsAsFactors = FALSE)
    required <- c("nu_fixed", "gamma_fixed", "I0_fixed", "est_ReproNum_pop", "est_S0_pop",
                  "ll_importanceSampling.BIC")
    if (!all(required %in% names(d))) stop("Missing result columns in ", p, ": ",
      paste(setdiff(required, names(d)), collapse = ", "))
    d$weekly_source_file <- basename(p)
    d$weekly_source_row <- seq_len(nrow(d))
    d
  })
  cols <- unique(unlist(lapply(frames, names)))
  frames <- lapply(frames, function(d) {
    for (n in setdiff(cols, names(d))) d[[n]] <- NA
    d[, cols, drop = FALSE]
  })
  message(province, " ", year, ": ", paste(basename(paths), collapse = ", "))
  do.call(rbind, frames)
}

# Raw file paths are unchanged. Select weeks/days explicitly if needed;
# auto uses calendar-date spacing (7 days per weekly time unit).
raw_time_units <- "auto"  # "auto", "weeks", or "days"
weekly_time <- function(t, dates = NULL) {
  t <- suppressWarnings(as.numeric(as.character(t)))
  unit <- raw_time_units
  if (unit == "auto") {
    unit <- "weeks"
    if (!is.null(dates)) {
      ok <- is.finite(t) & !is.na(dates)
      tt <- t[ok]; dd <- as.numeric(dates[ok]); o <- order(tt)
      dt <- diff(tt[o]); ds <- diff(dd[o]); valid <- dt > 0 & ds > 0
      if (any(valid)) {
        ratio <- median(ds[valid] / dt[valid])
        if (abs(ratio - 7) < 0.15) unit <- "weeks"
        else if (abs(ratio - 1) < 0.03) unit <- "days"
        else stop("Raw time/date spacing is neither days nor weeks; check the input time origin/units.")
      } else warning("Insufficient dates to infer time units; assuming weeks.")
    } else warning("No calendar dates available; assuming raw times are weeks.")
  }
  if (!unit %in% c("weeks", "days")) stop("raw_time_units must be auto, weeks, or days")
  if (unit == "days") t <- t / 7
  t
}

numify <- function(x) suppressWarnings(as.numeric(as.character(x)))

first_existing <- function(paths, root, label) {
  paths <- file.path(root, paths)
  found <- paths[file.exists(paths)]
  if (length(found) == 0) stop("Missing ", label, ". Checked: ", paste(paths, collapse = "; "))
  found[1]
}

find_column <- function(x, choices) {
  hit <- choices[choices %in% names(x)]
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

read_raw <- function(cfg) {
  file <- first_existing(cfg$raw_file, base_path, paste0(cfg$province, " raw data"))
  x <- read.csv(file, check.names = FALSE, stringsAsFactors = FALSE)
  time_col <- find_column(x, c(paste0("Time (weeks) ", cfg$province), "Time (weeks)", "Time_weeks", paste0("Time (days) ", cfg$province), "Time (days)", "Time_days"))
  if (is.na(time_col)) time_col <- names(x)[grep("time|week|day", names(x), ignore.case = TRUE)[1]]
  y_col <- find_column(x, c(paste0("Y incidence ", cfg$province), "Y_True_Total", "Incidence"))
  if (is.na(y_col)) y_col <- names(x)[grep("incidence|case|y_true", names(x), ignore.case = TRUE)[1]]
  date_col <- find_column(x, c("week_start_date", paste0("week_start_date_", cfg$province),
                               "week_start_date_Ontario", "week_start_date_Alberta",
                               "week_start_date_Manitoba"))
  if (is.na(date_col)) date_col <- names(x)[grep("date", names(x), ignore.case = TRUE)[1]]
  date <- as.Date(as.character(x[[date_col]]), tryFormats = c("%m/%d/%Y", "%Y-%m-%d"))
  out <- data.frame(time = weekly_time(x[[time_col]], date), observed = numify(x[[y_col]]), date = date)
  out <- out[is.finite(out$time) & is.finite(out$observed), , drop = FALSE]
  out <- out[order(out$time), , drop = FALSE]
  if (sum(!is.na(out$date)) < 2) stop("Could not parse dates from ", file)
  out
}

get_se <- function(row, parameter_name, fallback_index) {
  pcols <- grep("^se_.*parameter[0-9]+$", names(row), value = TRUE)
  for (pc in pcols) {
    if (identical(as.character(row[[pc]][1]), parameter_name)) {
      index <- sub("^.*parameter([0-9]+)$", "\\1", pc)
      scol <- grep(paste0("^se_.*se", index, "$"), names(row), value = TRUE)
      if (length(scol)) return(numify(row[[scol[1]]][1]))
    }
  }
  scol <- grep(paste0("^se_.*se", fallback_index, "$"), names(row), value = TRUE)
  if (length(scol)) numify(row[[scol[1]]][1]) else NA_real_
}

read_fit <- function(id, cfg) {
  x <- weekly_scan(fitData_dir, cfg$province, cfg$wave)
  file <- paste(unique(x$weekly_source_file), collapse = ";")
  if ("success" %in% names(x)) x <- x[tolower(trimws(as.character(x$success))) %in% c("true", "1"), , drop = FALSE]
  bic_col <- find_column(x, c("ll_importanceSampling.BIC", "ll_linearization.BIC", "BIC"))
  if (is.na(bic_col)) stop("No BIC column in ", file)
  x$nu <- numify(x$nu_fixed)
  x$gamma <- numify(x$gamma_fixed)
  x$I0 <- numify(x$I0_fixed)
  x$BIC <- numify(x[[bic_col]])
  x$Reff0 <- numify(x$est_ReproNum_pop)
  x$S0 <- numify(x$est_S0_pop)
  x <- x[is.finite(x$nu) & is.finite(x$gamma) & abs(x$nu - nu_fixed) <= rate_tol & abs(x$gamma - gamma_fixed) <= rate_tol &
           is.finite(x$BIC) & is.finite(x$I0) & is.finite(x$Reff0) & is.finite(x$S0), , drop = FALSE]
  if (!nrow(x)) stop("No rows at nu=0.7 and gamma=0.875 week^-1 in ", file)
  x <- x[order(x$I0, x$BIC), , drop = FALSE]
  profile <- do.call(rbind, lapply(split(x, x$I0), function(z) z[which.min(z$BIC), , drop = FALSE]))
  profile <- profile[order(profile$I0), , drop = FALSE]
  profile$delta_BIC <- profile$BIC - min(profile$BIC)
  profile$delta_BIC_le_2 <- profile$delta_BIC <= delta_BIC_threshold
  best <- profile[which.min(profile$BIC), , drop = FALSE]
  supported <- profile[profile$delta_BIC_le_2, , drop = FALSE]
  summary <- data.frame(
    fit_id = id, province = cfg$province, wave = cfg$wave, source_file = best$weekly_source_file[1],
    nu_fixed = nu_fixed, gamma_fixed = gamma_fixed,
    best_I0 = best$I0[1],
    I0_deltaBIC_le_2_min = min(supported$I0),
    I0_deltaBIC_le_2_max = max(supported$I0),
    Reff0 = best$Reff0[1], Reff0_SE = get_se(best, "ReproNum_pop", 1),
    S0 = best$S0[1], S0_SE = get_se(best, "S0_pop", 2),
    best_BIC = best$BIC[1]
  )
  if (!is.finite(summary$Reff0_SE) || !is.finite(summary$S0_SE)) stop("Missing Reff0 or S0 SE in ", file)
  list(profile = profile, summary = summary)
}

simulate_k1 <- function(times, Reff0, S0, I0) {
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
  data.frame(time = sol$time, fitted_incidence = pmax(0, beta * sol$S * sol$I))
}

draw_positive <- function(n, mean, se) {
  x <- rnorm(n, mean, se)
  while (any(!is.finite(x) | x <= 0)) {
    bad <- which(!is.finite(x) | x <= 0)
    x[bad] <- rnorm(length(bad), mean, se)
  }
  x
}

make_wave <- function(id, cfg, fit, raw) {
  max_time <- max(raw$time)
  times <- seq(0, max_time, by = dense_dt_weeks)
  if (tail(times, 1) < max_time) times <- c(times, max_time)
  uncertainty_times <- seq(0, max_time, by = uncertainty_dt_weeks)
  if (tail(uncertainty_times, 1) < max_time) uncertainty_times <- c(uncertainty_times, max_time)
  set.seed(uncertainty_seed + match(id, names(fits)))
  dense <- simulate_k1(times, fit$summary$Reff0, fit$summary$S0, fit$summary$best_I0)
  cache_file <- file.path(cache_dir, paste0(id, "_uncertainty_bands.rds"))
  cache_key <- list(
    Reff0 = fit$summary$Reff0, Reff0_SE = fit$summary$Reff0_SE,
    S0 = fit$summary$S0, S0_SE = fit$summary$S0_SE, I0 = fit$summary$best_I0,
    nu = nu_fixed, gamma = gamma_fixed, n_sims = n_uncertainty_sims,
    dt = uncertainty_dt_weeks, max_time = max_time,
    seed = uncertainty_seed + match(id, names(fits))
  )
  cached <- if (!recompute_uncertainty && file.exists(cache_file)) readRDS(cache_file) else NULL
  if (!is.null(cached) && identical(cached$key, cache_key)) {
    uncertainty <- cached$uncertainty
  } else {
    r_draw <- draw_positive(n_uncertainty_sims, fit$summary$Reff0, fit$summary$Reff0_SE)
    s_draw <- draw_positive(n_uncertainty_sims, fit$summary$S0, fit$summary$S0_SE)
    sims <- matrix(NA_real_, nrow = length(uncertainty_times), ncol = n_uncertainty_sims)
    for (j in seq_len(n_uncertainty_sims)) {
      sims[, j] <- simulate_k1(uncertainty_times, r_draw[j], s_draw[j], fit$summary$best_I0)$fitted_incidence
    }
    ci <- t(apply(sims, 1, quantile, probs = c(0.025, 0.975), na.rm = TRUE))
    uncertainty <- data.frame(time = uncertainty_times, PI_95_low = pmax(0, ci[, 1]), PI_95_high = pmax(0, ci[, 2]))
    saveRDS(list(key = cache_key, uncertainty = uncertainty), cache_file)
  }
  origin <- median(as.numeric(raw$date) - 7 * raw$time, na.rm = TRUE)
  dense$date <- as.Date(origin + 7 * dense$time, origin = "1970-01-01")
  dense$PI_95_low <- approx(uncertainty$time, uncertainty$PI_95_low, xout = dense$time, rule = 2)$y
  dense$PI_95_high <- approx(uncertainty$time, uncertainty$PI_95_high, xout = dense$time, rule = 2)$y
  dense$fit_id <- id
  dense$province <- cfg$province
  dense$wave <- cfg$wave
  raw$fitted_incidence <- approx(dense$time, dense$fitted_incidence, xout = raw$time, rule = 2)$y
  raw$fit_id <- id
  raw$province <- cfg$province
  raw$wave <- cfg$wave
  list(dense = dense, overlay = raw)
}

### Display-only connector: first-wave values are retained until the next
### wave begins, then the same province-coloured trajectory continues.
stitch_dense <- function(first, second) {
  handoff <- min(second$date)
  out <- rbind(first[first$date < handoff, , drop = FALSE], second)
  out[order(out$date, out$time), , drop = FALSE]
}
stitch_observed <- function(first, second) {
  out <- rbind(first, second)
  out <- out[order(out$date, out$time), , drop = FALSE]
  out[!duplicated(out$date, fromLast = TRUE), , drop = FALSE]
}
even_ticks <- function(x, n = 9L) {
  z <- range(as.numeric(x), na.rm = TRUE)
  as.Date(seq(z[1], z[2], length.out = n), origin = "1970-01-01")
}

plot_figure1a <- function(dat, file, device = c("png", "pdf")) {
  device <- match.arg(device)
  if (device == "png") png(file, width = 3300, height = 1900, res = 220) else pdf(file, width = 13.2, height = 8.8)
  on.exit(dev.off(), add = TRUE)
  par(mar = c(6.0, 6.8, 1.3, 1.6), mgp = c(3.4, 1.0, 0), tcl = -0.35, cex.axis = 1.45, cex.lab = 1.95, family = "sans")
  all_dates <- do.call(c, lapply(dat, function(z) c(z$dense$date, z$overlay$date)))
  all_y <- do.call(c, lapply(dat, function(z) c(z$overlay$observed, z$dense$fitted_incidence, z$dense$PI_95_high)))
  xlim <- range(all_dates, na.rm = TRUE)
  ylim <- c(0, max(all_y, na.rm = TRUE) * 1.08)
  ticks <- even_ticks(all_dates)
  plot(xlim, ylim, type = "n", xlab = "", ylab = "Weekly measles case incidence", xaxt = "n", bty = "n")
  #axis(1, at = ticks, labels = format(ticks, "%m-%d-%y"), las = 2)
  axis(1, at = ticks, labels = FALSE)
  text(
    x = ticks,
    y = par("usr")[3] - 0.06 * diff(par("usr")[3:4]),
    labels = format(ticks, "%m-%d-%y"),
    srt = 35,
    adj = 1,
    xpd = NA,
    cex = 1.15
  )
  axis(2)
  abline(h = pretty(ylim, 6), v = ticks, col = "#D0D0D0", lty = 1, lwd = 0.9)
  for (province in names(dat)) {
    z <- dat[[province]]$dense
    col <- province_colours[[province]]
    polygon(c(z$date, rev(z$date)), c(z$PI_95_low, rev(z$PI_95_high)), border = NA, col = adjustcolor(col$fill, alpha.f = 0.36))
    lines(z$date, z$PI_95_low, col = adjustcolor(col$line, alpha.f = 0.70), lty = 2, lwd = 1.4)
    lines(z$date, z$PI_95_high, col = adjustcolor(col$line, alpha.f = 0.70), lty = 2, lwd = 1.4)
    lines(z$date, z$fitted_incidence, col = col$line, lwd = 3.6)
    p <- province_colours[[province]]
    points(dat[[province]]$overlay$date, dat[[province]]$overlay$observed,
           pch = p$pch, cex = 1.25, lwd = p$lwd,
           col = p$line)
  }
  provinces <- names(dat)
  legend(
    "topright",
    bty = "n",
    cex = 1.45,
    inset = c(0.01, 0.01),
    legend = as.vector(rbind(
      paste(provinces, "observed incidence"),
      paste(provinces, "best-fit trajectory")
    )),
    pch = as.vector(rbind(
      sapply(provinces, function(p) province_colours[[p]]$pch),
      rep(NA, length(provinces))
    )),
    lty = as.vector(rbind(
      rep(NA, length(provinces)),
      rep(1, length(provinces))
    )),
    lwd = as.vector(rbind(
      rep(NA, length(provinces)),
      rep(3.5, length(provinces))
    )),
    col = as.vector(rbind(
      sapply(provinces, function(p) province_colours[[p]]$line),
      sapply(provinces, function(p) province_colours[[p]]$line)
    ))
  )
#  legend("topright", bty = "n", cex = 1.45, inset = c(0.01, 0.01),
#         legend = as.vector(rbind(paste(provinces, "observed incidence"),
#                                  paste(provinces, "best-fit trajectory"),
##                                  paste(provinces, "95% uncertainty interval"))),
#         pch = as.vector(rbind(sapply(provinces, function(p) province_colours[[p]]$pch), rep(NA, 3), rep(15, 3))),
#         lty = as.vector(rbind(rep(NA, 3), rep(1, 3), rep(NA, 3))),
#         lwd = as.vector(rbind(rep(NA, 3), rep(3.5, 3), rep(NA, 3))),
#         col = as.vector(rbind(rep("#5A5A5A", 3), sapply(provinces, function(p) province_colours[[p]]$line),
#                               sapply(provinces, function(p) adjustcolor(province_colours[[p]]$fill, alpha.f = 0.85)))))
 box(bty = "o", lwd = 1.1)
}

### MAIN
fit_results <- list()
wave_data <- list()
for (id in names(fits)) {
  fit_results[[id]] <- read_fit(id, fits[[id]])
  wave_data[[id]] <- make_wave(id, fits[[id]], fit_results[[id]], read_raw(fits[[id]]))
  message(id, ": I0 = ", fit_results[[id]]$summary$best_I0,
          "; Delta BIC <= 2 range = ", fit_results[[id]]$summary$I0_deltaBIC_le_2_min,
          "-", fit_results[[id]]$summary$I0_deltaBIC_le_2_max)
}

summary_table <- do.call(rbind, lapply(fit_results, function(z) z$summary))
profile_table <- do.call(rbind, lapply(names(fit_results), function(id) {
  z <- fit_results[[id]]$profile[, c("I0", "BIC", "delta_BIC", "delta_BIC_le_2", "Reff0", "S0", "weekly_source_file"), drop = FALSE]
  z$fit_id <- id
  z
}))
write.csv(summary_table, file.path(table_dir, "K1_Figure1A_selected_fit_summary.csv"), row.names = FALSE)
write.csv(profile_table, file.path(table_dir, "K1_Figure1A_fixedRate_I0_profiles.csv"), row.names = FALSE)

province_data <- list(
  Ontario = wave_data$Ontario_2025,
  Alberta = list(dense = stitch_dense(wave_data$Alberta_2025$dense, wave_data$Alberta_2026$dense),
                 overlay = stitch_observed(wave_data$Alberta_2025$overlay, wave_data$Alberta_2026$overlay)),
  Manitoba = list(dense = stitch_dense(wave_data$Manitoba_2025$dense, wave_data$Manitoba_2026$dense),
                  overlay = stitch_observed(wave_data$Manitoba_2025$overlay, wave_data$Manitoba_2026$overlay))
)

combined_dense <- do.call(rbind, lapply(names(province_data), function(p) {
  z <- province_data[[p]]$dense
  z$display_province <- p
  z
}))
combined_overlay <- do.call(rbind, lapply(names(province_data), function(p) {
  z <- province_data[[p]]$overlay
  z$display_province <- p
  z
}))
write.csv(combined_dense, file.path(table_dir, "K1_Figure1A_continuous_dense_uncertainty_values.csv"), row.names = FALSE)
write.csv(combined_overlay, file.path(table_dir, "K1_Figure1A_continuous_observed_overlay_values.csv"), row.names = FALSE)

png_file <- file.path(plot_dir, "K1_Figure1A_Ontario_Alberta_Manitoba_continuous.png")
pdf_file <- file.path(plot_dir, "K1_Figure1A_Ontario_Alberta_Manitoba_continuous.pdf")
plot_figure1a(province_data, png_file, "png")
plot_figure1a(province_data, pdf_file, "pdf")
message("Saved Figure 1A: ", png_file)
message("Saved Figure 1A: ", pdf_file)

# All model times/rates use weeks. Output beta*S*I is the fitted instantaneous
# rate (cases/week); this update does not change it to integrated weekly counts.
# Supplementary figure: sensitivity to latent and infectious periods.
# Run from the project root containing fitData/, Ontario/, Alberta/, Manitoba/.
# Uses saved scan RESULTS (not grid-definition files); no Monolix runs needed.
# For each wave and period pair, choose the minimum-BIC row over its I0 scan.
# Display independent fitted waves: no stitching, shifting, or rescaling.

if (!requireNamespace("deSolve", quietly = TRUE)) stop("Install package 'deSolve' first.")
base_path <- getwd()
fit_dir <- file.path(base_path, "fitData")
out_dir <- file.path(fit_dir, "K1_supplement_period_sensitivity_wkly")
include_reference <- TRUE  # FALSE gives exactly the four corner curves.
dense_dt <- 0.1 / 7  # weeks
rate_tol <- 1e-8

periods <- data.frame(latent_days = c(9, 9, 11, 11),
                      infectious_days = c(7, 9, 7, 9))
if (include_reference) periods <- rbind(periods, c(10, 8))
periods$reference <- periods$latent_days == 10 & periods$infectious_days == 8
periods$label <- paste0("L = ", periods$latent_days, ", D = ",
                       periods$infectious_days, " days")
periods$label[periods$reference] <- paste0(periods$label[periods$reference], " (reference)")

waves <- list(
  Ontario_2025 = list(province = "Ontario", year = "2025", colour = "#2F6FAE", pch = 16,
    fit = c("K1_incidence_Ontario_fitR0S0csv_wkly.csv", "K1_incidence_Ontario_fitR0S0_wkly.csv"),
    raw = "Ontario/MeaslesOntData_Province_2025.csv"),
  Alberta_2025 = list(province = "Alberta", year = "2025", colour = "#D98C2B", pch = 0,
    fit = "K1_incidence_Alberta_fitR0S0_2025_wkly.csv",
    raw = "Alberta/MeaslesAlbData_Province_2025.csv"),
  Alberta_2026 = list(province = "Alberta", year = "2026", colour = "#D98C2B", pch = 0,
    fit = "K1_incidence_Alberta_fitR0S0_2026_wkly.csv",
    raw = "Alberta/MeaslesAlbData_Province_2026_revised.csv"),
  Manitoba_2025 = list(province = "Manitoba", year = "2025", colour = "#3F8F5A", pch = 17,
    fit = "K1_incidence_Manitoba_fitR0S0_2025_wkly.csv",
    raw = "Manitoba/MeaslesManData_Province_2025.csv"),
  Manitoba_2026 = list(province = "Manitoba", year = "2026", colour = "#3F8F5A", pch = 17,
    fit = "K1_incidence_Manitoba_fitR0S0_2026_wkly.csv",
    raw = "Manitoba/MeaslesManData_Province_2026.csv")
)
# To use another saved scan, edit the corresponding 'fit' above.
# weekly_scan() pools matching weekly coarse and fine results for each wave.


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

num <- function(x) suppressWarnings(as.numeric(as.character(x)))
first_file <- function(root, candidates) {
  paths <- file.path(root, candidates)
  found <- paths[file.exists(paths)]
  if (!length(found)) stop("Missing input. Checked: ", paste(paths, collapse = "; "))
  normalizePath(found[1], mustWork = TRUE)
}
column <- function(dat, exact, pattern = NULL) {
  hit <- exact[exact %in% names(dat)]
  if (!length(hit) && !is.null(pattern)) hit <- grep(pattern, names(dat), ignore.case = TRUE, value = TRUE)
  if (!length(hit)) stop("Cannot identify column: ", paste(exact, collapse = ", "))
  hit[1]
}

read_observations <- function(cfg) {
  path <- first_file(base_path, cfg$raw)
  dat <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  tc <- column(dat, c(paste0("Time (weeks) ", cfg$province), "Time (weeks)", "Time_weeks", paste0("Time (days) ", cfg$province), "Time (days)", "Time_days"), "time|week|day")
  yc <- column(dat, c(paste0("Y incidence ", cfg$province), "Y_True_Total", "Incidence"), "incidence|case|y_true")
  dc <- grep("date", names(dat), ignore.case = TRUE, value = TRUE)
  dates <- if (length(dc)) as.Date(as.character(dat[[dc[1]]]), tryFormats = c("%m/%d/%Y", "%Y-%m-%d")) else NULL
  out <- data.frame(time = weekly_time(dat[[tc]], dates), observed = num(dat[[yc]]))
  out <- out[is.finite(out$time) & is.finite(out$observed), , drop = FALSE]
  if (!nrow(out) || any(out$time < 0) || max(out$time) <= 0 || any(out$observed < 0)) {
    stop("Invalid observation times/counts in ", path)
  }
  out[order(out$time), , drop = FALSE]
}

read_scan <- function(cfg) {
  dat <- weekly_scan(fit_dir, cfg$province, cfg$year)
  path <- paste(unique(dat$weekly_source_file), collapse = ";")
  required <- c("success", "nu_fixed", "gamma_fixed", "I0_fixed", "est_ReproNum_pop", "est_S0_pop")
  if (!all(required %in% names(dat))) stop("Missing scan-result columns in ", path)
  bic_col <- column(dat, c("ll_importanceSampling.BIC", "ll_linearization.BIC", "BIC"))
  out <- data.frame(nu = num(dat$nu_fixed), gamma = num(dat$gamma_fixed),
                    I0 = num(dat$I0_fixed), Reff0 = num(dat$est_ReproNum_pop),
                    S0 = num(dat$est_S0_pop), BIC = num(dat[[bic_col]]),
                    source_row = seq_len(nrow(dat)))
  keep <- tolower(trimws(as.character(dat$success))) %in% c("true", "1") &
    apply(out, 1, function(z) all(is.finite(z))) &
    out$nu > 0 & out$gamma > 0 & out$I0 > 0 & out$Reff0 > 0 & out$S0 > 0
  out <- out[keep, , drop = FALSE]
  out$source_file <- dat$weekly_source_file[keep]
  out
}

simulate_curve <- function(best, times) {
  beta <- best$Reff0 * best$gamma / best$S0
  rhs <- function(t, y, p) {
    incidence <- beta * y[1] * y[3]
    list(c(-incidence, incidence - best$nu * y[2],
           best$nu * y[2] - best$gamma * y[3], best$gamma * y[3]))
  }
  sol <- deSolve::ode(y = c(S = best$S0, E = 0, I = best$I0, R = 0),
                      times = times, func = rhs, parms = NULL, method = "lsoda",
                      rtol = 1e-8, atol = 1e-8)
  if (nrow(sol) != length(times) || any(!is.finite(sol))) stop("Incomplete/nonfinite solution.")
  pmax(0, beta * sol[, "S"] * sol[, "I"])
}

panels <- list()
selected <- list()
for (id in names(waves)) {
  cfg <- waves[[id]]
  obs <- read_observations(cfg)
  scan <- read_scan(cfg)
  times <- sort(unique(c(0, seq(0, max(obs$time), by = dense_dt), max(obs$time))))
  curves <- matrix(NA_real_, nrow = length(times), ncol = nrow(periods))
  for (j in seq_len(nrow(periods))) {
    slice <- scan[abs(scan$nu - 7 / periods$latent_days[j]) < rate_tol &
                    abs(scan$gamma - 7 / periods$infectious_days[j]) < rate_tol, , drop = FALSE]
    if (!nrow(slice)) stop(id, ": no successful fits for ", periods$label[j])
    # Resolve exact BIC ties deterministically by smaller I0, then source row.
    slice <- slice[order(slice$BIC, slice$I0, slice$source_row), , drop = FALSE]
    best <- slice[1, , drop = FALSE]
    curves[, j] <- simulate_curve(best, times)
    selected[[length(selected) + 1L]] <- cbind(
      data.frame(wave = id, latent_days = periods$latent_days[j],
                 infectious_days = periods$infectious_days[j], reference = periods$reference[j],
                 examined_I0_min = min(slice$I0), examined_I0_max = max(slice$I0),
                 distinct_I0_values = length(unique(slice$I0))), best)
    message(id, "; ", periods$label[j], "; best I0 = ", best$I0, "; BIC = ", signif(best$BIC, 7))
  }
  # Match the corner-pair order from light to dark within each province.
  base_rgb <- as.numeric(col2rgb(cfg$colour)) / 255
  light <- rgb(t(base_rgb * 0.48 + 0.52))
  dark <- rgb(t(base_rgb * 0.52))
  cols <- colorRampPalette(c(light, cfg$colour, dark))(4)
  if (include_reference) cols <- c(cols, cfg$colour)
  panels[[id]] <- list(cfg = cfg, obs = obs, time = times, curves = curves, colours = cols)
}

draw_figure <- function() {
  # Five panels: Ontario spans the top row; Alberta and Manitoba follow.
  layout(matrix(c(1, 1, 2, 3, 4, 5), nrow = 3, byrow = TRUE), heights = c(1, 1, 1))
  par(mar = c(4.4, 4.8, 2.4, 1), mgp = c(2.8, 0.75, 0),
      las = 1, tcl = 0.3, bg = "white", cex.axis = 1.05, cex.lab = 1.1)
  for (i in seq_along(panels)) {
    dat <- panels[[i]]
    # Independent y scales preserve visibility of the smaller waves.
    ymax <- max(dat$obs$observed, dat$curves) * 1.22
    plot(NA, xlim = c(0, max(dat$time)), ylim = c(0, ymax),
         xlab = "Weeks from fitted wave origin", ylab = "Measles case incidence",
         main = paste0(LETTERS[i], "  ", dat$cfg$province, " ", dat$cfg$year),
         cex.main = 1.2, bty = "o")
    grid(col = "grey88", lty = 1)
    for (j in seq_len(nrow(periods))) {
      lines(dat$time, dat$curves[, j], col = dat$colours[j],
            lwd = if (periods$reference[j]) 2.7 else 2,
            lty = if (periods$reference[j]) 2 else 1)
    }
    points(dat$obs$time, dat$obs$observed, pch = dat$cfg$pch,
           col = "grey35", cex = 0.85, lwd = 1.3)
    legend("topright", legend = c("Observed incidence", periods$label),
           col = c("grey35", dat$colours), pch = c(dat$cfg$pch, rep(NA, nrow(periods))),
           lty = c(NA, ifelse(periods$reference, 2, 1)),
           lwd = c(NA, ifelse(periods$reference, 2.7, 2)),
           bty = "n", cex = 0.78*1.5, inset = 0.01)
    box(lwd = 1.2)
  }
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(out_dir, "FigureS_K1_period_sensitivity_fiveWaves")
write.csv(do.call(rbind, selected), paste0(stem, "_selected_fits.csv"), row.names = FALSE)
png(paste0(stem, ".png"), width = 3600, height = 3600, res = 300, bg = "white")
tryCatch(draw_figure(), finally = dev.off())
pdf(paste0(stem, ".pdf"), width = 12, height = 12, useDingbats = FALSE, bg = "white")
tryCatch(draw_figure(), finally = dev.off())
message("Saved supplementary figure and selected-fit table to: ", out_dir)

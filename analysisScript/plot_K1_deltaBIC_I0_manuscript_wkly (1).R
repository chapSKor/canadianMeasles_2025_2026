# Weekly K1 Delta-BIC profiles: Ontario 2025 and individual Alberta/Manitoba waves.
# Run from the project root containing fitData. Fine runs are optional.
# Old daily files and grid-definition CSVs are never used.
base_path <- getwd()
fitData_dir <- file.path(base_path, "fitData")
out_dir <- file.path(fitData_dir, "K1_manuscript_deltaBIC_I0_plots_wkly")
plot_dir <- file.path(out_dir, "plots")
table_dir <- file.path(out_dir, "tables")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
nu_fixed <- 7 / 10
gamma_fixed <- 7 / 8
rate_tolerance <- 1e-8
delta_BIC_threshold <- 2
png_width <- 1800
png_height <- 1800
png_res <- 220
pdf_width <- 7.2
pdf_height <- 7.2
axis_cex <- 1.35*1.3
label_cex <- 1.65
legend_cex <- 1.3
grid_col <- "#D0D0D0"
point_col <- "#5A5A5A"
star_col <- "#111111"
selected_slices <- list(
  Ontario = list(colour_line = "#2F6FAE", point_pch = 1, point_cex = 1.15*1.3, point_lwd = 1.6),
  Alberta = list(colour_line = "#D98C2B", point_pch = 0, point_cex = 1.18*1.3, point_lwd = 1.8),
  Manitoba = list(colour_line = "#3F8F5A", point_pch = 2, point_cex = 1.18*1.3, point_lwd = 1.8)
)
waves <- data.frame(province = c("Ontario", "Alberta", "Alberta", "Manitoba", "Manitoba"),
                    year = c(2025, 2025, 2026, 2025, 2026))
read_fit_csv <- function(path) {
  d <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("success", "nu_fixed", "gamma_fixed", "I0_fixed", "ll_importanceSampling.BIC")
  if (!all(required %in% names(d))) stop("Missing columns in ", path, ": ",
                                          paste(setdiff(required, names(d)), collapse = ", "))
  for (n in required[-1]) d[[n]] <- suppressWarnings(as.numeric(d[[n]]))
  keep <- as.logical(d$success) %in% TRUE &
    is.finite(d$nu_fixed) & is.finite(d$gamma_fixed) &
    abs(d$nu_fixed - nu_fixed) < rate_tolerance &
    abs(d$gamma_fixed - gamma_fixed) < rate_tolerance &
    is.finite(d$I0_fixed) & is.finite(d$ll_importanceSampling.BIC)
  d <- d[keep, , drop = FALSE]
  if (!nrow(d)) stop("No successful weekly 10-day/8-day fits in ", path)
  data.frame(I0_fixed = d$I0_fixed, BIC_for_plot = d$ll_importanceSampling.BIC)
}
make_profile <- function(d) {
  p <- aggregate(BIC_for_plot ~ I0_fixed, d, min)
  p <- p[order(p$I0_fixed), , drop = FALSE]
  p$delta_BIC <- p$BIC_for_plot - min(p$BIC_for_plot)
  p$equivalent_deltaBIC_le_2 <- p$delta_BIC <= delta_BIC_threshold
  p
}
nice_ylim <- function(y) {
  ymax <- max(c(y, delta_BIC_threshold), na.rm = TRUE)
  ymax <- ymax * 1.08
  c(0, ymax)
}

draw_deltaBIC_I0_plot <- function(profile, best, province, plot_title,
                                  colour_line, point_pch, point_cex, point_lwd,
                                  out_stem,
                                  ylim_override = NULL,
                                  show_legend = TRUE,
                                  legend_position = "topright",
                                  solid_white_panel = FALSE,
                                  show_grid = TRUE) {

  x <- profile$I0_fixed
  y <- profile$delta_BIC

  ylim <- if (is.null(ylim_override)) nice_ylim(y) else ylim_override
  xlim <- range(x, na.rm = TRUE)
  x_pad <- 0.04 * diff(xlim)
  if (!is.finite(x_pad) || x_pad == 0) x_pad <- 1
  xlim <- c(xlim[1] - x_pad, xlim[2] + x_pad)

  plot_one <- function() {
    par(
      mar = c(5.5, 6.2, 1.4, 1.4),
      mgp = c(3.4, 1.0, 0),
      tcl = -0.35,
      cex.axis = axis_cex,
      cex.lab = label_cex,
      family = "sans",
      bg = "white"
    )

    plot(
      x, y,
      type = "n",
      xlim = xlim,
      ylim = ylim,
      xlab = expression(I[0]),
      ylab = expression(Delta~BIC),
      xaxt = "n",
      yaxt = "n",
      bty = "n"
    )

    ### For inset use: paint the plotting region solid white so that
    ### grid lines from the underlying/outset plot do not show through
    ### when the PDF/PNG is layered in Keynote or Illustrator.
    if (isTRUE(solid_white_panel)) {
      usr <- par("usr")
      rect(usr[1], usr[3], usr[2], usr[4], col = "white", border = NA)
    }

    x_ticks <- pretty(xlim, n = 6)
    y_ticks <- pretty(ylim, n = 6)
    axis(1, at = x_ticks)
    axis(2, at = y_ticks, las = 1)

    if (isTRUE(show_grid)) {
      abline(h = y_ticks, v = x_ticks, col = grid_col, lty = 1, lwd = 1.0)
    }

    ### Delta BIC = 2 equivalence threshold.
    ### This is drawn only when it falls inside the plotting limits.
    if (delta_BIC_threshold >= ylim[1] && delta_BIC_threshold <= ylim[2]) {
      abline(h = delta_BIC_threshold, col = "#2E2E2E", lwd = 1.8, lty = 1)
    }

    ### Connecting line through the selected gamma/nu slice.
    lines(
      x, y,
      col = adjustcolor(colour_line, alpha.f = 0.45),
      lwd = 4.1,
      lty = 1
    )

    ### All examined I0 values.
    points(
      x, y,
      pch = point_pch,
      cex = point_cex,
      lwd = point_lwd,
      col = adjustcolor(point_col, alpha.f = 0.95)
    )

    ### Highlight Delta BIC <= 2 region with a slightly larger point.
    equiv <- profile$equivalent_deltaBIC_le_2
    if (any(equiv, na.rm = TRUE)) {
      points(
        x[equiv], y[equiv],
        pch = point_pch,
        cex = point_cex * 1.28,
        lwd = point_lwd + 0.4,
        col = adjustcolor(colour_line, alpha.f = 0.90)
      )
    }

    ### Star at the minimal BIC I0.
    points(
      best$I0_fixed[1],
      best$delta_BIC[1],
      pch = 8,
      cex = 2.2*1.3,
      lwd = 2.0,
      col = star_col
    )

    if (isTRUE(show_legend)) {
      legend(
        legend_position,
        bty = "n",
        cex = legend_cex,
        inset = c(0.015, 0.015),
        legend = c("Examined I0", expression(Delta~BIC <= 2), "Minimum BIC", expression(Delta~BIC == 2)),
        pch = c(point_pch, point_pch, 8, NA),
        lty = c(NA, NA, NA, 1),
        lwd = c(NA, NA, NA, 1.8),
        pt.lwd = c(point_lwd, point_lwd + 0.4, 2.0, NA),
        col = c(point_col, colour_line, star_col, "#2E2E2E"),
        pt.cex = c(point_cex, point_cex * 1.28, 1.8, NA)
      )
    }

    box(bty = "o", lwd = 1.1)

    invisible(NULL)
  }

  png_file <- file.path(plot_dir, paste0(out_stem, ".png"))
  pdf_file <- file.path(plot_dir, paste0(out_stem, ".pdf"))

  png(png_file, width = png_width, height = png_height, res = png_res, bg = "white")
  plot_one()
  dev.off()

  pdf(pdf_file, width = pdf_width, height = pdf_height, bg = "white")
  plot_one()
  dev.off()

  message("Saved: ", png_file)
  message("Saved: ", pdf_file)
}

save_profile_tables <- function(profile, best, out_stem) {
  write.csv(
    profile,
    file.path(table_dir, paste0(out_stem, "_profile.csv")),
    row.names = FALSE
  )

  write.csv(
    best,
    file.path(table_dir, paste0(out_stem, "_best.csv")),
    row.names = FALSE
  )
}


summary_rows <- list()
for (i in seq_len(nrow(waves))) {
  province <- waves$province[i]
  year <- waves$year[i]
  stem <- paste0("K1_incidence_", province, "_fitR0S0_", year)
  candidates <- paste0(stem, "_wkly.csv")
  if (province == "Ontario") candidates <- c(candidates,
    "K1_incidence_Ontario_fitR0S0_wkly.csv", "K1_incidence_Ontario_fitR0S0csv_wkly.csv")
  paths <- file.path(fitData_dir, candidates)
  paths <- paths[file.exists(paths)]
  if (!length(paths)) {
    warning("Skipping ", province, " ", year, ": no weekly fit file; expected ", candidates[1])
    next
  }
  coarse <- read_fit_csv(paths[1])
  roots <- unique(c(stem, sub("_wkly\\.csv$", "", basename(paths[1]))))
  fine_names <- unlist(lapply(roots, function(x) c(paste0(x, "_FineTuneI0_wkly.csv"),
                                                 paste0(x, "_wkly_FineTuneI0.csv"))))
  # Case-insensitive filename matching permits FineTuneI0 or fineTuneI0.
  available <- list.files(fitData_dir, full.names = TRUE)
  fine_paths <- available[tolower(basename(available)) %in% tolower(fine_names)]
  pooled <- coarse
  for (f in fine_paths) pooled <- rbind(pooled, read_fit_csv(f))
  full <- make_profile(pooled)
  prof <- make_profile(coarse)
  prof$delta_BIC <- prof$BIC_for_plot - min(full$BIC_for_plot)
  prof$equivalent_deltaBIC_le_2 <- prof$delta_BIC <= delta_BIC_threshold
  best <- full[which.min(full$BIC_for_plot), , drop = FALSE]
  cfg <- selected_slices[[province]]
  out_stem <- paste0(province, "_", year, "_deltaBIC_I0_wkly")
  draw_deltaBIC_I0_plot(prof, best, paste(province, year), "",
    cfg$colour_line, cfg$point_pch, cfg$point_cex, cfg$point_lwd,
    paste0(out_stem, "_fullRange"), legend_position = if (province == "Ontario") "topleft" else "topright")
  if (length(fine_paths)) {
    fine_i0 <- unique(unlist(lapply(fine_paths, function(f) read_fit_csv(f)$I0_fixed)))
    inset_profile <- full[full$I0_fixed %in% union(fine_i0, best$I0_fixed), , drop = FALSE]
    draw_deltaBIC_I0_plot(inset_profile, best, paste(province, year), "",
      cfg$colour_line, cfg$point_pch, cfg$point_cex, cfg$point_lwd,
      paste0(out_stem, "_fineInset"), ylim_override = c(0, 1),
      show_legend = FALSE, solid_white_panel = TRUE, show_grid = FALSE)
  } else message(province, " ", year, ": fine inset skipped until weekly fine results exist.")
  save_profile_tables(full, best, out_stem)
  supported <- full$I0_fixed[full$equivalent_deltaBIC_le_2]
  summary_rows[[out_stem]] <- data.frame(province = province, year = year,
    source_file = basename(paths[1]), nu_week_inverse = nu_fixed, gamma_week_inverse = gamma_fixed,
    best_I0 = paste(full$I0_fixed[full$BIC_for_plot == min(full$BIC_for_plot)], collapse = ";"),
    minimum_BIC = min(full$BIC_for_plot), support_min = min(supported), support_max = max(supported),
    supported_I0_values = paste(supported, collapse = ";"), fine_available = length(fine_paths) > 0,
    support_reaches_grid_boundary = any(range(full$I0_fixed) %in% supported))
}
if (!length(summary_rows)) stop("No weekly fit files found in ", fitData_dir)
write.csv(do.call(rbind, summary_rows), file.path(table_dir, "K1_deltaBIC_I0_wkly_summary.csv"), row.names = FALSE)
message("Plots: ", plot_dir, "\nTables: ", table_dir)

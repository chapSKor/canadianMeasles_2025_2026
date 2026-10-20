# STEP 2: Plot cached Figure 1 B-E; NO ODE simulations or scan-file reads.
# Run from the project root containing fitData/.
# source("plot_K1_Figure1_BCDE_from_cache.R")

run_plot <- function() {
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2: install.packages('ggplot2')")
library(ggplot2)
base_path <- getwd()  # EDIT if needed; use the same root as the analysis script.
out_dir <- file.path(base_path, "fitData", "K1_Figure1_BCDE_weekly_cached")
cache_file <- file.path(out_dir, "cache", "Figure1_BCDE_plot_cache.rds")
plot_dir <- file.path(out_dir, "plots")

# Appearance only: changing these never triggers analysis.
figure_width_inches <- 14.0
figure_height_inches <- 8.0
panel_height_to_width <- 0.50  # Each plotting rectangle is twice as wide as tall.
png_resolution <- 300
axis_number_size <- 18.9
axis_category_size <- 14.5
strip_text_size <- 14
point_size <- 4.1
errorbar_width <- 0.13
errorbar_linewidth <- 1.15
export_individual_panels <- TRUE
show_diagnostic_caption <- TRUE

if (!file.exists(cache_file)) stop("Missing completed analysis cache. Run analyse_K1_Figure1_BCDE_weekly.R first.")
state_path <- file.path(out_dir, "cache", "analysis_run_state.rds")
if (file.exists(state_path) && readRDS(state_path)$status != "complete")
  stop("Latest analysis is incomplete. Rerun the analysis script to resume its checkpoints.")
b <- readRDS(cache_file)
if (!identical(b$schema_version, 2L)) stop("Unsupported cache version. Rerun the matching analysis script.")
d <- b$plot_data
if (nrow(d) != 20L || anyDuplicated(paste(d$series_id, d$panel))) stop("Incomplete/duplicated plotting data.")
d$label <- factor(d$label, levels = b$label_order)
d$province <- factor(d$province, levels = names(b$province_style))
d$metric <- factor(d$metric, levels = b$metric_order)
colours <- vapply(b$province_style, `[[`, "", "colour")
shapes <- vapply(b$province_style, `[[`, numeric(1), "shape")
headers <- setNames(paste(c("B", "C", "D", "E"), b$metric_order), b$metric_order)
missing_bar <- !is.finite(d$ymin) | !is.finite(d$ymax)
if (any(missing_bar)) warning("Some intervals cannot be plotted; see diagnostics and the figure caption.", call. = FALSE)
notes <- character()
s <- b$summary
for (i in which(s$FWHM_undefined_draws > 0)) {
  conditional <- grepl("conditional", s$FWHM_interval_type[i])
  notes <- c(notes, paste0(s$province[i], " ", s$year[i], ": FWHM ",
    if (conditional) "interval conditional on " else "interval omitted; ",
    s$FWHM_defined_draws[i], "/", s$n_draws[i], " defined draws."))
}
if (any(missing_bar & d$panel != "E")) notes <- c(notes,
  "Non-finite intervals omitted; see cached summary.")
caption <- if (show_diagnostic_caption && length(notes))
  paste(strwrap(paste(notes, collapse = " "), width = 140), collapse = "\n") else NULL

make_panel <- function(data, columns = 2L, caption_text = NULL) {
  ggplot(data, aes(x = label, y = value, colour = province, shape = province)) +
    geom_errorbar(data = data[is.finite(data$ymin) & is.finite(data$ymax), , drop = FALSE],
      aes(ymin = ymin, ymax = ymax), width = errorbar_width,
      linewidth = errorbar_linewidth, alpha = 0.82) +
    geom_point(size = point_size, stroke = 1.55, na.rm = TRUE) +
    facet_wrap(~metric, scales = "free_y", ncol = columns,
      labeller = labeller(metric = headers)) +
    scale_colour_manual(values = colours) + scale_shape_manual(values = shapes) +
    labs(x = NULL, y = NULL, caption = caption_text) +
    theme_bw(base_size = 14) +
    theme(panel.grid.major = element_line(colour = "#D0D0D0", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.75),
      axis.ticks = element_line(colour = "black", linewidth = 0.5),
      axis.ticks.length = grid::unit(-0.13, "cm"),
      axis.text.x = element_text(size = axis_category_size, colour = "#1F1F1F", margin = margin(t = 8)),
      axis.text.y = element_text(size = axis_number_size, colour = "#1F1F1F", margin = margin(r = 8)),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.65),
      strip.text = element_text(size = strip_text_size, colour = "#1F1F1F", face = "bold"),
      panel.spacing.x = grid::unit(1.90, "lines"),
      panel.spacing.y = grid::unit(1.25, "lines"),
      aspect.ratio = panel_height_to_width,
      legend.position = "none", plot.margin = margin(14, 18, 14, 18),
      plot.caption = element_text(size = 10, hjust = 0)) + coord_cartesian(clip = "off")
}
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
save_plot <- function(p, stem, width, height) {
  for (extension in c("png", "pdf")) ggsave(file.path(plot_dir, paste0(stem, ".", extension)),
    p, width = width, height = height, units = "in", dpi = png_resolution, bg = "white")
}
p <- make_panel(d, caption_text = caption)
save_plot(p, "Figure1_BCDE_weekly_95CI_2x2", figure_width_inches, figure_height_inches)
if (export_individual_panels) {
  for (letter in c("B", "C", "D", "E")) {
    sub <- d[d$panel == letter, , drop = FALSE]
    note <- if (letter == "E") caption else NULL
    save_plot(make_panel(sub, columns = 1L, caption_text = note),
      paste0("Figure1", letter, "_weekly_95CI"), figure_width_inches/2, figure_height_inches/2+0.3)
  }
}
print(p)
message("Plots saved to: ", plot_dir)
invisible(p)
}
run_plot()

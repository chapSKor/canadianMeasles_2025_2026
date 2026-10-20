README: Weekly measles fits — Figure 1 and supplementary plots
Updated 13 October 2026
============================================================

QUICK START
-----------
Place the updated R files in the project root containing fitData/, Ontario/,
Alberta/, and Manitoba/. Set getwd() to that root before sourcing a script.
The completed Monolix scans are inputs; plotting does not refit the models.

Main Figure 1A:
  source("plot_K1_Figure1A_fixedRates_continuous_CacheFiles.R")

Main Figure 1B–D (parameter dot plots):
  source("plot_K1_inferred_parameter_dotplot_fixedRates_fiveWaves_Fig1BCD.R")

Ontario I0 profile and separate provincial profiles (Figure 1E source):
  source("plot_K1_deltaBIC_I0_manuscript_wkly.R")

Combined Alberta/Manitoba I0 profiles and separate inset (Figure 1F source):
  source("plot_K1_deltaBIC_I0_Alberta_Manitoba_wkly.R")

Supplementary period-sensitivity figure:
  source("plot_K1_supplement_period_sensitivity_fiveWaves_wkly.R")

Best-fit SEIR final size (manuscript Table 2; NOT a plotting script):
  Use the folder containing K1_fixedRates_fiveWaves_audit.csv.
  Keep compute_SEIR_final_size_bestfits.R in that SAME folder.
  Set getwd() to that folder, then run:
    source("compute_SEIR_final_size_bestfits.R")
  Output in that folder: SEIR_bestfit_final_sizes.csv

These plotting scripts are independent. The exploratory analysis script is
NOT a prerequisite for Figure 1A or the dot plots. Assemble the final panels
and position the separately exported inset manually in Keynote/Illustrator.

FINAL-SIZE CALCULATION: LOCATION, SCRIPT, INPUT, OUTPUT
--------------------------------------------------------
Purpose: calculate the asymptotic SEIR final size for the minimum-BIC best
fit of each of the five waves, for inclusion in manuscript Table 2.
This is a downstream calculation: it does NOT refit Monolix models, load the
raw Ontario/Alberta/Manitoba incidence files, or produce a figure.

WHERE TO RUN IT
  Work in the SAME folder that holds the five-wave audit CSV. The final-size
  R script has been left in that folder alongside the audit file. This is
  the calculation folder; it may differ from the project root used by the
  figure plotting scripts. Do not assume R's current working directory is
  automatically the folder where the script was saved.

  Required files together in that folder:
    K1_fixedRates_fiveWaves_audit.csv      <- input (five best-fit records)
    compute_SEIR_final_size_bestfits.R     <- calculation script

  Generated there when the script is run:
    SEIR_bestfit_final_sizes.csv           <- output (five waves)

HOW TO RUN
  Interactive R or RStudio on Mac:
    setwd("/full/path/to/folder/containing/the/audit/csv")
    source("compute_SEIR_final_size_bestfits.R")

  Alternatively, from Terminal after changing to that folder:
    Rscript compute_SEIR_final_size_bestfits.R

  Check location before running if necessary:
    getwd()
    file.exists("K1_fixedRates_fiveWaves_audit.csv")   # must be TRUE

  The original script uses relative input/output filenames, so placing the
  R script in the same folder as the audit CSV is not enough unless that
  folder is also the current working directory (or absolute paths are
  supplied to Rscript). Missing-input errors will otherwise cascade.

INPUT AND METHOD
  Input: K1_fixedRates_fiveWaves_audit.csv, containing one selected
  minimum-BIC fit per wave with unrounded S0, Reff0 and I0_minBIC.
  The script also carries forward the audit's weekly growth rate, doubling
  time and FWHM columns; it does not independently recompute these.
  For each wave, it numerically solves for S_infinity from:

    log(S0/S_infinity) = Reff0*(1 + I0/S0 - S_infinity/S0)

  under S(0)=S0, E(0)=R(0)=0 and I(0)=I0, then computes:

    model_final_size_new_infections = S0 - S_infinity
    model_total_including_initial_I0 = S0 - S_infinity + I0

  Final size counts infections occurring AFTER the start of the fitted
  wave; it does not include the initially infectious I0 individuals.
  The relation is exact for this deterministic, constant-parameter SEIR
  model, conditional on the fitted parameters. It does NOT equal or enforce
  agreement with observed cumulative reported cases, and the output is not
  an independent validation of the incidence fit. Final size itself is a
  count of individuals, not a weekly rate; the model is fitted in weeks.

OUTPUT
  Output: SEIR_bestfit_final_sizes.csv in the same working folder as the
  input. Important columns include S_infinity,
  model_final_size_new_infections, model_final_size_rounded,
  model_total_including_initial_I0 and model_attack_percent_of_S0.
  For Table 2 use model_final_size_rounded (one value per wave), with an
  explanation that this is the model-derived asymptotic final size.
  The script requires only base R. Re-running overwrites the output CSV.

SCOPE AND UNITS
---------------
Five independently fitted waves:
  Ontario 2025; Alberta 2025 and 2026; Manitoba 2025 and 2026.
There is no Ontario 2026 fit. Year labels identify fitted waves, not strict
January-to-December windows; the actual raw dates determine their placement.

Model time is in weeks. Main fixed rates:
  nu = 7/10 = 0.7 week^-1       (mean latent period 10 days)
  gamma = 7/8 = 0.875 week^-1   (mean infectious period 8 days)

K1 SEIR:
  dS/dt = -beta*S*I
  dE/dt = beta*S*I - nu*E
  dI/dt = nu*E - gamma*I
  dR/dt = gamma*I
  beta = Reff0*gamma/S0
  E(0) = R(0) = 0; S(0) = S0; I(0) = I0.

Reff0 is the initial effective reproduction number within the fitted
participating population. S0 is the effective susceptible population (also
labelled Neff in plots); I0 is the initial infectious seed size.

The plotted observation equation remains beta*S*I, an instantaneous rate in
cases/week. The scripts do not change it to an integrated weekly count.
Observed weekly counts are left unchanged. Corrected weekly fit estimates
must be used; do not combine them with previous daily-fit population estimates.

INPUT FILES
-----------
All fit-result CSVs are under fitData/.

Standard coarse files:
  K1_incidence_Ontario_fitR0S0_2025_wkly.csv
  K1_incidence_Alberta_fitR0S0_2025_wkly.csv
  K1_incidence_Alberta_fitR0S0_2026_wkly.csv
  K1_incidence_Manitoba_fitR0S0_2025_wkly.csv
  K1_incidence_Manitoba_fitR0S0_2026_wkly.csv

Ontario alternatives recognized by the current readers:
  K1_incidence_Ontario_fitR0S0_wkly.csv
  K1_incidence_Ontario_fitR0S0csv_wkly.csv

Fine-tuned results follow either pattern:
  K1_incidence_<Province>_fitR0S0_<Year>_wkly_FineTuneI0.csv
  K1_incidence_<Province>_fitR0S0_<Year>_FineTuneI0_wkly.csv
Matching fine suffixes are also recognized for the Ontario alternative names
by the five-wave readers and individual-profile script.

Use results files containing estimates and BIC. The *_grid.csv companions
are grid definitions and are not used for these plots. Old daily filenames
are not fallbacks. Keep only intended runs under the recognized filenames.

The five-wave readers pool all recognized coarse/fine result files. The
individual Delta-BIC script selects the first existing coarse alternative and
pools matching fine files. The combined Alberta/Manitoba Delta-BIC script
uses one coarse file and the first matching fine file per province: keep one
fine filename convention for each wave to avoid ambiguity.

Important fit columns:
  success, nu_fixed, gamma_fixed, I0_fixed,
  est_ReproNum_pop, est_S0_pop, ll_importanceSampling.BIC.
Trajectory uncertainty and dot plots also need Monolix standard-error columns.

Raw-data locations are unchanged:
  Ontario/MeaslesOntData_Province_2025.csv
  Alberta/MeaslesAlbData_Province_2025.csv
  Alberta/MeaslesAlbData_Province_2026_revised.csv
  Manitoba/MeaslesManData_Province_2025.csv
  Manitoba/MeaslesManData_Province_2026.csv

Trajectory scripts have raw_time_units = "auto". Calendar-date spacing is
used to detect daily versus weekly numeric time; daily coordinates are divided
by 7. If needed, explicitly choose "weeks" or "days". If dates cannot establish
the units, the helper assumes weeks. Preserve the fitted wave's original time
zero. Figure 1A reconstructs calendar dates using 7 days per model week.

SCRIPT GUIDE
============

1. plot_K1_Figure1A_fixedRates_continuous_CacheFiles.R
----------------------------------------------------
Purpose: Figure 1A, all five fitted waves on one calendar-date plot.

Reads weekly coarse/fine results directly, selects the minimum-BIC fit within
the fixed 10-day/8-day slice for each wave, simulates trajectories, and overlays
raw observations. Uses province-coloured continuous lines and shaded polygons.
Connections between successive waves are display-only: the models are fitted
independently and their fitted I0 values are retained.

Uncertainty: 500 independent positive normal draws of Reff0 and S0 using their
Monolix SEs; I0, nu and gamma remain fixed. The 2.5% and 97.5% pointwise quantiles
form the band. It excludes parameter covariance, residual observation noise,
and uncertainty across fixed I0/rate choices; it is not a prediction interval.

Output root:
  fitData/K1_figure1A_continuous_wkly/
  plots/: K1_Figure1A_Ontario_Alberta_Manitoba_continuous.png and .pdf
  tables/: selected-fit summary, I0 profiles, observed overlays, dense bands
  cache/: per-wave uncertainty .rds files

Settings: n_uncertainty_sims, uncertainty_seed, dense_dt_weeks,
uncertainty_dt_weeks, recompute_uncertainty, raw_time_units.
Default resolution is 0.1 day/7 for best-fit curves and 1 day/7 for uncertainty;
bands are interpolated onto the dense grid. Cache keys include fitted values,
SEs, rates, simulation count, uncertainty spacing, horizon and seed.
Set recompute_uncertainty = TRUE to force recalculation. Weekly cache/output
folders are separate from old daily outputs.

2. Fig1A_plot.R
--------------
Purpose: convenience launcher for script 1.
It sources plot_K1_Figure1A_fixedRates_continuous_CacheFiles.R from the project
root. Run either this launcher OR script 1, not both. The former simplified
52-week split implementation has been replaced by the actual-date workflow.

3. plot_K1_inferred_parameter_dotplot_fixedRates_fiveWaves_Fig1BCD.R
----------------------------------------------------------------
Purpose: main Figure 1B–D parameter summaries for all five waves.

B: Reff0; C: effective susceptible population; D: I0.
Coarse/fine results are pooled, and all point estimates come from the same
minimum-BIC fit per wave. Old hard-coded I0 overrides are removed.

Outputs include wide and vertically stacked PNG/PDF versions:
  Default: Reff0 and S0 error bars are +/-1 SE; I0 bars are the evaluated
           Delta BIC <= 2 support range.
  Alternative stacked: all three bars span values across the selected
           Delta BIC <= 2 I0 profile. These are model-support ranges.

Output root:
  fitData/K1_manuscript_parameter_summary_fixedRates_fiveWaves_wkly/
  plots/: Plot1BCD_K1_inferred_parameter_dotplot_fixedRates_fiveWaves.*
          corresponding *_square.* and *_deltaBIC_lt2_ranges_* versions
  tables/: selected estimates, plotting values, and per-wave I0 profiles

Some legacy filenames/column names retain "lt2"; the updated comparison is <=2.
S0/Neff is displayed in thousands. Settings: axis_text_size, strip_text_size,
point_size, errorbar_width, errorbar_lwd, se_multiplier, and output dimensions.
The Fig1BCD version sets y-axis text size directly to 18.9 in theme().

4. plot_K1_inferred_parameter_dotplot_fixedRates_fiveWaves_wkly.R
--------------------------------------------------------------
Purpose: updated simpler dot-plot variant from the supplied fiveWaves(1) file.
Uses weekly pooled fits and retains its supplied styling. Produces the standard
SE/support-range wide and stacked plots; it does not produce the additional
all-parameter profile-support version available in script 3.

Output root:
  fitData/K1_manuscript_parameter_summary_fixedRates_fiveWaves_basic_wkly/
Use script 3 for the principal manuscript panels; this is an alternative.

5. plot_K1_deltaBIC_I0_manuscript_wkly.R
--------------------------------------
Purpose: separate Ontario, Alberta and Manitoba I0 profiles. Use Ontario's
output for Figure 1E; the other separate profiles can be used for diagnostics
or alternative panel arrangements.

Year settings near the top:
  selected_year <- 2026  # change to 2025 and rerun for those AB/MB waves
  ontario_year <- 2025
  alberta_year <- selected_year
  manitoba_year <- selected_year

Exports full profiles and fine insets if fine files exist; otherwise produces
provisional coarse-point insets. It uses the fixed weekly 10/8-day slice.
Full profiles show coarse evaluated points, with Delta BIC measured relative
to the pooled coarse/fine minimum when fine results exist. Stars mark minima.

Output root:
  fitData/K1_manuscript_deltaBIC_I0_plots_wkly/
  plots/: <Province>_<Year>_deltaBIC_I0_wkly_fullRange.png/.pdf
          corresponding *_fineInset.* or *_provisionalInset.*
  tables/: profiles, minimum rows and year-selection-specific support summary

Settings: selected_year and province years; axis_cex for axis numbers,
label_cex for titles; selected_slices point_cex/point_lwd for markers.
The minimum star has a separate cex in draw_deltaBIC_I0_plot().

6. plot_K1_deltaBIC_I0_Alberta_Manitoba_wkly.R
-------------------------------------------
Purpose: combined Alberta/Manitoba full profiles and one separate white inset
with both provinces overlaid. Source for Figure 1F.

Set selected_year to 2025 or 2026 and rerun; alternatively set alberta_year and
manitoba_year independently. Inputs, legends and output names follow these
year settings. The main panel shows coarse points, a faint grid and a solid
Delta BIC = 2 reference line. The inset has no legend, uses y=0–1, and marks
minimum BIC values with black asterisks.

Fine runs are optional. When present, coarse/fine results are pooled for the
inset and the common BIC baseline within each province. When either fine file
is absent, the inset is named a provisional zoom and uses available evaluated
points. inset_xlim <- NULL automatically spans supported I0 values; replace
NULL with c(lower, upper) to set it manually.

Output root (different from the other scripts):
  plots/weekly_deltaBIC/
  PanelF_Alberta_<year>_Manitoba_<year>_coarse_deltaBIC_I0_wkly.png/.pdf
  corresponding *_fine_deltaBIC_I0_insets_wkly.* or
                *_provisional_zoom_deltaBIC_I0_wkly.*
  per-wave profiles and year-specific support-summary CSVs

Year-specific names keep 2025 and 2026 outputs separate. Rerunning the same
year/settings overwrites that selection's outputs. Individual per-wave outputs
can be regenerated by multiple runs; Ontario 2025 is the same selection.

7. plot_K1_supplement_period_sensitivity_fiveWaves_wkly.R
------------------------------------------------------
Purpose: supplementary five-panel figure of period sensitivity.
Panel A: Ontario 2025; B/C: Alberta 2025/2026; D/E: Manitoba 2025/2026.

For each wave, pool weekly coarse/fine results and choose the minimum-BIC I0
fit separately at each (latent days, infectious days) pair:
  (9,7), (9,9), (11,7), (11,9), plus the (10,8) reference.
Rates are 7/latent_days and 7/infectious_days per week. Light-to-dark provincial
curves show the four perturbations; the reference is dashed. Raw points are
grey. The x-axis is weeks from fitted wave origin. These waves are not stitched.

Output root:
  fitData/K1_supplement_period_sensitivity_wkly/
  FigureS_K1_period_sensitivity_fiveWaves.png/.pdf
  FigureS_K1_period_sensitivity_fiveWaves_selected_fits.csv

Settings: include_reference, dense_dt (weeks), raw_time_units.
Missing fits for a required period combination stop the script. The selected-fit
CSV describes the plotted curves; it is not the full nine-combination Table S1.

8. analyze_K1_profile_grid_with_overlay_v10_albertaColumns_wkly.R
---------------------------------------------------------------
Purpose: optional exploratory fit diagnostics and single-wave overlays.
It is not needed to create the main figures, which read fits directly.

Settings: province_vec, selected_year, fixed_rates_only.
Ontario uses 2025; Alberta and Manitoba use selected_year. With
fixed_rates_only=TRUE (default), analysis uses nu=0.7 and gamma=0.875.
FALSE explores all saved weekly rate combinations. Periods are reported in days
as 7/nu and 7/gamma. The plausible-window settings are 9–11 latent days and
7–9 infectious days. Best-fit overlays use 1000 uncertainty draws by default.

Output root:
  fitData/K1_analysis_wkly/<Province>_<Year>/fixedRates/
  or .../allRates/
  Each contains plots/ and tables/ for estimates, BIC/AIC diagnostics,
  parameter/SE comparisons and fitted-trajectory uncertainty.

This script has more likelihood-column requirements than the main plotters,
including AIC. It can be slower because it uses its own integration/uncertainty
workflow rather than the main Figure 1A cache.

MODEL SELECTION AND INTERVALS
----------------------------
For the main plots, select the fixed 10/8-day rate slice within each wave.
For repeated runs at the same I0, the main plotters keep the lowest successful
BIC. Delta BIC = BIC - minimum BIC within that wave and selected rate slice.
The support rule is Delta BIC <= 2. For sensitivity curves, select the minimum
separately within each period pair. Do not compare absolute BIC across waves
with different observations as a measure of relative outbreak fit quality.

I0 support endpoints span evaluated values, not interpolated crossings. Coarse
spacing or scan boundaries can truncate the range. Delta-BIC summaries include
boundary flags and supported-value lists to guide further fine scans.
Monolix SE bars, trajectory-uncertainty polygons and profile-support ranges
are different summaries; identify the displayed version in the caption.

COLOURS AND PRESENTATION
-----------------------
Ontario:  #2F6FAE; fill #BFD7F1; circles.
Alberta:  #D98C2B; fill #F2C38B; hollow squares.
Manitoba: #3F8F5A; fill #B7DDB8; triangles.
Figure 1A colours observations by province; supplementary observations are
grey. Individual Delta-BIC panels use grey examined points with coloured
support highlights. The combined Delta-BIC panel uses province-coloured points.
Tick directions, text settings and dimensions remain script-specific.

UPSTREAM SCANS AND OLDER SCRIPTS
-------------------------------
scanFits_Ontario_Alberta_K1K10_nu_gamma_I0.R is the historical scan-script name.
The completed user-run weekly scans supply the current CSVs. The recent plotting
updates did not rerun Monolix or replace the user's scan configuration.
Weekly scans require weekly observation times and nu=7/L_days, gamma=7/D_days;
filenames must identify province, wave year, weekly units and fine runs.

Superseded for current Figure 1:
  plot_K1_Ontario_Alberta_combined_v5_FINAL_cleanStyle.R
  plot_K1_deltaBIC_I0_manuscript_v5_insetNoGrid.R
  plot_K1_inferred_parameter_dotplot_v4_moreSpacing_biggerText.R
  earlier non-weekly five-wave variants

Do not follow the old Ontario/Alberta-only panel mapping or its 12/4-day and
10/9-day rate selections. The current mapping is A=trajectories, B–D=parameters,
E=Ontario I0 profile, F=combined Alberta/Manitoba I0 profiles. Year-specific
alternatives can be assembled as desired.

Older table-generation scripts for main-table growth rate, doubling time,
FWHM and the full nine-pair Table S1 were not part of the six-script weekly
update. They must be checked/updated before regenerating those manuscript
quantities from scratch. The SEPARATE, downstream script
compute_SEIR_final_size_bestfits.R calculates the final size for each
minimum-BIC fit using K1_fixedRates_fiveWaves_audit.csv, and copies the
existing growth-rate, doubling-time and FWHM values from that audit into
SEIR_bestfit_final_sizes.csv; see FINAL-SIZE CALCULATION above.
With weekly rates, r is per week, doubling time is in weeks, and simulated FWHM
is in weeks; convert r/7 to per day and times*7 to days when reporting daily units.

VALIDATION STATUS
-----------------
The six scripts in weekly_figure_scripts.zip passed syntax-tree parsing.
Full R execution/rendering against the new weekly CSVs was not performed here.
The Delta-BIC scripts received static checks in the earlier update.
Dependencies: deSolve for main/supplementary trajectories; ggplot2 for dot plots;
base R for Delta-BIC plots, the exploratory analysis workflow, and the
best-fit SEIR final-size script.

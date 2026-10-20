FIGURE 1 B-E: TWO-SCRIPT WEEKLY WORKFLOW

If RStudio still shows Browse[1] or Browse[2], enter Q to leave debugging.
Restart R (Session > Restart R) if recursive errors persist. Then use source()
to run each complete script; do not paste subsequent lines after an error.

Put both scripts in the project root containing fitData/. Set that folder as
the working directory. Install deSolve and ggplot2 if necessary:
  install.packages(c("deSolve", "ggplot2"))

1. Run analysis:
  source("analyse_K1_Figure1_BCDE_weekly.R")
2. Run plotting:
  source("plot_K1_Figure1_BCDE_from_cache.R")

For subsequent changes to colours, font sizes, dimensions or symbols, run
ONLY script 2. It does not read scan files, refit models or simulate trajectories.
Default output: wide 14 x 8 inch 2x2 PNG/PDF; each panel has height/width=0.5.
Individual panels are also exported. Change appearance at the top of script 2.

INPUTS
The same weekly coarse and FineTuneI0 results in fitData/ are pooled separately
for Ontario 2025, Alberta 2025/2026 and Manitoba 2025/2026. Grid-definition
files and daily fits are excluded. nu=0.7/week and gamma=0.875/week are fixed.
The minimum-BIC row is selected after retaining the best row for each I0.
Named Monolix standard errors are used. The previous audit CSV is no longer
required: growth rate, doubling time and FWHM are calculated in this analysis.

OUTPUT ROOT
fitData/K1_Figure1_BCDE_weekly_cached/
  cache/ : per-wave draw checkpoints and Figure1_BCDE_plot_cache.rds
  tables/: selected-fit summary, plotting values, profiles, every parameter
           draw and FWHM status counts/fractions
  plots/ : wide 2x2 and individual-panel PNG/PDF files

CHECKPOINTS
Each wave saves its parameter draws before simulation and checkpoints every
100 draws. An interrupted run resumes from the matching checkpoint. All draws
are generated reproducibly before simulation using a fixed per-wave seed.
Source-file hashes, estimates/SEs, rate constants, numerical settings and solver
version invalidate incompatible checkpoints. Different interval-reporting
settings can reuse trajectories. The final RDS contains all data needed to plot.
An incomplete analysis is flagged so plotting does not silently use an older run.
This version starts a fresh cache; old CSV draws without a matching provenance
signature are not silently reused.

UNCERTAINTY AND THE ORIGINAL FAILURE
B/C: approximate marginal 95% Wald confidence intervals, estimate +/-1.96 SE.
D/E: 2.5th/97.5th percentiles from 1000 independent positive-normal Reff0/S0
parameter draws. I0 and the latent/infectious periods remain fixed. Parameter
covariance, profile uncertainty, model uncertainty and observation noise are
not included. These are conditional parameter-uncertainty intervals, not full
prediction intervals.

The old code stopped the whole run when ANY draw had an undefined FWHM.
Now each draw records a diagnostic status: ok, peak_at_origin,
no_left_half_crossing, right_crossing_or_peak_beyond_horizon, or numerical error.
FWHM uses the global peak and the nearest half-height crossing on each side.
A trajectory can start too high to contain a rising half-height crossing; a
finite FWHM then cannot be measured over its nonnegative simulation times.
Such draws are NOT replaced, assigned zero, or silently redrawn.

When >=95% of draws have a defined width and no numerical failures occur, the
FWHM interval is computed from defined widths. If any widths are undefined,
it is explicitly labelled CONDITIONAL ON DEFINED FWHM. Counts/fractions appear
in the CSVs and a short plot caption. The 95% defined-draw threshold is only a
reporting safeguard: this is NOT a claim of unconditional 95% coverage.
If fewer than 95% are defined, or numerical failures occur, analysis still
finishes and the plot shows the best-fit point without a FWHM error bar, with
a diagnostic note. Examine diagnostics before publication. Do not remove the
conditional qualification from the manuscript caption if draws were excluded.
Non-growing draws retain infinite doubling times rather than being discarded;
an unbounded percentile interval is flagged and not drawn as a finite bar.

FWHM is calculated from beta*S*I in cases/week with E(0)=0, around each draw's
own peak. Crossing times are linearly interpolated on a 1/14-week output grid;
LSODA uses adaptive internal steps. The horizon extends if needed. All reported
time quantities use weeks. No additional Monolix runs are needed.

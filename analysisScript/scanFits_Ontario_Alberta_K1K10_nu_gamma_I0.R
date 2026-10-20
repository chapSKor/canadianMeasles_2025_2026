############################################################
### 3D PROFILE LIKELIHOOD GRID: nu x gamma x I0          ###
### Ontario + Alberta, K1 through K10                    ###
###                                                       ###
### At each grid point:                                  ###
###   - nu, gamma, and I0 are FIXED                      ###
###   - ReproNum, S0, a, and b are re-fit by MLE         ###
###   - Monolix outputs are streamed to CSV              ###
############################################################

### ========================================================
### 0. USER SETTINGS: EDIT THESE NEAR THE TOP
### ========================================================

### ---- Monolix connector path ----
pkg_path <- "/Applications/MonolixSuite2024R1.app/Contents/Resources/monolixSuite/connectors/lixoftConnectors.tar.gz"

# If needed, uncomment this line once:
# install.packages(pkg_path, repos = NULL, type = "source")

library(lixoftConnectors)
initializeLixoftConnectors(software = "monolix")
getLixoftConnectorsState()

### ---- Parent folder containing Ontario/ and Alberta/ ----
### This script is intended to live inside:
###   June23_ScanFits_Ontario_Alberta/
### with projects such as:
###   Ontario/K1_fitR0/K1Fit_R0.mlxtran
###   Ontario/K2_fitR0/K2Fit_R0.mlxtran
###   ...
###   Alberta/K10_fitR0/K10Fit_R0.mlxtran

base_path <- getwd()

### ---- Output folder ----
out_dir <- file.path(base_path, "fitData")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
######

# LINES 188 - 189 NEED TO BE EDITED IF YOU BROKE THE <K1_FitR0> FILE/DATA LABELLING PATTERN

######
### ---- Provinces and structural models to scan ----
#province_vec <- c("Ontario", "Alberta")
#province_vec <- c("Ontario")
#province_vec <- c("Alberta")
province_vec <- c("Manitoba")
K_vec        <- paste0("K", 1:1)
#K_vec        <- paste0("K", 1:10)
### ---- Data/model label used in output filenames ----
data_label <- "incidence"
#fit_label  <- "fitR0S0_2026_wkly"
fit_label  <- "fitR0S0_2025_wkly_FineTuneI0"

### ========================================================
### 1. SCAN RANGES: nu, gamma, and I0
### ========================================================

### nu = exposed -> infectious transition rate
### Latent period: 11, 10, 9 days; nu expressed in week^(-1)
nu_period_min_days <- 9
nu_period_max_days <- 11
nu_period_grid_days <- seq(
  from = nu_period_max_days,
  to   = nu_period_min_days,
  by   = -1
)
nu_period_grid_weeks <- nu_period_grid_days / 7
nu_grid <- 1 / nu_period_grid_weeks

### gamma = infectious -> recovered transition rate
### Infectious period: 9, 8, 7 days; gamma expressed in week^(-1)
gamma_period_min_days <- 7
gamma_period_max_days <- 9
gamma_period_grid_days <- seq(
  from = gamma_period_max_days,
  to   = gamma_period_min_days,
  by   = -1
)
gamma_period_grid_weeks <- gamma_period_grid_days / 7
gamma_grid <- 1 / gamma_period_grid_weeks

### I0 = number of infectious individuals seeding the outbreak
### EDIT THIS RANGE as needed. A linear integer scan is used here.
I0_min <- 1
I0_max <- 15
I0_by  <- 1
I0_grid <- seq(from = I0_min, to = I0_max, by = I0_by)
#overwrite to make smaller
#I0_grid <- c(1,10, 15, 20, 25, 30, 35, 40, )

stopifnot(all(is.finite(nu_grid)),    all(nu_grid > 0))
stopifnot(all(is.finite(gamma_grid)), all(gamma_grid > 0))
stopifnot(all(is.finite(I0_grid)),    all(I0_grid > 0))

message("nu scan: latent period ", min(nu_period_grid_days), " to ",
        max(nu_period_grid_days), " days; nu range ",
        signif(min(nu_grid), 4), " to ", signif(max(nu_grid), 4))

message("gamma scan: infectious period ", min(gamma_period_grid_days), " to ",
        max(gamma_period_grid_days), " days; gamma range ",
        signif(min(gamma_grid), 4), " to ", signif(max(gamma_grid), 4))

message("I0 scan: ", min(I0_grid), " to ", max(I0_grid),
        " infectious individuals, by ", I0_by)

### ========================================================
### 2. INITIAL VALUES FOR MLE PARAMETERS
### ========================================================

### These are only starting values for parameters that are re-fit.
### They do not fix the final values.
initial_by_province <- list(
  Ontario  = list(ReproNum = 1.8, S0 = 3400, a = 1, b = 0.2),
  Alberta  = list(ReproNum = 1.8, S0 = 3400, a = 1, b = 0.2),
  Manitoba = list(ReproNum = 1.5, S0 = 500, a = 1.56, b = 0.36)
)

### ========================================================
### 3. HELPER FUNCTIONS
### ========================================================

flatten_named <- function(x, prefix) {
  if (inherits(x, "error") || is.null(x)) {
    return(setNames(NA_character_, paste0(prefix, "_error")))
  }

  df <- NULL
  if (is.data.frame(x)) {
    df <- x
  } else {
    df <- tryCatch(as.data.frame(x, stringsAsFactors = FALSE),
                   error = function(e) NULL)
  }

  if (!is.null(df) && "name" %in% names(df)) {
    val_col <- intersect(c("value", "estimate", "val", "Value"), names(df))[1]
    if (is.na(val_col)) val_col <- setdiff(names(df), "name")[1]
    vals <- df[[val_col]]
    names(vals) <- paste0(prefix, "_", df$name)
    return(vals)
  }

  if (!is.null(df) && nrow(df) == 1) {
    vals <- unlist(df[1, ], use.names = TRUE)
    names(vals) <- paste0(prefix, "_", names(vals))
    return(vals)
  }

  vec <- unlist(x, use.names = TRUE)
  if (is.null(names(vec))) {
    names(vec) <- paste0(prefix, "_V", seq_along(vec))
  } else {
    names(vec) <- paste0(prefix, "_", names(vec))
  }
  vec
}

append_row_csv <- function(row_df, out_file) {
  if (!file.exists(out_file)) {
    write.csv(row_df, out_file, row.names = FALSE, quote = TRUE)
    return(invisible(TRUE))
  }

  header <- names(read.csv(out_file, nrows = 1, check.names = FALSE))
  row_names <- names(row_df)

  if (identical(header, row_names)) {
    write.table(row_df, out_file, sep = ",", row.names = FALSE,
                col.names = FALSE, append = TRUE, quote = TRUE)
    return(invisible(TRUE))
  }

  old <- read.csv(out_file, check.names = FALSE)
  all_cols <- union(names(old), row_names)

  for (cc in setdiff(all_cols, names(old))) old[[cc]] <- NA
  for (cc in setdiff(all_cols, row_names)) row_df[[cc]] <- NA

  old    <- old[, all_cols, drop = FALSE]
  row_df <- row_df[, all_cols, drop = FALSE]

  write.csv(rbind(old, row_df), out_file, row.names = FALSE, quote = TRUE)
  invisible(TRUE)
}

make_project_path <- function(base_path, province, modelString) {
  file.path(
    base_path,
    province,
    #NEEDS TO BE EDITED IF YOU BROKE THE PATTERN
    #paste0(modelString, "_fitR0"),
    paste0("k1_fitR0_2026"),
    #paste0(modelString, "Fit_R0.mlxtran")
    #paste0(modelString, "Fit_R0_Manitoba_2026.mlxtran")
    paste0("k1_fitR0_Alberta2026.mlxtran")
  )
}

make_result_file <- function(out_dir, modelString, province, data_label, fit_label) {
  file.path(
    out_dir,
    paste0(modelString, "_", data_label, "_", province, "_", fit_label, ".csv")
  )
}

make_grid_file <- function(out_dir, modelString, province, data_label, fit_label) {
  file.path(
    out_dir,
    paste0(modelString, "_", data_label, "_", province, "_", fit_label, "_grid.csv")
  )
}


run_and_save_3d <- function(province,
                            modelString,
                            nu_val,
                            gamma_val,
                            I0_val,
                            nu_index,
                            gamma_index,
                            I0_index,
                            iteration_id,
                            out_file) {

  init <- initial_by_province[[province]]
  
  if (is.null(init)) {
    stop("No initial values defined for province: ", province)
  }

  setPopulationParameterInformation(
    nu_pop       = list(initialValue = as.numeric(nu_val),        method = "FIXED"),
    gamma_pop    = list(initialValue = as.numeric(gamma_val),     method = "FIXED"),
    I0_pop       = list(initialValue = as.numeric(I0_val),        method = "FIXED"),
    ReproNum_pop = list(initialValue = as.numeric(init$ReproNum), method = "MLE"),
    S0_pop       = list(initialValue = as.numeric(init$S0),       method = "MLE"),
    a            = list(initialValue = as.numeric(init$a),        method = "MLE"),
    b            = list(initialValue = as.numeric(init$b),        method = "MLE")
  )

  success <- TRUE
  error_msg <- NA_character_

  tryCatch(
    runScenario(),
    error = function(e) {
      success <<- FALSE
      error_msg <<- conditionMessage(e)
    }
  )

  est <- if (success) tryCatch(getEstimatedPopulationParameters(), error = function(e) e) else NULL
  se  <- if (success) tryCatch(getEstimatedStandardErrors(),       error = function(e) e) else NULL
  ll  <- if (success) tryCatch(getEstimatedLogLikelihood(),        error = function(e) e) else NULL

  meta <- c(
    province    = province,
    model       = modelString,
    iteration   = iteration_id,
    nu_index    = nu_index,
    gamma_index = gamma_index,
    I0_index    = I0_index,
    nu_fixed    = nu_val,
    gamma_fixed = gamma_val,
    I0_fixed    = I0_val,
    latent_period_days     = 1 / nu_val,
    infectious_period_days = 1 / gamma_val,
    success     = success,
    error_msg   = error_msg
  )

  row <- c(
    meta,
    flatten_named(est, "est"),
    flatten_named(se,  "se"),
    flatten_named(ll,  "ll")
  )

  row_df <- as.data.frame(as.list(row), stringsAsFactors = FALSE)
  append_row_csv(row_df, out_file)

  invisible(success)
}

### ========================================================
### 4. MASTER LOOP: province x K-model x gamma x nu x I0
### ========================================================

for (province in province_vec) {
  for (modelString in K_vec) {

    message("\n============================================================")
    message("Starting ", province, " / ", modelString)
    message("============================================================")

    project_path <- path.expand(make_project_path(base_path, province, modelString))

    if (!file.exists(project_path)) {
      warning("Project file not found, skipping: ", project_path)
      next
    }

    out_file <- make_result_file(out_dir, modelString, province, data_label, fit_label)
    grid_file <- make_grid_file(out_dir, modelString, province, data_label, fit_label)

    if (file.exists(out_file)) {
      message("Appending to existing result file: ", out_file)
    } else {
      message("Creating result file: ", out_file)
    }

    grid_df <- expand.grid(
      gamma_index = seq_along(gamma_grid),
      nu_index    = seq_along(nu_grid),
      I0_index    = seq_along(I0_grid)
    )
    grid_df$province <- province
    grid_df$model <- modelString
    grid_df$gamma_fixed <- gamma_grid[grid_df$gamma_index]
    grid_df$nu_fixed    <- nu_grid[grid_df$nu_index]
    grid_df$I0_fixed    <- I0_grid[grid_df$I0_index]
    grid_df$latent_period_days     <- 1 / grid_df$nu_fixed
    grid_df$infectious_period_days <- 1 / grid_df$gamma_fixed

    write.csv(grid_df, grid_file, row.names = FALSE)
    message("Saved grid definition to: ", grid_file)

    loadProject(projectFile = project_path)
    
    setScenario(
      tasks = list(
        populationParameterEstimation = TRUE,
        conditionalDistributionSampling = FALSE,
        conditionalModeEstimation = FALSE,
        standardErrorEstimation = TRUE,
        logLikelihoodEstimation = TRUE,
        plots = FALSE
      ),
      linearization = FALSE
    )
    
    print(getScenario())

    total_iters <- nrow(grid_df)
    pb <- txtProgressBar(min = 0, max = total_iters, style = 3)
    iter <- 0

    for (j in seq_along(gamma_grid)) {
      for (i in seq_along(nu_grid)) {
        for (m in seq_along(I0_grid)) {

          iter <- iter + 1
          setTxtProgressBar(pb, iter)

          run_and_save_3d(
            province     = province,
            modelString  = modelString,
            nu_val       = nu_grid[i],
            gamma_val    = gamma_grid[j],
            I0_val       = I0_grid[m],
            nu_index     = i,
            gamma_index  = j,
            I0_index     = m,
            iteration_id = iter,
            out_file     = out_file
          )
        }
      }
    }

    close(pb)
    cat("\nCompleted ", province, " / ", modelString, "\n", sep = "")
    cat("Results written to: ", out_file, "\n", sep = "")
    cat("Grid definition written to: ", grid_file, "\n", sep = "")
  }
}

cat("\nAll requested scans complete.\n")
cat("All CSV files are in: ", out_dir, "\n", sep = "")

# Notes:
# 1. AIC/BIC/BICc are expected in the flattened output from getEstimatedLogLikelihood(),
#    as in your previous profile_likelihood_results_K1.csv file.
# 2. The RSE columns returned by Monolix are ordered according to the MLE parameters
#    specified in setPopulationParameterInformation(): ReproNum_pop, S0_pop, a, b.
# 3. nu, gamma, and I0 are fixed at each grid point and therefore should appear as
#    fixed estimates in the Monolix output.

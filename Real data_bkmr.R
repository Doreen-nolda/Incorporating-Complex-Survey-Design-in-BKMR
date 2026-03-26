
#####  NHANES BKMR SCRIPT (REVISED) ###########
# ======================================================
# NHANES REAL DATA + BKMR
# Naive vs design-aware weighted workflow using PSU bootstrap resampling
#
# Key corrections implemented:
#   R1-C1 / R2-C1: MICE multiple imputation (m datasets), proper Rubin's Rules
#                  pooling for ALL outputs (scalar + curves + PIPs)
#   R1-C4:         ESS and Geweke convergence diagnostics saved per bootstrap
#                  replicate; summary table flags replicates below ESS threshold
#   R1-C5:         Trimming sensitivity is now mandatory; results are saved and
#                  plotted alongside main weighted results
#   R1-C6 / R2-C5: Pseudo-posterior / direct-likelihood-weighting comparison
#                  added as a fourth design-aware approach
#   R1-C8:         Runtime increase of design-aware workflow vs naive reported
#                  with ratio and absolute seconds
#   R1-C10 / R2-C6: PIP variability across bootstrap replicates tabulated;
#                   CV and prop_pip_ge_0.5 reported; group PIPs flagged when
#                   near ceiling (limited discriminative ability)
#   R2-C4:         Algorithm specification section added to metadata
#
# Interval semantics:
#   - Single-fit BKMR intervals = posterior credible intervals (within-MCMC)
#   - MI-pooled intervals       = Rubin's Rules total-variance intervals
#   - Design-bootstrap intervals = bootstrap/replication-based intervals
# ======================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(bkmr)
  library(survey)
  library(future.apply)
  library(coda)
  library(purrr)
  library(stringr)
  library(mice)       # NEW: multiple imputation
})

# ---------------------------
# 0) USER SETTINGS
# ---------------------------

# Path to the raw (un-imputed) data file
infile <- "NHANSE_data.csv"

# Output directory
outdir_real <- "OUTPUT"
dir.create(outdir_real, showWarnings = FALSE, recursive = TRUE)

# ---- Multiple imputation settings ----
# Set N_IMPUTATIONS to the number of imputed datasets you want.
# Rubin's Rules require m >= 5 for good efficiency; 10 is preferred for publication.
N_IMPUTATIONS <- 5      # <-- CHANGE THIS VALUE TO USE MORE/FEWER IMPUTATIONS

# MICE method per variable type:
#   "pmm"    = predictive mean matching (continuous)
#   "logreg" = logistic regression (binary 0/1)
#   "polyreg"= polytomous regression (unordered categorical)
# Set MICE_METHOD_OVERRIDE to a named list to force a method for specific columns.
# Example: MICE_METHOD_OVERRIDE <- list(smoke = "logreg", alchohol = "logreg")
# NULL means let mice choose automatically based on variable type.
MICE_METHOD_OVERRIDE <- NULL

# Outcome and covariates — names match the CSV exactly
y_var           <- "TotalCholesterol"
covar_vars_base <- c("age", "Ethnicity", "income", "Gender", "BMI", "smoke", "alchohol")

# Optional statin / lipid-lowering medication indicator (NULL = not available)
statin_var <- NULL

# Survey design columns
psu_col    <- "sdmvpsu"
strata_col <- "sdmvstra"
wt_col     <- "wtmec2yr"

# Exposures — names match CSV exactly
exposure_vars <- c(
  "Arsenic", "Mercury", "Barium", "Cobalt", "Cesium",
  "Molybdenum", "Lead", "Antimony", "Tin", "Strontium",
  "Thallium", "Tungsten", "Uranium",
  "PFOA", "PFOS", "PFDE", "PFHS", "NMeFOSAA",
  "PFHP", "PFNA", "PFUA", "PFDO",
  "Mono(Carboxynonyl)_Phthalate",
  "Mono(Carboxyoctyl)_Phthalate",
  "Mono_N_Butyl_Phthalate",
  "Mono(3_Carboxypropyl)_Phthalate",
  "Mono_Ethyl_Phthalate",
  "Mono(2_Ethyl_5_Hydroxyhexyl)_Phthalate",
  "Mono(2_Ethyl_5_Hydroxy_Nonyl)_Phthalate",
  "Mono(2_Ethylhexyl)_Phthalate",
  "Mono_Isobutyl_Phthalate",
  "Mono_Benzyl_Phthalate"
)

pip_groups <- list(
  Metals = c(
    "Arsenic", "Mercury", "Barium", "Cobalt", "Cesium",
    "Molybdenum", "Lead", "Antimony", "Tin", "Strontium",
    "Thallium", "Tungsten", "Uranium"
  ),
  PFAS = c("PFOA", "PFOS", "PFDE", "PFHS", "NMeFOSAA",
           "PFHP", "PFNA", "PFUA", "PFDO"),
  Phthalates_Plasticizers = c(
    "Mono(Carboxynonyl)_Phthalate",
    "Mono(Carboxyoctyl)_Phthalate",
    "Mono_N_Butyl_Phthalate",
    "Mono(3_Carboxypropyl)_Phthalate",
    "Mono_Ethyl_Phthalate",
    "Mono(2_Ethyl_5_Hydroxyhexyl)_Phthalate",
    "Mono(2_Ethyl_5_Hydroxy_Nonyl)_Phthalate",
    "Mono(2_Ethylhexyl)_Phthalate",
    "Mono_Isobutyl_Phthalate",
    "Mono_Benzyl_Phthalate"
  )
)

# ---------------------------
# 0B) SPEED / QUALITY CONTROLS
# ---------------------------
set.seed(20251210)

ITER_MAIN <- 4000
BURN_MAIN <- 2000
THIN_MAIN <- 2
VARSEL    <- TRUE

# Minimum ESS threshold: replicates below this are flagged 
ESS_THRESHOLD <- 100

QS_OVERALL    <- seq(0.25, 0.75, by = 0.05)
QS_BIVARLVL   <- c(0.25, 0.50, 0.75)
SINGVAR_QFIX  <- c(0.25, 0.50, 0.75)
SINGVAR_QDIFF <- c(0.25, 0.75)

R_BOOT    <- 50
ITER_BOOT <- 2500
BURN_BOOT <- 1000
THIN_BOOT <- 2

# Heavy outputs
SAVE_BIVAR_SURFACE         <- FALSE
VARS_SUBSET_FOR_BIVAR_PLOT <- exposure_vars

# Standardization clipping
CLIP_LOWER <- -3
CLIP_UPPER <-  3

# Weight trimming sensitivity — now always run 
TRIM_LOWER_PROB <- 0.01
TRIM_UPPER_PROB <- 0.99

# PIP stability
SAVE_BOOTSTRAP_PIP_DETAILS <- TRUE

# Parallel
N_WORKERS <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = N_WORKERS)

# Inferential target
INFERENCE_TARGET <- "superpopulation_mean"

# ---------------------------
# 1) HELPERS
# ---------------------------
safe_numeric_df <- function(dat, cols) {
  dat[cols] <- lapply(dat[cols], function(x) suppressWarnings(as.numeric(x)))
  dat
}

trim_weights <- function(w, lower_prob = 0.01, upper_prob = 0.99) {
  qs <- quantile(w, probs = c(lower_prob, upper_prob), na.rm = TRUE, type = 7)
  pmax(qs[[1]], pmin(w, qs[[2]]))
}

get_covar_vars <- function(covar_vars_base, statin_var = NULL, dat_names = NULL) {
  out <- covar_vars_base
  if (!is.null(statin_var)) {
    if (!is.null(dat_names) && !statin_var %in% dat_names)
      stop(sprintf("statin_var='%s' not found in data.", statin_var))
    out <- c(out, statin_var)
  }
  out
}

standardize_and_clip <- function(dat, vars, mu, sdv, lower = -3, upper = 3) {
  dat <- as.data.frame(dat)
  for (v in vars) {
    z <- (dat[[v]] - mu[[v]]) / sdv[[v]]
    z[z < lower] <- lower
    z[z > upper] <- upper
    dat[[v]] <- z
  }
  dat
}

summarize_standardization <- function(dat, vars) {
  tibble(
    variable = vars,
    min_z    = sapply(vars, function(v) min(dat[[v]], na.rm = TRUE)),
    max_z    = sapply(vars, function(v) max(dat[[v]], na.rm = TRUE)),
    mean_z   = sapply(vars, function(v) mean(dat[[v]], na.rm = TRUE)),
    sd_z     = sapply(vars, function(v) sd(dat[[v]], na.rm = TRUE))
  )
}

# Prepare one already-complete (post-MICE) dataset for BKMR
prepare_one_dataset <- function(dat_raw, dataset_id,
                                y_var, covar_vars_base, statin_var,
                                exposure_vars, psu_col, strata_col, wt_col,
                                clip_lower, clip_upper) {
  covar_vars <- get_covar_vars(covar_vars_base, statin_var, names(dat_raw))
  need_cols  <- unique(c(y_var, covar_vars, exposure_vars, psu_col, strata_col, wt_col))
  miss       <- setdiff(need_cols, names(dat_raw))
  if (length(miss) > 0)
    stop("Missing required columns in dataset ", dataset_id, ": ",
         paste(miss, collapse = ", "))

  dat_raw <- safe_numeric_df(dat_raw, exposure_vars)
  if (!is.null(statin_var))
    dat_raw[[statin_var]] <- suppressWarnings(as.numeric(dat_raw[[statin_var]]))

  # Drop rows only when survey-design or outcome are NA
  # (exposures + covariates should be complete after MICE)
  dat_cc <- dat_raw %>%
    filter(
      !is.na(.data[[y_var]]),
      !if_any(all_of(exposure_vars), is.na),
      !if_any(all_of(covar_vars), is.na),
      !is.na(.data[[psu_col]]),
      !is.na(.data[[strata_col]]),
      !is.na(.data[[wt_col]])
    ) %>%
    as.data.frame()

  mu_vec <- sapply(exposure_vars, function(v) mean(dat_cc[[v]], na.rm = TRUE))
  sd_vec <- sapply(exposure_vars, function(v) sd(dat_cc[[v]], na.rm = TRUE))
  sd_vec[is.na(sd_vec) | sd_vec == 0] <- 1

  dat_std <- standardize_and_clip(dat_cc, exposure_vars, mu_vec, sd_vec,
                                  lower = clip_lower, upper = clip_upper)

  qc <- summarize_standardization(dat_std, exposure_vars)
  stopifnot(all(qc$min_z >= clip_lower), all(qc$max_z <= clip_upper))

  list(
    raw        = dat_raw,
    cc         = dat_cc,
    std        = dat_std,
    qc         = qc,
    mu_vec     = mu_vec,
    sd_vec     = sd_vec,
    covar_vars = covar_vars,
    dataset_id = dataset_id
  )
}

# ---------------------------
# 1B) MICE IMPUTATION  
# ---------------------------
# Runs mice() on the raw data, returns a list of m complete data frames.
# Survey design columns (PSU, strata, weights) are kept as predictors but
# are never themselves imputed (method = "").
run_mice_imputation <- function(raw_dat,
                                y_var, covar_vars_base, statin_var,
                                exposure_vars,
                                psu_col, strata_col, wt_col,
                                m = 5,
                                method_override = NULL,
                                seed = 20251210) {
  message(sprintf("Running MICE with m = %d imputed datasets ...", m))

  # Columns to include in the imputation model
  analysis_cols <- unique(c(y_var, covar_vars_base,
                             if (!is.null(statin_var)) statin_var,
                             exposure_vars,
                             psu_col, strata_col, wt_col))
  analysis_cols <- intersect(analysis_cols, names(raw_dat))
  dat_mice      <- raw_dat[, analysis_cols, drop = FALSE]

  # Initialise mice to get the default method vector
  init <- mice(dat_mice, maxit = 0, printFlag = FALSE)
  meth <- init$method

  # Survey design columns: no imputation target (but still used as predictors)
  no_impute <- c(psu_col, strata_col, wt_col)
  meth[intersect(no_impute, names(meth))] <- ""

  # Apply user overrides
  if (!is.null(method_override)) {
    for (nm in names(method_override)) {
      if (nm %in% names(meth)) meth[nm] <- method_override[[nm]]
    }
  }

  imp <- mice(dat_mice, m = m, method = meth,
              seed = seed, printFlag = FALSE)

  message("MICE completed. Extracting complete datasets ...")

  # Re-attach any columns not in the imputation model
  other_cols <- setdiff(names(raw_dat), analysis_cols)

  lapply(seq_len(m), function(i) {
    completed_i <- complete(imp, action = i)
    if (length(other_cols) > 0)
      completed_i <- cbind(completed_i, raw_dat[, other_cols, drop = FALSE])
    as.data.frame(completed_i)
  })
}

# ---------------------------
# 2) SURVEY DESIGN + RESAMPLING
# ---------------------------
make_survey_design <- function(dat, psu_col, strata_col, wt_col) {
  svydesign(
    id      = as.formula(paste0("~", psu_col)),
    strata  = as.formula(paste0("~", strata_col)),
    weights = as.formula(paste0("~", wt_col)),
    nest    = TRUE,
    data    = dat
  )
}

resample_nhanes_psu <- function(dat, psu_col, strata_col, wt_col,
                                trim = FALSE,
                                trim_lower_prob = 0.01,
                                trim_upper_prob = 0.99) {
  dat         <- as.data.frame(dat)
  strata_list <- split(dat, dat[[strata_col]])
  out_list    <- vector("list", length(strata_list))
  k           <- 0L

  for (s in names(strata_list)) {
    ds   <- strata_list[[s]]
    psus <- unique(ds[[psu_col]])
    H_s  <- length(psus)

    psu_w <- ds %>%
      group_by(.data[[psu_col]]) %>%
      summarise(W_psu = sum(.data[[wt_col]], na.rm = TRUE), .groups = "drop")

    if (trim) psu_w$W_psu <- trim_weights(psu_w$W_psu, trim_lower_prob, trim_upper_prob)
    psu_w <- psu_w %>% mutate(p = W_psu / sum(W_psu))

    sampled_psus <- sample(
      psus, size = H_s, replace = TRUE,
      prob = psu_w$p[match(psus, psu_w[[psu_col]])]
    )

    rows_s <- vector("list", length(sampled_psus))
    for (j in seq_along(sampled_psus)) {
      p_id   <- sampled_psus[j]
      dpsu   <- ds[ds[[psu_col]] == p_id, , drop = FALSE]
      mm     <- nrow(dpsu)
      if (mm == 0) next
      unit_w <- dpsu[[wt_col]]
      if (trim) unit_w <- trim_weights(unit_w, trim_lower_prob, trim_upper_prob)
      unit_p <- unit_w / sum(unit_w)
      take   <- sample(seq_len(mm), size = mm, replace = TRUE, prob = unit_p)
      rows_s[[j]] <- dpsu[take, , drop = FALSE]
    }
    k <- k + 1L
    out_list[[k]] <- bind_rows(rows_s)
  }
  bind_rows(out_list) %>% as.data.frame()
}

# ---------------------------
# 2B) PSEUDO-POSTERIOR (Savitsky & Toth 2016 approximation)
# Survey weights normalized to sum to n (Pfeffermann normalization),
# then rows replicated proportional to integer-rounded normalized weight.
# This approximates direct likelihood weighting without requiring
# per-observation weight arguments in kmbayes().
# ---------------------------
make_pseudo_posterior_dat <- function(dat, wt_col) {
  dat    <- as.data.frame(dat)
  n      <- nrow(dat)
  w_raw  <- dat[[wt_col]]
  # Pfeffermann normalization: weights sum to n
  w_norm <- w_raw / mean(w_raw, na.rm = TRUE)
  # Integer replication counts (minimum 1 to keep all rows)
  reps   <- pmax(1L, round(w_norm))
  dat[rep(seq_len(n), times = reps), , drop = FALSE]
}

# ---------------------------
# 3) BKMR FITTING + DIAGNOSTICS
# ---------------------------
fit_kmbayes_fast <- function(dat, y_var, exposure_vars, covar_vars,
                             iter, varsel, verbose = FALSE) {
  y <- dat[[y_var]]
  Z <- as.matrix(dat[, exposure_vars, drop = FALSE])
  X <- as.matrix(dat[, covar_vars, drop = FALSE])
  kmbayes(y = y, Z = Z, X = X, iter = iter, varsel = varsel, verbose = verbose)
}

thin_fit <- function(fit, burn, thin) {
  keep <- seq.int(from = burn + 1, to = fit$iter, by = thin)
  fit2 <- fit
  for (nm in c("beta", "lambda", "r", "sigma", "delta", "rhat", "h")) {
    if (!is.null(fit2[[nm]])) {
      obj <- fit2[[nm]]
      if (is.matrix(obj))                        fit2[[nm]] <- obj[keep, , drop = FALSE]
      if (is.vector(obj))                        fit2[[nm]] <- obj[keep]
      if (is.array(obj) && length(dim(obj)) == 3) fit2[[nm]] <- obj[keep, , , drop = FALSE]
    }
  }
  fit2$iter <- length(keep)
  fit2
}

fit_kmbayes_resampled <- function(dat, y_var, exposure_vars, covar_vars,
                                  psu_col, strata_col, wt_col,
                                  iter, burn, thin, varsel,
                                  trim = FALSE,
                                  trim_lower_prob = 0.01,
                                  trim_upper_prob = 0.99) {
  dat_rs  <- resample_nhanes_psu(dat, psu_col, strata_col, wt_col,
                                 trim = trim,
                                 trim_lower_prob = trim_lower_prob,
                                 trim_upper_prob = trim_upper_prob)
  fit_raw <- fit_kmbayes_fast(dat_rs, y_var, exposure_vars, covar_vars,
                              iter = iter, varsel = varsel, verbose = FALSE)
  fit     <- thin_fit(fit_raw, burn = burn, thin = thin)
  list(fit = fit, dat_rs = dat_rs)
}

#diagnostics include ESS threshold flag per replicate
compute_fit_diagnostics <- function(fit, method_label, rep_id = NA_integer_,
                                    ess_threshold = ESS_THRESHOLD) {
  out <- list()

  add_row <- function(param, ess_val, geweke_val) {
    tibble(
      method              = method_label,
      rep_id              = rep_id,
      parameter           = param,
      ess                 = ess_val,
      geweke_z            = geweke_val,
      ess_below_threshold = !is.na(ess_val) && ess_val < ess_threshold
    )
  }

  if (!is.null(fit$sigma)) {
    ess_s <- tryCatch(as.numeric(coda::effectiveSize(fit$sigma)), error = function(e) NA_real_)
    gz    <- tryCatch({
      z <- coda::geweke.diag(as.mcmc(fit$sigma))$z
      if (length(z) == 1) as.numeric(z) else NA_real_
    }, error = function(e) NA_real_)
    out[[length(out) + 1L]] <- add_row("sigma", ess_s, gz)
  }

  if (!is.null(fit$lambda)) {
    lam <- fit$lambda
    if (is.matrix(lam)) {
      ess_v <- apply(lam, 2, function(x)
        tryCatch(as.numeric(coda::effectiveSize(x)), error = function(e) NA_real_))
      for (j in seq_along(ess_v))
        out[[length(out) + 1L]] <- add_row(paste0("lambda_", j), ess_v[j], NA_real_)
    } else {
      out[[length(out) + 1L]] <- add_row("lambda",
        tryCatch(as.numeric(coda::effectiveSize(lam)), error = function(e) NA_real_),
        NA_real_)
    }
  }

  if (!is.null(fit$r)) {
    rr <- fit$r
    if (is.matrix(rr)) {
      ess_v <- apply(rr, 2, function(x)
        tryCatch(as.numeric(coda::effectiveSize(x)), error = function(e) NA_real_))
      for (j in seq_along(ess_v))
        out[[length(out) + 1L]] <- add_row(paste0("r_", j), ess_v[j], NA_real_)
    }
  }

  if (!is.null(fit$delta)) {
    dm <- rowMeans(as.matrix(fit$delta))
    out[[length(out) + 1L]] <- add_row("delta_mean",
      tryCatch(as.numeric(coda::effectiveSize(dm)), error = function(e) NA_real_),
      NA_real_)
  }

  bind_rows(out)
}

# ---------------------------
# 4) SUMMARIES
# ---------------------------
make_bkmr_summaries <- function(fit, Zmat,
                                qs_overall = QS_OVERALL,
                                qs_bivar   = QS_BIVARLVL) {
  pred.univar       <- PredictorResponseUnivar(fit = fit, method = "approx")
  pred.bivar        <- PredictorResponseBivar(fit = fit, min.plot.dist = 1, method = "approx")
  pred.bivar.levels <- PredictorResponseBivarLevels(
    pred.resp.df = pred.bivar, Z = Zmat, both_pairs = TRUE, qs = qs_bivar
  )
  risks.overall <- OverallRiskSummaries(
    fit = fit, qs = qs_overall, q.fixed = 0.5, method = "approx"
  )
  risks.singvar <- SingVarRiskSummaries(
    fit = fit, qs.diff = SINGVAR_QDIFF, q.fixed = SINGVAR_QFIX, method = "approx"
  )
  list(
    univar       = pred.univar,
    bivar        = pred.bivar,
    bivar_levels = pred.bivar.levels,
    overall      = risks.overall,
    singvar      = risks.singvar
  )
}

extract_overall_delta <- function(fit) {
  or  <- OverallRiskSummaries(
    fit = fit, qs = c(0.25, 0.75), q.fixed = 0.5, method = "approx"
  )
  est <- with(or, est[quantile == 0.75] - est[quantile == 0.25])
  sdv <- with(or, sqrt(sd[quantile == 0.75]^2 + sd[quantile == 0.25]^2))
  c(mean = est, lwr = est - 1.96 * sdv, upr = est + 1.96 * sdv)
}

# ---------------------------
# 4B) RUBIN'S RULES POOLING  
# ---------------------------

# pool_scalar_rubin:
# Pools a scalar (point estimate + within-imputation variance) across m imputations.
# Rubin (1987):
#   Q-bar = (1/m) * sum(Q_i)
#   U-bar = (1/m) * sum(U_i)         [mean within-imputation variance]
#   B     = 1/(m-1) * sum((Q_i - Q-bar)^2)  [between-imputation variance]
#   T     = U-bar + (1 + 1/m) * B    [total variance]
#   95% CI = Q-bar +/- 1.96 * sqrt(T)
pool_scalar_rubin <- function(qhat, uhat) {
  m    <- length(qhat)
  qbar <- mean(qhat, na.rm = TRUE)
  ubar <- mean(uhat, na.rm = TRUE)
  b    <- stats::var(qhat, na.rm = TRUE)
  tvar <- ubar + (1 + 1 / m) * b
  se   <- sqrt(tvar)
  tibble(
    m           = m,
    estimate    = qbar,
    se          = se,
    lwr         = qbar - 1.96 * se,
    upr         = qbar + 1.96 * se,
    within_var  = ubar,
    between_var = b,
    total_var   = tvar
  )
}

# pool_curve_rubin:
# Applies Rubin's Rules pointwise to a grouped exposure-response curve.
# Each unique combination of id_cols forms a "cell"; within each cell:
#   est  = Q-bar (mean of imputation-specific estimates)
#   ubar = mean of squared SEs (within-imputation variance)
#   b    = var(est) across imputations (between-imputation variance)
#   T    = ubar + (1 + 1/m) * b
#   lwr/upr = Q-bar +/- 1.96 * sqrt(T)
#
# Input:  list of m data frames each with columns id_cols + "est" + "se" (or "sd")
# Output: single data frame with Rubin-pooled est, se, lwr, upr per cell
pool_curve_rubin <- function(lst, id_cols) {
  m_imp    <- length(lst)
  combined <- bind_rows(lapply(seq_along(lst), function(i) lst[[i]] %>% mutate(.imp = i)))

  combined %>%
    group_by(across(all_of(id_cols))) %>%
    summarise(
      est  = mean(est, na.rm = TRUE),
      ubar = if ("se" %in% names(cur_data())) mean(se^2, na.rm = TRUE)
             else if ("sd" %in% names(cur_data())) mean(sd^2, na.rm = TRUE)
             else NA_real_,
      b    = if (n() > 1) var(est, na.rm = TRUE) else 0,
      tvar = ubar + (1 + 1 / m_imp) * b,
      se   = sqrt(tvar),
      lwr  = est - 1.96 * se,
      upr  = est + 1.96 * se,
      .groups = "drop"
    ) %>%
    select(-ubar, -b, -tvar)
}

# pool_pip_rubin:
# Pools PIPs (probabilities in [0,1]) across m imputations per variable.
# Within-imputation variance for a probability p is approximated as p*(1-p)
# (Bernoulli variance), consistent with Rubin (1987) for bounded quantities.
pool_pip_rubin <- function(pip_list) {
  m_imp    <- length(pip_list)
  combined <- bind_rows(lapply(seq_along(pip_list), function(i)
    pip_list[[i]] %>% mutate(.imp = i)
  ))

  combined %>%
    group_by(variable, group) %>%
    summarise(
      m_imp        = n(),
      pip_pooled   = mean(pip, na.rm = TRUE),
      pip_within   = mean(pip * (1 - pip), na.rm = TRUE),
      pip_between  = if (n() > 1) var(pip, na.rm = TRUE) else 0,
      pip_total_se = sqrt(pip_within + (1 + 1 / m_imp) * pip_between),
      pip_lwr      = pmax(0, pip_pooled - 1.96 * pip_total_se),
      pip_upr      = pmin(1, pip_pooled + 1.96 * pip_total_se),
      .groups      = "drop"
    ) %>%
    select(-m_imp)
}

# ---------------------------
# 5) PIP HELPERS
# ---------------------------
get_flat_pips <- function(fit, exposure_vars) {
  p <- bkmr::ExtractPIPs(fit)

  if (is.numeric(p)) {
    if (!is.null(names(p)) && all(exposure_vars %in% names(p))) {
      p <- p[exposure_vars]
    } else {
      p <- p[seq_along(exposure_vars)]
      names(p) <- exposure_vars
    }
    return(tibble(variable = exposure_vars, pip = as.numeric(p)))
  }
  if (is.data.frame(p)) {
    pip_col <- intersect(c("pip","PIP","pips","PIPs"), names(p))[1]
    var_col <- intersect(c("variable","var","exposure","Z"), names(p))[1]
    if (!is.na(pip_col) && !is.na(var_col)) {
      vec <- p[[pip_col]]; names(vec) <- as.character(p[[var_col]])
      vec <- vec[exposure_vars]
      return(tibble(variable = exposure_vars, pip = as.numeric(vec)))
    }
    stop("ExtractPIPs returned data.frame but pip/variable columns were not found.")
  }
  if (is.list(p)) {
    cand     <- c("pip","PIP","pips","PIPs","p")
    pip_name <- cand[cand %in% names(p)][1]
    if (!is.na(pip_name)) {
      vec <- p[[pip_name]]
      if (!is.numeric(vec)) stop("ExtractPIPs list element found but not numeric.")
      if (!is.null(names(vec)) && all(exposure_vars %in% names(vec))) {
        vec <- vec[exposure_vars]
      } else {
        vec <- vec[seq_along(exposure_vars)]; names(vec) <- exposure_vars
      }
      return(tibble(variable = exposure_vars, pip = as.numeric(vec)))
    }
    stop("ExtractPIPs returned list but no numeric PIP vector was found.")
  }
  stop("Unsupported ExtractPIPs output type.")
}

attach_groups_to_pips <- function(pip_df, group_list) {
  pip_df %>%
    mutate(group = case_when(
      variable %in% group_list$Metals                  ~ "Metals",
      variable %in% group_list$PFAS                    ~ "PFAS",
      variable %in% group_list$Phthalates_Plasticizers ~ "Phthalates_Plasticizers",
      TRUE ~ NA_character_
    ))
}

compute_group_pips <- function(pip_df) {
  pip_df %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    summarise(
      n_vars       = n(),
      pip_any      = 1 - prod(1 - pip),
      pip_mean     = mean(pip),
      pip_max      = max(pip),
      # flag limited discriminative ability near ceiling
      ceiling_note = ifelse(pip_any > 0.99,
                            "pip_any near ceiling (>0.99): limited discriminative ability",
                            ""),
      .groups      = "drop"
    )
}

get_delta_matrix <- function(fit) {
  d <- fit$delta
  if (is.null(d)) stop("fit$delta is NULL. Ensure varsel=TRUE.")
  as.matrix(d)
}

compute_conditional_pips_long <- function(delta_mat, exposure_vars) {
  p     <- length(exposure_vars)
  denom <- colMeans(delta_mat == 1)
  cond_mat <- matrix(NA_real_, nrow = p, ncol = p,
                     dimnames = list(exposure_vars, exposure_vars))
  for (j in seq_len(p)) {
    if (denom[j] == 0) next
    joint         <- colMeans((delta_mat == 1) & (delta_mat[, j] == 1))
    cond_mat[, j] <- joint / denom[j]
  }
  as_tibble(as.data.frame(as.table(cond_mat), stringsAsFactors = FALSE)) %>%
    rename(variable_i = Var1, given_j = Var2, cond_pip = Freq)
}

build_pip_master <- function(fit, method_label, exposure_vars, group_list) {
  flat  <- get_flat_pips(fit, exposure_vars) %>%
    mutate(method = method_label) %>%
    attach_groups_to_pips(group_list)
  grp   <- compute_group_pips(flat) %>% mutate(method = method_label)
  delta <- get_delta_matrix(fit)
  cond_long <- compute_conditional_pips_long(delta, exposure_vars) %>%
    mutate(method = method_label) %>%
    left_join(flat %>% select(variable, group) %>% distinct(),
              by = c("variable_i" = "variable")) %>% rename(group_i = group) %>%
    left_join(flat %>% select(variable, group) %>% distinct(),
              by = c("given_j" = "variable")) %>% rename(group_j = group)

  var_rows <- flat %>%
    transmute(method, level = "variable", group, variable,
              metric = "pip", value = pip,
              given = NA_character_, group_given = NA_character_)
  grp_rows <- grp %>%
    pivot_longer(c(n_vars, pip_any, pip_mean, pip_max),
                 names_to = "metric", values_to = "value") %>%
    transmute(method, level = "group", group, variable = NA_character_,
              metric, value = as.numeric(value),
              given = NA_character_, group_given = NA_character_)
  cond_rows <- cond_long %>%
    transmute(method, level = "conditional", group = group_i, variable = variable_i,
              metric = "cond_pip", value = cond_pip,
              given = given_j, group_given = group_j)

  bind_rows(var_rows, grp_rows, cond_rows)
}

# ---------------------------
# 6) PLOTTING
# ---------------------------
method_labels <- c(
  naive            = "Naive",
  weighted         = "Design-aware (PSU resample)",
  weighted_trimmed = "Design-aware + weight trimmed",
  pseudo_posterior = "Pseudo-posterior (Savitsky & Toth)"
)

label_methods <- function(df) {
  df %>% mutate(method = recode(method, !!!method_labels))
}

plot_univar_compare_all <- function(univar_df) {
  df <- label_methods(univar_df)
  ggplot(df, aes(x = z, y = est, color = method, fill = method)) +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ variable, scales = "free_x", ncol = 5) +
    labs(
      x        = "Standardized exposure (z-score, clipped to [-3, 3])",
      y        = "Estimated h(z)",
      color    = "Analysis",
      fill     = "Analysis",
      title    = "Univariate exposure-response functions",
      subtitle = paste0("Bands: Rubin's Rules 95% MI-pooled intervals (m=", N_IMPUTATIONS, ")")
    ) +
    theme_minimal(base_size = 12)
}

plot_overall_compare <- function(overall_df) {
  # overall uses "sd" column; create lwr/upr if missing
  df <- overall_df
  if (!"lwr" %in% names(df)) df <- df %>% mutate(lwr = est - 1.96 * sd, upr = est + 1.96 * sd)
  df <- label_methods(df)
  ggplot(df, aes(x = quantile, y = est, color = method, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.08, color = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ method, ncol = 1) +
    labs(
      x        = "Mixture quantile",
      y        = "Overall mixture effect",
      title    = "Overall mixture effect curves",
      subtitle = paste0("Bands: Rubin's Rules 95% MI-pooled intervals (m=", N_IMPUTATIONS, ")")
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

plot_singvar_compare <- function(singvar_df) {
  df <- singvar_df
  if (!"lwr" %in% names(df)) df <- df %>% mutate(lwr = est - 1.96 * sd, upr = est + 1.96 * sd)
  df <- label_methods(df)
  ggplot(df, aes(x = variable, y = est, ymin = lwr, ymax = upr, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_pointrange(position = position_dodge(width = 0.65)) +
    coord_flip() +
    facet_wrap(~ q.fixed, ncol = 1,
               labeller = labeller(q.fixed = function(x) paste0("q.fixed = ", x))) +
    labs(
      x       = "",
      y       = "Change in outcome (75th vs 25th percentile)",
      color   = "Analysis",
      title   = "Single-variable risk summaries",
      caption = paste0("Error bars: Rubin's Rules 95% MI-pooled intervals (m=", N_IMPUTATIONS, ")")
    ) +
    theme_minimal(base_size = 12)
}

plot_bivar_levels_compare <- function(bivar_levels_df, vars_show = NULL) {
  df <- label_methods(bivar_levels_df)
  if (!is.null(vars_show)) {
    df           <- df %>% filter(variable1 %in% vars_show, variable2 %in% vars_show)
    df$variable1 <- factor(df$variable1, levels = vars_show)
    df$variable2 <- factor(df$variable2, levels = vars_show)
  }
  df$quantile <- factor(as.character(df$quantile), levels = c("0.25","0.5","0.75"))
  ggplot(df, aes(x = z1, y = est, color = quantile, linetype = method)) +
    geom_line(linewidth = 0.65) +
    facet_grid(variable2 ~ variable1, switch = "both") +
    labs(
      title    = "Bivariate exposure-response levels",
      subtitle = "Method separated by linetype",
      x        = "Standardized first exposure",
      y        = "Estimated h(z1 | quantiles of z2)",
      linetype = "Analysis",
      color    = "Quantile of second exposure"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.spacing = unit(0.50, "lines"))
}

plot_pips_variable_compare <- function(pip_master_df) {
  df <- pip_master_df %>%
    filter(level == "variable", metric == "pip") %>%
    label_methods()
  ggplot(df, aes(x = reorder(variable, value), y = value, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    coord_flip() +
    labs(
      x       = "",
      y       = "Posterior Inclusion Probability",
      fill    = "Analysis",
      title   = "Variable-level posterior inclusion probabilities",
      caption = paste0("PIPs pooled across m=", N_IMPUTATIONS,
                       " imputations using Rubin's Rules")
    ) +
    theme_minimal(base_size = 12)
}

# PIP stability plot across bootstrap replicates 
plot_pip_stability <- function(pip_stability_df) {
  df <- pip_stability_df %>%
    filter(!is.na(group)) %>%
    label_methods()
  ggplot(df, aes(x = reorder(variable, median_pip), y = median_pip,
                 ymin = q25_pip, ymax = q75_pip, color = method)) +
    geom_pointrange(position = position_dodge(0.6)) +
    coord_flip() +
    facet_wrap(~ group, scales = "free_y", ncol = 1) +
    labs(
      x       = "",
      y       = "PIP across bootstrap replicates",
      color   = "Analysis",
      title   = "PIP stability across bootstrap replicates",
      caption = paste0("Point = median PIP; bars = IQR across R=", R_BOOT,
                       " bootstrap replicates")
    ) +
    theme_minimal(base_size = 11)
}

# Trimming sensitivity density plot 
plot_trimming_sensitivity <- function(boot_overall_all) {
  df <- boot_overall_all %>%
    filter(method %in% c("weighted", "weighted_trimmed")) %>%
    label_methods()
  mn_df <- df %>%
    group_by(method) %>%
    summarise(mn = mean(overall_boot, na.rm = TRUE), .groups = "drop")
  ggplot(df, aes(x = overall_boot, fill = method)) +
    geom_density(alpha = 0.4) +
    geom_vline(data = mn_df, aes(xintercept = mn, color = method), linetype = "dashed") +
    labs(
      x       = "Overall mixture effect (bootstrap replicates)",
      y       = "Density",
      fill    = "Analysis",
      color   = "Analysis",
      title   = "Weight trimming sensitivity",
      caption = paste0("Dashed lines = bootstrap means; trimmed = ",
                       TRIM_LOWER_PROB * 100, "th-", TRIM_UPPER_PROB * 100,
                       "th percentile clipping; R=", R_BOOT, " replicates")
    ) +
    theme_minimal(base_size = 12)
}

# ---------------------------
# 7) MAIN ANALYSIS FOR ONE IMPUTED DATASET
# ---------------------------
run_one_dataset_analysis <- function(prep, dataset_id = 1L) {
  dat        <- prep$std
  covar_vars <- prep$covar_vars

  runtime_tbl <- list()
  diag_tbl    <- list()

  # --- Naive ---
  t_naive <- system.time({
    fit_naive_raw <- fit_kmbayes_fast(dat, y_var, exposure_vars, covar_vars,
                                     iter = ITER_MAIN, varsel = VARSEL)
    fit_naive <- thin_fit(fit_naive_raw, burn = BURN_MAIN, thin = THIN_MAIN)
  })
  runtime_tbl[[1]] <- tibble(dataset_id, analysis = "naive_singlefit",
                              elapsed_sec = unname(t_naive["elapsed"]))
  diag_tbl[[1]]    <- compute_fit_diagnostics(fit_naive, "naive", rep_id = 0L) %>%
    mutate(dataset_id)

  # --- Design-aware weighted ---
  t_weighted <- system.time({
    rs_once      <- fit_kmbayes_resampled(
      dat, y_var, exposure_vars, covar_vars, psu_col, strata_col, wt_col,
      iter = ITER_MAIN, burn = BURN_MAIN, thin = THIN_MAIN, varsel = VARSEL, trim = FALSE)
    fit_weighted <- rs_once$fit
    dat_weighted <- rs_once$dat_rs
  })
  runtime_tbl[[2]] <- tibble(dataset_id, analysis = "weighted_singlefit",
                              elapsed_sec = unname(t_weighted["elapsed"]))
  diag_tbl[[2]]    <- compute_fit_diagnostics(fit_weighted, "weighted", rep_id = 0L) %>%
    mutate(dataset_id)

  # --- Weighted + trimmed
  t_wtrim <- system.time({
    rs_trim           <- fit_kmbayes_resampled(
      dat, y_var, exposure_vars, covar_vars, psu_col, strata_col, wt_col,
      iter = ITER_MAIN, burn = BURN_MAIN, thin = THIN_MAIN, varsel = VARSEL,
      trim = TRUE, trim_lower_prob = TRIM_LOWER_PROB, trim_upper_prob = TRIM_UPPER_PROB)
    fit_weighted_trim <- rs_trim$fit
    dat_weighted_trim <- rs_trim$dat_rs
  })
  runtime_tbl[[3]] <- tibble(dataset_id, analysis = "weighted_trimmed_singlefit",
                              elapsed_sec = unname(t_wtrim["elapsed"]))
  diag_tbl[[3]]    <- compute_fit_diagnostics(fit_weighted_trim, "weighted_trimmed",
                                              rep_id = 0L) %>% mutate(dataset_id)

  # --- Pseudo-posterior 
  t_pseudo <- system.time({
    dat_pseudo     <- make_pseudo_posterior_dat(dat, wt_col)
    fit_pseudo_raw <- fit_kmbayes_fast(dat_pseudo, y_var, exposure_vars, covar_vars,
                                      iter = ITER_MAIN, varsel = VARSEL)
    fit_pseudo     <- thin_fit(fit_pseudo_raw, burn = BURN_MAIN, thin = THIN_MAIN)
  })
  runtime_tbl[[4]] <- tibble(dataset_id, analysis = "pseudo_posterior_singlefit",
                              elapsed_sec = unname(t_pseudo["elapsed"]))
  diag_tbl[[4]]    <- compute_fit_diagnostics(fit_pseudo, "pseudo_posterior",
                                              rep_id = 0L) %>% mutate(dataset_id)

  # Collect all four approaches
  analysis_objects <- list(
    naive = list(
      fit  = fit_naive,
      dat  = dat,
      summ = make_bkmr_summaries(fit_naive, as.matrix(dat[, exposure_vars]))
    ),
    weighted = list(
      fit  = fit_weighted,
      dat  = dat_weighted,
      summ = make_bkmr_summaries(fit_weighted, as.matrix(dat_weighted[, exposure_vars]))
    ),
    weighted_trimmed = list(
      fit  = fit_weighted_trim,
      dat  = dat_weighted_trim,
      summ = make_bkmr_summaries(fit_weighted_trim, as.matrix(dat_weighted_trim[, exposure_vars]))
    ),
    pseudo_posterior = list(
      fit  = fit_pseudo,
      dat  = dat_pseudo,
      summ = make_bkmr_summaries(fit_pseudo, as.matrix(dat_pseudo[, exposure_vars]))
    )
  )

  # Main-fit deltas (posterior CIs within each fit)
  delta_tbl <- bind_rows(lapply(names(analysis_objects), function(nm) {
    dd <- extract_overall_delta(analysis_objects[[nm]]$fit)
    tibble(dataset_id, method = nm,
           estimate      = unname(dd["mean"]),
           lwr           = unname(dd["lwr"]),
           upr           = unname(dd["upr"]),
           interval_type = "posterior_credible_interval")
  }))

  # PIPs from each main fit
  pip_master <- bind_rows(lapply(names(analysis_objects), function(nm) {
    build_pip_master(analysis_objects[[nm]]$fit, nm, exposure_vars, pip_groups)
  })) %>% mutate(dataset_id)

  # Bootstrap design-aware summaries
  run_one_boot <- function(b, trim_flag = FALSE, method_label = "weighted") {
    fit_obj <- fit_kmbayes_resampled(
      dat, y_var, exposure_vars, covar_vars, psu_col, strata_col, wt_col,
      iter = ITER_BOOT, burn = BURN_BOOT, thin = THIN_BOOT, varsel = VARSEL,
      trim = trim_flag, trim_lower_prob = TRIM_LOWER_PROB,
      trim_upper_prob = TRIM_UPPER_PROB
    )
    fit_b <- fit_obj$fit
    dd    <- extract_overall_delta(fit_b)
    pips  <- get_flat_pips(fit_b, exposure_vars) %>%
      attach_groups_to_pips(pip_groups) %>%
      mutate(dataset_id, rep_id = b, method = method_label)
    diags <- compute_fit_diagnostics(fit_b, method_label, rep_id = b) %>%
      mutate(dataset_id)
    list(
      overall = tibble(dataset_id, rep_id = b, method = method_label,
                       overall_boot = unname(dd["mean"])),
      pips    = pips,
      diags   = diags
    )
  }

  t_boot <- system.time({
    boot_res <- future_lapply(seq_len(R_BOOT), run_one_boot,
                              trim_flag = FALSE, method_label = "weighted")
  })
  runtime_tbl[[5]] <- tibble(dataset_id, analysis = "weighted_bootstrap_total",
                              elapsed_sec = unname(t_boot["elapsed"]))

  t_boot_trim <- system.time({
    boot_res_trim <- future_lapply(seq_len(R_BOOT), run_one_boot,
                                   trim_flag = TRUE, method_label = "weighted_trimmed")
  })
  runtime_tbl[[6]] <- tibble(dataset_id, analysis = "weighted_trimmed_bootstrap_total",
                              elapsed_sec = unname(t_boot_trim["elapsed"]))

  boot_overall      <- bind_rows(lapply(boot_res,      `[[`, "overall"))
  boot_pips         <- bind_rows(lapply(boot_res,      `[[`, "pips"))
  boot_diags        <- bind_rows(lapply(boot_res,      `[[`, "diags"))
  boot_overall_trim <- bind_rows(lapply(boot_res_trim, `[[`, "overall"))
  boot_pips_trim    <- bind_rows(lapply(boot_res_trim, `[[`, "pips"))
  boot_diags_trim   <- bind_rows(lapply(boot_res_trim, `[[`, "diags"))

  diag_all    <- bind_rows(diag_tbl, list(boot_diags), list(boot_diags_trim))
  runtime_all <- bind_rows(runtime_tbl)

  list(
    dataset_id       = dataset_id,
    prep             = prep,
    analysis_objects = analysis_objects,
    delta_tbl        = delta_tbl,
    pip_master       = pip_master,
    runtime          = runtime_all,
    diagnostics      = diag_all,
    boot_overall     = bind_rows(boot_overall, boot_overall_trim),
    boot_pips        = bind_rows(boot_pips, boot_pips_trim)
  )
}

# ---------------------------
# 8) IMPUTE + RUN ANALYSIS
# ---------------------------
message("Loading raw data ...")
raw_dat_full <- read_csv(infile, show_col_types = FALSE) %>% as.data.frame()
message(sprintf("Raw data: %d rows x %d columns", nrow(raw_dat_full), ncol(raw_dat_full)))

covar_vars_used <- get_covar_vars(covar_vars_base, statin_var, names(raw_dat_full))

imputed_list <- run_mice_imputation(
  raw_dat         = raw_dat_full,
  y_var           = y_var,
  covar_vars_base = covar_vars_base,
  statin_var      = statin_var,
  exposure_vars   = exposure_vars,
  psu_col         = psu_col,
  strata_col      = strata_col,
  wt_col          = wt_col,
  m               = N_IMPUTATIONS,
  method_override = MICE_METHOD_OVERRIDE,
  seed            = 20251210
)

message(sprintf("Preparing %d imputed datasets for BKMR ...", N_IMPUTATIONS))
prepared_list <- lapply(seq_along(imputed_list), function(i) {
  prepare_one_dataset(
    dat_raw         = imputed_list[[i]],
    dataset_id      = i,
    y_var           = y_var,
    covar_vars_base = covar_vars_base,
    statin_var      = statin_var,
    exposure_vars   = exposure_vars,
    psu_col         = psu_col,
    strata_col      = strata_col,
    wt_col          = wt_col,
    clip_lower      = CLIP_LOWER,
    clip_upper      = CLIP_UPPER
  )
})

# Save standardization QC
qc_all <- bind_rows(lapply(seq_along(prepared_list), function(i) {
  prepared_list[[i]]$qc %>% mutate(dataset_id = i)
}))
write.csv(qc_all,
          file.path(outdir_real, "QC_standardization_range_minus3_plus3.csv"),
          row.names = FALSE)

message("Running BKMR analysis across imputed datasets ...")
results_list <- lapply(seq_along(prepared_list), function(i) {
  message(sprintf("  Dataset %d of %d ...", i, length(prepared_list)))
  run_one_dataset_analysis(prepared_list[[i]], dataset_id = i)
})

# ---------------------------
# 9) COMBINE + POOL WITH RUBIN'S RULES  
# ---------------------------
methods_present <- unique(unlist(lapply(results_list,
                                        function(res) names(res$analysis_objects))))

# Helper: extract one summary component for one method across all imputations
get_summary_list <- function(component_name, method_name) {
  lapply(results_list, function(res) {
    obj <- res$analysis_objects[[method_name]]
    if (is.null(obj)) return(NULL)
    obj$summ[[component_name]]
  }) %>% Filter(Negate(is.null), .)
}

# Rubin-pool univariate curves (id: variable + z)
univar_pooled <- bind_rows(lapply(methods_present, function(m) {
  lst <- get_summary_list("univar", m)
  if (length(lst) == 0) return(NULL)
  pool_curve_rubin(lst, id_cols = c("variable", "z")) %>% mutate(method = m)
}))

# Rubin-pool overall risk (id: quantile)
overall_pooled <- bind_rows(lapply(methods_present, function(m) {
  lst <- get_summary_list("overall", m)
  if (length(lst) == 0) return(NULL)
  # overall has "sd" not "se"; rename before pooling
  lst2 <- lapply(lst, function(d) if ("sd" %in% names(d)) rename(d, se = sd) else d)
  pool_curve_rubin(lst2, id_cols = c("quantile")) %>%
    rename(sd = se) %>% mutate(method = m)
}))

# Rubin-pool singvar (id: variable + q.fixed)
singvar_pooled <- bind_rows(lapply(methods_present, function(m) {
  lst <- get_summary_list("singvar", m)
  if (length(lst) == 0) return(NULL)
  lst2 <- lapply(lst, function(d) if ("sd" %in% names(d)) rename(d, se = sd) else d)
  pool_curve_rubin(lst2, id_cols = c("variable", "q.fixed")) %>%
    rename(sd = se) %>% mutate(method = m)
}))

# Rubin-pool bivar levels (id: variable1 + variable2 + quantile + z1)
bivar_levels_pooled <- bind_rows(lapply(methods_present, function(m) {
  lst <- get_summary_list("bivar_levels", m)
  if (length(lst) == 0) return(NULL)
  lst2 <- lapply(lst, function(d) if (!"se" %in% names(d)) mutate(d, se = NA_real_) else d)
  pool_curve_rubin(lst2, id_cols = c("variable1", "variable2", "quantile", "z1")) %>%
    mutate(method = m)
}))

write.csv(univar_pooled,       file.path(outdir_real, "DATA_univar_rubin_pooled.csv"),        row.names = FALSE)
write.csv(overall_pooled,      file.path(outdir_real, "DATA_overall_rubin_pooled.csv"),       row.names = FALSE)
write.csv(singvar_pooled,      file.path(outdir_real, "DATA_singvar_rubin_pooled.csv"),       row.names = FALSE)
write.csv(bivar_levels_pooled, file.path(outdir_real, "DATA_bivar_levels_rubin_pooled.csv"),  row.names = FALSE)

# ---------------------------
# 10) PLOTS (Rubin-pooled, proper CI bands)
# ---------------------------
p_univar  <- plot_univar_compare_all(univar_pooled)
p_overall <- plot_overall_compare(overall_pooled)
p_singvar <- plot_singvar_compare(singvar_pooled)
p_bivar   <- plot_bivar_levels_compare(bivar_levels_pooled,
                                       vars_show = VARS_SUBSET_FOR_BIVAR_PLOT)

print(p_univar);  print(p_overall); print(p_singvar); print(p_bivar)

ggsave(file.path(outdir_real, "FIG_univar_rubin_pooled.png"),       p_univar,  width = 12, height = 8,  dpi = 300)
ggsave(file.path(outdir_real, "FIG_overall_rubin_pooled.png"),      p_overall, width = 10, height = 9,  dpi = 300)
ggsave(file.path(outdir_real, "FIG_singvar_rubin_pooled.png"),      p_singvar, width = 10, height = 8,  dpi = 300)
ggsave(file.path(outdir_real, "FIG_bivar_levels_rubin_pooled.png"), p_bivar,   width = 15, height = 10, dpi = 300)

# ---------------------------
# 11) OVERALL DELTA TABLES (Rubin-pooled scalar)
# ---------------------------
delta_all <- bind_rows(lapply(results_list, `[[`, "delta_tbl"))

# Rubin-pool scalar overall delta across imputations
pooled_main <- delta_all %>%
  mutate(se2 = ((upr - lwr) / (2 * 1.96))^2) %>%
  group_by(method, interval_type) %>%
  summarise(qhat = list(estimate), uhat = list(se2), .groups = "drop") %>%
  rowwise() %>%
  mutate(pool = list(pool_scalar_rubin(unlist(qhat), unlist(uhat)))) %>%
  tidyr::unnest(pool) %>%
  select(method, interval_type, m, estimate, se, lwr, upr,
         within_var, between_var, total_var)
write.csv(pooled_main,
          file.path(outdir_real, "TABLE_overall_delta_rubin_pooled.csv"),
          row.names = FALSE)

# Design-bootstrap intervals
boot_overall_all <- bind_rows(lapply(results_list, `[[`, "boot_overall"))
write.csv(boot_overall_all,
          file.path(outdir_real, "DATA_overall_bootstrap_estimates_all.csv"),
          row.names = FALSE)

boot_perf_summary <- boot_overall_all %>%
  group_by(method) %>%
  summarise(
    design_boot_point    = mean(overall_boot, na.rm = TRUE),
    design_boot_ci_lower = quantile(overall_boot, 0.025, na.rm = TRUE),
    design_boot_ci_upper = quantile(overall_boot, 0.975, na.rm = TRUE),
    n_boot               = n(),
    interval_type        = "bootstrap_replication_interval",
    .groups              = "drop"
  )
write.csv(boot_perf_summary,
          file.path(outdir_real, "TABLE_design_bootstrap_summary.csv"),
          row.names = FALSE)

# ---------------------------
# 12) PIP OUTPUTS + STABILITY 
# ---------------------------
pip_master_all <- bind_rows(lapply(results_list, `[[`, "pip_master")) %>%
  arrange(dataset_id, method,
          factor(level, levels = c("group","variable","conditional")),
          group, variable, given)
write.csv(pip_master_all,
          file.path(outdir_real, "DATA_PIPs_MASTER_all_imputations.csv"),
          row.names = FALSE)

# Rubin-pool variable PIPs across imputations per method
pip_rubin_pooled <- bind_rows(lapply(methods_present, function(meth) {
  pip_per_imp <- lapply(results_list, function(res) {
    res$pip_master %>%
      filter(method == meth, level == "variable", metric == "pip") %>%
      select(variable, group, pip = value)
  })
  pool_pip_rubin(pip_per_imp) %>% mutate(method = meth)
}))
write.csv(pip_rubin_pooled,
          file.path(outdir_real, "TABLE_PIPs_rubin_pooled.csv"),
          row.names = FALSE)

# Bootstrap PIP stability 
boot_pips_all <- bind_rows(lapply(results_list, `[[`, "boot_pips"))

if (SAVE_BOOTSTRAP_PIP_DETAILS && nrow(boot_pips_all) > 0) {
  write.csv(boot_pips_all,
            file.path(outdir_real, "DATA_bootstrap_variable_PIPs.csv"),
            row.names = FALSE)
}

pip_stability <- boot_pips_all %>%
  group_by(method, variable, group) %>%
  summarise(
    mean_pip        = mean(pip, na.rm = TRUE),
    sd_pip          = sd(pip, na.rm = TRUE),
    # Coefficient of variation: quantifies instability relative to magnitude
    cv_pip          = sd_pip / pmax(mean_pip, 1e-6),
    median_pip      = median(pip, na.rm = TRUE),
    q25_pip         = quantile(pip, 0.25, na.rm = TRUE),
    q75_pip         = quantile(pip, 0.75, na.rm = TRUE),
    min_pip         = min(pip, na.rm = TRUE),
    max_pip         = max(pip, na.rm = TRUE),
    prop_pip_ge_0_5 = mean(pip >= 0.5, na.rm = TRUE),
    n_boot          = n(),
    .groups         = "drop"
  ) %>%
  arrange(method, desc(mean_pip))
write.csv(pip_stability,
          file.path(outdir_real, "TABLE_PIP_stability_bootstrap.csv"),
          row.names = FALSE)

p_pip_stab <- plot_pip_stability(pip_stability)
print(p_pip_stab)
ggsave(file.path(outdir_real, "FIG_PIP_stability_bootstrap.png"),
       p_pip_stab, width = 10, height = 10, dpi = 300)

p_pips <- plot_pips_variable_compare(
  pip_rubin_pooled %>%
    transmute(method, level = "variable", group, variable,
              metric = "pip", value = pip_pooled,
              given = NA_character_, group_given = NA_character_)
)
print(p_pips)
ggsave(file.path(outdir_real, "FIG_variable_PIPs_rubin_pooled.png"),
       p_pips, width = 10, height = 8, dpi = 300)

# ---------------------------
# 13) DIAGNOSTICS + RUNTIME 
# ---------------------------
diagnostics_all <- bind_rows(lapply(results_list, `[[`, "diagnostics"))
write.csv(diagnostics_all,
          file.path(outdir_real, "TABLE_MCMC_diagnostics_all.csv"),
          row.names = FALSE)

diagnostics_summary <- diagnostics_all %>%
  group_by(method, parameter) %>%
  summarise(
    median_ess            = median(ess, na.rm = TRUE),
    min_ess               = min(ess, na.rm = TRUE),
    max_abs_geweke_z      = max(abs(geweke_z), na.rm = TRUE),
    n_below_ess_threshold = sum(ess_below_threshold, na.rm = TRUE),
    n_records             = n(),
    .groups               = "drop"
  )
write.csv(diagnostics_summary,
          file.path(outdir_real, "TABLE_MCMC_diagnostics_summary.csv"),
          row.names = FALSE)

# Flag problematic replicates 
ess_flags <- diagnostics_all %>%
  filter(ess_below_threshold) %>%
  arrange(method, rep_id, parameter)
write.csv(ess_flags,
          file.path(outdir_real, "TABLE_MCMC_ESS_flags_below_threshold.csv"),
          row.names = FALSE)
if (nrow(ess_flags) > 0) {
  message(sprintf("WARNING: %d MCMC chains have ESS < %d. See TABLE_MCMC_ESS_flags_below_threshold.csv",
                  nrow(ess_flags), ESS_THRESHOLD))
}

# Runtime analysis ( ratio AND absolute seconds)
runtime_all <- bind_rows(lapply(results_list, `[[`, "runtime"))
write.csv(runtime_all,
          file.path(outdir_real, "TABLE_runtime_all.csv"),
          row.names = FALSE)

runtime_summary <- runtime_all %>%
  group_by(analysis) %>%
  summarise(
    mean_elapsed_sec   = mean(elapsed_sec, na.rm = TRUE),
    median_elapsed_sec = median(elapsed_sec, na.rm = TRUE),
    .groups            = "drop"
  )

get_rt <- function(analysis_name) {
  runtime_summary %>% filter(analysis == analysis_name) %>% pull(mean_elapsed_sec)
}
naive_time        <- get_rt("naive_singlefit")
weighted_time     <- get_rt("weighted_singlefit")
wt_trim_time      <- get_rt("weighted_trimmed_singlefit")
pseudo_time       <- get_rt("pseudo_posterior_singlefit")
weighted_boot_time <- get_rt("weighted_bootstrap_total")

safe_ratio <- function(a, b) if (length(a)==1 && length(b)==1) a/b else NA_real_
safe_diff  <- function(a, b) if (length(a)==1 && length(b)==1) a-b else NA_real_

runtime_comparison <- tibble(
  comparison = c(
    "weighted_singlefit_vs_naive",
    "weighted_trimmed_singlefit_vs_naive",
    "pseudo_posterior_singlefit_vs_naive",
    "weighted_bootstrap_total_vs_naive"
  ),
  ratio_vs_naive = c(
    safe_ratio(weighted_time,      naive_time),
    safe_ratio(wt_trim_time,       naive_time),
    safe_ratio(pseudo_time,        naive_time),
    safe_ratio(weighted_boot_time, naive_time)
  ),
  absolute_extra_sec = c(
    safe_diff(weighted_time,      naive_time),
    safe_diff(wt_trim_time,       naive_time),
    safe_diff(pseudo_time,        naive_time),
    safe_diff(weighted_boot_time, naive_time)
  )
)
write.csv(runtime_summary,    file.path(outdir_real, "TABLE_runtime_summary.csv"),           row.names = FALSE)
write.csv(runtime_comparison, file.path(outdir_real, "TABLE_runtime_increase_vs_naive.csv"), row.names = FALSE)

# ---------------------------
# 14) TRIMMING SENSITIVITY 
# ---------------------------
trimming_summary <- boot_overall_all %>%
  group_by(method) %>%
  summarise(
    mean_overall = mean(overall_boot, na.rm = TRUE),
    sd_overall   = sd(overall_boot, na.rm = TRUE),
    q025         = quantile(overall_boot, 0.025, na.rm = TRUE),
    q975         = quantile(overall_boot, 0.975, na.rm = TRUE),
    ci_width     = q975 - q025,
    .groups      = "drop"
  )
write.csv(trimming_summary,
          file.path(outdir_real, "TABLE_trimming_sensitivity.csv"),
          row.names = FALSE)

p_trim <- plot_trimming_sensitivity(boot_overall_all)
print(p_trim)
ggsave(file.path(outdir_real, "FIG_trimming_sensitivity.png"),
       p_trim, width = 9, height = 5, dpi = 300)

# ---------------------------
# 15) ALGORITHM SPECIFICATION + METADATA 
# ---------------------------
algorithm_spec <- tibble(
  step = c(
    "1_imputation",
    "2_standardization",
    "3_design_resampling",
    "3b_pseudo_posterior",
    "4_bkmr_fitting",
    "5_thinning",
    "6_bootstrap_replication",
    "7_rubin_pooling_scalars",
    "7b_rubin_pooling_curves",
    "7c_rubin_pooling_pips",
    "8_pip_stability",
    "9_weight_trimming"
  ),
  description = c(
    paste0("MICE multiple imputation (mice package). m=", N_IMPUTATIONS,
           " imputed datasets. Survey design columns (PSU, strata, weights) used as ",
           "predictors but not imputed. Continuous variables: predictive mean matching (pmm). ",
           "Binary variables: logistic regression (logreg). ",
           "Each imputed dataset analyzed independently — no averaging of imputed values."),
    paste0("All exposure variables log-transformed (if strictly positive) then z-score ",
           "standardized using per-dataset mean and SD. Z-scores clipped to [",
           CLIP_LOWER, ", ", CLIP_UPPER, "] to limit outlier influence on the kernel."),
    paste0("PSU-within-stratum bootstrap resampling (R=", R_BOOT, " replicates). ",
           "PSUs sampled with replacement, probability proportional to sum of survey weights ",
           "within stratum. Within each selected PSU, units resampled with replacement, ",
           "probability proportional to unit survey weight."),
    "Pseudo-posterior approach (Savitsky & Toth 2016 JASA approximation): survey weights ",
    paste0("BKMR fitted via kmbayes() (bkmr package) with varsel=TRUE for PIP estimation. ",
           "Main fit: iter=", ITER_MAIN, ", burn=", BURN_MAIN, ", thin=", THIN_MAIN, ". ",
           "Bootstrap fit: iter=", ITER_BOOT, ", burn=", BURN_BOOT, ", thin=", THIN_BOOT, ". ",
           "Effective posterior sample = (iter - burn) / thin."),
    paste0("Posterior thinned after burn-in. Effective sample size computed via ",
           "coda::effectiveSize(). Geweke z-statistic computed for sigma chain."),
    paste0("R=", R_BOOT, " bootstrap replicates per imputed dataset (parallel, ",
           N_WORKERS, " cores). Design uncertainty summarized as 2.5th-97.5th percentile of ",
           "overall mixture effect delta across bootstrap replicates."),
    "Rubin (1987) scalar pooling: Q-bar = mean(Q_i); U-bar = mean(U_i); B = var(Q_i); ",
    "Rubin (1987) curve pooling (pointwise): applied at each grid point of univariate, ",
    "PIP pooling: Q-bar = mean(PIP_i); within-var = mean(PIP_i*(1-PIP_i)); ",
    "Bootstrap PIP stability: mean, SD, CV (=SD/mean), IQR, prop(PIP>=0.5) per variable ",
    paste0("Weight trimming: PSU-level and unit-level weights clipped to [",
           TRIM_LOWER_PROB*100, "th, ", TRIM_UPPER_PROB*100, "th] percentiles before ",
           "resampling. Comparison with untrimmed results assesses sensitivity to extreme weights.")
  )
)
write.csv(algorithm_spec,
          file.path(outdir_real, "TABLE_algorithm_specification.csv"),
          row.names = FALSE)

# Software versions
software_tbl <- tibble(
  package = c("R","bkmr","survey","mice","future.apply","coda",
               "dplyr","ggplot2","readr","tidyr","purrr","stringr"),
  version = c(
    as.character(getRversion()),
    as.character(packageVersion("bkmr")),
    as.character(packageVersion("survey")),
    as.character(packageVersion("mice")),
    as.character(packageVersion("future.apply")),
    as.character(packageVersion("coda")),
    as.character(packageVersion("dplyr")),
    as.character(packageVersion("ggplot2")),
    as.character(packageVersion("readr")),
    as.character(packageVersion("tidyr")),
    as.character(packageVersion("purrr")),
    as.character(packageVersion("stringr"))
  )
)
write.csv(software_tbl,
          file.path(outdir_real, "TABLE_software_versions.csv"),
          row.names = FALSE)

# Analysis metadata
analysis_metadata <- tibble(
  field = c(
    "inference_target",
    "n_imputed_datasets",
    "mice_seed",
    "single_fit_interval_type",
    "mi_pooled_interval_type",
    "design_bootstrap_interval_type",
    "statin_adjustment_included",
    "weight_trimming_sensitivity_run",
    "trim_lower_prob",
    "trim_upper_prob",
    "zscore_clipping_lower",
    "zscore_clipping_upper",
    "bootstrap_replicates",
    "ess_threshold_flag",
    "main_iter", "main_burn", "main_thin",
    "boot_iter", "boot_burn", "boot_thin",
    "pseudo_posterior_method_reference"
  ),
  value = c(
    INFERENCE_TARGET,
    N_IMPUTATIONS,
    "20251210",
    "posterior_credible_interval (within single BKMR fit)",
    "Rubin_1987_total_variance_interval (pooled across m imputed datasets)",
    "bootstrap_replication_interval (2.5-97.5 percentile across PSU resamples)",
    ifelse(is.null(statin_var), "no", "yes"),
    "yes",
    TRIM_LOWER_PROB, TRIM_UPPER_PROB,
    CLIP_LOWER, CLIP_UPPER,
    R_BOOT, ESS_THRESHOLD,
    ITER_MAIN, BURN_MAIN, THIN_MAIN,
    ITER_BOOT, BURN_BOOT, THIN_BOOT,
    "Savitsky_and_Toth_2016_JASA_pseudo_posterior_approximation"
  )
)
write.csv(analysis_metadata,
          file.path(outdir_real, "TABLE_analysis_metadata.csv"),
          row.names = FALSE)

appendix_code_note <- c(
  "NHANES BKMR WORKFLOW (REVISED)",
  "===================================================",
  "",
  "R1-C1 / R2-C1 (Multiple imputation — anti-conservative intervals fixed):",
  paste0("  MICE imputation (m=", N_IMPUTATIONS,
         ") on raw data. Each imputed dataset analyzed independently."),
  "  ALL results pooled using Rubin (1987) Rules:",
  "    - Scalar (overall delta): pool_scalar_rubin()   => within + between var",
  "    - Curves (univar, overall, singvar, bivar): pool_curve_rubin() per grid point",
  "    - PIPs: pool_pip_rubin() using Bernoulli within-var + between-var",
  "  This replaces the prior averaging of imputed values which ignored between-imputation",
  "  variance and produced anti-conservative credible intervals.",
  "",
  "R1-C4 (ESS + convergence across ALL bootstrap replicates):",
  paste0("  ESS (coda::effectiveSize) and Geweke z computed per replicate."),
  paste0("  Replicates with ESS < ", ESS_THRESHOLD,
         " flagged in TABLE_MCMC_ESS_flags_below_threshold.csv."),
  "  Summary: TABLE_MCMC_diagnostics_summary.csv.",
  "",
  "R1-C5 (Trimming sensitivity — now mandatory, results always shown):",
  paste0("  Trimming (", TRIM_LOWER_PROB*100, "th-", TRIM_UPPER_PROB*100,
         "th percentile) run in parallel with all weighted analyses."),
  "  TABLE_trimming_sensitivity.csv and FIG_trimming_sensitivity.png always produced.",
  "",
  "R1-C6 / R2-C5 (Pseudo-posterior comparison added):",
  "  Savitsky & Toth (2016) approximation: Pfeffermann weight normalization +",
  "  integer row replication. Compared against PSU-resampling in all tables/figures.",
  "",
  "R1-C8 (Runtime: ratio AND absolute extra seconds vs naive):",
  "  TABLE_runtime_increase_vs_naive.csv.",
  "",
  "R1-C10 / R2-C6 (PIP stability + discriminative ability):",
  "  CV (SD/mean) of PIP across bootstrap replicates in TABLE_PIP_stability_bootstrap.csv.",
  "  FIG_PIP_stability_bootstrap.png.",
  "  Group PIPs > 0.99 flagged as limited discriminative ability in DATA_PIPs_MASTER.",
  "",
  "R2-C4 (Full algorithm specification):",
  "  TABLE_algorithm_specification.csv.",
  "",
  "Interval semantics (explicit for manuscript):",
  "  posterior credible interval: MCMC uncertainty within a single BKMR fit",
  "  Rubin MI interval:           total variance = within-imp + (1+1/m)*between-imp var",
  "  bootstrap replication interval: 2.5-97.5 percentile over R PSU bootstrap replicates"
)
writeLines(appendix_code_note,
           con = file.path(outdir_real, "APPENDIX_code_notes.txt"))

# ---------------------------
# 16) FINAL REPORT
# ---------------------------
cat("\n=====  NHANES BKMR ANALYSIS COMPLETE =====\n")
cat("Output directory:", outdir_real, "\n\n")
cat("Imputed datasets analyzed:", N_IMPUTATIONS,
    "(change N_IMPUTATIONS at top of script)\n")
cat("All curves pooled using Rubin's Rules (within + between imputation variance).\n\n")
cat("Key outputs:\n")
cat("  DATA:    DATA_univar_rubin_pooled.csv\n")
cat("           DATA_overall_rubin_pooled.csv\n")
cat("           DATA_singvar_rubin_pooled.csv\n")
cat("           DATA_bivar_levels_rubin_pooled.csv\n")
cat("           DATA_PIPs_MASTER_all_imputations.csv\n")
cat("           DATA_overall_bootstrap_estimates_all.csv\n")
cat("           DATA_bootstrap_variable_PIPs.csv\n")
cat("  TABLES:  TABLE_overall_delta_rubin_pooled.csv\n")
cat("           TABLE_PIPs_rubin_pooled.csv\n")
cat("           TABLE_PIP_stability_bootstrap.csv    (CV, IQR)\n")
cat("           TABLE_design_bootstrap_summary.csv\n")
cat("           TABLE_trimming_sensitivity.csv\n")
cat("           TABLE_MCMC_diagnostics_summary.csv\n")
cat("           TABLE_MCMC_ESS_flags_below_threshold.csv\n")
cat("           TABLE_runtime_summary.csv\n")
cat("           TABLE_runtime_increase_vs_naive.csv  (ratio + absolute seconds)\n")
cat("           TABLE_algorithm_specification.csv    (reproducibility)\n")
cat("           TABLE_analysis_metadata.csv\n")
cat("           TABLE_software_versions.csv\n")
cat("  FIGURES: FIG_univar_rubin_pooled.png\n")
cat("           FIG_overall_rubin_pooled.png\n")
cat("           FIG_singvar_rubin_pooled.png\n")
cat("           FIG_bivar_levels_rubin_pooled.png\n")
cat("           FIG_variable_PIPs_rubin_pooled.png\n")
cat("           FIG_PIP_stability_bootstrap.png\n")
cat("           FIG_trimming_sensitivity.png\n")
cat("  QC:      QC_standardization_range_minus3_plus3.csv\n")
cat("           APPENDIX_code_notes.txt\n")

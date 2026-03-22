
##### REVIEWER-RESPONSIVE NHANES BKMR SCRIPT ###########
# ======================================================
# NHANES REAL DATA + BKMR
# Naive vs design-aware weighted workflow using PSU bootstrap resampling
#
# Key reviewer-driven corrections implemented:
#   1) Distinguishes BKMR posterior credible intervals from design-bootstrap intervals
#   2) Supports optional statin/lipid-lowering medication adjustment
#   3) Supports optional multiple-imputation workflow (fit each imputed dataset separately)
#   4) Reports convergence / ESS diagnostics
#   5) Adds weight-trimming sensitivity analysis
#   6) Quantifies runtime increase of design-aware workflow
#   7) States inferential target explicitly in machine-readable output
#   8) Assesses PIP variability across bootstrap replicates
#   9) Writes software / reproducibility metadata
#
# Notes:
# - BKMR itself is still fit with kmbayes() from the bkmr package.
# - "Weighted" in this script refers to a design-aware workflow that uses
#   PSU-within-stratum resampling with probabilities informed by survey weights.
# - Single-fit BKMR intervals are posterior credible intervals.
# - Design-aware intervals are bootstrap / replication-based intervals across resamples.
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
})

# ---------------------------
# 0) USER SETTINGS
# ---------------------------
# EITHER provide one analysis file:
infile <- "/Users/doreennoldajehu-appiah/Desktop/PAPER_3/PIP3/RealData_Sim/TC.csv"

# OR provide multiple imputed files; leave NULL to use infile above.
# Example:
# imputed_files <- list.files("/path/to/imputed/", pattern = "\\.csv$", full.names = TRUE)
imputed_files <- NULL

outdir_real <- "/Users/doreennoldajehu-appiah/Desktop/PAPER_3/PIP3/RealData_Sim/REVISION"
dir.create(outdir_real, showWarnings = FALSE, recursive = TRUE)

# Outcome and covariates
y_var <- "TotalCholesterol"
covar_vars_base <- c("age", "Ethnicity", "income", "Gender", "BMI", "smq020", "alq101")

# Optional statin / lipid-lowering medication indicator
# Set to NULL if not available yet.
statin_var <- NULL
# Example:
# statin_var <- "statin_use"

# Survey design columns
psu_col    <- "sdmvpsu"
strata_col <- "sdmvstra"
wt_col     <- "wtmec2yr"

# Exposures
exposure_vars <- c(
  "Arsenic","Mercury","Barium","Cobalt","Cesium",
  "Molybdenum","Lead","Antimony","Tin","Strontium",
  "Thallium","Tungsten","Uranium",
  "PFOA","PFOS","PFDE","PFHS","NMeFOSAA",
  "PFHP","PFNA","PFUA","PFDO",
  "Mono(carboxynonyl) phthalate",
  "Mono(carboxyoctyl) phthalate",
  "Mono-n-butyl phthalate",
  "Mono(3-carboxypropyl) phthalate",
  "urdmc1lc",
  "Mono-ethyl phthalate",
  "Mono(2-ethyl-5-hydroxyhexyl) phthalate",
  "Mono(2-ethyl-5-hydroxy-nonyl) phthalate",
  "Mono(2-ethylhexyl) phthalate",
  "Mono-isobutyl phthalate",
  "Mono-benzyl phthalate"
)

pip_groups <- list(
  Metals = c(
    "Arsenic","Mercury","Barium","Cobalt","Cesium",
    "Molybdenum","Lead","Antimony","Tin","Strontium",
    "Thallium","Tungsten","Uranium"
  ),
  PFAS = c("PFOA","PFOS","PFDE","PFHS","NMeFOSAA","PFHP","PFNA","PFUA","PFDO"),
  Phthalates_Plasticizers = c(
    "Mono(carboxynonyl) phthalate",
    "Mono(carboxyoctyl) phthalate",
    "Mono-n-butyl phthalate",
    "Mono(3-carboxypropyl) phthalate",
    "urdmc1lc",
    "Mono-ethyl phthalate",
    "Mono(2-ethyl-5-hydroxyhexyl) phthalate",
    "Mono(2-ethyl-5-hydroxy-nonyl) phthalate",
    "Mono(2-ethylhexyl) phthalate",
    "Mono-isobutyl phthalate",
    "Mono-benzyl phthalate"
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

QS_OVERALL   <- seq(0.25, 0.75, by = 0.05)
QS_BIVARLVL  <- c(0.25, 0.50, 0.75)
SINGVAR_QFIX <- c(0.25, 0.50, 0.75)
SINGVAR_QDIFF <- c(0.25, 0.75)

R_BOOT    <- 50
ITER_BOOT <- 2500
BURN_BOOT <- 1000
THIN_BOOT <- 2

# Heavy outputs
SAVE_BIVAR_SURFACE <- FALSE
VARS_SUBSET_FOR_BIVAR_PLOT <- exposure_vars

# Standardization clipping
CLIP_LOWER <- -3
CLIP_UPPER <-  3

# Optional trimming sensitivity
RUN_TRIMMING_SENSITIVITY <- TRUE
TRIM_LOWER_PROB <- 0.01
TRIM_UPPER_PROB <- 0.99

# PIP stability
SAVE_BOOTSTRAP_PIP_DETAILS <- TRUE

# Parallel
N_WORKERS <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = N_WORKERS)

# Inferential target
# Choose one explicitly for the paper.
INFERENCE_TARGET <- "superpopulation_mean"

# ---------------------------
# 1) HELPERS
# ---------------------------
get_analysis_files <- function(infile, imputed_files = NULL) {
  if (!is.null(imputed_files) && length(imputed_files) > 0) {
    return(imputed_files)
  }
  infile
}

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
    if (!is.null(dat_names) && !statin_var %in% dat_names) {
      stop(sprintf("statin_var='%s' not found in data.", statin_var))
    }
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
    min_z = sapply(vars, function(v) min(dat[[v]], na.rm = TRUE)),
    max_z = sapply(vars, function(v) max(dat[[v]], na.rm = TRUE)),
    mean_z = sapply(vars, function(v) mean(dat[[v]], na.rm = TRUE)),
    sd_z = sapply(vars, function(v) sd(dat[[v]], na.rm = TRUE))
  )
}

prepare_one_dataset <- function(path, y_var, covar_vars_base, statin_var,
                                exposure_vars, psu_col, strata_col, wt_col,
                                clip_lower, clip_upper) {
  dat_raw <- read_csv(path, show_col_types = FALSE) %>% as.data.frame()
  covar_vars <- get_covar_vars(covar_vars_base, statin_var, names(dat_raw))
  need_cols <- unique(c(y_var, covar_vars, exposure_vars, psu_col, strata_col, wt_col))
  miss <- setdiff(need_cols, names(dat_raw))
  if (length(miss) > 0) stop("Missing required columns in ", basename(path), ": ", paste(miss, collapse = ", "))

  dat_raw <- safe_numeric_df(dat_raw, exposure_vars)
  if (!is.null(statin_var)) {
    dat_raw[[statin_var]] <- suppressWarnings(as.numeric(dat_raw[[statin_var]]))
  }

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
    raw = dat_raw,
    cc = dat_cc,
    std = dat_std,
    qc = qc,
    mu_vec = mu_vec,
    sd_vec = sd_vec,
    covar_vars = covar_vars,
    file = path
  )
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
  dat <- as.data.frame(dat)
  strata_list <- split(dat, dat[[strata_col]])

  out_list <- vector("list", length(strata_list))
  k <- 0L

  for (s in names(strata_list)) {
    ds <- strata_list[[s]]
    psus <- unique(ds[[psu_col]])
    H_s <- length(psus)

    psu_w <- ds %>%
      group_by(.data[[psu_col]]) %>%
      summarise(W_psu = sum(.data[[wt_col]], na.rm = TRUE), .groups = "drop")

    if (trim) {
      psu_w$W_psu <- trim_weights(psu_w$W_psu, trim_lower_prob, trim_upper_prob)
    }

    psu_w <- psu_w %>% mutate(p = W_psu / sum(W_psu))

    sampled_psus <- sample(
      psus, size = H_s, replace = TRUE,
      prob = psu_w$p[match(psus, psu_w[[psu_col]])]
    )

    rows_s <- vector("list", length(sampled_psus))
    for (j in seq_along(sampled_psus)) {
      p_id <- sampled_psus[j]
      dpsu <- ds[ds[[psu_col]] == p_id, , drop = FALSE]
      m <- nrow(dpsu)
      if (m == 0) next

      unit_w <- dpsu[[wt_col]]
      if (trim) {
        unit_w <- trim_weights(unit_w, trim_lower_prob, trim_upper_prob)
      }
      unit_p <- unit_w / sum(unit_w)

      take <- sample(seq_len(m), size = m, replace = TRUE, prob = unit_p)
      rows_s[[j]] <- dpsu[take, , drop = FALSE]
    }

    k <- k + 1L
    out_list[[k]] <- bind_rows(rows_s)
  }

  bind_rows(out_list) %>% as.data.frame()
}

# ---------------------------
# 3) BKMR FITTING + DIAGNOSTICS
# ---------------------------
fit_kmbayes_fast <- function(dat, y_var, exposure_vars, covar_vars, iter, varsel, verbose = FALSE) {
  y <- dat[[y_var]]
  Z <- as.matrix(dat[, exposure_vars, drop = FALSE])
  X <- as.matrix(dat[, covar_vars, drop = FALSE])

  kmbayes(
    y = y, Z = Z, X = X,
    iter = iter,
    varsel = varsel,
    verbose = verbose
  )
}

thin_fit <- function(fit, burn, thin) {
  keep <- seq.int(from = burn + 1, to = fit$iter, by = thin)
  fit2 <- fit

  for (nm in c("beta", "lambda", "r", "sigma", "delta", "rhat", "h")) {
    if (!is.null(fit2[[nm]])) {
      obj <- fit2[[nm]]
      if (is.matrix(obj)) fit2[[nm]] <- obj[keep, , drop = FALSE]
      if (is.vector(obj)) fit2[[nm]] <- obj[keep]
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
  dat_rs <- resample_nhanes_psu(
    dat, psu_col, strata_col, wt_col,
    trim = trim,
    trim_lower_prob = trim_lower_prob,
    trim_upper_prob = trim_upper_prob
  )

  fit_raw <- fit_kmbayes_fast(dat_rs, y_var, exposure_vars, covar_vars,
                              iter = iter, varsel = varsel, verbose = FALSE)
  fit <- thin_fit(fit_raw, burn = burn, thin = thin)
  list(fit = fit, dat_rs = dat_rs)
}

compute_fit_diagnostics <- function(fit, method_label, rep_id = NA_integer_) {
  out <- list()

  if (!is.null(fit$sigma)) {
    ess_sigma <- tryCatch(as.numeric(coda::effectiveSize(fit$sigma)), error = function(e) NA_real_)
    geweke_sigma <- tryCatch({
      gz <- coda::geweke.diag(as.mcmc(fit$sigma))$z
      if (length(gz) == 1) as.numeric(gz) else NA_real_
    }, error = function(e) NA_real_)
    out[[length(out) + 1L]] <- tibble(
      method = method_label, rep_id = rep_id, parameter = "sigma",
      ess = ess_sigma, geweke_z = geweke_sigma
    )
  }

  if (!is.null(fit$lambda)) {
    lam <- fit$lambda
    if (is.matrix(lam)) {
      ess_vals <- apply(lam, 2, function(x) tryCatch(as.numeric(coda::effectiveSize(x)), error = function(e) NA_real_))
      out[[length(out) + 1L]] <- tibble(
        method = method_label, rep_id = rep_id,
        parameter = paste0("lambda_", seq_along(ess_vals)),
        ess = as.numeric(ess_vals), geweke_z = NA_real_
      )
    } else {
      out[[length(out) + 1L]] <- tibble(
        method = method_label, rep_id = rep_id, parameter = "lambda",
        ess = tryCatch(as.numeric(coda::effectiveSize(lam)), error = function(e) NA_real_),
        geweke_z = NA_real_
      )
    }
  }

  if (!is.null(fit$r)) {
    rr <- fit$r
    if (is.matrix(rr)) {
      ess_vals <- apply(rr, 2, function(x) tryCatch(as.numeric(coda::effectiveSize(x)), error = function(e) NA_real_))
      out[[length(out) + 1L]] <- tibble(
        method = method_label, rep_id = rep_id,
        parameter = paste0("r_", seq_along(ess_vals)),
        ess = as.numeric(ess_vals), geweke_z = NA_real_
      )
    }
  }

  if (!is.null(fit$delta)) {
    delta_mean <- rowMeans(as.matrix(fit$delta))
    out[[length(out) + 1L]] <- tibble(
      method = method_label, rep_id = rep_id,
      parameter = "delta_mean",
      ess = tryCatch(as.numeric(coda::effectiveSize(delta_mean)), error = function(e) NA_real_),
      geweke_z = NA_real_
    )
  }

  bind_rows(out)
}

# ---------------------------
# 4) SUMMARIES
# ---------------------------
make_bkmr_summaries <- function(fit, Zmat,
                                qs_overall = QS_OVERALL,
                                qs_bivar = QS_BIVARLVL) {
  pred.univar <- PredictorResponseUnivar(fit = fit, method = "approx")
  pred.bivar  <- PredictorResponseBivar(fit = fit, min.plot.dist = 1, method = "approx")

  pred.bivar.levels <- PredictorResponseBivarLevels(
    pred.resp.df = pred.bivar, Z = Zmat,
    both_pairs = TRUE, qs = qs_bivar
  )

  risks.overall <- OverallRiskSummaries(
    fit = fit, qs = qs_overall, q.fixed = 0.5, method = "approx"
  )

  risks.singvar <- SingVarRiskSummaries(
    fit = fit, qs.diff = SINGVAR_QDIFF,
    q.fixed = SINGVAR_QFIX,
    method = "approx"
  )

  list(
    univar = pred.univar,
    bivar = pred.bivar,
    bivar_levels = pred.bivar.levels,
    overall = risks.overall,
    singvar = risks.singvar
  )
}

extract_overall_delta <- function(fit) {
  or <- OverallRiskSummaries(
    fit = fit, qs = c(0.25, 0.75), q.fixed = 0.5, method = "approx"
  )
  est <- with(or, est[quantile == 0.75] - est[quantile == 0.25])
  sdv <- with(or, sqrt(sd[quantile == 0.75]^2 + sd[quantile == 0.25]^2))
  c(mean = est, lwr = est - 1.96 * sdv, upr = est + 1.96 * sdv)
}

# Rubin-style pooling for scalar summaries only
pool_scalar_rubin <- function(qhat, uhat) {
  m <- length(qhat)
  qbar <- mean(qhat, na.rm = TRUE)
  ubar <- mean(uhat, na.rm = TRUE)
  b <- stats::var(qhat, na.rm = TRUE)
  tvar <- ubar + (1 + 1/m) * b
  se <- sqrt(tvar)
  tibble(
    m = m,
    estimate = qbar,
    se = se,
    lwr = qbar - 1.96 * se,
    upr = qbar + 1.96 * se,
    within_var = ubar,
    between_var = b,
    total_var = tvar
  )
}

pool_curve_mean <- function(lst, id_cols) {
  bind_rows(lst, .id = "imputation_id") %>%
    group_by(across(all_of(id_cols))) %>%
    summarise(
      est = mean(est, na.rm = TRUE),
      se = if ("se" %in% names(cur_data())) sqrt(mean(se^2, na.rm = TRUE)) else NA_real_,
      sd = if ("sd" %in% names(cur_data())) sqrt(mean(sd^2, na.rm = TRUE)) else NA_real_,
      .groups = "drop"
    )
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
    pip_col <- intersect(c("pip", "PIP", "pips", "PIPs"), names(p))[1]
    var_col <- intersect(c("variable", "var", "exposure", "Z"), names(p))[1]
    if (!is.na(pip_col) && !is.na(var_col)) {
      vec <- p[[pip_col]]
      names(vec) <- as.character(p[[var_col]])
      vec <- vec[exposure_vars]
      return(tibble(variable = exposure_vars, pip = as.numeric(vec)))
    }
    stop("ExtractPIPs returned data.frame but pip/variable columns were not found.")
  }

  if (is.list(p)) {
    cand <- c("pip", "PIP", "pips", "PIPs", "p")
    pip_name <- cand[cand %in% names(p)][1]
    if (!is.na(pip_name)) {
      vec <- p[[pip_name]]
      if (!is.numeric(vec)) stop("ExtractPIPs list element found but not numeric.")
      if (!is.null(names(vec)) && all(exposure_vars %in% names(vec))) {
        vec <- vec[exposure_vars]
      } else {
        vec <- vec[seq_along(exposure_vars)]
        names(vec) <- exposure_vars
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
      variable %in% group_list$Metals ~ "Metals",
      variable %in% group_list$PFAS ~ "PFAS",
      variable %in% group_list$Phthalates_Plasticizers ~ "Phthalates_Plasticizers",
      TRUE ~ NA_character_
    ))
}

compute_group_pips <- function(pip_df) {
  pip_df %>%
    filter(!is.na(group)) %>%
    group_by(group) %>%
    summarise(
      n_vars   = n(),
      pip_any  = 1 - prod(1 - pip),
      pip_mean = mean(pip),
      pip_max  = max(pip),
      .groups  = "drop"
    )
}

get_delta_matrix <- function(fit) {
  d <- fit$delta
  if (is.null(d)) stop("fit$delta is NULL. Ensure varsel=TRUE.")
  as.matrix(d)
}

compute_conditional_pips_long <- function(delta_mat, exposure_vars) {
  p <- length(exposure_vars)
  denom <- colMeans(delta_mat == 1)

  cond_mat <- matrix(NA_real_, nrow = p, ncol = p,
                     dimnames = list(exposure_vars, exposure_vars))

  for (j in seq_len(p)) {
    if (denom[j] == 0) next
    joint <- colMeans((delta_mat == 1) & (delta_mat[, j] == 1))
    cond_mat[, j] <- joint / denom[j]
  }

  as_tibble(as.data.frame(as.table(cond_mat), stringsAsFactors = FALSE)) %>%
    rename(variable_i = Var1, given_j = Var2, cond_pip = Freq)
}

build_pip_master <- function(fit, method_label, exposure_vars, group_list) {
  flat <- get_flat_pips(fit, exposure_vars) %>%
    mutate(method = method_label) %>%
    attach_groups_to_pips(group_list)

  grp <- compute_group_pips(flat) %>% mutate(method = method_label)

  delta <- get_delta_matrix(fit)

  cond_long <- compute_conditional_pips_long(delta, exposure_vars) %>%
    mutate(method = method_label) %>%
    left_join(flat %>% select(variable, group) %>% distinct(),
              by = c("variable_i" = "variable")) %>%
    rename(group_i = group) %>%
    left_join(flat %>% select(variable, group) %>% distinct(),
              by = c("given_j" = "variable")) %>%
    rename(group_j = group)

  var_rows <- flat %>%
    transmute(
      method, level = "variable", group = group, variable = variable,
      metric = "pip", value = pip, given = NA_character_, group_given = NA_character_
    )

  grp_rows <- grp %>%
    pivot_longer(c(n_vars, pip_any, pip_mean, pip_max),
                 names_to = "metric", values_to = "value") %>%
    transmute(
      method, level = "group", group = group, variable = NA_character_,
      metric = metric, value = as.numeric(value),
      given = NA_character_, group_given = NA_character_
    )

  cond_rows <- cond_long %>%
    transmute(
      method, level = "conditional", group = group_i, variable = variable_i,
      metric = "cond_pip", value = cond_pip, given = given_j, group_given = group_j
    )

  bind_rows(var_rows, grp_rows, cond_rows)
}

# ---------------------------
# 6) PLOTTING
# ---------------------------
plot_univar_compare_all <- function(univar_df) {
  df <- univar_df %>%
    mutate(method = recode(method, naive = "Naive", weighted = "Design-aware weighted", weighted_trimmed = "Weighted + trimmed"))

  ggplot(df, aes(x = z, y = est, color = method, fill = method)) +
    geom_ribbon(aes(ymin = est - 1.96 * se, ymax = est + 1.96 * se), alpha = 0.10, color = NA) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ variable, scales = "free_x", ncol = 5) +
    labs(
      x = "Standardized exposure (z-score, clipped to [-3, 3])",
      y = "Estimated h(z)",
      color = "Analysis",
      fill = "Analysis",
      title = "Univariate exposure-response functions",
      subtitle = "Shaded bands are within-fit 95% posterior credible intervals"
    ) +
    theme_minimal(base_size = 12)
}

plot_overall_compare <- function(overall_df) {
  df <- overall_df %>%
    mutate(method = recode(method, naive = "Naive", weighted = "Design-aware weighted", weighted_trimmed = "Weighted + trimmed"))

  ggplot(df, aes(x = quantile, y = est, color = method, fill = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(aes(ymin = est - 1.96 * sd, ymax = est + 1.96 * sd), alpha = 0.08, color = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ method, ncol = 1) +
    labs(
      x = "Mixture quantile",
      y = "Overall mixture effect",
      title = "Overall mixture effect curves",
      subtitle = "Bands are within-fit 95% posterior credible intervals; design-bootstrap intervals are reported separately in tables"
    ) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
}

plot_singvar_compare <- function(singvar_df) {
  df <- singvar_df %>%
    mutate(method = recode(method, naive = "Naive", weighted = "Design-aware weighted", weighted_trimmed = "Weighted + trimmed"))

  ggplot(df, aes(x = variable, y = est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_pointrange(position = position_dodge(width = 0.65)) +
    coord_flip() +
    facet_wrap(~ q.fixed, ncol = 1,
               labeller = labeller(q.fixed = function(x) paste0("q.fixed = ", x))) +
    labs(
      x = "",
      y = "Change in outcome (75th vs 25th percentile)",
      color = "Analysis",
      title = "Single-variable risk summaries"
    ) +
    theme_minimal(base_size = 12)
}

plot_bivar_levels_compare <- function(bivar_levels_df, vars_show = NULL) {
  df <- bivar_levels_df %>%
    mutate(method = recode(method, naive = "Naive", weighted = "Design-aware weighted", weighted_trimmed = "Weighted + trimmed"))

  if (!is.null(vars_show)) {
    df <- df %>% filter(variable1 %in% vars_show, variable2 %in% vars_show)
    df$variable1 <- factor(df$variable1, levels = vars_show)
    df$variable2 <- factor(df$variable2, levels = vars_show)
  }

  df$quantile <- factor(as.character(df$quantile), levels = c("0.25","0.5","0.75"))

  ggplot(df, aes(x = z1, y = est, color = quantile, linetype = method)) +
    geom_line(linewidth = 0.65) +
    facet_grid(variable2 ~ variable1, switch = "both") +
    labs(
      title = "Bivariate exposure-response levels",
      subtitle = "Method separated by linetype to reduce overlap",
      x = "Standardized first exposure",
      y = "Estimated h(z1 | quantiles of z2)",
      linetype = "Analysis",
      color = "Quantile of second exposure"
    ) +
    theme_minimal(base_size = 11) +
    theme(panel.spacing = unit(0.50, "lines"))
}

plot_pips_variable_compare <- function(pip_master_df) {
  df <- pip_master_df %>%
    filter(level == "variable", metric == "pip") %>%
    mutate(method = recode(method, naive = "Naive", weighted = "Design-aware weighted", weighted_trimmed = "Weighted + trimmed"))

  ggplot(df, aes(x = reorder(variable, value), y = value, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    coord_flip() +
    labs(
      x = "",
      y = "Posterior Inclusion Probability",
      fill = "Analysis",
      title = "Variable-level posterior inclusion probabilities"
    ) +
    theme_minimal(base_size = 12)
}

# ---------------------------
# 7) MAIN ANALYSIS FOR ONE DATASET
# ---------------------------
run_one_dataset_analysis <- function(prep, dataset_id = 1L) {
  dat <- prep$std
  covar_vars <- prep$covar_vars

  design_obj <- make_survey_design(dat, psu_col, strata_col, wt_col)

  runtime_tbl <- list()
  diag_tbl <- list()

  t_naive <- system.time({
    fit_naive_raw <- fit_kmbayes_fast(dat, y_var, exposure_vars, covar_vars,
                                      iter = ITER_MAIN, varsel = VARSEL)
    fit_naive <- thin_fit(fit_naive_raw, burn = BURN_MAIN, thin = THIN_MAIN)
  })
  runtime_tbl[[length(runtime_tbl) + 1L]] <- tibble(
    dataset_id = dataset_id, analysis = "naive_singlefit",
    elapsed_sec = unname(t_naive["elapsed"])
  )
  diag_tbl[[length(diag_tbl) + 1L]] <- compute_fit_diagnostics(fit_naive, "naive", rep_id = 0L) %>% mutate(dataset_id = dataset_id)

  t_weighted <- system.time({
    rs_once <- fit_kmbayes_resampled(
      dat, y_var, exposure_vars, covar_vars, psu_col, strata_col, wt_col,
      iter = ITER_MAIN, burn = BURN_MAIN, thin = THIN_MAIN, varsel = VARSEL,
      trim = FALSE
    )
    fit_weighted <- rs_once$fit
    dat_weighted <- rs_once$dat_rs
  })
  runtime_tbl[[length(runtime_tbl) + 1L]] <- tibble(
    dataset_id = dataset_id, analysis = "weighted_singlefit",
    elapsed_sec = unname(t_weighted["elapsed"])
  )
  diag_tbl[[length(diag_tbl) + 1L]] <- compute_fit_diagnostics(fit_weighted, "weighted", rep_id = 0L) %>% mutate(dataset_id = dataset_id)

  if (RUN_TRIMMING_SENSITIVITY) {
    t_weighted_trim <- system.time({
      rs_trim <- fit_kmbayes_resampled(
        dat, y_var, exposure_vars, covar_vars, psu_col, strata_col, wt_col,
        iter = ITER_MAIN, burn = BURN_MAIN, thin = THIN_MAIN, varsel = VARSEL,
        trim = TRUE,
        trim_lower_prob = TRIM_LOWER_PROB,
        trim_upper_prob = TRIM_UPPER_PROB
      )
      fit_weighted_trim <- rs_trim$fit
      dat_weighted_trim <- rs_trim$dat_rs
    })
    runtime_tbl[[length(runtime_tbl) + 1L]] <- tibble(
      dataset_id = dataset_id, analysis = "weighted_trimmed_singlefit",
      elapsed_sec = unname(t_weighted_trim["elapsed"])
    )
    diag_tbl[[length(diag_tbl) + 1L]] <- compute_fit_diagnostics(fit_weighted_trim, "weighted_trimmed", rep_id = 0L) %>% mutate(dataset_id = dataset_id)
  } else {
    fit_weighted_trim <- NULL
    dat_weighted_trim <- NULL
  }

  summ_naive <- make_bkmr_summaries(fit_naive, as.matrix(dat[, exposure_vars, drop = FALSE]))
  summ_weighted <- make_bkmr_summaries(fit_weighted, as.matrix(dat_weighted[, exposure_vars, drop = FALSE]))
  if (!is.null(fit_weighted_trim)) {
    summ_weighted_trim <- make_bkmr_summaries(fit_weighted_trim, as.matrix(dat_weighted_trim[, exposure_vars, drop = FALSE]))
  } else {
    summ_weighted_trim <- NULL
  }

  analysis_objects <- list(
    naive = list(fit = fit_naive, dat = dat, summ = summ_naive),
    weighted = list(fit = fit_weighted, dat = dat_weighted, summ = summ_weighted)
  )
  if (!is.null(fit_weighted_trim)) {
    analysis_objects$weighted_trimmed <- list(fit = fit_weighted_trim, dat = dat_weighted_trim, summ = summ_weighted_trim)
  }

  # Main-fit deltas
  delta_tbl <- bind_rows(lapply(names(analysis_objects), function(nm) {
    dd <- extract_overall_delta(analysis_objects[[nm]]$fit)
    tibble(
      dataset_id = dataset_id,
      method = nm,
      estimate = unname(dd["mean"]),
      lwr = unname(dd["lwr"]),
      upr = unname(dd["upr"]),
      interval_type = "posterior_credible_interval"
    )
  }))

  # PIPs
  pip_master <- bind_rows(lapply(names(analysis_objects), function(nm) {
    build_pip_master(analysis_objects[[nm]]$fit, nm, exposure_vars, pip_groups)
  })) %>%
    mutate(dataset_id = dataset_id)

  # Bootstrap design-aware summaries
  run_one_boot <- function(b, trim_flag = FALSE, method_label = "weighted") {
    fit_obj <- fit_kmbayes_resampled(
      dat, y_var, exposure_vars, covar_vars,
      psu_col, strata_col, wt_col,
      iter = ITER_BOOT, burn = BURN_BOOT, thin = THIN_BOOT, varsel = VARSEL,
      trim = trim_flag,
      trim_lower_prob = TRIM_LOWER_PROB,
      trim_upper_prob = TRIM_UPPER_PROB
    )
    fit_b <- fit_obj$fit
    dd <- extract_overall_delta(fit_b)
    pips <- get_flat_pips(fit_b, exposure_vars) %>%
      attach_groups_to_pips(pip_groups) %>%
      mutate(dataset_id = dataset_id, rep_id = b, method = method_label)
    diags <- compute_fit_diagnostics(fit_b, method_label, rep_id = b) %>% mutate(dataset_id = dataset_id)

    list(
      overall = tibble(
        dataset_id = dataset_id,
        rep_id = b,
        method = method_label,
        overall_boot = unname(dd["mean"])
      ),
      pips = pips,
      diags = diags
    )
  }

  t_boot <- system.time({
    boot_res <- future_lapply(seq_len(R_BOOT), run_one_boot,
                              trim_flag = FALSE, method_label = "weighted")
  })
  runtime_tbl[[length(runtime_tbl) + 1L]] <- tibble(
    dataset_id = dataset_id,
    analysis = "weighted_bootstrap_total",
    elapsed_sec = unname(t_boot["elapsed"])
  )

  boot_overall <- bind_rows(lapply(boot_res, `[[`, "overall"))
  boot_pips <- bind_rows(lapply(boot_res, `[[`, "pips"))
  boot_diags <- bind_rows(lapply(boot_res, `[[`, "diags"))

  if (RUN_TRIMMING_SENSITIVITY) {
    t_boot_trim <- system.time({
      boot_res_trim <- future_lapply(seq_len(R_BOOT), run_one_boot,
                                     trim_flag = TRUE, method_label = "weighted_trimmed")
    })
    runtime_tbl[[length(runtime_tbl) + 1L]] <- tibble(
      dataset_id = dataset_id,
      analysis = "weighted_trimmed_bootstrap_total",
      elapsed_sec = unname(t_boot_trim["elapsed"])
    )

    boot_overall_trim <- bind_rows(lapply(boot_res_trim, `[[`, "overall"))
    boot_pips_trim <- bind_rows(lapply(boot_res_trim, `[[`, "pips"))
    boot_diags_trim <- bind_rows(lapply(boot_res_trim, `[[`, "diags"))
  } else {
    boot_overall_trim <- tibble()
    boot_pips_trim <- tibble()
    boot_diags_trim <- tibble()
  }

  diag_all <- bind_rows(diag_tbl, list(boot_diags), list(boot_diags_trim))
  runtime_all <- bind_rows(runtime_tbl)

  list(
    dataset_id = dataset_id,
    prep = prep,
    design_obj = design_obj,
    analysis_objects = analysis_objects,
    delta_tbl = delta_tbl,
    pip_master = pip_master,
    runtime = runtime_all,
    diagnostics = diag_all,
    boot_overall = bind_rows(boot_overall, boot_overall_trim),
    boot_pips = bind_rows(boot_pips, boot_pips_trim)
  )
}

# ---------------------------
# 8) RUN ANALYSIS ACROSS DATASETS / IMPUTATIONS
# ---------------------------
analysis_files <- get_analysis_files(infile, imputed_files)

prepared_list <- lapply(analysis_files, function(f) {
  prepare_one_dataset(
    path = f,
    y_var = y_var,
    covar_vars_base = covar_vars_base,
    statin_var = statin_var,
    exposure_vars = exposure_vars,
    psu_col = psu_col,
    strata_col = strata_col,
    wt_col = wt_col,
    clip_lower = CLIP_LOWER,
    clip_upper = CLIP_UPPER
  )
})

qc_all <- bind_rows(lapply(seq_along(prepared_list), function(i) {
  prepared_list[[i]]$qc %>% mutate(dataset_id = i, file = basename(prepared_list[[i]]$file))
}))
write.csv(qc_all,
          file.path(outdir_real, "QC_standardization_range_minus3_plus3.csv"),
          row.names = FALSE)

results_list <- lapply(seq_along(prepared_list), function(i) {
  message(sprintf("Running dataset %d of %d: %s", i, length(prepared_list), basename(prepared_list[[i]]$file)))
  run_one_dataset_analysis(prepared_list[[i]], dataset_id = i)
})

# ---------------------------
# 9) COMBINE FIT-LEVEL SUMMARIES
# ---------------------------
get_summary_df <- function(results_list, component_name, method_name) {
  bind_rows(lapply(results_list, function(res) {
    obj <- res$analysis_objects[[method_name]]
    if (is.null(obj)) return(NULL)
    obj$summ[[component_name]] %>% mutate(dataset_id = res$dataset_id, method = method_name)
  }))
}

methods_present <- unique(unlist(lapply(results_list, function(res) names(res$analysis_objects))))

univar_all_df <- bind_rows(lapply(methods_present, function(m) get_summary_df(results_list, "univar", m)))
overall_all_df <- bind_rows(lapply(methods_present, function(m) get_summary_df(results_list, "overall", m)))
singvar_all_df <- bind_rows(lapply(methods_present, function(m) get_summary_df(results_list, "singvar", m)))
bivar_levels_all_df <- bind_rows(lapply(methods_present, function(m) get_summary_df(results_list, "bivar_levels", m)))
if (SAVE_BIVAR_SURFACE) {
  bivar_surface_all_df <- bind_rows(lapply(methods_present, function(m) get_summary_df(results_list, "bivar", m)))
}

# Pooled or averaged figure datasets
if (length(results_list) > 1) {
  univar_save <- univar_all_df %>%
    group_by(method, variable, z) %>%
    summarise(est = mean(est, na.rm = TRUE),
              se = sqrt(mean(se^2, na.rm = TRUE)),
              .groups = "drop")

  overall_save <- overall_all_df %>%
    group_by(method, quantile) %>%
    summarise(est = mean(est, na.rm = TRUE),
              sd = sqrt(mean(sd^2, na.rm = TRUE)),
              .groups = "drop")

  singvar_save <- singvar_all_df %>%
    group_by(method, variable, q.fixed) %>%
    summarise(est = mean(est, na.rm = TRUE),
              sd = sqrt(mean(sd^2, na.rm = TRUE)),
              .groups = "drop")

  bivar_levels_save <- bivar_levels_all_df %>%
    group_by(method, variable1, variable2, quantile, z1) %>%
    summarise(est = mean(est, na.rm = TRUE), .groups = "drop")
} else {
  univar_save <- univar_all_df %>% select(-dataset_id)
  overall_save <- overall_all_df %>% select(-dataset_id)
  singvar_save <- singvar_all_df %>% select(-dataset_id)
  bivar_levels_save <- bivar_levels_all_df %>% select(-dataset_id)
}

write.csv(univar_save, file.path(outdir_real, "DATA_univar_methods_STD_clip_m3_p3.csv"), row.names = FALSE)
write.csv(overall_save, file.path(outdir_real, "DATA_overallrisk_methods_STD_clip_m3_p3.csv"), row.names = FALSE)
write.csv(singvar_save, file.path(outdir_real, "DATA_singvarrisk_methods_STD_clip_m3_p3.csv"), row.names = FALSE)
write.csv(bivar_levels_save, file.path(outdir_real, "DATA_bivar_levels_methods_STD_clip_m3_p3.csv"), row.names = FALSE)
if (SAVE_BIVAR_SURFACE) {
  write.csv(bivar_surface_all_df, file.path(outdir_real, "DATA_bivar_surface_methods_STD_clip_m3_p3.csv"), row.names = FALSE)
}

# ---------------------------
# 10) PLOTS
# ---------------------------
p_univar <- plot_univar_compare_all(univar_save)
p_overall <- plot_overall_compare(overall_save)
p_singvar <- plot_singvar_compare(singvar_save)
p_bivar_levels <- plot_bivar_levels_compare(bivar_levels_save, vars_show = VARS_SUBSET_FOR_BIVAR_PLOT)

print(p_univar)
print(p_overall)
print(p_singvar)
print(p_bivar_levels)

ggsave(file.path(outdir_real, "FIG_univar_methods_STD_clip_m3_p3.png"), p_univar, width = 12, height = 8, dpi = 300)
ggsave(file.path(outdir_real, "FIG_overall_methods_STD_clip_m3_p3.png"), p_overall, width = 10, height = 9, dpi = 300)
ggsave(file.path(outdir_real, "FIG_singvar_methods_STD_clip_m3_p3.png"), p_singvar, width = 10, height = 8, dpi = 300)
ggsave(file.path(outdir_real, "FIG_bivar_levels_methods_STD_clip_m3_p3.png"), p_bivar_levels, width = 15, height = 10, dpi = 300)

# ---------------------------
# 11) OVERALL DELTA TABLES
# ---------------------------
delta_all <- bind_rows(lapply(results_list, `[[`, "delta_tbl"))

# Design-bootstrap intervals
boot_overall_all <- bind_rows(lapply(results_list, `[[`, "boot_overall"))
write.csv(boot_overall_all,
          file.path(outdir_real, "DATA_overall_bootstrap_estimates_methods_STD_clip_m3_p3.csv"),
          row.names = FALSE)

boot_perf_summary <- boot_overall_all %>%
  group_by(method) %>%
  summarise(
    design_boot_point = mean(overall_boot, na.rm = TRUE),
    design_boot_ci_lower = quantile(overall_boot, 0.025, na.rm = TRUE),
    design_boot_ci_upper = quantile(overall_boot, 0.975, na.rm = TRUE),
    n_boot = n(),
    interval_type = "bootstrap_replication_interval",
    .groups = "drop"
  )

if (nrow(delta_all) > 0) {
  if (length(results_list) > 1) {
    pooled_main <- delta_all %>%
      mutate(se2 = ((upr - lwr) / (2 * 1.96))^2) %>%
      group_by(method, interval_type) %>%
      summarise(
        qhat = list(estimate),
        uhat = list(se2),
        .groups = "drop"
      ) %>%
      rowwise() %>%
      mutate(pool = list(pool_scalar_rubin(unlist(qhat), unlist(uhat)))) %>%
      tidyr::unnest(pool) %>%
      select(method, interval_type, m, estimate, se, lwr, upr, within_var, between_var, total_var)

    write.csv(pooled_main,
              file.path(outdir_real, "TABLE_overall_delta_mainfit_pooled_MI.csv"),
              row.names = FALSE)
  } else {
    write.csv(delta_all,
              file.path(outdir_real, "TABLE_overall_delta_mainfit_single_dataset.csv"),
              row.names = FALSE)
  }
}

write.csv(boot_perf_summary,
          file.path(outdir_real, "TABLE_performance_summary_designaware_bootstrap.csv"),
          row.names = FALSE)

# ---------------------------
# 12) PIP OUTPUTS + STABILITY
# ---------------------------
pip_master_all <- bind_rows(lapply(results_list, `[[`, "pip_master")) %>%
  arrange(dataset_id, method, factor(level, levels = c("group","variable","conditional")), group, variable, given)

write.csv(pip_master_all,
          file.path(outdir_real, "DATA_PIPs_MASTER_group_variable_conditional_methods_STD_clip_m3_p3.csv"),
          row.names = FALSE)

boot_pips_all <- bind_rows(lapply(results_list, `[[`, "boot_pips"))

if (SAVE_BOOTSTRAP_PIP_DETAILS && nrow(boot_pips_all) > 0) {
  write.csv(boot_pips_all,
            file.path(outdir_real, "DATA_bootstrap_variable_PIPs_methods.csv"),
            row.names = FALSE)
}

pip_stability <- boot_pips_all %>%
  group_by(method, variable, group) %>%
  summarise(
    mean_pip = mean(pip, na.rm = TRUE),
    sd_pip = sd(pip, na.rm = TRUE),
    median_pip = median(pip, na.rm = TRUE),
    q25_pip = quantile(pip, 0.25, na.rm = TRUE),
    q75_pip = quantile(pip, 0.75, na.rm = TRUE),
    min_pip = min(pip, na.rm = TRUE),
    max_pip = max(pip, na.rm = TRUE),
    prop_pip_ge_0_5 = mean(pip >= 0.5, na.rm = TRUE),
    n_boot = n(),
    .groups = "drop"
  ) %>%
  arrange(method, desc(mean_pip))

write.csv(pip_stability,
          file.path(outdir_real, "TABLE_PIP_stability_across_bootstrap_replicates.csv"),
          row.names = FALSE)

p_pips <- plot_pips_variable_compare(
  pip_master_all %>% filter(dataset_id == min(dataset_id))
)
print(p_pips)
ggsave(file.path(outdir_real, "FIG_variable_PIPs_methods.png"), p_pips, width = 10, height = 8, dpi = 300)

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
    median_ess = median(ess, na.rm = TRUE),
    min_ess = min(ess, na.rm = TRUE),
    max_abs_geweke_z = max(abs(geweke_z), na.rm = TRUE),
    n_records = n(),
    .groups = "drop"
  )
write.csv(diagnostics_summary,
          file.path(outdir_real, "TABLE_MCMC_diagnostics_summary.csv"),
          row.names = FALSE)

runtime_all <- bind_rows(lapply(results_list, `[[`, "runtime"))
write.csv(runtime_all,
          file.path(outdir_real, "TABLE_runtime_all.csv"),
          row.names = FALSE)

runtime_summary <- runtime_all %>%
  group_by(analysis) %>%
  summarise(
    mean_elapsed_sec = mean(elapsed_sec, na.rm = TRUE),
    median_elapsed_sec = median(elapsed_sec, na.rm = TRUE),
    .groups = "drop"
  )

naive_time <- runtime_summary %>% filter(analysis == "naive_singlefit") %>% pull(mean_elapsed_sec)
weighted_time <- runtime_summary %>% filter(analysis == "weighted_singlefit") %>% pull(mean_elapsed_sec)
weighted_boot_time <- runtime_summary %>% filter(analysis == "weighted_bootstrap_total") %>% pull(mean_elapsed_sec)

runtime_comparison <- tibble(
  comparison = c("weighted_singlefit_vs_naive", "weighted_bootstrap_total_vs_naive"),
  ratio = c(
    ifelse(length(naive_time) == 1 && length(weighted_time) == 1, weighted_time / naive_time, NA_real_),
    ifelse(length(naive_time) == 1 && length(weighted_boot_time) == 1, weighted_boot_time / naive_time, NA_real_)
  )
)

write.csv(runtime_summary,
          file.path(outdir_real, "TABLE_runtime_summary.csv"),
          row.names = FALSE)
write.csv(runtime_comparison,
          file.path(outdir_real, "TABLE_runtime_increase_vs_naive.csv"),
          row.names = FALSE)

# ---------------------------
# 14) WEIGHT TRIMMING SENSITIVITY
# ---------------------------
trimming_summary <- boot_overall_all %>%
  group_by(method) %>%
  summarise(
    mean_overall = mean(overall_boot, na.rm = TRUE),
    sd_overall = sd(overall_boot, na.rm = TRUE),
    q025 = quantile(overall_boot, 0.025, na.rm = TRUE),
    q975 = quantile(overall_boot, 0.975, na.rm = TRUE),
    .groups = "drop"
  )
write.csv(trimming_summary,
          file.path(outdir_real, "TABLE_weight_trimming_sensitivity.csv"),
          row.names = FALSE)

# ---------------------------
# 15) INFERENCE TARGET + SOFTWARE / REPRODUCIBILITY METADATA
# ---------------------------
software_tbl <- tibble(
  package = c("R", "bkmr", "survey", "future.apply", "coda", "dplyr", "ggplot2", "readr", "tidyr", "purrr", "stringr"),
  version = c(
    as.character(getRversion()),
    as.character(packageVersion("bkmr")),
    as.character(packageVersion("survey")),
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

analysis_metadata <- tibble(
  field = c(
    "inference_target",
    "single_fit_interval_interpretation",
    "design_bootstrap_interval_interpretation",
    "uses_multiple_imputation",
    "n_imputed_datasets",
    "statin_adjustment_included",
    "statin_variable",
    "weight_trimming_sensitivity_run",
    "trim_lower_prob",
    "trim_upper_prob",
    "zscore_clipping_lower",
    "zscore_clipping_upper",
    "bootstrap_replicates",
    "main_iter",
    "main_burn",
    "main_thin",
    "boot_iter",
    "boot_burn",
    "boot_thin"
  ),
  value = c(
    INFERENCE_TARGET,
    "posterior_credible_interval",
    "bootstrap_replication_interval",
    ifelse(length(analysis_files) > 1, "yes", "no"),
    length(analysis_files),
    ifelse(is.null(statin_var), "no", "yes"),
    ifelse(is.null(statin_var), "", statin_var),
    ifelse(RUN_TRIMMING_SENSITIVITY, "yes", "no"),
    TRIM_LOWER_PROB,
    TRIM_UPPER_PROB,
    CLIP_LOWER,
    CLIP_UPPER,
    R_BOOT,
    ITER_MAIN,
    BURN_MAIN,
    THIN_MAIN,
    ITER_BOOT,
    BURN_BOOT,
    THIN_BOOT
  )
)
write.csv(analysis_metadata,
          file.path(outdir_real, "TABLE_analysis_metadata.csv"),
          row.names = FALSE)

# Appendix-ready code note
appendix_code_note <- c(
  "This script implements the reviewer-responsive NHANES BKMR workflow.",
  "Single-fit BKMR intervals are posterior credible intervals.",
  "Design-aware uncertainty is summarized with bootstrap / replication-based intervals across PSU resamples.",
  "Optional statin adjustment is controlled through statin_var.",
  "Optional multiple-imputation support is controlled through imputed_files."
)
writeLines(appendix_code_note,
           con = file.path(outdir_real, "APPENDIX_code_notes.txt"))

# ---------------------------
# 16) FINAL REPORT
# ---------------------------
cat("\n===== REVIEWER-RESPONSIVE NHANES BKMR ANALYSIS COMPLETE =====\n")
cat("Saved to:\n", outdir_real, "\n\n")
cat("Key outputs:\n")
cat("- QC_standardization_range_minus3_plus3.csv\n")
cat("- DATA_univar_methods_STD_clip_m3_p3.csv\n")
cat("- DATA_overallrisk_methods_STD_clip_m3_p3.csv\n")
cat("- DATA_singvarrisk_methods_STD_clip_m3_p3.csv\n")
cat("- DATA_bivar_levels_methods_STD_clip_m3_p3.csv\n")
cat("- DATA_PIPs_MASTER_group_variable_conditional_methods_STD_clip_m3_p3.csv\n")
cat("- DATA_overall_bootstrap_estimates_methods_STD_clip_m3_p3.csv\n")
cat("- TABLE_performance_summary_designaware_bootstrap.csv\n")
cat("- TABLE_MCMC_diagnostics_all.csv\n")
cat("- TABLE_MCMC_diagnostics_summary.csv\n")
cat("- TABLE_runtime_summary.csv\n")
cat("- TABLE_runtime_increase_vs_naive.csv\n")
cat("- TABLE_weight_trimming_sensitivity.csv\n")
cat("- TABLE_PIP_stability_across_bootstrap_replicates.csv\n")
cat("- TABLE_analysis_metadata.csv\n")
cat("- TABLE_software_versions.csv\n")
cat("- APPENDIX_code_notes.txt\n")

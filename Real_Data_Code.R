#####STANDARDIZED FIGURES (CLIPPED -3 to +3)###########
# ======================================================
# NHANES REAL DATA (TC.csv) + BKMR
# Naive vs Weighted (PSU bootstrap resampling) + STANDARDIZED exposures
# IMPORTANT: Standardization = z-score, then CLIP to [-3, +3]

# Saves (CSV + PNG):
#   - Figure datasets: univar, overall, singvar, bivar_levels (+ optional bivar_surface)
#   - Figures: univar, overall, singvar, bivar_levels
#   - Comparison table: naive vs weighted single-fit overall delta
#   - Design-aware overall bootstrap estimates + perf summary
#   - MASTER PIP dataset (group + variable + conditional) for naive & weighted
#   - QC: standardization range check confirming [-3, +3]
#
# Speed controls:
#   - ITER/BURN/THIN for main fits + bootstrap fits
#   - Parallel PSU bootstrap via future.apply
#   - Optional bivar surface OFF by default (very heavy)
# ======================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(bkmr)
  library(survey)
  library(future.apply)
})

# ---------------------------
# 0) USER SETTINGS
# ---------------------------
infile      <- "/Users/doreennoldajehu-appiah/Desktop/PAPER_3/PIP3/RealData_Sim/TC.csv"
outdir_real <- "/Users/doreennoldajehu-appiah/Desktop/PAPER_3/PIP3/RealData_Sim"
dir.create(outdir_real, showWarnings = FALSE, recursive = TRUE)

y_var <- "TotalCholesterol"
covar_vars <- c("age", "Ethnicity", "income", "Gender", "BMI", "smq020", "alq101")

# Design columns (must match TC.csv)
psu_col    <- "sdmvpsu"
strata_col <- "sdmvstra"
wt_col     <- "wtmec2yr"

# Exposures (must match TC.csv exactly)
exposure_vars <- c(
  # Metals
  "Arsenic","Mercury","Barium","Cobalt","Cesium",
  "Molybdenum","Lead","Antimony","Tin","Strontium",
  "Thallium","Tungsten","Uranium",
  
  # PFAS
  "PFOA","PFOS","PFDE","PFHS","NMeFOSAA",
  "PFHP","PFNA","PFUA","PFDO",
  
  # Phthalates / Plasticizers
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

# Group definitions for group PIPs
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

ITER_MAIN <- 2500
BURN_MAIN <- 1000
THIN_MAIN <- 2
VARSEL    <- TRUE

QS_OVERALL  <- seq(0.25, 0.75, by = 0.05)
QS_BIVARLVL <- c(0.25, 0.50, 0.75)

R_BOOT    <- 30
ITER_BOOT <- 1500
BURN_BOOT <- 600
THIN_BOOT <- 2

# Parallel bootstrap workers
N_WORKERS <- max(1, parallel::detectCores() - 1)
future::plan(future::multisession, workers = N_WORKERS)

# Heavy outputs
SAVE_BIVAR_SURFACE <- FALSE

# Do NOT limit variables shown in bivar plot: keep as exposure_vars
VARS_SUBSET_FOR_BIVAR_PLOT <- exposure_vars

# Standardization clipping bounds
CLIP_LOWER <- -3
CLIP_UPPER <-  3

# ---------------------------
# 1) LOAD + VALIDATE REQUIRED COLUMNS
# ---------------------------
nhanes_raw <- read_csv(infile, show_col_types = FALSE) %>% as.data.frame()

need_cols <- unique(c(y_var, covar_vars, exposure_vars, psu_col, strata_col, wt_col))
miss <- setdiff(need_cols, names(nhanes_raw))
if (length(miss) > 0) stop("Missing required columns in TC.csv: ", paste(miss, collapse = ", "))

# Ensure exposures are numeric (prevents silent character issues)
nhanes_raw[exposure_vars] <- lapply(nhanes_raw[exposure_vars], function(x) suppressWarnings(as.numeric(x)))

# ---------------------------
# 2) FILTER COMPLETE CASES
# ---------------------------
nhanes_cc <- nhanes_raw %>%
  dplyr::filter(
    !is.na(.data[[y_var]]),
    !if_any(all_of(exposure_vars), is.na),
    !if_any(all_of(covar_vars), is.na),
    !is.na(.data[[psu_col]]),
    !is.na(.data[[strata_col]]),
    !is.na(.data[[wt_col]])
  ) %>% as.data.frame()

# ---------------------------
# 3) STANDARDIZE EXPOSURES (z-score) + CLIP TO [-3, +3]
#    Same mu/sd used for naive and resampled analyses
# ---------------------------
mu_vec <- sapply(exposure_vars, function(v) mean(nhanes_cc[[v]], na.rm = TRUE))
sd_vec <- sapply(exposure_vars, function(v) sd(nhanes_cc[[v]], na.rm = TRUE))
sd_vec[is.na(sd_vec) | sd_vec == 0] <- 1

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

nhanes_std <- standardize_and_clip(
  nhanes_cc, exposure_vars, mu_vec, sd_vec,
  lower = CLIP_LOWER, upper = CLIP_UPPER
)

# QC: confirm hard bounds [-3, +3]
range_check <- data.frame(
  variable = exposure_vars,
  min_z = sapply(exposure_vars, \(v) min(nhanes_std[[v]], na.rm = TRUE)),
  max_z = sapply(exposure_vars, \(v) max(nhanes_std[[v]], na.rm = TRUE))
)
print(range_check)
stopifnot(all(range_check$min_z >= CLIP_LOWER), all(range_check$max_z <= CLIP_UPPER))

write.csv(range_check,
          file.path(outdir_real, "QC_standardization_range_minus3_plus3.csv"),
          row.names = FALSE)

# ---------------------------
# 4) SURVEY DESIGN OBJECT (diagnostics only; BKMR itself is unweighted)
# ---------------------------
nhanes_design <- svydesign(
  id      = as.formula(paste0("~", psu_col)),
  strata  = as.formula(paste0("~", strata_col)),
  weights = as.formula(paste0("~", wt_col)),
  nest    = TRUE,
  data    = nhanes_std
)

# ---------------------------
# 5) PSU BOOTSTRAP RESAMPLING (NHANES-like)
# ---------------------------
resample_nhanes_psu <- function(dat, psu_col, strata_col, wt_col) {
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
      summarise(W_psu = sum(.data[[wt_col]], na.rm = TRUE), .groups = "drop") %>%
      mutate(p = W_psu / sum(W_psu))
    
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
      
      p_unit <- dpsu[[wt_col]] / sum(dpsu[[wt_col]])
      take <- sample(seq_len(m), size = m, replace = TRUE, prob = p_unit)
      rows_s[[j]] <- dpsu[take, , drop = FALSE]
    }
    
    k <- k + 1L
    out_list[[k]] <- bind_rows(rows_s)
  }
  
  bind_rows(out_list) %>% as.data.frame()
}

# ---------------------------
# 6) BKMR FIT HELPERS (burn-in + thinning)
# ---------------------------
fit_kmbayes_fast <- function(dat, y_var, exposure_vars, covar_vars,
                             iter, varsel, verbose = FALSE) {
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
                                  iter, burn, thin, varsel) {
  dat_rs <- resample_nhanes_psu(dat, psu_col, strata_col, wt_col)
  
  # IMPORTANT: dat is already standardized+clipped. Resampling preserves that.
  fit_rs_raw <- fit_kmbayes_fast(dat_rs, y_var, exposure_vars, covar_vars, iter = iter, varsel = varsel, verbose = FALSE)
  fit_rs <- thin_fit(fit_rs_raw, burn = burn, thin = thin)
  list(fit = fit_rs, dat_rs = dat_rs)
}

# ---------------------------
# 7) BKMR SUMMARIES (DATA FOR FIGURES)
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
    fit = fit, qs.diff = c(0.25, 0.75),
    q.fixed = c(0.25, 0.50, 0.75),
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
  or <- OverallRiskSummaries(fit = fit, qs = c(0.25, 0.75), q.fixed = 0.5, method = "approx")
  est <- with(or, est[quantile == 0.75] - est[quantile == 0.25])
  sdv <- with(or, sqrt(sd[quantile == 0.75]^2 + sd[quantile == 0.25]^2))
  c(mean = est, lwr = est - 1.96 * sdv, upr = est + 1.96 * sdv)
}

# ---------------------------
# 8) PIP HELPERS (variable + group + conditional) -> MASTER
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
    stop("ExtractPIPs returned data.frame but couldn't detect pip/variable columns. Try: str(ExtractPIPs(fit))")
  }
  
  if (is.list(p)) {
    cand <- c("pip", "PIP", "pips", "PIPs", "p")
    pip_name <- cand[cand %in% names(p)][1]
    if (!is.na(pip_name)) {
      vec <- p[[pip_name]]
      if (!is.numeric(vec)) stop("ExtractPIPs list element found but not numeric. Try: str(ExtractPIPs(fit))")
      if (!is.null(names(vec)) && all(exposure_vars %in% names(vec))) {
        vec <- vec[exposure_vars]
      } else {
        vec <- vec[seq_along(exposure_vars)]
        names(vec) <- exposure_vars
      }
      return(tibble(variable = exposure_vars, pip = as.numeric(vec)))
    }
    stop("ExtractPIPs returned list but no PIP vector found. Try: str(ExtractPIPs(fit))")
  }
  
  stop("Unsupported ExtractPIPs output type. Try: str(ExtractPIPs(fit))")
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
  denom <- colMeans(delta_mat == 1)  # P(j in)
  
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
  if (ncol(delta) != length(exposure_vars)) {
    warning("delta matrix columns != length(exposure_vars). Assuming ordering matches Z columns.")
  }
  
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
      method,
      level = "variable",
      group = group,
      variable = variable,
      metric = "pip",
      value = pip,
      given = NA_character_,
      group_given = NA_character_
    )
  
  grp_rows <- grp %>%
    pivot_longer(c(n_vars, pip_any, pip_mean, pip_max),
                 names_to = "metric", values_to = "value") %>%
    transmute(
      method,
      level = "group",
      group = group,
      variable = NA_character_,
      metric = metric,
      value = as.numeric(value),
      given = NA_character_,
      group_given = NA_character_
    )
  
  cond_rows <- cond_long %>%
    transmute(
      method,
      level = "conditional",
      group = group_i,
      variable = variable_i,
      metric = "cond_pip",
      value = cond_pip,
      given = given_j,
      group_given = group_j
    )
  
  bind_rows(var_rows, grp_rows, cond_rows)
}

# ---------------------------
# 9) FIT NAIVE + SINGLE WEIGHTED
# ---------------------------
message("Fitting naive BKMR ...")
fit_naive_raw <- fit_kmbayes_fast(nhanes_std, y_var, exposure_vars, covar_vars, iter = ITER_MAIN, varsel = VARSEL)
fit_naive <- thin_fit(fit_naive_raw, burn = BURN_MAIN, thin = THIN_MAIN)

message("Fitting weighted BKMR (single PSU bootstrap resample) ...")
rs_once <- fit_kmbayes_resampled(
  nhanes_std, y_var, exposure_vars, covar_vars,
  psu_col, strata_col, wt_col,
  iter = ITER_MAIN, burn = BURN_MAIN, thin = THIN_MAIN, varsel = VARSEL
)
fit_weighted <- rs_once$fit
dat_weighted <- rs_once$dat_rs

# ---------------------------
# 10) SUMMARIES (DATA FOR FIGURES) + SAVE CSV
# ---------------------------
message("Computing summaries (naive) ...")
summ_naive <- make_bkmr_summaries(
  fit_naive,
  Zmat = as.matrix(nhanes_std[, exposure_vars, drop = FALSE])
)

message("Computing summaries (weighted) ...")
summ_weighted <- make_bkmr_summaries(
  fit_weighted,
  Zmat = as.matrix(dat_weighted[, exposure_vars, drop = FALSE])
)

univar_all_df <- bind_rows(
  summ_naive$univar %>% mutate(method = "naive"),
  summ_weighted$univar %>% mutate(method = "weighted")
)
write.csv(univar_all_df, file.path(outdir_real, "DATA_univar_naive_weighted_STD_clip_m3_p3.csv"), row.names = FALSE)

overall_all_df <- bind_rows(
  summ_naive$overall %>% mutate(method = "naive"),
  summ_weighted$overall %>% mutate(method = "weighted")
)
write.csv(overall_all_df, file.path(outdir_real, "DATA_overallrisk_naive_weighted_STD_clip_m3_p3.csv"), row.names = FALSE)

singvar_all_df <- bind_rows(
  summ_naive$singvar %>% mutate(method = "naive"),
  summ_weighted$singvar %>% mutate(method = "weighted")
)
write.csv(singvar_all_df, file.path(outdir_real, "DATA_singvarrisk_naive_weighted_STD_clip_m3_p3.csv"), row.names = FALSE)

bivar_levels_all_df <- bind_rows(
  summ_naive$bivar_levels %>% mutate(method = "naive"),
  summ_weighted$bivar_levels %>% mutate(method = "weighted")
)
write.csv(bivar_levels_all_df, file.path(outdir_real, "DATA_bivar_levels_naive_weighted_STD_clip_m3_p3.csv"), row.names = FALSE)

if (SAVE_BIVAR_SURFACE) {
  bivar_surface_all_df <- bind_rows(
    summ_naive$bivar %>% mutate(method = "naive"),
    summ_weighted$bivar %>% mutate(method = "weighted")
  )
  write.csv(bivar_surface_all_df, file.path(outdir_real, "DATA_bivar_surface_naive_weighted_STD_clip_m3_p3.csv"), row.names = FALSE)
}

# ---------------------------
# 11) PLOTS (SHOW + SAVE PNG) — NO VARIABLE LIMITING
# ---------------------------
plot_univar_compare_all <- function(univar_df) {
  df <- univar_df %>% mutate(method = recode(method, naive = "Naive", weighted = "Weighted"))
  
  ggplot(df, aes(x = z, y = est, ymin = est - 1.96 * se, ymax = est + 1.96 * se, color = method)) +
    geom_ribbon(alpha = 0.12, aes(fill = method, color = NULL)) +
    geom_line(linewidth = 0.7) +
    facet_wrap(~ variable, scales = "free_x", ncol = 5) +
    labs(x = "z (standardized + clipped to [-3,3])",
         y = "h(z)",
         color = "Method", fill = "Method",
         title = "Univariate exposure–response (Naive vs Weighted)") +
    theme_minimal(base_size = 12)
}

plot_overall_compare <- function(overall_df) {
  df <- overall_df %>% mutate(method = recode(method, naive = "Naive", weighted = "Weighted"))
  
  ggplot(df, aes(x = quantile, y = est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_ribbon(alpha = 0.10, aes(fill = method, color = NULL)) +
    geom_line(linewidth = 0.8) +
    labs(x = "Mixture quantile",
         y = "Overall mixture effect",
         color = "Method", fill = "Method",
         title = "Overall mixture effect curve (Naive vs Weighted)") +
    theme_minimal(base_size = 12)
}

plot_singvar_compare <- function(singvar_df) {
  df <- singvar_df %>% mutate(method = recode(method, naive = "Naive", weighted = "Weighted"))
  
  ggplot(df, aes(x = variable, y = est, ymin = est - 1.96 * sd, ymax = est + 1.96 * sd, color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_pointrange(position = position_dodge(width = 0.6)) +
    coord_flip() +
    facet_wrap(~ q.fixed, ncol = 1,
               labeller = labeller(q.fixed = function(x) paste0("q_fixed = ", x))) +
    labs(x = "",
         y = "Δ (75th vs 25th percentile)",
         color = "Method",
         title = "Single-variable risk summaries (Naive vs Weighted)") +
    theme_minimal(base_size = 12)
}

plot_bivar_levels_compare <- function(bivar_levels_df, vars_show = NULL) {
  df <- bivar_levels_df %>% mutate(Method = recode(method, naive = "Naive", weighted = "Weighted"))
  
  if (!is.null(vars_show)) {
    df <- df %>% filter(variable1 %in% vars_show, variable2 %in% vars_show)
    df$variable1 <- factor(df$variable1, levels = vars_show)
    df$variable2 <- factor(df$variable2, levels = vars_show)
  }
  
  df$quantile <- factor(df$quantile, levels = c("0.25","0.5","0.75"))
  
  ggplot(df, aes(x = z1, y = est, linetype = Method, color = quantile)) +
    geom_line(linewidth = 0.7) +
    facet_grid(variable2 ~ variable1, switch = "both") +
    labs(title = "Bivariate exposure–response levels (Naive vs Weighted)",
         x = "z1 (standardized + clipped)",
         y = "h(z1 | quantiles of z2)",
         linetype = "Method", color = "Quantile of 2nd exposure") +
    theme_minimal(base_size = 11) +
    theme(panel.spacing = unit(0.45, "lines"))
}

p_univar       <- plot_univar_compare_all(univar_all_df)
p_overall      <- plot_overall_compare(overall_all_df)
p_singvar      <- plot_singvar_compare(singvar_all_df)
p_bivar_levels <- plot_bivar_levels_compare(bivar_levels_all_df, vars_show = VARS_SUBSET_FOR_BIVAR_PLOT)

# Print to screen (guaranteed)
print(p_univar)
print(p_overall)
print(p_singvar)
print(p_bivar_levels)

# Save PNG
ggsave(file.path(outdir_real, "FIG_univar_naive_vs_weighted_STD_clip_m3_p3.png"), p_univar, width = 12, height = 8, dpi = 300)
ggsave(file.path(outdir_real, "FIG_overall_naive_vs_weighted_STD_clip_m3_p3.png"), p_overall, width = 9, height = 6, dpi = 300)
ggsave(file.path(outdir_real, "FIG_singvar_naive_vs_weighted_STD_clip_m3_p3.png"), p_singvar, width = 10, height = 8, dpi = 300)
ggsave(file.path(outdir_real, "FIG_bivar_levels_naive_vs_weighted_STD_clip_m3_p3.png"), p_bivar_levels, width = 15, height = 10, dpi = 300)

# ---------------------------
# 12) COMPARISON SUMMARY TABLE (single-fit naive vs weighted)
# ---------------------------
overall_naive_delta    <- extract_overall_delta(fit_naive)
overall_weighted_delta <- extract_overall_delta(fit_weighted)

comparison_summary <- data.frame(
  method = c("naive_singlefit", "weighted_singlefit"),
  est = c(overall_naive_delta["mean"], overall_weighted_delta["mean"]),
  lwr = c(overall_naive_delta["lwr"],  overall_weighted_delta["lwr"]),
  upr = c(overall_naive_delta["upr"],  overall_weighted_delta["upr"])
)
write.csv(comparison_summary,
          file.path(outdir_real, "TABLE_comparison_summary_naive_vs_weighted_STD_clip_m3_p3.csv"),
          row.names = FALSE)


plot_pips_variable_compare <- function(pip_master_df) {
  df <- pip_master_df %>%
    filter(level == "variable", metric == "pip") %>%
    mutate(
      method = recode(method, naive = "Naive", weighted = "Weighted")
    )
  
  ggplot(df, aes(x = reorder(variable, value), y = value, fill = method)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    coord_flip() +
    labs(
      x = "",
      y = "Posterior Inclusion Probability (PIP)",
      fill = "Method",
      title = "Variable-level PIPs (Naive vs Weighted)"
    ) +
    theme_minimal(base_size = 12)
}


# ---------------------------
# 13) MASTER PIP DATASET (variable + group + conditional) for BOTH methods
# ---------------------------
message("Building PIP master (naive + weighted) ...")
pip_master <- bind_rows(
  build_pip_master(fit_naive, "naive", exposure_vars, pip_groups),
  build_pip_master(fit_weighted, "weighted", exposure_vars, pip_groups)
) %>%
  arrange(method, factor(level, levels = c("group","variable","conditional")), group, variable, given)

write.csv(pip_master,
          file.path(outdir_real, "DATA_PIPs_MASTER_group_variable_conditional_naive_weighted_STD_clip_m3_p3.csv"),
          row.names = FALSE)

# ---------------------------
# 14) DESIGN-AWARE OVERALL MIXTURE EFFECT (many PSU bootstraps) — PARALLEL
# ---------------------------
run_one_boot <- function(b, dat) {
  rs <- fit_kmbayes_resampled(
    dat, y_var, exposure_vars, covar_vars,
    psu_col, strata_col, wt_col,
    iter = ITER_BOOT, burn = BURN_BOOT, thin = THIN_BOOT, varsel = VARSEL
  )
  extract_overall_delta(rs$fit)["mean"]
}

message(sprintf("Running design-aware PSU bootstrap (R=%d) with %d workers ...", R_BOOT, N_WORKERS))
boot_est <- future_sapply(seq_len(R_BOOT), run_one_boot, dat = nhanes_std)

boot_df <- data.frame(overall_boot = as.numeric(boot_est))
write.csv(boot_df,
          file.path(outdir_real, "DATA_overall_bootstrap_estimates_STD_clip_m3_p3.csv"),
          row.names = FALSE)

perf_summary <- data.frame(
  metric = c("design_aware_point", "design_aware_ci_lower", "design_aware_ci_upper"),
  value  = c(mean(boot_est, na.rm = TRUE),
             as.numeric(quantile(boot_est, 0.025, na.rm = TRUE)),
             as.numeric(quantile(boot_est, 0.975, na.rm = TRUE)))
)
write.csv(perf_summary,
          file.path(outdir_real, "TABLE_performance_summary_designaware_STD_clip_m3_p3.csv"),
          row.names = FALSE)

# ---------------------------
# 15) FINAL REPORT
# ---------------------------
cat("\n✅ SAVED TO: ", outdir_real, "\n\n",
    "FIGURES (PNG):\n",
    "- FIG_univar_naive_vs_weighted_STD_clip_m3_p3.png\n",
    "- FIG_overall_naive_vs_weighted_STD_clip_m3_p3.png\n",
    "- FIG_singvar_naive_vs_weighted_STD_clip_m3_p3.png\n",
    "- FIG_bivar_levels_naive_vs_weighted_STD_clip_m3_p3.png\n\n",
    "FIGURE DATASETS (CSV):\n",
    "- DATA_univar_naive_weighted_STD_clip_m3_p3.csv\n",
    "- DATA_overallrisk_naive_weighted_STD_clip_m3_p3.csv\n",
    "- DATA_singvarrisk_naive_weighted_STD_clip_m3_p3.csv\n",
    "- DATA_bivar_levels_naive_weighted_STD_clip_m3_p3.csv\n",
    if (SAVE_BIVAR_SURFACE) "- DATA_bivar_surface_naive_weighted_STD_clip_m3_p3.csv\n" else "",
    "\nTABLES (CSV):\n",
    "- TABLE_comparison_summary_naive_vs_weighted_STD_clip_m3_p3.csv\n",
    "- TABLE_performance_summary_designaware_STD_clip_m3_p3.csv\n\n",
    "BOOTSTRAP (CSV):\n",
    "- DATA_overall_bootstrap_estimates_STD_clip_m3_p3.csv\n\n",
    "PIPs (CSV):\n",
    "- DATA_PIPs_MASTER_group_variable_conditional_naive_weighted_STD_clip_m3_p3.csv\n\n",
    "QC (CSV):\n",
    "- QC_standardization_range_minus3_plus3.csv\n",
    sep = ""
)

cat("\n===== Overall mixture effect (25th→75th) =====\n")
cat("Naive single fit:\n"); print(overall_naive_delta)
cat("\nWeighted single fit:\n"); print(overall_weighted_delta)
cat("\nDesign-aware (PSU bootstrap) point + CI:\n"); print(perf_summary)
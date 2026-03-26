# ===============================
# Simulation: BKMR under Complex Survey Sampling (p = 10)
# Compare: Naïve BKMR vs Survey-aware BKMR via weighted resampling (PSU bootstrap)
# Now ALSO extracting PIPs (per exposure) for naive vs weighted fits + PIP plots
# ===============================

suppressPackageStartupMessages({
  library(MASS)          # mvrnorm
  library(dplyr)
  library(ggplot2)
  library(bkmr)          # kmbayes + BKMR helpers
  library(future.apply)  # parallel
  library(tidyr)
  library(scales)
  library(readr)
})

set.seed(20250921)  # reproducible launcher seed

# ===============================
# 0) CONFIG
# ===============================
CFG <- list(
  M_exposures = 3,               # *** p = 10 exposures ***
  N_pop       = 20000,            # finite population size
  n_clusters  = 120,              # PSUs in population
  n_strata    = 6,                # strata count
  clusters_per_stratum = rep(4, 6), # sampled clusters/stratum
  nsim        = 100,              # Monte Carlo reps per scenario
  warmup      = 2000,             # (if you later use sel to drop burn-in)
  iter        = 4000,             # total iterations per BKMR fit
  varsel      = TRUE,             # PIP: must be TRUE to populate delta / PIPs
  psu_bootstrap = TRUE            # TRUE = PSU bootstrap resampling
)

EXPOS <- paste0("Z", 1:CFG$M_exposures)  # Z1,...,Z10

# PIP: one group per exposure (simulation has no natural groups)
groups_sim <- 1:CFG$M_exposures

# Helper: mark which exposures truly drive the outcome
# In this DGP: Z1, Z2, Z3 are strong; Z4–Z10 are weak
truth_signal_label <- function(exposure_names) {
  sapply(exposure_names, function(z) {
    idx <- as.integer(sub("Z", "", z))
    if (!is.na(idx) && idx <= 3) "Strong" else "Weak"
  })
}

# ===============================
# 1) POPULATION GENERATION
# ===============================
# h_fun now takes a vector Z (length >= 3), used for both truth and outcome.
h_fun <- function(Z) {
  # Strong structure on first 3 exposures (same as before)
  val <- 1.0 * Z[1]*Z[2] +
    0.5 * Z[1]^2  -
    0.6 * Z[2]^2  +
    0.30 * Z[3]   +
    0.25 * Z[3]^2
  
  # Optional weaker linear effects on the remaining exposures (Z4–Z10)
  if (length(Z) > 3) {
    k <- length(Z) - 3
    coeffs <- rep(0.1, k)   # all 0.1, small but nonzero effect
    val <- val + sum(coeffs * Z[4:length(Z)])
  }
  val
}

generate_population <- function(N_pop, M_exposures, rho, n_clusters, ICC, sigma_error = 1) {
  # exposures: exchangeable covariance (ρ off-diagonal, 1 on diagonal)
  Sigma <- matrix(rho, nrow = M_exposures, ncol = M_exposures); diag(Sigma) <- 1
  Z <- MASS::mvrnorm(n = N_pop, mu = rep(0, M_exposures), Sigma = Sigma)
  colnames(Z) <- paste0("Z", seq_len(M_exposures))
  
  # cluster assignment (PSUs) ~ uniform
  cluster <- sample(seq_len(n_clusters), N_pop, replace = TRUE)
  
  # cluster random effects (induces ICC in Y)
  sigma_cluster <- sqrt(ICC)
  u_cluster <- rnorm(n_clusters, mean = 0, sd = sigma_cluster)
  u <- u_cluster[cluster]
  
  data.frame(ID = seq_len(N_pop), cluster = cluster, u = u, Z, check.names = FALSE)
}

simulate_outcome <- function(pop, sigma_error = 1) {
  # h_true uses ALL exposures via h_fun(Z)
  pop$h_true <- apply(pop[, EXPOS, drop = FALSE], 1, h_fun)
  pop$Y      <- pop$h_true + pop$u + rnorm(nrow(pop), 0, sigma_error)
  pop
}

# ===============================
# 2) STRATIFICATION & SAMPLING
# ===============================
assign_strata <- function(cluster, n_strata) {
  clus <- sort(unique(cluster))
  strata <- rep(seq_len(n_strata), length.out = length(clus))
  m <- setNames(strata, clus)
  unname(m[as.character(cluster)])
}

draw_sample <- function(pop, n_strata, clusters_per_stratum, n_per_cluster,
                        informative_weights = TRUE) {
  stopifnot(is.data.frame(pop),
            "cluster" %in% names(pop),
            length(clusters_per_stratum) == n_strata)
  
  pop$stratum <- assign_strata(pop$cluster, n_strata)
  sel_ids <- integer(0)
  C_all <- sort(unique(pop$cluster))
  P_cl  <- setNames(rep(0, length(C_all)), C_all)
  
  for (h in seq_len(n_strata)) {
    clus_h <- unique(pop$cluster[pop$stratum == h])
    if (length(clus_h) == 0) next
    
    k <- min(clusters_per_stratum[h], length(clus_h))
    if (k < clusters_per_stratum[h]) {
      warning(sprintf("Stratum %d: requested %d clusters, only %d available; using %d.",
                      h, clusters_per_stratum[h], length(clus_h), k))
    }
    set_h <- sample(clus_h, size = k, replace = FALSE)
    P_h <- k / length(clus_h)
    P_cl[as.character(set_h)] <- P_h
    
    for (cl in set_h) {
      idx <- which(pop$cluster == cl)
      Nc <- length(idx)
      mc <- min(n_per_cluster, Nc)
      if (mc <= 0) next
      
      if (informative_weights && "Z1" %in% names(pop)) {
        z1 <- pop$Z1[idx]
        pr <- z1 - min(z1, na.rm = TRUE) + 1e-6
        pr[!is.finite(pr) | is.na(pr)] <- 0
        s  <- sum(pr)
        if (s > 0) {
          pr <- pr / s
          take <- sample(idx, size = mc, replace = FALSE, prob = pr)
        } else {
          take <- sample(idx, size = mc, replace = FALSE)
        }
      } else {
        take <- sample(idx, size = mc, replace = FALSE)
      }
      sel_ids <- c(sel_ids, take)
    }
  }
  
  if (length(sel_ids) == 0) stop("Sampling produced 0 rows; check inputs.")
  samp <- pop[sel_ids, , drop = FALSE]
  rownames(samp) <- NULL
  
  size_by_cl <- table(pop$cluster)
  n_by_cl    <- table(samp$cluster)
  P_cl_used  <- P_cl[as.character(samp$cluster)]
  Nc_vec     <- as.numeric(size_by_cl[as.character(samp$cluster)])
  mc_vec     <- as.numeric(n_by_cl[as.character(samp$cluster)])
  
  if (informative_weights && "Z1" %in% names(samp)) {
    p_within <- numeric(nrow(samp))
    for (cl in names(n_by_cl)) {
      idc <- which(samp$cluster == as.integer(cl))
      z1c <- samp$Z1[idc]
      prc <- z1c - min(z1c, na.rm = TRUE) + 1e-6
      prc[!is.finite(prc) | is.na(prc)] <- 0
      sc  <- sum(prc)
      if (sc > 0) prc <- prc / sc else prc <- rep(1/length(idc), length(idc))
      p_within[idc] <- pmin(1, length(idc) * prc)
    }
    pi_within <- p_within
  } else {
    pi_within <- mc_vec / Nc_vec
  }
  
  pi_cluster <- P_cl_used
  pi_hat     <- pmax(pi_cluster * pi_within, 1e-10)
  w          <- 1 / pi_hat
  samp$weight <- as.numeric(w / mean(w))
  samp
}

# ===============================
# 3) GROUND-TRUTH ESTIMAND
# ===============================
true_overall_effect <- function(pop) {
  qs <- apply(pop[, EXPOS, drop = FALSE], 2, quantile, probs = c(0.25, 0.75))
  Z_lo <- qs[1, ]
  Z_hi <- qs[2, ]
  Y_hi <- h_fun(Z_hi)
  Y_lo <- h_fun(Z_lo)
  unname(Y_hi - Y_lo)
}

# ===============================
# 4) BKMR FITS + SUMMARIES + PIPs
# ===============================
fit_kmbayes_unweighted <- function(samp, expos_names,
                                   iter = CFG$iter,
                                   varsel = CFG$varsel) {
  Z <- as.matrix(samp[, expos_names, drop = FALSE])
  y <- samp$Y
  kmbayes(
    y = y,
    Z = Z,
    iter   = iter,
    varsel = varsel,
    groups = groups_sim,   # PIP: one group per exposure
    verbose = FALSE
  )
}

resample_by_weight_unit <- function(df, size = NULL) {
  w <- df$weight
  p <- w / sum(w)
  if (is.null(size)) size <- nrow(df)
  take <- sample(seq_len(nrow(df)), size = size, replace = TRUE, prob = p)
  df[take, , drop = FALSE]
}

resample_by_weight_psu <- function(df, size_psu = NULL) {
  clus <- sort(unique(df$cluster))
  w_cl <- sapply(clus, function(c) mean(df$weight[df$cluster == c]))
  p_cl <- w_cl / sum(w_cl)
  if (is.null(size_psu)) size_psu <- length(clus)
  sel_psu <- sample(clus, size = size_psu, replace = TRUE, prob = p_cl)
  
  out <- list()
  for (c in sel_psu) {
    sub <- df[df$cluster == c, , drop = FALSE]
    m   <- nrow(sub)
    out[[length(out)+1]] <- resample_by_weight_unit(sub, size = m)
  }
  do.call(rbind, out)
}

fit_kmbayes_resampled <- function(samp, expos_names,
                                  iter = CFG$iter,
                                  varsel = CFG$varsel,
                                  psu_bootstrap = CFG$psu_bootstrap) {
  if (psu_bootstrap) {
    samp_rs <- resample_by_weight_psu(samp)
  } else {
    samp_rs <- resample_by_weight_unit(samp)
  }
  Z <- as.matrix(samp_rs[, expos_names, drop = FALSE])
  y <- samp_rs$Y
  fit <- kmbayes(
    y = y,
    Z = Z,
    iter   = iter,
    varsel = varsel,
    groups = groups_sim,  # PIP: one group per exposure
    verbose = FALSE
  )
  list(fit = fit, samp_rs = samp_rs)
}

make_bkmr_summaries <- function(fit, Zmat) {
  pred.resp.univar <- PredictorResponseUnivar(fit = fit, method = "approx")
  pred.resp.bivar  <- PredictorResponseBivar(fit = fit, min.plot.dist = 1, method = "approx")
  pred.resp.bivar.levels <- PredictorResponseBivarLevels(
    pred.resp.df = pred.resp.bivar, Z = Zmat, both_pairs = TRUE, qs = c(0.25, 0.5, 0.75)
  )
  risks.overall <- OverallRiskSummaries(fit = fit, qs = seq(0.25, 0.75, by = 0.05),
                                        q.fixed = 0.5, method = "approx")
  risks.singvar <- SingVarRiskSummaries(fit = fit, qs.diff = c(0.25, 0.75),
                                        q.fixed = c(0.25, 0.50, 0.75), method = "approx")
  risks.int <- SingVarIntSummaries(fit = fit, qs.diff = c(0.25, 0.75),
                                   qs.fixed = c(0.25, 0.75))
  list(univar = pred.resp.univar,
       bivar = pred.resp.bivar,
       bivar_levels = pred.resp.bivar.levels,
       overall = risks.overall,
       singvar = risks.singvar,
       int = risks.int)
}

extract_overall_delta <- function(fit, method = "approx") {
  or <- OverallRiskSummaries(fit = fit, qs = c(0.25, 0.75), q.fixed = 0.5, method = method)
  est <- with(or, est[quantile == 0.75] - est[quantile == 0.25])
  sdv <- with(or, sqrt(sd[quantile == 0.75]^2 + sd[quantile == 0.25]^2))
  lwr <- est - 1.96 * sdv
  upr <- est + 1.96 * sdv
  c(mean = est, lwr = lwr, upr = upr)
}

# ---------- PIP helper (mirrors NHANES code) ----------
compute_pips <- function(fit, expos_names, groups_vec = NULL) {
  # 1) PIPs from summary()
  s <- summary(fit)
  pip_summary <- as.data.frame(s$pip)
  
  # 2) Exposure-level PIPs from delta
  delta <- fit$delta
  
  if (is.null(colnames(delta)) || length(colnames(delta)) == 0) {
    colnames(delta) <- expos_names
  }
  
  conditional_pips <- colMeans(delta)  # per exposure
  
  # 3) Group PIPs (here each exposure is its own group, but keep structure)
  if (!is.null(groups_vec)) {
    group_ids <- sort(unique(groups_vec))
    group_pips <- sapply(group_ids, function(g) {
      mean(rowSums(delta[, groups_vec == g, drop = FALSE]) > 0)
    })
  } else {
    group_pips <- NULL
  }
  
  # 4) Combined table
  pip_df <- data.frame(
    Exposure        = colnames(delta),
    Conditional_PIP = conditional_pips,
    stringsAsFactors = FALSE
  )
  if (!is.null(groups_vec)) {
    pip_df$Group     <- groups_vec
    pip_df$Group_PIP <- group_pips[pip_df$Group]
  }
  
  list(
    summary_pip     = pip_summary,
    conditional_pip = conditional_pips,
    group_pip       = group_pips,
    pip_table       = pip_df
  )
}

# ===============================
# 5) ONE REPLICATE
# ===============================
one_replicate <- function(rho, n, ICC, informative) {
  pop <- generate_population(CFG$N_pop, CFG$M_exposures, rho, CFG$n_clusters, ICC)
  pop <- simulate_outcome(pop, sigma_error = 1)
  truth <- true_overall_effect(pop)
  
  n_per_cluster <- ceiling(n / sum(CFG$clusters_per_stratum))
  samp <- draw_sample(pop,
                      n_strata = CFG$n_strata,
                      clusters_per_stratum = CFG$clusters_per_stratum,
                      n_per_cluster = n_per_cluster,
                      informative_weights = informative)
  
  # Naïve BKMR
  fit_unw <- fit_kmbayes_unweighted(samp, EXPOS, iter = CFG$iter, varsel = CFG$varsel)
  est_unw <- extract_overall_delta(fit_unw)
  summaries_unw <- make_bkmr_summaries(fit_unw, Zmat = as.matrix(samp[, EXPOS, drop = FALSE]))
  pips_unw      <- compute_pips(fit_unw, expos_names = EXPOS, groups_vec = groups_sim)
  
  # Survey-aware BKMR (weighted resample)
  rs <- fit_kmbayes_resampled(samp, EXPOS, iter = CFG$iter, varsel = CFG$varsel,
                              psu_bootstrap = CFG$psu_bootstrap)
  fit_rs <- rs$fit
  summaries_wrs <- make_bkmr_summaries(fit_rs, Zmat = as.matrix(rs$samp_rs[, EXPOS, drop = FALSE]))
  est_wrs <- extract_overall_delta(fit_rs)
  pips_wrs <- compute_pips(fit_rs, expos_names = EXPOS, groups_vec = groups_sim)
  
  data.frame(
    rho = rho, n = n, ICC = ICC, informative = informative,
    truth = truth,
    est_unw = est_unw["mean"], lwr_unw = est_unw["lwr"], upr_unw = est_unw["upr"],
    est_wrs = est_wrs["mean"], lwr_wrs = est_wrs["lwr"], upr_wrs = est_wrs["upr"],
    summaries_naive    = I(list(summaries_unw)),
    summaries_weighted = I(list(summaries_wrs)),
    pips_naive         = I(list(pips_unw)),    # PIP: store naive PIPs
    pips_weighted      = I(list(pips_wrs))     # PIP: store weighted PIPs
  )
}

# ===============================
# 6) FACTORIAL MONTE CARLO
# ===============================
future::plan(future::multisession)

run_factorial <- function(nsim = CFG$nsim,
                          corr_levels = c(0.0, 0.8),
                          n_levels    = c(300, 800),
                          ICC_levels  = c(0.0, 0.15),
                          inf_levels  = c(FALSE, TRUE),
                          keep_first_plots_only = TRUE) {
  
  scenarios <- expand.grid(rho = corr_levels, n = n_levels, ICC = ICC_levels, informative = inf_levels)
  results_list <- vector("list", nrow(scenarios))
  
  for (i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i, ]
    cat(sprintf("Scenario %d/%d: rho=%.1f, n=%d, ICC=%.2f, inf=%s\n",
                i, nrow(scenarios), sc$rho, sc$n, sc$ICC, sc$informative))
    
    reps <- future.apply::future_lapply(
      seq_len(nsim),
      function(r) {
        set.seed(1000 + i*100 + r)
        one_replicate(sc$rho, sc$n, sc$ICC, sc$informative)
      },
      future.seed = TRUE
    )
    
    df <- do.call(rbind, reps)
    
    if (keep_first_plots_only) {
      if (nrow(df) > 1) {
        df$summaries_naive[-1]    <- list(NULL)
        df$summaries_weighted[-1] <- list(NULL)
        df$pips_naive[-1]         <- list(NULL)  # PIP: keep only first replicate PIPs
        df$pips_weighted[-1]      <- list(NULL)
      }
    }
    results_list[[i]] <- df
  }
  do.call(rbind, results_list)
}

# ===============================
# 7) PLOT HELPERS (COMPARING NAÏVE vs WEIGHTED)
# ===============================
plot_univar_compare <- function(univar_naive, univar_wt, vars_to_show = NULL) {
  univar_naive$method <- "Naïve"
  univar_wt$method    <- "Weighted"
  df <- bind_rows(univar_naive, univar_wt)
  if (!is.null(vars_to_show)) {
    df <- df %>% filter(variable %in% vars_to_show)
  }
  ggplot(df, aes(z, est, ymin = est - 1.96*se, ymax = est + 1.96*se, color = method)) +
    geom_ribbon(alpha = 0.1, aes(fill = method, color = NULL)) +
    geom_line(size = 0.7) +
    facet_wrap(~ variable, scales = "free_x") +
    ylab("h(z)") +
    theme_minimal(base_size = 12)
}

plot_overall_compare <- function(overall_naive, overall_wt) {
  overall_naive$method <- "Naïve"
  overall_wt$method    <- "Weighted"
  df <- bind_rows(overall_naive, overall_wt)
  ggplot(df, aes(quantile, est, ymin = est - 1.96*sd, ymax = est + 1.96*sd,
                 color = method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
    geom_ribbon(aes(fill = method, color = NULL), alpha = 0.1) +
    geom_line(size = 0.8) +
    scale_y_continuous(name = "Overall mixture effect") +
    theme_minimal(base_size = 12)
}

# ===============================
# 8) PILOT RUN (ADJUST nsim LATER)
# ===============================
cat("NOTE: BKMR MCMC is compute-heavy; start small to test plumbing.\n")
pilot <- run_factorial(
  nsim = 5,
  corr_levels = c(0.0, 0.8),
  n_levels    = c(300, 800),
  ICC_levels  = c(0.0, 0.15),
  inf_levels  = c(FALSE, TRUE),
  keep_first_plots_only = TRUE
)

# ===============================
# 9) PERFORMANCE SUMMARY (NAÏVE vs WEIGHTED)
# ===============================
perf <- pilot %>%
  mutate(covered_unw = as.numeric(truth >= lwr_unw & truth <= upr_unw),
         covered_wrs = as.numeric(truth >= lwr_wrs & truth <= upr_wrs),
         width_unw   = upr_unw - lwr_unw,
         width_wrs   = upr_wrs - lwr_wrs) %>%
  group_by(rho, n, ICC, informative) %>%
  summarise(
    bias_unw  = mean(est_unw - truth, na.rm = TRUE),
    bias_wrs  = mean(est_wrs - truth, na.rm = TRUE),
    cover_unw = mean(covered_unw, na.rm = TRUE),
    cover_wrs = mean(covered_wrs, na.rm = TRUE),
    width_unw = mean(width_unw, na.rm = TRUE),
    width_wrs = mean(width_wrs, na.rm = TRUE),
    .groups   = "drop"
  )

print(perf)

outdir <- "Output"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write.csv(perf, file.path(outdir, "performance_summary_p3.csv"), row.names = FALSE)

# Scenario-wise deltas
cmp <- perf %>%
  mutate(
    d_bias   = bias_wrs - bias_unw,
    d_cover  = cover_wrs - cover_unw,
    d_width  = width_wrs - width_unw
  ) %>%
  arrange(desc(informative), desc(ICC), rho, n)

cmp_view <- cmp %>%
  mutate(
    winner = case_when(
      abs(bias_wrs) < abs(bias_unw) & cover_wrs >= cover_unw ~ "weighted",
      abs(bias_wrs) < abs(bias_unw) & cover_wrs <  cover_unw ~ "weighted (bias)",
      abs(bias_wrs) >= abs(bias_unw) & cover_wrs > cover_unw ~ "weighted (coverage)",
      TRUE ~ "naive or tie"
    )
  )

ts <- format(Sys.time(), "%Y%m%d-%H%M%S")
cmp_file <- file.path(outdir, paste0("scenario_comparison_deltas_", ts, ".csv"))
cmp_view_file <- file.path(outdir, paste0("scenario_comparison_winner_", ts, ".csv"))
write.csv(cmp, cmp_file, row.names = FALSE)
write.csv(cmp_view, cmp_view_file, row.names = FALSE)

message("Wrote:\n  - ", cmp_file, "\n  - ", cmp_view_file)

# ===============================
# 10) METRIC FIGURES (NAÏVE vs WEIGHTED)
# ===============================
lab_scenario <- function(df) {
  df %>%
    mutate(
      inf_lab = if_else(informative, "Informative", "Non-informative"),
      ICC_lab = paste0("ICC=", formatC(ICC, format = "f", digits = 2)),
      rho_lab = paste0("rho=", formatC(rho, format = "f", digits = 1)),
      n_lab   = paste0("n=", n),
      scenario = paste(rho_lab, ICC_lab, n_lab, inf_lab, sep = " | ")
    )
}

cmp      <- lab_scenario(cmp)
cmp_view <- lab_scenario(cmp_view)

p_winner <- ggplot(cmp_view, aes(x = factor(n), y = interaction(rho, ICC, informative, sep = " | "))) +
  geom_tile(aes(fill = winner), color = "white") +
  scale_fill_manual(values = c("weighted" = "#2c7fb8",
                               "weighted (bias)" = "#74a9cf",
                               "weighted (coverage)" = "#a6bddb",
                               "naive or tie" = "#fdbb84")) +
  labs(x = "Sample size (n)",
       y = "Scenario (rho | ICC | informative)",
       fill = "Winner",
       title = "Naïve vs Weighted BKMR — winner by scenario (p = 3)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "cmp_winner_heatmap_p10.png"), p_winner, width = 10, height = 6, dpi = 300)

cmp_long_bias <- cmp_view %>%
  select(rho, n, ICC, informative, bias_unw, bias_wrs, scenario) %>%
  pivot_longer(c(bias_unw, bias_wrs), names_to = "method", values_to = "bias") %>%
  mutate(method = recode(method, bias_unw = "Naïve", bias_wrs = "Weighted"))

p_bias <- ggplot(cmp_long_bias,
                 aes(x = bias, y = factor(paste0("rho=", rho, ", n=", n)))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(aes(group = interaction(rho, n)), linewidth = 0.6, color = "grey60") +
  geom_point(aes(color = method), size = 2.6) +
  facet_grid(informative ~ ICC, labeller = label_both) +
  labs(x = "Bias (overall mixture effect)",
       y = "Scenario row",
       color = "Method",
       title = "Bias by scenario: Naïve vs Weighted BKMR (p = 3)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "cmp_bias_dumbbell_p10.png"), p_bias, width = 11, height = 7, dpi = 300)

p_cover <- ggplot(cmp_view, aes(x = cover_unw, y = cover_wrs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "grey50") +
  geom_hline(yintercept = 0.95, linetype = "dashed", color = "grey50") +
  geom_vline(xintercept = 0.95, linetype = "dashed", color = "grey50") +
  geom_point(aes(size = n, color = informative, shape = factor(ICC)), alpha = 0.8) +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_color_manual(values = c(`TRUE` = "#1b9e77", `FALSE` = "#7570b3"),
                     labels = c(`TRUE` = "Informative", `FALSE` = "Non-informative")) +
  labs(x = "Coverage (Naïve)",
       y = "Coverage (Weighted)",
       color = "Design",
       shape = "ICC",
       size  = "n",
       title = "Coverage comparison — closer to upper-left is better for weighted (p = 3)") +
  theme_minimal(base_size = 12)
ggsave(file.path(outdir, "cmp_coverage_scatter_p3.png"), p_cover, width = 8, height = 6, dpi = 300)

# ===============================
# 11) BKMR SHAPE FIGURES (NAÏVE vs WEIGHTED) + PIPs FOR ONE SCENARIO
# ===============================
first_with_plots <- pilot %>%
  filter(!sapply(summaries_naive, is.null)) %>%
  slice(1)

summaries_naive    <- first_with_plots$summaries_naive[[1]]
summaries_weighted <- first_with_plots$summaries_weighted[[1]]

# Overall risk curve: Naïve vs Weighted
p_overall <- plot_overall_compare(
  summaries_naive$overall,
  summaries_weighted$overall
)
ggsave(file.path(outdir, "bkmr_overall_naive_vs_weighted_p3.png"),
       p_overall, width = 9, height = 6, dpi = 300)

# Univariate curves for first 3 exposures as example (Z1–Z3)
p_univar <- plot_univar_compare(
  summaries_naive$univar,
  summaries_weighted$univar,
  vars_to_show = c("Z1", "Z2", "Z3")
)
ggsave(file.path(outdir, "bkmr_univar_naive_vs_weighted_p3.png"),
       p_univar, width = 9, height = 6, dpi = 300)

# ----- PIPs for the same scenario (Naïve vs Weighted) -----
first_with_pips <- pilot %>%
  filter(!sapply(pips_naive, is.null)) %>%
  slice(1)

pips_naive_1    <- first_with_pips$pips_naive[[1]]
pips_weighted_1 <- first_with_pips$pips_weighted[[1]]

# Save PIP tables (per exposure; group == exposure index)
write.csv(pips_naive_1$pip_table,
          file.path(outdir, "bkmr_pips_naive_firstscenario_p3.csv"),
          row.names = FALSE)

write.csv(pips_weighted_1$pip_table,
          file.path(outdir, "bkmr_pips_weighted_firstscenario_p3.csv"),
          row.names = FALSE)

# ===============================
# 12) PIP PLOTS (Naïve vs Weighted + Truth vs Strong/Weak)
# ===============================

pip_tab_naive    <- pips_naive_1$pip_table
pip_tab_weighted <- pips_weighted_1$pip_table

pip_tab_naive$Method    <- "Naïve"
pip_tab_weighted$Method <- "Weighted"

pip_long <- bind_rows(pip_tab_naive, pip_tab_weighted)

# Make Exposure an ordered factor: Z1, Z2, ..., Z10
pip_long$Exposure <- factor(pip_long$Exposure,
                            levels = paste0("Z", 1:CFG$M_exposures))

# ----- Plot 1: Naïve vs Weighted PIPs across exposures -----
p_pip_compare <- ggplot(pip_long,
                        aes(x = Exposure,
                            y = Conditional_PIP,
                            color = Method,
                            group = Method)) +
  geom_point(position = position_dodge(width = 0.4), size = 2.5) +
  geom_line(position = position_dodge(width = 0.4), linewidth = 0.7) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Posterior Inclusion Probabilities (Naïve vs Weighted)",
    x     = "Exposure",
    y     = "PIP",
    color = "Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

print(p_pip_compare)
ggsave(file.path(outdir, "bkmr_pips_naive_vs_weighted_p10.png"),
       p_pip_compare, width = 8, height = 5, dpi = 300)

# ----- Plot 2: Truth vs PIP (Strong vs Weak exposures) -----
pip_long$Truth <- truth_signal_label(as.character(pip_long$Exposure))
pip_long$Truth <- factor(pip_long$Truth, levels = c("Strong", "Weak"))

p_pip_truth <- ggplot(pip_long,
                      aes(x = Exposure,
                          y = Conditional_PIP,
                          color = Truth,
                          shape = Method)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70") +
  geom_point(position = position_dodge(width = 0.4), size = 2.8) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "PIPs vs True Signal Status (Strong vs Weak)",
    x     = "Exposure",
    y     = "PIP",
    color = "True signal?",
    shape = "Method"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )

print(p_pip_truth)
ggsave(file.path(outdir, "bkmr_pips_truth_vs_method_p3.png"),
       p_pip_truth, width = 8, height = 5, dpi = 300)

cat("Done. p = 10 simulation with naïve vs weighted results, figures, and PIP diagnostics generated.\n")






# ===============================
# 13) BKMR SHAPE FIGURES (Naïve vs Weighted) — FULL SUITE
# ===============================

# Overall risk curve (Naive vs Weighted)
plot_overall_both <- function(overall_unw, overall_wrs) {
  df_unw <- overall_unw %>% mutate(Method = "Naive")
  df_wrs <- overall_wrs %>% mutate(Method = "Weighted")
  df_all <- bind_rows(df_unw, df_wrs)
  
  ggplot(df_all,
         aes(x = quantile, y = est,
             ymin = est - 1.96 * sd,
             ymax = est + 1.96 * sd,
             color = Method, fill = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_ribbon(alpha = 0.1, color = NA) +
    geom_line(linewidth = 0.8) +
    labs(x = "Quantile of joint mixture",
         y = "Overall mixture effect",
         color = "Method",
         fill  = "Method",
         title = "Overall risk curve (Naive vs Weighted)") +
    theme_minimal(base_size = 12)
}

# Univariate exposure–response curves
plot_univar_both <- function(univar_unw, univar_wrs, vars_show = NULL) {
  df_unw <- univar_unw %>% mutate(Method = "Naive")
  df_wrs <- univar_wrs %>% mutate(Method = "Weighted")
  df_all <- bind_rows(df_unw, df_wrs)
  
  if (!is.null(vars_show)) {
    df_all <- df_all %>% filter(variable %in% vars_show)
  }
  
  ggplot(df_all,
         aes(x = z, y = est,
             ymin = est - 1.96 * se,
             ymax = est + 1.96 * se,
             color = Method, fill = Method)) +
    geom_ribbon(alpha = 0.08, color = NA) +
    geom_line(linewidth = 0.8) +
    facet_wrap(~ variable, scales = "free_x") +
    labs(x = "Exposure (z, standardized)",
         y = "h(z)",
         color = "Method",
         fill  = "Method",
         title = "Univariate exposure–response (Naive vs Weighted)") +
    theme_minimal(base_size = 12)
}

# Single–variable risk summaries
plot_singvar_both <- function(sv_unw, sv_wrs) {
  df_unw <- sv_unw %>% mutate(Method = "Naive")
  df_wrs <- sv_wrs %>% mutate(Method = "Weighted")
  df_all <- bind_rows(df_unw, df_wrs)
  
  ggplot(df_all,
         aes(x = variable, y = est,
             ymin = est - 1.96 * sd,
             ymax = est + 1.96 * sd,
             color = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_pointrange(
      position = position_dodge(width = 0.6),
      linewidth = 0.5
    ) +
    coord_flip() +
    facet_wrap(~ q.fixed, ncol = 1,
               labeller = label_bquote(q[fixed] == .(q.fixed))) +
    labs(x = "",
         y = expression(Delta~"(75th vs 25th percentile)"),
         color = "Method",
         title = "Single-variable risk summaries (Naive vs Weighted)") +
    theme_minimal(base_size = 12)
}

# Pairwise interaction summaries
plot_int_both <- function(int_unw, int_wrs) {
  df_unw <- int_unw %>% mutate(Method = "Naive")
  df_wrs <- int_wrs %>% mutate(Method = "Weighted")
  df_all <- bind_rows(df_unw, df_wrs)
  
  ggplot(df_all,
         aes(x = variable, y = est,
             ymin = est - 1.96 * sd,
             ymax = est + 1.96 * sd,
             color = Method)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_pointrange(
      position = position_dodge(width = 0.6),
      linewidth = 0.5
    ) +
    coord_flip() +
    labs(x = "",
         y = "Interaction effect",
         color = "Method",
         title = "Pairwise interaction summaries (Naive vs Weighted)") +
    theme_minimal(base_size = 12)
}

# Bivariate surface (faceted by Method)
plot_bivar_both <- function(biv_unw, biv_wrs) {
  df_unw <- biv_unw %>% mutate(Method = "Naive")
  df_wrs <- biv_wrs %>% mutate(Method = "Weighted")
  df_all <- bind_rows(df_unw, df_wrs)
  
  ggplot(df_all, aes(x = z1, y = z2, fill = est)) +
    geom_raster() +
    facet_grid(variable2 ~ interaction(variable1, Method)) +
    scale_fill_gradient2() +
    labs(x = "z1", y = "z2", fill = "h(z1, z2)",
         title = "Bivariate exposure–response surfaces\n(Naive vs Weighted)") +
    theme_minimal(base_size = 10)
}

# Bivariate exposure–response by quantiles of second exposure
plot_bivar_levels_both <- function(bl_unw, bl_wrs) {
  df_unw <- bl_unw %>% mutate(Method = "Naive")
  df_wrs <- bl_wrs %>% mutate(Method = "Weighted")
  df_all <- bind_rows(df_unw, df_wrs)
  
  ggplot(df_all,
         aes(x = z1, y = est,
             color = factor(quantile),
             linetype = Method)) +
    geom_line(linewidth = 0.7) +
    facet_grid(variable2 ~ variable1) +
    labs(x = "z1",
         y = "h(z1 | quantiles of z2)",
         color = "Quantile of 2nd exposure",
         linetype = "Method",
         title = "Bivariate exposure–response (Naive vs Weighted)") +
    theme_minimal(base_size = 10)
}

# ---- Build and save combined figures ----

# IMPORTANT: use summaries_naive and summaries_weighted (NOT 'summaries')
first_both <- pilot %>%
  filter(!sapply(summaries_naive,    is.null) &
           !sapply(summaries_weighted, is.null)) %>%
  slice(1)

plots_n <- first_both$summaries_naive[[1]]     # naive summaries
plots_w <- first_both$summaries_weighted[[1]]  # weighted summaries

# Build plots
p_overall   <- plot_overall_both       (plots_n$overall,       plots_w$overall)
p_univar    <- plot_univar_both        (plots_n$univar,        plots_w$univar)
p_singvar   <- plot_singvar_both       (plots_n$singvar,       plots_w$singvar)
p_int       <- plot_int_both           (plots_n$int,           plots_w$int)
p_bivar     <- plot_bivar_both         (plots_n$bivar,         plots_w$bivar)
p_bivar_lv  <- plot_bivar_levels_both  (plots_n$bivar_levels,  plots_w$bivar_levels)

# Print to screen (optional)
print(p_overall);  print(p_univar)
print(p_singvar);  print(p_int)
print(p_bivar);    print(p_bivar_lv)

# Save to disk
ggsave(file.path(outdir, "OVERALL_naive_weighted.png"),
       p_overall,  width = 9,  height = 6, dpi = 300)
ggsave(file.path(outdir, "UNIVAR_naive_weighted.png"),
       p_univar,   width = 11, height = 7, dpi = 300)
ggsave(file.path(outdir, "SINGVAR_naive_weighted.png"),
       p_singvar,  width = 11, height = 8, dpi = 300)
ggsave(file.path(outdir, "INT_naive_weighted.png"),
       p_int,      width = 9,  height = 6, dpi = 300)
ggsave(file.path(outdir, "BIVAR_surface_naive_weighted.png"),
       p_bivar,    width = 11, height = 7, dpi = 300)
ggsave(file.path(outdir, "BIVAR_levels_naive_weighted.png"),
       p_bivar_lv, width = 11, height = 7, dpi = 300)

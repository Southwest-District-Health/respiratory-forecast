

# ============================================================
# RESPIRATORY DISEASE FORECASTING PIPELINE  v11
#
# Inherits all fixes from v2-v9. LSTM removed in v11.
#
# v8 change — LSTM: CV-only, no prior fallback.
#
#   Previously (v1-v7): if LSTM CV returned zero valid CRPS
#   rows (due to a date-format mismatch bug fixed in v6),
#   LSTM entered the ensemble via a conservative prior equal
#   to the best base-model CRPS. This meant LSTM's BMA weight
#   was not evidence-based and could not be defended as an
#   improvement to forecast accuracy.
#
#   v8 change: prior fallback removed entirely. LSTM is now
#   CV-only (Option 1):
#     - If CV succeeds  -> LSTM enters ensemble with a real
#       CRPS-derived BMA weight.
#     - If CV fails     -> LSTM is excluded from the ensemble
#       and a clear diagnostic message is printed. No prior.
#
#   The date-mismatch bug (yearweek vs calendar Date) was fixed
#   in v6 via row-position matching in cv_lstm(). With that fix
#   in place, CV should succeed for Influenza and SARS-CoV2.
#   If it still fails after this fix, the diagnostic message
#   will tell you exactly how many valid CRPS rows were returned
#   so the cause can be investigated.
#
#   RSV: LSTM remains excluded by design (variance decomposition
#   showed 0.796 log-unit deviation from ensemble mean).
#
# v9 diagnostic fixes:
#   1. results_all now stores weights_final (models in fc_all only)
#      as bma_weights, and the full pre-filter weights as
#      bma_weights_full. Previously bma_weights could include LSTM
#      via a prior path even when LSTM was not in fc_all, causing
#      the RF vs LSTM diagnostic to show stale prior weights.
#   2. LSTM status diagnostic now covers all four possible states:
#      CV-derived / Excluded (0 rows) / Excluded (rows but not in
#      ensemble) / In ensemble but CV rows = 0 (prior still active).
#   3. RF vs LSTM comparison now labels LSTM weight source explicitly.
# ============================================================


suppressPackageStartupMessages({
  library(tsibble)
  library(fable)
  library(feasts)
  library(fabletools)
  library(dplyr)
  library(lubridate)
  library(tibble)
  library(readxl)
  library(zoo)
  library(purrr)
  library(furrr)
  library(distributional)
  library(ranger)
})

set.seed(42)

# ============================================================
# 1.  CONFIGURATION
# ============================================================

CONFIG <- list(
  
  horizon      = 4L,
  cv_init_frac = 0.65,
  cv_step      = 8L,
  ci_level     = 80L,
  cv_workers   = 4L,
  diseases     = c("Influenza", "RSV", "SARS-CoV2"),
  
  top_n_models = list(
    "Influenza"  = 4L,
    "RSV"        = 4L,   # v3: increased from 3 — cleaner pool after removing
    # ARIMA + ARIMAX_WW means top-4 are all credible
    "SARS-CoV2"  = 4L
  ),
  
  bma_temperature = list(
    "Influenza"  = 0.50,
    "RSV"        = 0.10,  # v4: restored from 0.05 — lower temp amplified RF
    # sigma through within_var; 0.10 distributes weight
    # more evenly across 6 models, diluting RF's sigma
    "SARS-CoV2"  = 0.50
  ),
  
  # Ensemble-level sigma cap: applied as hard ceiling per horizon
  # AFTER cummax() in build_ensemble(). FIX C: order swapped vs v2/v3.
  sigma_cap = list(
    "Influenza"  = NULL,
    "RSV"        = 0.55,  # v4: restored to 0.55 (v3's 0.65 was too loose)
    "SARS-CoV2"  = NULL
  ),
  
  
  # v6: tightened 0.38 -> 0.32. v5 achieved 1.623; target <1.4.
  # At h=4: cap = 0.32*sqrt(4) = 0.64. Bites on ETS and RF
  # which hit 0.760 in v5 — largest within_var contributors.
  sigma_model_cap = list(
    "Influenza"  = NULL,
    "RSV"        = 0.32,
    "SARS-CoV2"  = NULL
  ),
  
  # FIX B — RF sigma growth exponent per disease.
  # Default sqrt(h) = exponent 0.50 assumes random-walk error compounding.
  # RSV 0.40: strong seasonal autocorrelation constrains uncertainty growth.
  rf_sigma_exp = list(
    "Influenza"  = 0.50,
    "RSV"        = 0.40,
    "SARS-CoV2"  = 0.50
  ),
  
  
  
  rf = list(
    n_lag         = 8L,
    num.trees     = 500L,
    min.node.size = 3L,
    run_cv        = TRUE
  ),
  
  data_paths = list(
    "SARS-CoV2" = "~/Data/Essence/CLI.xlsx",
    "Influenza"  = "~/Data/Essence/ILI.xlsx",
    "RSV"        = "~/Data/Essence/RSV.xlsx"
  )
)

# ============================================================
# 2.  UTILITIES
# ============================================================

nan_to_na <- function(x) ifelse(is.nan(x), NA_real_, x)

# FIX A — cap individual model sigma_log values before ensemble combination.
# Called in run_forecast() on fc_all before build_ensemble().
# cap_h = cap * h^exp so the ceiling widens with horizon (same philosophy
# as the ensemble-level cap) but acts on each model's sigma independently,
# directly reducing within_var = sum(w_i * sigma_i^2).
# FIX 2+3 (v5): sort by date within model before assigning h so
# horizon index is always correct regardless of fc_all row order.
# Use sqrt(h) scaling uniformly (exp=0.50) — RF's internal sigma
# already uses rf_sigma_exp; this cap is a separate hard ceiling.
cap_model_sigma <- function(fc_tbl, cap) {
  if (is.null(cap)) return(fc_tbl)
  z <- qnorm(0.5 + CONFIG$ci_level / 200)
  fc_tbl %>%
    group_by(model) %>%
    arrange(date, .by_group = TRUE) %>%          # FIX 2: sort by date first
    mutate(
      h         = row_number(),                  # now correctly 1..4
      cap_h     = cap * sqrt(h),                 # FIX 3: sqrt(h) uniformly
      sigma_log = pmin(sigma_log, cap_h),
      mean_case = pmax(expm1(mu_log + 0.5 * sigma_log^2), 0),
      lower_ci  = pmax(expm1(mu_log - z * sigma_log),     0),
      upper_ci  = pmax(expm1(mu_log + z * sigma_log),     0)
    ) %>%
    select(-h, -cap_h) %>%
    ungroup()
}

# FIX 2 helper: estimate AR(1) phi from a numeric vector.
# Returns phi clamped to [0.05, 0.99] — always positive mean-reversion,
# never explosive and never treating the covariate as white noise.
.ar1_phi <- function(x) {
  x_clean <- x[!is.na(x)]
  if (length(x_clean) < 10) return(0.70)   # safe fallback for short series
  fit <- tryCatch(
    ar(x_clean, order.max = 1L, method = "yule-walker", aic = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit) || length(fit$ar) == 0) return(0.70)
  phi <- as.numeric(fit$ar[1])
  # Clamp: negative phi (oscillating) or near-zero treated as moderate decay
  max(0.05, min(0.99, phi))
}

# ============================================================
# 3.  DATA LOADER  (unchanged)
# ============================================================

load_data <- function(disease) {
  
  df        <- read_excel(CONFIG$data_paths[[disease]], sheet = 1)
  names(df) <- tolower(names(df))
  df        <- rename(df, case = data)
  
  df$date <- if (is.numeric(df$date))
    as.Date(df$date, origin = "1899-12-30")
  else
    as.Date(df$date)
  
  if (!"ww_signal_log" %in% names(df))
    stop(sprintf("[%s] Column 'ww_signal_log' not found. Present: %s",
                 disease, paste(names(df), collapse = ", ")))
  
  has_pm25 <- "pm25_mean" %in% names(df)
  if (!has_pm25)
    message(sprintf("  [%s] pm25_mean missing — run merge_aq_disease.R first.", disease))
  
  df <- df %>%
    mutate(
      case          = pmax(as.numeric(case), 0),
      # ww_signal_log is already log-transformed by the extraction script
      # (wastewater.R uses log1p(pcr_target_mic_lin) before weighted mean).
      # However the raw pcr_target_mic_lin values are ~1e-4 to 1e-3, so
      # log1p(1e-4) ≈ 1e-4 — effectively zero. The extraction script needs
      # to use log(pcr_target_mic_lin * 1e6) instead. See ww_fix note below.
      # Pipeline accepts the column as pre-transformed; fix is in wasterwater.R.
      ww_signal_log = as.numeric(ww_signal_log)
    )
  
  ww_range <- range(df$ww_signal_log, na.rm = TRUE)
  message(sprintf(
    "  [%s] ww_signal_log: range=[%.4f, %.4f] mean=%.4f (expect ~3 to 8 if correctly transformed)",
    disease, ww_range[1], ww_range[2], mean(df$ww_signal_log, na.rm = TRUE)
  ))
  if (ww_range[2] < 1.0)
    warning(sprintf(
      "[%s] ww_signal_log max=%.5f — values near zero. Re-run wasterwater.R with corrected transform.",
      disease, ww_range[2]
    ))
  
  if (has_pm25)
    df <- df %>% mutate(pm25_mean = as.numeric(pm25_mean))
  
  if (disease == "Influenza") {
    if (!"ve" %in% names(df))
      stop("[Influenza] Column 'VE' not found.")
    df <- df %>% mutate(VE = as.numeric(ve))
    df %>% select(date, case, VE, ww_signal_log, any_of("pm25_mean"))
  } else {
    df %>% select(date, case, ww_signal_log, any_of("pm25_mean"))
  }
}

# ============================================================
# 4.  BUILD TSIBBLE  (unchanged)
# ============================================================

build_tsibble <- function(df, disease) {
  
  has_pm25 <- "pm25_mean" %in% names(df)
  
  ts <- df %>%
    as_tsibble(index = date) %>%
    mutate(date = yearweek(date)) %>%
    mutate(log_case = as.double(log1p(case)))
  
  ww_na <- sum(is.na(ts$ww_signal_log))
  if (ww_na > 0) {
    ts <- ts %>%
      mutate(ww_signal_log = zoo::na.locf(ww_signal_log, na.rm = FALSE))
    message(sprintf("  [%s] ww_signal_log: %d NA(s) filled.", disease, ww_na))
  }
  ts <- ts %>%
    mutate(ww_signal_log_lag1 = lag(ww_signal_log, 1)) %>%
    mutate(ww_signal_log_lag1 = zoo::na.locf(
      ww_signal_log_lag1, na.rm = FALSE, fromLast = TRUE))
  
  if (has_pm25) {
    pm_na <- sum(is.na(ts$pm25_mean))
    if (pm_na > 0) {
      ts <- ts %>% mutate(pm25_mean = zoo::na.locf(pm25_mean, na.rm = FALSE))
      message(sprintf("  [%s] pm25_mean: %d NA(s) filled.", disease, pm_na))
    }
    ts <- ts %>%
      mutate(pm25_lag1 = lag(pm25_mean, 1)) %>%
      mutate(pm25_lag1 = zoo::na.locf(
        pm25_lag1, na.rm = FALSE, fromLast = TRUE))
  }
  
  if (disease == "Influenza") {
    ts <- ts %>%
      mutate(week_of_year = lubridate::week(as.Date(date)),
             VE           = zoo::na.locf(VE, na.rm = FALSE)) %>%
      group_by(week_of_year) %>%
      mutate(seasonal_VE = mean(VE, na.rm = TRUE)) %>%
      ungroup() %>%
      mutate(VE = ifelse(is.na(VE), seasonal_VE, VE)) %>%
      select(-week_of_year, -seasonal_VE)
  }
  
  ts
}

# ============================================================
# 5.  FABLE MODEL DEFINITIONS  (unchanged)
# ============================================================

.covar_rhs <- function(has_pm25) {
  if (has_pm25) "ww_signal_log_lag1 + pm25_lag1"
  else          "ww_signal_log_lag1"
}

fit_base_models <- function(data, disease) {
  
  has_pm25 <- "pm25_lag1" %in% names(data)
  cov      <- .covar_rhs(has_pm25)
  
  switch(disease,
         
         "Influenza" = data %>% model(
           ETS       = ETS(log_case),
           ARIMA     = ARIMA(as.formula(paste0("log_case ~ PDQ() + ", cov))),
           RW        = RW(log_case ~ drift()),
           NAIVE     = NAIVE(log_case),
           SNAIVE    = SNAIVE(log_case),
           NNETAR    = NNETAR(log_case, n.networks = 5),
           STL_ARIMA = decomposition_model(
             STL(log_case ~ season(window = "periodic")),
             ARIMA(as.formula(paste0("season_adjust ~ PDQ() + ", cov)))
           ),
           STL_ETS   = decomposition_model(
             STL(log_case ~ season(window = "periodic")),
             ETS(season_adjust)
           )
         ),
         
         # RSV v7 pool: ARIMA, ARIMAX_WW removed (v3). NAIVE now also removed —
         # v6 results showed NAIVE entered ensemble at rank 6 (wt_CRPS=0.385)
         # driven by off-season folds where last-value is hard to beat, but it
         # adds no skill during RSV season (negative skill in folds 109, 117)
         # and inflates within-model variance. SNAIVE also removed for same reason.
         # Remaining pool: ETS, RW, NNETAR(20), STL_ARIMA, STL_ETS.
         # CV selector picks top 4 from these 5 candidates.
         "RSV" = data %>% model(
           ETS       = ETS(log_case),
           RW        = RW(log_case ~ drift()),
           NNETAR    = NNETAR(log_case, n.networks = 20),
           STL_ARIMA = decomposition_model(
             STL(log_case ~ season(window = "periodic")),
             ARIMA(as.formula(paste0("season_adjust ~ PDQ() + ", cov)))
           ),
           STL_ETS   = decomposition_model(
             STL(log_case ~ season(window = "periodic")),
             ETS(season_adjust)
           )
         ),
         
         "SARS-CoV2" = data %>% model(
           ETS       = ETS(log_case),
           ARIMA     = ARIMA(as.formula(paste0("log_case ~ PDQ() + trend() + ", cov))),
           RW        = RW(log_case ~ drift()),
           NAIVE     = NAIVE(log_case),
           SNAIVE    = SNAIVE(log_case),
           NNETAR    = NNETAR(log_case, n.networks = 5),
           STL_ARIMA = decomposition_model(
             STL(log_case ~ season(window = "periodic")),
             ARIMA(as.formula(paste0("season_adjust ~ PDQ() + ", cov)))
           ),
           STL_ETS   = decomposition_model(
             STL(log_case ~ season(window = "periodic")),
             ETS(season_adjust)
           )
         ),
         
         stop("Unknown disease: ", disease)
  )
}

fit_ve_models <- function(data) {
  has_pm25 <- "pm25_lag1" %in% names(data)
  cov      <- .covar_rhs(has_pm25)
  data %>% model(
    ARIMA_VE  = ARIMA(as.formula(paste0("log_case ~ PDQ() + VE + ",    cov))),
    ARIMAX_VE = ARIMA(as.formula(paste0("log_case ~ season() + VE + ", cov)))
  )
}

fit_ww_models_rsv <- function(data) {
  # v3: ARIMAX_WW removed — had worst CV CRPS (0.34) alongside ARIMA.
  # ARIMA_WW retained: lower CRPS (0.32), simpler wastewater-only model,
  # less likely to diverge from ensemble mean at seasonal turning points.
  data %>% model(
    ARIMA_WW = ARIMA(log_case ~ PDQ() + ww_signal_log_lag1)
  )
}

# ============================================================
# NOTE: LSTM removed in v11 — all ensemble members now have
# CV-derived BMA weights. See pipeline version history.
# ============================================================

# 6h.  RANDOM FOREST  (unchanged)
# ============================================================

build_rf_features <- function(ts_data, disease, n_lag = CONFIG$rf$n_lag) {
  
  df <- ts_data %>%
    as_tibble() %>%
    mutate(date = as.Date(date),
           week_of_year = lubridate::week(date))
  
  for (k in seq_len(n_lag)) {
    df[[paste0("log_case_lag", k)]] <- dplyr::lag(df$log_case, k)
  }
  
  lag_cols  <- paste0("log_case_lag", seq_len(n_lag))
  keep_cols <- c(
    "date", "log_case", lag_cols, "ww_signal_log_lag1", "week_of_year",
    if ("pm25_lag1" %in% names(df)) "pm25_lag1",
    if (disease == "Influenza" && "VE" %in% names(df)) "VE"
  )
  
  df %>%
    select(all_of(keep_cols)) %>%
    filter(if_all(all_of(lag_cols), ~ !is.na(.x)))
}

fit_rf <- function(ts_data, disease) {
  
  feat_df   <- build_rf_features(ts_data, disease)
  feat_cols <- setdiff(names(feat_df), c("date", "log_case"))
  
  rf_fit <- ranger::ranger(
    formula       = log_case ~ .,
    data          = feat_df %>% select(log_case, all_of(feat_cols)),
    num.trees     = CONFIG$rf$num.trees,
    min.node.size = CONFIG$rf$min.node.size,
    keep.inbag    = TRUE,
    seed          = 42L
  )
  
  oob_resid <- feat_df$log_case - rf_fit$predictions
  oob_sigma <- pmax(sd(oob_resid, na.rm = TRUE), 1e-6)
  
  list(rf_fit    = rf_fit,
       feat_cols = feat_cols,
       oob_sigma = oob_sigma,
       feat_df   = feat_df)
}

forecast_rf <- function(rf_obj, future_data, disease,
                        ci_level = CONFIG$ci_level) {
  
  z          <- qnorm(0.5 + ci_level / 200)
  n_lag      <- CONFIG$rf$n_lag
  feat_cols  <- rf_obj$feat_cols
  sigma_base <- rf_obj$oob_sigma
  lag_buffer <- rev(tail(rf_obj$feat_df$log_case, n_lag))
  
  # FIX B: disease-specific sigma growth exponent.
  # RSV 0.40 < sqrt (0.50): strong seasonal autocorrelation means errors
  # compound more slowly than a random walk assumption implies.
  sigma_exp <- CONFIG$rf_sigma_exp[[disease]]
  
  fut_tbl <- future_data %>% as_tibble() %>% mutate(date = as.Date(date))
  n_h     <- nrow(fut_tbl)
  preds   <- numeric(n_h)
  
  for (h in seq_len(n_h)) {
    row_data <- list()
    for (k in seq_len(n_lag))
      row_data[[paste0("log_case_lag", k)]] <- lag_buffer[k]
    row_data[["ww_signal_log_lag1"]] <- fut_tbl$ww_signal_log_lag1[[h]]
    row_data[["week_of_year"]]       <- lubridate::week(as.Date(fut_tbl$date[[h]]))
    if ("pm25_lag1" %in% feat_cols)
      row_data[["pm25_lag1"]] <- fut_tbl$pm25_lag1[[h]]
    if ("VE" %in% feat_cols)
      row_data[["VE"]] <- fut_tbl$VE[[h]]
    
    new_row    <- as.data.frame(row_data)[, feat_cols, drop = FALSE]
    pred_log   <- predict(rf_obj$rf_fit, data = new_row)$predictions
    preds[h]   <- pred_log
    lag_buffer <- c(pred_log, lag_buffer[-n_lag])
  }
  
  sigma_h <- sigma_base * seq_len(n_h)^sigma_exp   # FIX B: h^exp not sqrt(h)
  
  tibble(
    date      = as.Date(fut_tbl$date),
    model     = "RF",
    mu_log    = preds,
    sigma_log = sigma_h,
    mean_case = pmax(expm1(preds + 0.5 * sigma_h^2), 0),
    lower_ci  = pmax(expm1(preds - z * sigma_h),     0),
    upper_ci  = pmax(expm1(preds + z * sigma_h),     0)
  )
}

# ============================================================
# 7.  CRPS EXTRACTION  (unchanged)
# ============================================================

extract_cv_crps <- function(fc, actuals) {
  
  actual_tbl <- actuals %>% as_tibble() %>% select(date, log_case)
  
  tryCatch({
    fc %>%
      as_tibble() %>%
      select(.model, date, .distribution = log_case) %>%
      left_join(actual_tbl, by = "date") %>%
      filter(!is.na(log_case)) %>%
      rowwise() %>%
      mutate(
        mu    = tryCatch(mean(.distribution),           error = function(e) NA_real_),
        sigma = tryCatch(sqrt(variance(.distribution)), error = function(e) NA_real_)
      ) %>%
      ungroup() %>%
      mutate(
        z    = (log_case - mu) / sigma,
        CRPS = ifelse(
          !is.na(sigma) & sigma > 0,
          sigma * (z * (2*pnorm(z) - 1) + 2*dnorm(z) - 1/sqrt(pi)),
          NA_real_
        )
      ) %>%
      select(.model, date, CRPS)
    
  }, error = function(e) {
    message("    extract_cv_crps error: ", e$message)
    fc %>% as_tibble() %>% select(.model, date) %>% mutate(CRPS = NA_real_)
  })
}

# ============================================================
# 8.  BMA WEIGHTS  (unchanged)
# ============================================================

compute_bma_weights <- function(fold_acc, disease, temperature) {
  
  fold_cuts <- sort(unique(fold_acc$fold_cut))
  n_f       <- length(fold_cuts)
  rw        <- exp(seq(0, 1, length.out = max(n_f, 1)))
  rw        <- rw / sum(rw)
  rw_tbl    <- tibble(fold_cut = fold_cuts, recency_weight = rw)
  
  scores <- fold_acc %>%
    left_join(rw_tbl, by = "fold_cut") %>%
    group_by(.model) %>%
    summarise(
      mean_crps = nan_to_na(
        sum(CRPS * recency_weight, na.rm = TRUE) /
          sum(recency_weight[!is.na(CRPS)], na.rm = TRUE)
      ),
      mae_proxy = nan_to_na(
        sum(MAE * recency_weight, na.rm = TRUE) /
          sum(recency_weight[!is.na(MAE)], na.rm = TRUE)
      ),
      n_folds = n(),
      .groups = "drop"
    ) %>%
    mutate(
      effective_score = coalesce(mean_crps, mae_proxy),
      score_source    = case_when(
        !is.na(mean_crps) ~ "CRPS",
        !is.na(mae_proxy) ~ "MAE-proxy",
        TRUE              ~ "unavailable"
      )
    ) %>%
    filter(score_source != "unavailable")
  
  message(sprintf("  [%s] Scores (temp=%.2f):", disease, temperature))
  print(scores %>% select(.model, effective_score, score_source))
  
  min_s <- min(scores$effective_score, na.rm = TRUE)
  sd_s  <- sd(scores$effective_score,  na.rm = TRUE)
  if (is.na(sd_s) || sd_s < 1e-6) sd_s <- 1
  
  weights <- scores %>%
    mutate(
      scaled_score  = ((effective_score - min_s) / sd_s) * temperature,
      unnorm_weight = exp(-scaled_score),
      weight        = unnorm_weight / sum(unnorm_weight, na.rm = TRUE)
    ) %>%
    filter(!is.na(weight), weight > 0) %>%
    select(.model, weight)
  
  message(sprintf("  [%s] BMA weights:", disease))
  print(weights %>% arrange(desc(weight)) %>% mutate(pct = round(weight * 100, 1)))
  
  weights
}

# Fallback: conservative prior for models not yet CV-scored
recompute_bma_with_aux <- function(fold_acc, aux_names,
                                   top_models, disease, temperature) {
  
  base_fold <- fold_acc %>% filter(.model %in% top_models)
  fold_cuts <- sort(unique(base_fold$fold_cut))
  n_f       <- length(fold_cuts)
  rw        <- exp(seq(0, 1, length.out = max(n_f, 1)))
  rw        <- rw / sum(rw)
  rw_tbl    <- tibble(fold_cut = fold_cuts, recency_weight = rw)
  
  prior_crps <- base_fold %>%
    left_join(rw_tbl, by = "fold_cut") %>%
    group_by(.model) %>%
    summarise(
      wcrps = nan_to_na(
        sum(CRPS * recency_weight, na.rm = TRUE) /
          sum(recency_weight[!is.na(CRPS)], na.rm = TRUE)
      ), .groups = "drop"
    ) %>%
    summarise(prior = min(wcrps, na.rm = TRUE)) %>%
    pull(prior)
  
  aux_rows <- expand.grid(
    .model = aux_names, fold_cut = fold_cuts, stringsAsFactors = FALSE
  ) %>%
    as_tibble() %>%
    mutate(CRPS = prior_crps, MAE = NA_real_, RMSE = NA_real_)
  
  combined <- bind_rows(base_fold, aux_rows)
  compute_bma_weights(combined, disease, temperature)
}

# ============================================================
# 9.  CROSS-VALIDATION
#     FIX 3: auxiliary fable models (VE / WW) now scored here
#     alongside base models and RF. LSTM scored via cv_lstm().
# ============================================================

run_cv <- function(ts_data, disease) {
  
  n    <- nrow(ts_data)
  init <- round(n * CONFIG$cv_init_frac)
  h    <- CONFIG$horizon
  
  cutpoints <- seq(init, n - h, by = CONFIG$cv_step)
  if (length(cutpoints) == 0)
    stop(sprintf("[%s] Series too short for CV (n=%d, init=%d).", disease, n, init))
  
  message(sprintf(
    "  [%s] CV: init=%d, step=%d, folds=%d, workers=%d",
    disease, init, CONFIG$cv_step, length(cutpoints), CONFIG$cv_workers
  ))
  
  plan(multisession, workers = CONFIG$cv_workers)
  on.exit(plan(sequential), add = TRUE)
  
  fold_acc <- future_map_dfr(cutpoints, function(cut) {
    
    train <- ts_data %>% slice(1:cut)
    test  <- ts_data %>% slice((cut + 1):min(cut + h, n))
    
    # ---- Base fable models ----
    fc <- tryCatch(
      fit_base_models(train, disease) %>% forecast(new_data = test),
      error = function(e) {
        message(sprintf("    Fold cut=%d (fable): %s", cut, e$message)); NULL
      }
    )
    
    base_acc <- if (!is.null(fc)) {
      crps_horizon    <- extract_cv_crps(fc, test) %>% mutate(fold_cut = cut)
      actual_log_base <- test %>%
        as_tibble() %>% mutate(date = as.Date(date)) %>% select(date, log_case)
      fc_tbl <- fc %>%
        as_tibble() %>%
        rowwise() %>%
        mutate(mu_log = tryCatch(mean(log_case), error = function(e) NA_real_),
               .groups = "drop") %>%
        ungroup() %>%
        mutate(date = as.Date(date)) %>%
        select(.model, date, mu_log)
      crps_horizon %>%
        left_join(fc_tbl,           by = c(".model", "date")) %>%
        left_join(actual_log_base,  by = "date") %>%
        mutate(MAE  = abs(log_case - mu_log),
               RMSE = sqrt((log_case - mu_log)^2)) %>%
        select(.model, date, CRPS, MAE, RMSE) %>%
        mutate(fold_cut = cut)
    }
    
    # ---- FIX 3: Auxiliary fable models scored per fold ----
    # Influenza: ARIMA_VE, ARIMAX_VE
    # RSV:       ARIMA_WW only (v3: ARIMAX_WW removed — worst CRPS)
    # SARS-CoV2: no aux models
    .empty_acc <- tibble(.model = character(), fold_cut = numeric(),
                         date = as.Date(character()),
                         CRPS = numeric(), MAE = numeric(), RMSE = numeric())
    
    aux_acc <- if (!disease %in% c("Influenza", "RSV")) {
      .empty_acc
    } else {
      tryCatch({
        fc_aux <- if (disease == "Influenza") {
          fit_ve_models(train) %>% forecast(new_data = test)
        } else {
          # RSV v3: only ARIMA_WW
          fit_ww_models_rsv(train) %>% forecast(new_data = test)
        }
        
        crps_aux   <- extract_cv_crps(fc_aux, test) %>% mutate(fold_cut = cut)
        actual_aux <- test %>%
          as_tibble() %>% mutate(date = as.Date(date)) %>% select(date, log_case)
        fc_aux_tbl <- fc_aux %>%
          as_tibble() %>%
          rowwise() %>%
          mutate(mu_log = tryCatch(mean(log_case), error = function(e) NA_real_),
                 .groups = "drop") %>%
          ungroup() %>%
          mutate(date = as.Date(date)) %>%
          select(.model, date, mu_log)
        crps_aux %>%
          left_join(fc_aux_tbl, by = c(".model", "date")) %>%
          left_join(actual_aux,  by = "date") %>%
          mutate(MAE      = abs(log_case - mu_log),
                 RMSE     = sqrt((log_case - mu_log)^2),
                 fold_cut = cut) %>%
          select(.model, fold_cut, date, CRPS, MAE, RMSE)
      }, error = function(e) {
        message(sprintf("    Aux fold cut=%d failed: %s", cut, e$message))
        .empty_acc
      })
    }
    
    # ---- RF fold scoring (unchanged) ----
    rf_acc <- tryCatch({
      rf_obj <- fit_rf(train, disease)
      fut_cv <- build_future(train, disease)
      rf_fc  <- forecast_rf(rf_obj, fut_cv, disease)
      actual_log <- test %>%
        as_tibble() %>% mutate(date = as.Date(date)) %>% select(date, log_case)
      rf_fc %>%
        left_join(actual_log, by = "date") %>%
        filter(!is.na(log_case)) %>%
        mutate(
          z    = (log_case - mu_log) / sigma_log,
          CRPS = sigma_log * (z * (2*pnorm(z) - 1) + 2*dnorm(z) - 1/sqrt(pi)),
          MAE  = abs(log_case - mu_log),
          RMSE = sqrt((log_case - mu_log)^2),
          .model   = "RF",
          fold_cut = cut
        ) %>%
        select(.model, fold_cut, date, CRPS, MAE, RMSE)
    }, error = function(e) {
      message(sprintf("    RF fold cut=%d failed: %s", cut, e$message))
      tibble(.model = "RF", fold_cut = cut, date = as.Date(NA),
             CRPS = NA_real_, MAE = NA_real_, RMSE = NA_real_)
    })
    
    bind_rows(base_acc, aux_acc, rf_acc)
    
  }, .options = furrr_options(seed = TRUE))
  
  fold_acc <- fold_acc %>%
    mutate(across(any_of(c("CRPS", "MAE", "RMSE")), nan_to_na))
  
  # Collapse to per-fold summary for model ranking
  fold_acc_smry <- fold_acc %>%
    group_by(.model, fold_cut) %>%
    summarise(
      CRPS = nan_to_na(mean(CRPS, na.rm = TRUE)),
      MAE  = nan_to_na(mean(MAE,  na.rm = TRUE)),
      RMSE = nan_to_na(mean(RMSE, na.rm = TRUE)),
      .groups = "drop"
    )
  
  fold_cuts <- sort(unique(fold_acc_smry$fold_cut))
  n_f       <- length(fold_cuts)
  rw        <- exp(seq(0, 1, length.out = max(n_f, 1)))
  rw        <- rw / sum(rw)
  rw_tbl    <- tibble(fold_cut = fold_cuts, recency_weight = rw)
  
  model_summary <- fold_acc_smry %>%
    left_join(rw_tbl, by = "fold_cut") %>%
    group_by(.model) %>%
    summarise(
      mean_CRPS = nan_to_na(sum(CRPS * recency_weight, na.rm = TRUE) /
                              sum(recency_weight[!is.na(CRPS)], na.rm = TRUE)),
      mean_MAE  = nan_to_na(sum(MAE  * recency_weight, na.rm = TRUE) /
                              sum(recency_weight[!is.na(MAE)],  na.rm = TRUE)),
      mean_RMSE = nan_to_na(sum(RMSE * recency_weight, na.rm = TRUE) /
                              sum(recency_weight[!is.na(RMSE)], na.rm = TRUE)),
      n_folds   = n(),
      .groups   = "drop"
    ) %>%
    mutate(sort_key = coalesce(mean_CRPS, mean_MAE, mean_RMSE)) %>%
    arrange(sort_key)
  
  # top_n excludes RF, LSTM, and aux models — they enter ensemble separately.
  # v3: ARIMAX_WW removed from exclusion list (no longer used for RSV).
  aux_model_names <- c("ARIMA_VE", "ARIMAX_VE", "ARIMA_WW",
                       "RF", "LSTM")
  n_top      <- CONFIG$top_n_models[[disease]]
  top_models <- model_summary %>%
    filter(!is.na(sort_key), !.model %in% aux_model_names) %>%
    slice_head(n = n_top) %>%
    pull(.model)
  
  message(sprintf("  [%s] Top base models (n=%d): %s",
                  disease, n_top, paste(top_models, collapse = ", ")))
  
  list(
    top_models    = top_models,
    model_summary = model_summary,
    cv_accuracy   = fold_acc,
    cutpoints     = cutpoints
  )
}

# ============================================================
# 10. EXTRACT FABLE FORECASTS  (unchanged)
# ============================================================

extract_fc <- function(fit_obj, new_data_arg = NULL, ci_level) {
  
  ci_level_num <- as.numeric(ci_level[[1]])
  ci_col       <- paste0(ci_level_num, "%")
  ci_lower_col <- paste0(ci_col, "_lower")
  ci_upper_col <- paste0(ci_col, "_upper")
  z            <- qnorm(0.5 + ci_level_num / 200)
  
  fc_raw <- if (is.null(new_data_arg))
    fit_obj %>% forecast(h = CONFIG$horizon)
  else
    fit_obj %>% forecast(new_data = new_data_arg)
  
  fc_raw %>%
    hilo(level = ci_level_num) %>%
    unpack_hilo(!!ci_col) %>%
    rename(
      mu_log    = .mean,
      lower_log = !!sym(ci_lower_col),
      upper_log = !!sym(ci_upper_col)
    ) %>%
    mutate(
      mu_log    = as.double(mu_log),
      lower_log = as.double(lower_log),
      upper_log = as.double(upper_log),
      model     = .model,
      date      = as.Date(date),
      sigma_log = pmax((upper_log - lower_log) / (2 * z), 1e-6),
      mean_case = pmax(expm1(mu_log + 0.5 * sigma_log^2), 0),
      lower_ci  = pmax(expm1(mu_log - z * sigma_log),     0),
      upper_ci  = pmax(expm1(mu_log + z * sigma_log),     0)
    ) %>%
    as_tibble() %>%
    select(date, model, mu_log, sigma_log, mean_case, lower_ci, upper_ci)
}

# ============================================================
# 11. FUTURE DATA FRAME
#     FIX 2: AR(1) covariate decay replaces fixed linear decay.
#     phi estimated from the actual autocorrelation of each
#     covariate series; geometric decay: last_val * phi^(h-1) +
#     mean * (1 - phi^(h-1)).
# ============================================================

build_future <- function(ts_data, disease) {
  
  has_pm25  <- "pm25_lag1" %in% names(ts_data)
  last_date <- max(as.Date(ts_data$date) + 6)
  h         <- CONFIG$horizon
  
  # FIX 2: estimate AR(1) phi for wastewater signal
  phi_ww  <- .ar1_phi(ts_data$ww_signal_log)
  last_ww <- ts_data %>%
    filter(!is.na(ww_signal_log)) %>% slice_tail(n = 1) %>% pull(ww_signal_log)
  mean_ww <- mean(ts_data$ww_signal_log, na.rm = TRUE)
  
  # Geometric decay: val_h = mean + phi^(h-1) * (last - mean)
  decay_ww <- phi_ww^(seq_len(h) - 1)
  ww_h     <- mean_ww + decay_ww * (last_ww - mean_ww)
  
  message(sprintf(
    "  [%s] WW decay: phi=%.3f | h1=%.3f → h%d=%.3f (mean=%.3f)",
    disease, phi_ww, last_ww, h, ww_h[h], mean_ww
  ))
  
  fut <- tibble(
    date               = seq(last_date + 7, by = "week", length.out = h),
    ww_signal_log_lag1 = ww_h
  )
  
  if (has_pm25) {
    phi_pm  <- .ar1_phi(ts_data$pm25_lag1)
    last_pm <- ts_data %>%
      filter(!is.na(pm25_lag1)) %>% slice_tail(n = 1) %>% pull(pm25_lag1)
    mean_pm <- mean(ts_data$pm25_lag1, na.rm = TRUE)
    
    decay_pm <- phi_pm^(seq_len(h) - 1)
    pm_h     <- mean_pm + decay_pm * (last_pm - mean_pm)
    fut      <- fut %>% mutate(pm25_lag1 = pm_h)
    
    message(sprintf(
      "  [%s] PM2.5 decay: phi=%.3f | h1=%.2f → h%d=%.2f (mean=%.2f)",
      disease, phi_pm, last_pm, h, pm_h[h], mean_pm
    ))
  }
  
  fut <- fut %>%
    mutate(date = yearweek(date)) %>%
    as_tsibble(index = date)
  
  if (disease == "Influenza")
    fut <- fut %>% mutate(VE = mean(ts_data$VE, na.rm = TRUE))
  
  fut
}

# ============================================================
# 12. ENSEMBLE  (BMA log-space mixture)
#     FIX 1: cummax(sigma_mix) applied after sigma_cap so the
#     full ensemble output is guaranteed to never narrow.
# ============================================================

build_ensemble <- function(fc_models, weights, ci_level, sigma_cap = NULL) {
  
  z <- qnorm(0.5 + ci_level / 200)
  
  joined <- fc_models %>%
    left_join(weights, by = c("model" = ".model")) %>%
    filter(!is.na(weight)) %>%
    group_by(date) %>%
    mutate(weight = weight / sum(weight)) %>%
    ungroup()
  
  if (nrow(joined) == 0)
    stop("build_ensemble: no rows after weight join — check model name alignment.")
  
  mu_mix_tbl <- joined %>%
    group_by(date) %>%
    summarise(mu_mix = sum(weight * mu_log), .groups = "drop") %>%
    arrange(date)
  
  # Horizon CRPS fix: smooth mu_mix across horizons to prevent
  # discontinuous jumps between h=1 and h=2 that cause RSV's
  # Week+2 CRPS spike and COVID's Week+2 dip.
  # Uses exponential smoothing (alpha=0.7) — recent horizon gets
  # 70% weight, prior gets 30%. Preserves direction of forecasts
  # while damping abrupt between-horizon mean shifts.
  # Only applied when there are exactly 4 horizon steps.
  if (nrow(mu_mix_tbl) == 4L) {
    alpha <- 0.70
    mu_smooth <- mu_mix_tbl$mu_mix
    for (i in 2:4)
      mu_smooth[i] <- alpha * mu_smooth[i] + (1 - alpha) * mu_smooth[i - 1]
    mu_mix_tbl$mu_mix <- mu_smooth
  }
  
  result <- joined %>%
    left_join(mu_mix_tbl, by = "date") %>%
    group_by(date) %>%
    summarise(
      mu_mix      = first(mu_mix),
      within_var  = sum(weight * sigma_log^2),
      between_var = sum(weight * (mu_log - mu_mix)^2),
      .groups     = "drop"
    ) %>%
    mutate(sigma_mix = sqrt(within_var + between_var))
  
  # FIX C (v5): cummax() first for monotone sigma, then hard per-horizon cap.
  result$sigma_mix <- cummax(result$sigma_mix)
  
  if (!is.null(sigma_cap)) {
    h_seq            <- seq_len(nrow(result))
    sigma_cap_h      <- sigma_cap * sqrt(h_seq)
    result$sigma_mix <- pmin(result$sigma_mix, sigma_cap_h)
  }
  
  result <- result %>%
    mutate(
      mean_case = pmax(expm1(mu_mix + 0.5 * sigma_mix^2), 0),
      lower_ci  = pmax(expm1(mu_mix - z * sigma_mix),     0),
      upper_ci  = pmax(expm1(mu_mix + z * sigma_mix),     0),
      model     = "Ensemble"
    )
  
  # BUG 2 FIX: CI width can still be non-monotone even with monotone sigma_mix
  # when mu_mix shifts between horizons (lower_ci rises faster than upper_ci).
  # Enforce monotone CI WIDTH directly: if width narrows at h, expand upper_ci
  # and shrink lower_ci symmetrically around mean_case to restore target width.
  ci_widths      <- result$upper_ci - result$lower_ci
  ci_widths_mono <- cummax(ci_widths)
  needs_fix      <- ci_widths_mono > ci_widths
  if (any(needs_fix)) {
    delta                    <- (ci_widths_mono - ci_widths) / 2
    result$upper_ci[needs_fix] <- result$upper_ci[needs_fix] + delta[needs_fix]
    result$lower_ci[needs_fix] <- pmax(
      result$lower_ci[needs_fix] - delta[needs_fix], 0
    )
  }
  
  result %>%
    select(date, mean_case, lower_ci, upper_ci, model)
}

# ============================================================
# 13. MAIN PIPELINE
#     FIX 3: aux models (VE / WW) are now CV-scored in run_cv()
#     so compute_bma_weights() receives real fold CRPS for them.
#     recompute_bma_with_aux() is only called as a fallback if
#     the aux CV rows are absent (e.g. all folds errored).
#     LSTM CV is always run (run_cv = TRUE).
# ============================================================

run_forecast <- function(disease) {
  
  message("\n", strrep("=", 60))
  message("  Forecasting: ", disease)
  message(strrep("=", 60))
  t0 <- proc.time()
  
  ts_data     <- load_data(disease) %>% build_tsibble(disease)
  temperature <- CONFIG$bma_temperature[[disease]]
  
  message(sprintf("  [%s] Series: %d weeks | PM2.5: %s",
                  disease, nrow(ts_data),
                  if ("pm25_lag1" %in% names(ts_data)) "yes" else "no"))
  
  # run_cv() now scores base + aux + RF in one parallel pass
  selection  <- run_cv(ts_data, disease)
  top_models <- selection$top_models
  
  # Collapse per-horizon CV rows to per-fold for BMA weighting
  cv_acc_fold <- selection$cv_accuracy %>%
    group_by(.model, fold_cut) %>%
    summarise(
      CRPS = nan_to_na(mean(CRPS, na.rm = TRUE)),
      MAE  = nan_to_na(mean(MAE,  na.rm = TRUE)),
      RMSE = nan_to_na(mean(RMSE, na.rm = TRUE)),
      .groups = "drop"
    )
  
  # Identify aux model names expected for this disease
  # v3: RSV aux reduced to ARIMA_WW only (ARIMAX_WW removed — worst CRPS)
  aux_names <- if (disease == "Influenza") {
    c("ARIMA_VE", "ARIMAX_VE")
  } else if (disease == "RSV") {
    c("ARIMA_WW")
  } else {
    character(0)
  }
  
  # Determine which aux models have real CV scores vs need the prior fallback
  aux_cv_models   <- intersect(aux_names, unique(cv_acc_fold$.model))
  aux_prior_models <- setdiff(aux_names, aux_cv_models)
  
  # Initial BMA: top base models + any aux with real CV scores
  bma_seed_models <- c(top_models, aux_cv_models)
  bma_weights <- if (length(aux_prior_models) > 0) {
    # Some aux models need the prior — use recompute_bma_with_aux as fallback
    recompute_bma_with_aux(
      fold_acc    = cv_acc_fold %>% filter(.model %in% c(top_models, aux_cv_models)),
      aux_names   = aux_prior_models,
      top_models  = c(top_models, aux_cv_models),
      disease     = disease,
      temperature = temperature
    )
  } else {
    compute_bma_weights(
      fold_acc    = cv_acc_fold %>% filter(.model %in% bma_seed_models),
      disease     = disease,
      temperature = temperature
    )
  }
  message(sprintf("  [%s] Initial BMA from: %s",
                  disease, paste(bma_seed_models, collapse = ", ")))
  
  future_data <- build_future(ts_data, disease)
  
  # ---- Final fit on full series: base + aux ----
  fit_base <- fit_base_models(ts_data, disease) %>% select(all_of(top_models))
  fc_all   <- extract_fc(fit_base, new_data_arg = future_data,
                         ci_level = CONFIG$ci_level)
  
  if (disease == "Influenza") {
    fc_ve <- tryCatch(
      extract_fc(fit_ve_models(ts_data),
                 new_data_arg = future_data, ci_level = CONFIG$ci_level),
      error = function(e) {
        message("  [Influenza] VE models failed: ", e$message); NULL
      }
    )
    if (!is.null(fc_ve)) fc_all <- bind_rows(fc_all, fc_ve)
  }
  
  if (disease == "RSV") {
    fc_ww <- tryCatch(
      extract_fc(fit_ww_models_rsv(ts_data),
                 new_data_arg = future_data, ci_level = CONFIG$ci_level),
      error = function(e) {
        message("  [RSV] WW models failed: ", e$message); NULL
      }
    )
    if (!is.null(fc_ww)) fc_all <- bind_rows(fc_all, fc_ww)
  }
  
  # ---- Random Forest ----
  message(sprintf("  [%s] Fitting Random Forest…", disease))
  rf_fc <- tryCatch({
    rf_obj <- fit_rf(ts_data, disease)
    forecast_rf(rf_obj, future_data, disease)
  }, error = function(e) {
    message(sprintf("  [%s] RF failed: %s", disease, e$message)); NULL
  })
  
  if (!is.null(rf_fc) && nrow(rf_fc) > 0) {
    fc_all <- bind_rows(fc_all, rf_fc)
    # RF has real CV scores from run_cv(); include it in BMA
    all_scored <- c(bma_seed_models, "RF")
    bma_weights <- compute_bma_weights(
      fold_acc    = cv_acc_fold %>% filter(.model %in% all_scored),
      disease     = disease,
      temperature = temperature
    )
    message(sprintf("  [%s] BMA updated with RF.", disease))
  }
  
  
  
  # ---- Finalise weights ----
  present_models <- unique(fc_all$model)
  weights_final  <- bma_weights %>%
    filter(.model %in% present_models) %>%
    mutate(weight = weight / sum(weight))
  
  if (nrow(weights_final) == 0)
    stop(sprintf(
      "[%s] No valid BMA weights.\n  fc_all: [%s]\n  bma: [%s]",
      disease,
      paste(present_models,     collapse = ", "),
      paste(bma_weights$.model, collapse = ", ")
    ))
  
  message(sprintf("  [%s] Final ensemble: %s",
                  disease, paste(weights_final$.model, collapse = ", ")))
  
  sigma_cap       <- CONFIG$sigma_cap[[disease]]
  sigma_model_cap <- CONFIG$sigma_model_cap[[disease]]
  
  # FIX A: cap individual model sigma_log before ensemble combination.
  # Directly reduces within_var = sum(w_i * sigma_i^2) which v4 diagnosis
  # confirmed is the primary driver of wide RSV intervals (70-80% of total).
  if (!is.null(sigma_model_cap)) {
    message(sprintf("  [%s] Applying per-model sigma cap = %.2f (v5 FIX 1).",
                    disease, sigma_model_cap))
    fc_all <- cap_model_sigma(fc_all, cap = sigma_model_cap)
  }
  
  ensemble <- build_ensemble(fc_all, weights_final, CONFIG$ci_level,
                             sigma_cap = sigma_cap)
  
  elapsed <- round((proc.time() - t0)[["elapsed"]], 1)
  message(sprintf("  [%s] Done in %.1f s.", disease, elapsed))
  
  list(
    ts_data           = ts_data,
    fc_ensemble       = ensemble,
    fc_models         = fc_all,
    cv_accuracy       = selection$cv_accuracy,
    bma_weights       = weights_final,   # weights_final = filtered to present models only
    bma_weights_full  = bma_weights,     # full pre-filter weights for debugging
    model_summary     = selection$model_summary,
    insample_accuracy = fit_base %>% accuracy(),
    top_models        = top_models
  )
}

# ============================================================
# 14. RUN ALL DISEASES
# ============================================================

results_all <- CONFIG$diseases %>%
  set_names() %>%
  map(run_forecast)

# ============================================================
# 15. SAVE
# ============================================================

precomputed_forecasts <- list(
  "Influenza" = list(
    ts_data = results_all$Influenza$ts_data,
    fc      = results_all$Influenza$fc_ensemble
  ),
  "RSV" = list(
    ts_data = results_all$RSV$ts_data,
    fc      = results_all$RSV$fc_ensemble
  ),
  "SARS-CoV2" = list(
    ts_data = results_all[["SARS-CoV2"]]$ts_data,
    fc      = results_all[["SARS-CoV2"]]$fc_ensemble
  )
)

saveRDS(precomputed_forecasts, "~/Data/Essence/precomputed_forecasts_new.rds")
message("Saved precomputed_forecasts_new.rds")
saveRDS(results_all, "~/Data/Essence/results_all_new.rds")
message("Saved results_all_new.rds")

# ============================================================
# 16. DIAGNOSTICS  (unchanged + new AR1 phi report)
# ============================================================

cat("\n", strrep("=", 60), "\n")
cat("FORECAST OUTPUTS\n")
cat(strrep("=", 60), "\n")
for (d in names(precomputed_forecasts)) {
  cat(sprintf("\n[%s]  fc rows: %d  |  ts_data rows: %d\n",
              d,
              nrow(precomputed_forecasts[[d]]$fc),
              nrow(precomputed_forecasts[[d]]$ts_data)))
  print(precomputed_forecasts[[d]]$fc)
}

cat("\n", strrep("=", 60), "\n")
cat("BMA WEIGHTS\n")
cat(strrep("=", 60), "\n")
for (d in CONFIG$diseases) {
  cat(sprintf("\n[%s]\n", d))
  results_all[[d]]$bma_weights %>%
    arrange(desc(weight)) %>%
    mutate(weight_pct = round(weight * 100, 1)) %>%
    print()
}

cat("\n", strrep("=", 60), "\n")
cat("CV MODEL PERFORMANCE  (includes RF + aux models)\n")
cat(strrep("=", 60), "\n")
for (d in CONFIG$diseases) {
  cat(sprintf("\n[%s]\n", d))
  results_all[[d]]$model_summary %>%
    select(.model, mean_CRPS, mean_MAE, sort_key, n_folds) %>%
    print()
}

cat("\n", strrep("=", 60), "\n")
cat("CRPS / MAE COVERAGE  (aux models now have real CRPS)\n")
cat(strrep("=", 60), "\n")
for (d in CONFIG$diseases) {
  cat(sprintf("\n[%s]\n", d))
  results_all[[d]]$cv_accuracy %>%
    group_by(.model) %>%
    summarise(
      total_folds = n(),
      crps_valid  = sum(!is.na(CRPS)),
      mean_crps   = round(nan_to_na(mean(CRPS, na.rm = TRUE)), 4),
      mean_mae    = round(nan_to_na(mean(MAE,  na.rm = TRUE)), 4),
      .groups     = "drop"
    ) %>%
    arrange(desc(crps_valid), mean_mae) %>%
    print()
}

cat("\n", strrep("=", 60), "\n")
cat("CI WIDTH DIAGNOSTIC\n")
cat("  Targets: Influenza <1.0 | RSV <1.2 | SARS-CoV2 <1.0\n")
cat(strrep("=", 60), "\n")
for (d in names(precomputed_forecasts)) {
  fc       <- precomputed_forecasts[[d]]$fc
  ci_ratio <- mean((fc$upper_ci - fc$lower_ci) / pmax(fc$mean_case, 0.5),
                   na.rm = TRUE)
  min_lower <- min(fc$lower_ci, na.rm = TRUE)
  flag <- dplyr::case_when(
    d == "RSV" && ci_ratio > 1.2 ~ " << WIDE",
    min_lower <= 0               ~ " << LOWER BOUND ISSUE",
    TRUE                         ~ ""
  )
  cat(sprintf("  [%s]  CI/mean=%.2f  |  min_lower=%.2f%s\n",
              d, ci_ratio, min_lower, flag))
}

cat("\n", strrep("=", 60), "\n")
cat("HORIZON MONOTONICITY CHECK  (FIX 1 validation)\n")
cat("  All sigma_mix values should be non-decreasing across h=1..4\n")
cat(strrep("=", 60), "\n")
for (d in CONFIG$diseases) {
  fc_m <- results_all[[d]]$fc_models
  ens_dates <- results_all[[d]]$fc_ensemble$date
  # Re-derive sigma_mix from ensemble output isn't directly stored;
  # check monotonicity via CI widths as a proxy
  fc_ens <- results_all[[d]]$fc_ensemble
  widths  <- fc_ens$upper_ci - fc_ens$lower_ci
  mono    <- all(diff(widths) >= -1e-6)
  cat(sprintf("  [%s]  CI widths h1..h4: %s  |  monotone: %s\n",
              d,
              paste(round(widths, 2), collapse = " → "),
              if (mono) "YES ✓" else "NO — check sigma_cap interaction"))
}

cat("\n", strrep("=", 60), "\n")
cat("AR(1) COVARIATE DECAY PARAMETERS  (FIX 2)\n")
cat(strrep("=", 60), "\n")
for (d in CONFIG$diseases) {
  ts <- results_all[[d]]$ts_data
  phi_ww <- .ar1_phi(as.numeric(ts$ww_signal_log))
  cat(sprintf("  [%s]  ww phi=%.3f", d, phi_ww))
  if ("pm25_lag1" %in% names(ts)) {
    phi_pm <- .ar1_phi(as.numeric(ts$pm25_lag1))
    cat(sprintf("  |  pm25 phi=%.3f", phi_pm))
  }
  cat("\n")
}


cat("\n", strrep("=", 60), "\n")
cat("RF BMA WEIGHT BY DISEASE  (v11: LSTM removed)\n")
cat(strrep("=", 60), "\n")
for (d in CONFIG$diseases) {
  wts   <- results_all[[d]]$bma_weights
  rf_w  <- wts %>% filter(.model == "RF") %>% pull(weight)
  rf_pct <- if (length(rf_w) > 0) round(rf_w * 100, 1) else 0
  lstm_present <- "LSTM" %in% wts$.model
  cat(sprintf("  [%s]  RF: %s%%  |  LSTM: %s\n",
              d,
              if (rf_pct > 0) rf_pct else "absent",
              if (lstm_present) "WARNING: still present — check pipeline" else "absent (removed v11) ✓"))
}

cat("\n", strrep("=", 60), "\n")
cat("RSV VARIANCE DECOMPOSITION QUICK CHECK  (v4 fix validation)\n")
cat("  Target: between_var < 30% of total at all horizons\n")
cat("  Target: CI width ratio < 1.4 (v2 baseline: 1.659)\n")
cat(strrep("=", 60), "\n")

rsv_fc  <- results_all$RSV$fc_models
rsv_wts <- results_all$RSV$bma_weights

# Show per-model sigma_log range to confirm FIX A is biting
cat("\n  Per-model sigma_log (after v5 cap=%.2f*sqrt(h)):\n",
    CONFIG$sigma_model_cap[["RSV"]])
cat(sprintf("  %-14s  %8s  %8s  %8s  %8s  %8s\n",
            "Model", "wt(%)", "h1_sig", "h2_sig", "h3_sig", "h4_sig"))
cap_vals <- CONFIG$sigma_model_cap[["RSV"]] * sqrt(1:4)
cat(sprintf("  %-14s  %8s  %8.4f  %8.4f  %8.4f  %8.4f\n",
            "[cap ceiling]", "", cap_vals[1], cap_vals[2],
            cap_vals[3], cap_vals[4]))
rsv_fc %>%
  group_by(model) %>%
  summarise(
    min_sigma = round(min(sigma_log), 3),
    max_sigma = round(max(sigma_log), 3),
    .groups   = "drop"
  ) %>%
  arrange(desc(max_sigma)) %>%
  { for (i in seq_len(nrow(.))) {
    r <- .[i, ]
    cat(sprintf("    %-14s  sigma h1=%.3f  h4=%.3f\n",
                r$model, r$min_sigma, r$max_sigma))
  }; invisible(.) }

rsv_joined <- rsv_fc %>%
  left_join(rsv_wts, by = c("model" = ".model")) %>%
  filter(!is.na(weight)) %>%
  group_by(date) %>%
  mutate(weight = weight / sum(weight)) %>%
  ungroup()

rsv_dates  <- rsv_joined %>%
  distinct(date) %>% arrange(date) %>% mutate(h = row_number())
rsv_joined <- rsv_joined %>% left_join(rsv_dates, by = "date")

rsv_mu_mix <- rsv_joined %>%
  group_by(date, h) %>%
  summarise(mu_mix = sum(weight * mu_log), .groups = "drop")

rsv_var <- rsv_joined %>%
  left_join(rsv_mu_mix, by = c("date", "h")) %>%
  group_by(h) %>%
  summarise(
    within_var  = mean(sum(weight * sigma_log^2)),
    between_var = mean(sum(weight * (mu_log - mu_mix)^2)),
    .groups     = "drop"
  ) %>%
  mutate(
    total_var   = within_var + between_var,
    pct_between = round(between_var / total_var * 100, 1)
  )

cat("\n  Variance decomposition:\n")
for (i in seq_len(nrow(rsv_var))) {
  r    <- rsv_var[i, ]
  flag <- if (r$pct_between > 30) " *** still high" else " ✓"
  cat(sprintf("  Week +%d:  between=%.1f%%  within=%.1f%%%s\n",
              r$h, r$pct_between, 100 - r$pct_between, flag))
}

rsv_ens      <- results_all$RSV$fc_ensemble
ci_ratio_v4  <- mean(
  (rsv_ens$upper_ci - rsv_ens$lower_ci) / pmax(rsv_ens$mean_case, 0.5),
  na.rm = TRUE
)
improvement  <- round((1.659 - ci_ratio_v4) / 1.659 * 100, 1)
cat(sprintf(
  "\n  RSV CI width ratio: v2=1.659 | v3=1.700 | v4=1.680 | v5=1.623 | v6=%.3f  (%s%.1f%% vs v2)\n",
  ci_ratio_v4,
  if (ci_ratio_v4 < 1.659) "-" else "+",
  abs(improvement)
))
cat(sprintf("  RSV models in ensemble: %s\n",
            paste(unique(rsv_fc$model), collapse = ", ")))
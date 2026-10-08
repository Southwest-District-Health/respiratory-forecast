
# =============================================================================
# NWSS wastewater pull via data.cdc.gov (Socrata SODA API)
# Replaces the three read.csv() calls with API queries + weekly MMWR aggregation
#
# Datasets (all updated weekly on Fridays):
#   Influenza A  ymmh-divb
#   RSV          45cq-cw4i
#   SARS-CoV-2   j9g8-acpt
#
# Run order:
#   1. nwss_schema()  - see what columns the API actually serves
#   2. fetch_nwss()   - pull filtered rows, auto-mapping columns
#   3. to_weekly()    - your existing aggregation logic
#   4. compare_to_legacy() - reconcile against current CSV-based output
# =============================================================================

library(dplyr)
library(lubridate)
library(httr2)
library(readr)
library(purrr)

# ── Config ───────────────────────────────────────────────────────────────────
NWSS_IDS   <- c(flu = "ymmh-divb", rsv = "45cq-cw4i", covid = "j9g8-acpt")
HD3_SITES  <- c(370, 372, 399, 1739)
START_DATE <- "2020-01-16"
PAGE_SIZE  <- 50000L
APP_TOKEN  <- Sys.getenv("SOCRATA_APP_TOKEN")   # optional, set in ~/.Renviron

# Candidate column names, in priority order. The first one present wins.
# Add to these lists if the schema probe turns up something different.
COL_CANDIDATES <- list(
  site = c("site", "wwtp_id", "sample_site_id", "key_plot_id"),
  date = c("sample_collect_date", "collection_date", "date_start",
           "sample_collect_datetime"),
  conc = c("pcr_target_mic_lin", "pcr_target_avg_conc", "pcr_conc_lin",
           "normalized_conc_lin", "pcr_target_flowpop_lin"),
  pop  = c("population_served", "pop_served"),
  det  = c("pcr_target_detect", "detect_flag")   # optional
)

# ── Low-level request helper ─────────────────────────────────────────────────
soda <- function(id, ...) {
  req <- request(sprintf("https://data.cdc.gov/resource/%s.csv", id)) |>
    req_url_query(...) |>
    req_error(body = function(resp) resp_body_string(resp))   # surface SoQL errors
  
  if (nzchar(APP_TOKEN)) req <- req_headers(req, `X-App-Token` = APP_TOKEN)
  
  req |>
    req_retry(max_tries = 3) |>
    req_perform() |>
    resp_body_string() |>
    read_csv(show_col_types = FALSE, progress = FALSE)
}

# ── 1. Schema probe ──────────────────────────────────────────────────────────
nwss_schema <- function(id, n = 5) {
  probe <- soda(id, `$limit` = n)
  message(sprintf("[%s] %d columns:", id, ncol(probe)))
  message(paste(" ", paste(names(probe), collapse = ", ")))
  probe
}

pick_col <- function(nms, key) {
  hit <- intersect(COL_CANDIDATES[[key]], nms)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

# ── 2. Fetch with server-side filtering + pagination ─────────────────────────
fetch_nwss <- function(id, sites = HD3_SITES, since = START_DATE,
                       verbose = TRUE) {
  
  probe <- soda(id, `$limit` = 1)
  nms   <- names(probe)
  
  map <- vapply(names(COL_CANDIDATES), \(k) pick_col(nms, k), character(1))
  
  missing <- names(map)[is.na(map) & names(map) != "det"]
  if (length(missing) > 0) {
    stop(sprintf(
      "Could not map required column(s): %s\nAvailable columns:\n  %s\nAdd the right name to COL_CANDIDATES and re-run.",
      paste(missing, collapse = ", "), paste(nms, collapse = ", ")
    ), call. = FALSE)
  }
  
  # Site IDs may be stored as number or text; try unquoted, fall back to quoted.
  site_clause_num <- sprintf("%s in (%s)", map[["site"]],
                             paste(sites, collapse = ","))
  site_clause_txt <- sprintf("%s in (%s)", map[["site"]],
                             paste0("'", sites, "'", collapse = ","))
  date_clause <- sprintf("%s >= '%sT00:00:00.000'", map[["date"]], since)
  
  keep <- unname(map[!is.na(map)])
  sel  <- paste(keep, collapse = ", ")
  
  pull_page <- function(where, offset) {
    soda(id,
         `$select` = sel,
         `$where`  = where,
         `$order`  = sprintf("%s, %s", map[["date"]], map[["site"]]),
         `$limit`  = PAGE_SIZE,
         `$offset` = offset)
  }
  
  where <- paste(site_clause_num, "AND", date_clause)
  first <- tryCatch(pull_page(where, 0), error = function(e) {
    if (verbose) message("Numeric site filter rejected; retrying as text.")
    where <<- paste(site_clause_txt, "AND", date_clause)
    pull_page(where, 0)
  })
  
  pages  <- list(first)
  offset <- PAGE_SIZE
  while (nrow(pages[[length(pages)]]) == PAGE_SIZE) {
    if (verbose) message(sprintf("  ...fetching rows from offset %d", offset))
    pages[[length(pages) + 1]] <- pull_page(where, offset)
    offset <- offset + PAGE_SIZE
  }
  
  out <- bind_rows(pages)
  
  # Rename to the canonical names the rest of the pipeline expects
  rename_vec <- c(site = map[["site"]], sample_collect_date = map[["date"]],
                  pcr_target_mic_lin = map[["conc"]],
                  population_served = map[["pop"]])
  if (!is.na(map[["det"]])) {
    rename_vec <- c(rename_vec, pcr_target_detect = map[["det"]])
  }
  out <- rename(out, all_of(rename_vec))
  
  if (!"pcr_target_detect" %in% names(out)) out$pcr_target_detect <- NA_character_
  
  if (verbose) {
    message(sprintf(
      "[%s] %d rows | sites: %s | %s to %s",
      id, nrow(out), paste(sort(unique(out$site)), collapse = ", "),
      min(as.Date(out$sample_collect_date)), max(as.Date(out$sample_collect_date))
    ))
  }
  out
}

# ── 3. Weekly aggregation (unchanged logic, applied once) ────────────────────
to_weekly <- function(df) {
  df |>
    mutate(
      sample_collect_date = as.Date(sample_collect_date),
      population_served   = as.numeric(gsub(",", "", as.character(population_served))),
      pcr_target_mic_lin  = as.numeric(pcr_target_mic_lin)
    ) |>
    filter(!is.na(pcr_target_mic_lin)) |>
    mutate(
      mmwr_year       = epiyear(sample_collect_date),
      mmwr_week       = epiweek(sample_collect_date),
      mmwr_week_start = floor_date(sample_collect_date, "week", week_start = 7),
      mmwr_week_end   = mmwr_week_start + 6L
    ) |>
    group_by(mmwr_year, mmwr_week, mmwr_week_start, mmwr_week_end) |>
    summarise(
      ww_signal     = weighted.mean(pcr_target_mic_lin, population_served,
                                    na.rm = TRUE),
      ww_signal_log = weighted.mean(log(pmax(pcr_target_mic_lin, 1e-10) * 1e6),
                                    population_served, na.rm = TRUE),
      n_samples     = n(),
      sites_present = paste(sort(unique(site)), collapse = ", "),
      frac_detected = mean(pcr_target_detect == "yes", na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(mmwr_year, mmwr_week)
}

# ── 4. Run ───────────────────────────────────────────────────────────────────
# Step 1 first, on its own, so you can eyeball the column names:
#   walk(NWSS_IDS, nwss_schema)

ww_raw    <- map(NWSS_IDS, fetch_nwss)
ww_weekly <- map(ww_raw, to_weekly)

ww_flu_weekly   <- ww_weekly$flu
ww_rsv_weekly   <- ww_weekly$rsv
ww_covid_weekly <- ww_weekly$covid

# ── Sanity checks ────────────────────────────────────────────────────────────
iwalk(ww_weekly, function(df, nm) {
  rng <- range(df$ww_signal_log, na.rm = TRUE)
  message(sprintf(
    "[%s] weeks=%d | %s to %s | samples=%d | ww_signal_log range=[%.3f, %.3f] mean=%.3f %s",
    nm, nrow(df), min(df$mmwr_week_start), max(df$mmwr_week_start),
    sum(df$n_samples), rng[1], rng[2], mean(df$ww_signal_log, na.rm = TRUE),
    if (rng[2] > 1.0) "OK" else "WARNING: near zero - check units"
  ))
})

# Weeks where the site roster differs from the expected three-site combination
iwalk(ww_weekly, function(df, nm) {
  odd <- filter(df, sites_present != "370, 372, 399")
  if (nrow(odd) > 0) {
    message(sprintf("[%s] %d weeks with a different site roster:", nm, nrow(odd)))
    print(select(odd, mmwr_year, mmwr_week, mmwr_week_start, sites_present, n_samples))
  }
})

# ── 5. Reconcile against the CSV-based output before cutting over ────────────
# Run your existing CSV script first so ww_flu_weekly_old etc. are in memory,
# or point this at saved copies.
compare_to_legacy <- function(new, old, label = "") {
  j <- inner_join(
    select(new, mmwr_year, mmwr_week, new_signal = ww_signal,
           new_log = ww_signal_log, new_n = n_samples),
    select(old, mmwr_year, mmwr_week, old_signal = ww_signal,
           old_log = ww_signal_log, old_n = n_samples),
    by = c("mmwr_year", "mmwr_week")
  ) |>
    mutate(
      pct_diff_signal = 100 * (new_signal - old_signal) / old_signal,
      diff_log        = new_log - old_log,
      n_diff          = new_n - old_n
    )
  
  message(sprintf(
    "[%s] %d overlapping weeks (new=%d, old=%d) | median |pct diff| in ww_signal: %.2f%% | max |diff| in ww_signal_log: %.3f | weeks with differing sample counts: %d",
    label, nrow(j), nrow(new), nrow(old),
    median(abs(j$pct_diff_signal), na.rm = TRUE),
    max(abs(j$diff_log), na.rm = TRUE),
    sum(j$n_diff != 0, na.rm = TRUE)
  ))
  j
}

# cmp_flu <- compare_to_legacy(ww_flu_weekly, ww_flu_weekly_old, "flu")
# build_forecast.R ----------------------------------------------------------
# Turns the ensemble forecast results into index.html: a standalone page
# with the same content as the forecast Shiny app, and no server behind it.
#
# Author: Lekshmi Rita-Venugopal, Southwest District Health
#
# Weekly order:
#   1. wasterwater_API.R
#   2. merge_ww_essence.R
#   3. Resp_forecast_ensemble_final.R   (writes the forecast .rds)
#   4. source("build_forecast.R", echo = TRUE)
#
# Needs in the project folder:
#   forecast_template.html
# ---------------------------------------------------------------------------

PROJECT_DIR <- "~/R/ESSENCE"
setwd(PROJECT_DIR)

suppressPackageStartupMessages({
  library(tidyverse)
  library(tsibble)   # the forecast .rds holds tsibbles
  library(jsonlite)
})

# --- Settings --------------------------------------------------------------
FORECAST_RDS <- "~/Data/Essence/precomputed_forecasts_new.rds"
TEMPLATE     <- "forecast_template.html"
OUT_HTML     <- "index.html"
OUT_FC_LOG   <- "published_forecast_log.csv"   # every forecast ever published

DEFAULT_DISEASE <- "flu"

# The "current activity" sentence is written from the activity levels page's
# log, so the two pages always agree. If the log is missing, or its latest
# week doesn't match the forecast data, the hand-typed values below are used.
LEVELS_LOG <- "~/Respiratory activity map/published_levels_log.csv"

# Fallback only. One of "minimal", "low", "moderate", "high", "very high",
# or NA to hide the sentence for that disease.
CURRENT_ACTIVITY <- list(flu = "minimal", rsv = "minimal", covid = "minimal")

DISEASES <- list(
  # name: used in "visits for the flu"; act_name: used in "Current flu activity"
  flu   = list(short = "Flu",      key = "Influenza", color = "#1A6FA8",
               name = "the flu",   act_name = "flu",
               formal = "influenza-like illness (ILI)"),
  rsv   = list(short = "RSV",      key = "RSV",       color = "#0F6E56",
               name = "RSV",       act_name = "RSV",
               formal = "RSV-like illness"),
  covid = list(short = "COVID-19", key = "SARS-CoV2", color = "#534AB7",
               name = "COVID-19",  act_name = "COVID-19",
               formal = "COVID-like illness (CLI)")
)

for (f in c(TEMPLATE, FORECAST_RDS)) {
  if (!file.exists(f)) stop(f, " not found.", call. = FALSE)
}

# --- 1. Read the forecast results ------------------------------------------
rds_age <- as.numeric(difftime(Sys.time(), file.mtime(FORECAST_RDS), units = "days"))
if (rds_age > 6) {
  warning("The forecast file is ", round(rds_age), " days old. ",
          "Rerun Resp_forecast_ensemble_final.R before publishing.",
          call. = FALSE)
}

precomputed <- readRDS(FORECAST_RDS)

# --- 2. Summaries ---------
change_word <- function(pct, flat_band = 5) {
  if (is.na(pct)) "no data"
  else if (pct >  flat_band) "up"
  else if (pct < -flat_band) "down"
  else "about the same"
}

summarise_forecast <- function(res) {
  ts <- res$ts_data %>% as_tibble() %>%
    mutate(date = as.Date(date), week_end = date + 6)
  fc <- res$fc %>% as_tibble() %>%
    mutate(across(c(mean_case, lower_ci, upper_ci),
                  ~ ifelse(is.finite(.), ., NA_real_)),
           date = as.Date(date), week_end = date + 6)

  last_obs  <- as.numeric(tail(ts$case, 1))[1]
  last_end  <- max(ts$week_end, na.rm = TRUE)
  recent    <- mean(tail(ts$case, 4),          na.rm = TRUE)
  prior4    <- mean(head(tail(ts$case, 8), 4), na.rm = TRUE)
  hist_mean <- mean(ts$case, na.rm = TRUE)
  hist_sd   <- sd(ts$case,   na.rm = TRUE)
  fc_vals   <- fc$mean_case[is.finite(fc$mean_case)]
  fc_mean   <- if (length(fc_vals)) mean(fc_vals) else NA_real_

  yr_ago_date <- tail(ts$date, 1) - 364
  yr_ago      <- ts %>% filter(date == yr_ago_date) %>% pull(case)
  yr_ago      <- if (length(yr_ago)) as.numeric(yr_ago[1]) else NA_real_

  chg_vs_yr_ago <- if (!is.na(yr_ago) && yr_ago > 0)
    (last_obs - yr_ago) / yr_ago * 100 else NA_real_

  # Stability factors damp marginal drift calls at low counts
  f1 <- if (!is.na(recent) && !is.na(prior4) && prior4 > 0) {
    p <- (recent - prior4) / prior4; a <- abs(recent - prior4)
    if (p > 0.30 && a > 3) 2L else if (p > 0.10 && a > 3) 1L else 0L
  } else 0L
  f2 <- if (!is.na(recent) && !is.na(hist_sd) && hist_sd > 0) {
    z <- (recent - hist_mean) / hist_sd
    if (z > 1.5) 2L else if (z > 0.5) 1L else 0L
  } else 0L
  f3 <- if (length(fc_vals) >= 4 && fc_vals[1] > 0) {
    a <- (fc_vals[4] - fc_vals[1]) / fc_vals[1]; ab <- abs(fc_vals[4] - fc_vals[1])
    if (a > 0.20 && ab > 3) 2L else if (a > 0.05 && ab > 3) 1L else 0L
  } else 0L

  trend <- tryCatch({
    if (length(fc_vals) > 0 && isTRUE(is.finite(last_obs)) && isTRUE(last_obs > 0)) {
      pct <- (fc_mean - last_obs) / last_obs
      if      (isTRUE(pct >  0.15)) "Increasing"
      else if (isTRUE(pct < -0.15)) "Decreasing"
      else if (isTRUE(pct >  0.05)) "Stable-up"
      else if (isTRUE(pct < -0.05)) "Stable-down"
      else                          "Stable"
    } else "Stable"
  }, error = function(e) "Stable")

  if (f1 == 0 && f2 == 0 && f3 == 0 &&
      trend %in% c("Stable-up", "Stable-down")) trend <- "Stable"

  fc_pct <- if (isTRUE(is.finite(last_obs)) && last_obs > 0)
    (fc_mean - last_obs) / last_obs * 100 else NA_real_

  list(
    obs = list(d = I(format(ts$week_end)),
               y = I(round(as.numeric(ts$case), 1))),
    fc  = list(d    = I(format(fc$week_end)),
               mean = I(round(fc$mean_case, 1)),
               lo   = I(round(pmax(fc$lower_ci, 0), 1)),
               hi   = I(round(fc$upper_ci, 1))),
    lastObs  = last_obs,
    lastEnd  = format(last_end),
    yrAgo    = yr_ago,
    chgYrAgo = round(chg_vs_yr_ago, 1),
    yrWord   = change_word(chg_vs_yr_ago),
    fcMean   = round(fc_mean, 1),
    fcWord   = change_word(fc_pct),
    trend    = trend
  )
}

forecasts <- map(DISEASES, function(d) {
  res <- precomputed[[d$key]]
  if (is.null(res)) stop("'", d$key, "' not found in ", FORECAST_RDS, call. = FALSE)
  summarise_forecast(res)
})

data_week <- max(as.Date(map_chr(forecasts, "lastEnd")))

# --- 3. Current activity sentence ------------------------------------------
join_names <- function(x) {
  if (length(x) <= 1) x
  else if (length(x) == 2) paste(x, collapse = " and ")
  else paste0(paste(head(x, -1), collapse = ", "), ", and ", tail(x, 1))
}

activity_from_log <- function() {
  if (!file.exists(LEVELS_LOG)) {
    message("No activity levels log found; using CURRENT_ACTIVITY.")
    return(NULL)
  }
  lg <- read_csv(LEVELS_LOG, col_types = cols(.default = "c")) %>%
    mutate(Week = as.Date(WeekEnding, format = "%m/%d/%y"),
           Level = as.integer(Level))
  latest <- max(lg$Week, na.rm = TRUE)
  if (latest != data_week) {
    warning("Activity levels log is for week ending ", latest,
            " but the forecast data runs through ", data_week,
            ". Using CURRENT_ACTIVITY instead. Rebuild the activity levels page ",
            "first if you want the two to match.", call. = FALSE)
    return(NULL)
  }
  wk <- filter(lg, Week == latest)
  map(names(DISEASES), function(id) {
    rows <- filter(wk, Disease == id)
    if (!nrow(rows)) return(NA_character_)
    lo <- min(rows$Level); hi <- max(rows$Level)
    nm <- DISEASES[[id]]$act_name
    lv <- function(i) tolower(rows$LevelName[match(i, rows$Level)])
    if (lo == hi) {
      paste0("Current ", nm, " activity is <b>", lv(hi), "</b> in all six counties.")
    } else {
      paste0("Current ", nm, " activity ranges from <b>", lv(lo), "</b> to <b>",
             lv(hi), "</b> across the six counties, highest in ",
             join_names(rows$County[rows$Level == hi]), ".")
    }
  }) %>% set_names(names(DISEASES))
}

activity <- activity_from_log()
if (is.null(activity)) {
  activity <- imap(CURRENT_ACTIVITY, function(a, id) {
    if (is.na(a) || !nzchar(a)) NA_character_
    else paste0("Current ", DISEASES[[id]]$act_name, " activity in the region is <b>", a, "</b>.")
  })
}

# --- 4. Keep every published forecast, so they can be scored later ---------
week_label <- format(data_week, "%m/%d/%y")
new_rows <- imap_dfr(forecasts, function(f, id) {
  tibble(PublishedFor = week_label, Disease = id,
         TargetWeek = f$fc$d, Expected = f$fc$mean,
         Lower = f$fc$lo, Upper = f$fc$hi, Trend = f$trend)
}) %>% mutate(across(everything(), as.character))

if (file.exists(OUT_FC_LOG)) {
  prior <- read_csv(OUT_FC_LOG, col_types = cols(.default = "c")) %>%
    filter(PublishedFor != week_label)   # rerunning a week replaces it
  new_rows <- bind_rows(prior, new_rows)
}
write_csv(new_rows, OUT_FC_LOG)

# --- 5. Build the page ------------------------------------------------------
payload <- list(
  defaultDisease = DEFAULT_DISEASE,
  diseases = imap(DISEASES, ~ list(id = .y, short = .x$short, color = .x$color,
                                   name = .x$name, formal = .x$formal)) %>% unname(),
  forecast = forecasts,
  activity = activity
)

json <- toJSON(payload, auto_unbox = TRUE, na = "null", null = "null", digits = NA)

tpl   <- paste(readLines(TEMPLATE, warn = FALSE), collapse = "\n")
parts <- strsplit(tpl, "/*__DATA__*/", fixed = TRUE)[[1]]
if (length(parts) != 2) {
  stop("Could not find the /*__DATA__*/ placeholder in ", TEMPLATE, call. = FALSE)
}
writeLines(paste0(parts[1], json, parts[2]), OUT_HTML, useBytes = TRUE)

# --- 6. Report --------------------------------------------------------------
message("\nWrote ", normalizePath(OUT_HTML), " (",
        round(file.size(OUT_HTML) / 1024), " KB)")
message("Data through week ending ", format(data_week, "%m/%d/%Y"),
        " | forecast file ", round(rds_age, 1), " days old\n")

imap_dfr(forecasts, ~ tibble(Disease = .y, DataThrough = .x$lastEnd,
                             LastWeek = .x$lastObs, Next4wkAvg = .x$fcMean,
                             Trend = .x$trend)) %>%
  print()

message("\nActivity sentences:")
walk2(names(activity), activity,
      ~ message("  ", .x, ": ", gsub("</?b>", "", ifelse(is.na(.y), "(hidden)", .y))))

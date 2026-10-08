library(dplyr)
library(lubridate)
library(writexl)
library(readxl)
library(ggplot2)

# ============================================================
# IMPORTANT: Run wasterwater.R FIRST to regenerate
# ww_flu_weekly, ww_rsv_weekly, ww_covid_weekly 
# ============================================================

# Helper: scale WW to overlay on disease count axis.
# Finds a multiplier so WW signal peak ≈ 80% of disease peak.
# This gives a clean visual without hardcoding 1e6.
ww_scale <- function(ww, disease_counts) {
  ww_peak      <- max(ww,             na.rm = TRUE)
  disease_peak <- max(disease_counts, na.rm = TRUE)
  if (ww_peak == 0) return(1)
  (disease_peak * 0.80) / ww_peak
}

# ============================================================
# ILI
# ============================================================

ILI_ed <- read_excel("~/Data/Essence/ILI_ed.xlsx")

ili_clean <- ILI_ed %>%
  mutate(
    Date      = as.Date(Date, format = "%d%b%y"),
    mmwr_year = epiyear(Date),
    mmwr_week = epiweek(Date)
  )

ili_ww <- ili_clean %>%
  left_join(
    ww_flu_weekly %>%
      select(mmwr_year, mmwr_week, ww_signal_log, ww_signal,
             n_samples, sites_present, frac_detected),
    by = c("mmwr_year", "mmwr_week")
  ) %>%
  arrange(Date)

cat(sprintf(
  "ILI rows: %d | Matched WW: %d | Unmatched: %d\n",
  nrow(ili_ww),
  sum(!is.na(ili_ww$ww_signal_log)),
  sum(is.na(ili_ww$ww_signal_log))
))

ILI <- ili_ww %>%
  filter(!is.na(ww_signal_log)) %>%
  select(Date, Data, VE, ww_signal_log)

write_xlsx(ILI, "~/Data/Essence/ILI.xlsx")
message(sprintf("ILI saved: %d rows | ww_signal_log range [%.2f, %.2f]",
                nrow(ILI),
                min(ILI$ww_signal_log, na.rm = TRUE),
                max(ILI$ww_signal_log, na.rm = TRUE)))

# Plot — scale WW to fit disease axis naturally
sc_ili <- ww_scale(ILI$ww_signal_log, ILI$Data)

ggplot(ILI, aes(x = Date)) +
  geom_line(aes(y = Data,                         colour = "ILI cases"),
            linewidth = 0.8) +
  geom_line(aes(y = ww_signal_log * sc_ili,       colour = "WW signal (log)"),
            linewidth = 0.8, linetype = "dashed") +
  scale_y_continuous(
    name     = "ILI cases",
    sec.axis = sec_axis(~ . / sc_ili,
                        name   = "WW signal log(conc × 1e6)",
                        breaks = scales::pretty_breaks(5))
  ) +
  scale_colour_manual(values = c("ILI cases" = "#0969DA",
                                 "WW signal (log)" = "#CF222E")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(title  = "ILI vs Wastewater Signal (Influenza)",
       colour = NULL,
       x      = NULL)

# ============================================================
# RSV
# ============================================================

RSV_ed <- read_excel("~/Data/Essence/RSV_ed.xlsx")

rsv_clean <- RSV_ed %>%
  mutate(
    Date      = as.Date(Date, format = "%d%b%y"),
    mmwr_year = epiyear(Date),
    mmwr_week = epiweek(Date)
  )

rsv_ww <- rsv_clean %>%
  left_join(
    ww_rsv_weekly %>%
      select(mmwr_year, mmwr_week, ww_signal_log, ww_signal,
             n_samples, sites_present, frac_detected),
    by = c("mmwr_year", "mmwr_week")
  ) %>%
  arrange(Date)

cat(sprintf(
  "RSV rows: %d | Matched WW: %d | Unmatched: %d\n",
  nrow(rsv_ww),
  sum(!is.na(rsv_ww$ww_signal_log)),
  sum(is.na(rsv_ww$ww_signal_log))
))

RSV <- rsv_ww %>%
  filter(!is.na(ww_signal_log)) %>%
  select(Date, Data, ww_signal_log)

write_xlsx(RSV, "~/Data/Essence/RSV.xlsx")
message(sprintf("RSV saved: %d rows | ww_signal_log range [%.2f, %.2f]",
                nrow(RSV),
                min(RSV$ww_signal_log, na.rm = TRUE),
                max(RSV$ww_signal_log, na.rm = TRUE)))

sc_rsv <- ww_scale(RSV$ww_signal_log, RSV$Data)

ggplot(RSV, aes(x = Date)) +
  geom_line(aes(y = Data,                         colour = "RSV cases"),
            linewidth = 0.8) +
  geom_line(aes(y = ww_signal_log * sc_rsv,       colour = "WW signal (log)"),
            linewidth = 0.8, linetype = "dashed") +
  scale_y_continuous(
    name     = "RSV cases",
    sec.axis = sec_axis(~ . / sc_rsv,
                        name   = "WW signal log(conc × 1e6)",
                        breaks = scales::pretty_breaks(5))
  ) +
  scale_colour_manual(values = c("RSV cases"      = "#0969DA",
                                 "WW signal (log)" = "#CF222E")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(title  = "RSV vs Wastewater Signal",
       colour = NULL,
       x      = NULL)

# ============================================================
# CLI (SARS-CoV-2)
# ============================================================

CLI_ed <- read_excel("~/Data/Essence/CLI_ed.xlsx")

CLI_clean <- CLI_ed %>%
  mutate(
    Date      = as.Date(Date, format = "%d%b%y"),
    mmwr_year = epiyear(Date),
    mmwr_week = epiweek(Date)
  )

CLI_ww <- CLI_clean %>%
  left_join(
    ww_covid_weekly %>%
      select(mmwr_year, mmwr_week, ww_signal_log, ww_signal,
             n_samples, sites_present, frac_detected),
    by = c("mmwr_year", "mmwr_week")
  ) %>%
  arrange(Date)

cat(sprintf(
  "CLI rows: %d | Matched WW: %d | Unmatched: %d\n",
  nrow(CLI_ww),
  sum(!is.na(CLI_ww$ww_signal_log)),
  sum(is.na(CLI_ww$ww_signal_log))
))

CLI <- CLI_ww %>%
  filter(!is.na(ww_signal_log)) %>%
  select(Date, Data, ww_signal_log)

write_xlsx(CLI, "~/Data/Essence/CLI.xlsx")
message(sprintf("CLI saved: %d rows | ww_signal_log range [%.2f, %.2f]",
                nrow(CLI),
                min(CLI$ww_signal_log, na.rm = TRUE),
                max(CLI$ww_signal_log, na.rm = TRUE)))

sc_cli <- ww_scale(CLI$ww_signal_log, CLI$Data)

ggplot(CLI, aes(x = Date)) +
  geom_line(aes(y = Data,                         colour = "CLI cases"),
            linewidth = 0.8) +
  geom_line(aes(y = ww_signal_log * sc_cli,       colour = "WW signal (log)"),
            linewidth = 0.8, linetype = "dashed") +
  scale_y_continuous(
    name     = "CLI cases",
    sec.axis = sec_axis(~ . / sc_cli,
                        name   = "WW signal log(conc × 1e6)",
                        breaks = scales::pretty_breaks(5))
  ) +
  scale_colour_manual(values = c("CLI cases"       = "#0969DA",
                                 "WW signal (log)" = "#CF222E")) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom") +
  labs(title  = "CLI vs Wastewater Signal (SARS-CoV-2)",
       colour = NULL,
       x      = NULL)

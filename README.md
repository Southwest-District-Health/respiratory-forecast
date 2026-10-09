# Respiratory Illness Forecast

A four-week forecast of emergency room visits for flu, RSV, and COVID-19
across the six counties of Health District 3: Adams, Canyon, Gem, Owyhee,
Payette, and Washington.

**Live page:** https://southwest-district-health.github.io/respiratory-forecast/

Updated every Monday during respiratory virus season, September through May.

---

## What the page shows

For each illness, the page shows last week's emergency room visits, the same
week a year ago, and what we expect over the next four weeks. The chart puts
what actually happened next to the forecast, with a shaded band showing the
range the forecast could reasonably fall in. Below the chart, a short summary
explains in plain language whether visits are expected to rise, fall, or stay
about the same.

## How the forecast works

Several statistical and machine learning models each make their own forecast
from past emergency room visits, seasonal patterns, and wastewater data. The
forecasts are then combined, with each model weighted by how well it has
predicted recent weeks. The "could be" range comes from that combined model.

The most recent week can be revised as more reports come in, and forecasts are
estimates, not certainties.

## What's in this repo

| File | What it is |
|---|---|
| `index.html` | The published page. Everything it needs is inside this one file. |
| `wasterwater_API.R` | Pulls wastewater data for the district's sites from the CDC and averages it by week. |
| `merge_ww_essence.R` | Joins the weekly wastewater signal to the ESSENCE emergency department data. |
| `Resp_forecast_ensemble_final.R` | Runs the forecasting models and combines them into one forecast. |
| `build_forecast.R` | Reads the forecast results and writes `index.html`. |
| `forecast_template.html` | Page layout, styling, and wording. |
| `REFRESH.md` | Step-by-step weekly update instructions. |

The emergency department and wastewater data themselves are kept on district
machines and are not part of this repo.

## Input data

**Emergency department data.** Three weekly exports from ESSENCE, one per
illness:

| File | Illness |
|---|---|
| `ILI_ed.xlsx` | Flu (influenza-like illness) |
| `RSV_ed.xlsx` | RSV |
| `CLI_ed.xlsx` | COVID-19 (COVID-like illness) |

Each file needs these columns, spelled exactly as shown:

| Column | What it holds | Needed in |
|---|---|---|
| `Date` | The week, written like `05Oct26` (day, three-letter month, two-digit year) | All three |
| `Data` | Weekly emergency department visits for that illness | All three |
| `VE` | An extra value used by the flu models | `ILI_ed.xlsx` only |

**Wastewater data.** `wasterwater_API.R` downloads this directly from the
CDC's National Wastewater Surveillance System for the district's sites, so
there's no file to prepare. Setting `SOCRATA_APP_TOKEN` in `~/.Renviron` is
optional but helps avoid rate limits.

**Merged files.** `merge_ww_essence.R` writes `ILI.xlsx`, `RSV.xlsx`, and
`CLI.xlsx`, which the forecasting script reads. They keep `Date`, `Data`, and
`VE` (flu only), and add `ww_signal_log`, the weekly wastewater signal on a
log scale. Weeks with no wastewater data are dropped.

A few things to keep in mind:

- Run `wasterwater_API.R` and `merge_ww_essence.R` in the same R session.
  The merge script uses the weekly wastewater tables the first script leaves
  in memory, not a saved file.
- If `ww_signal_log` is missing from a merged file, the forecasting script
  stops and lists the columns it did find. If `VE` is missing from the flu
  file, it stops too.

## Related pages

- [Respiratory activity levels by county](https://southwest-district-health.github.io/respiratory-activity-levels/)

## Contact

Maintained by Lekshmi Rita-Venugopal, Epidemiologist Program Manager,
Southwest District Health.

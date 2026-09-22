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
| `build_forecast.R` | Reads the forecast results and writes `index.html`. |
| `forecast_template.html` | Page layout, styling, and wording. |
| `REFRESH.md` | Step-by-step weekly update instructions. |

The forecasting models themselves, and the emergency room and wastewater data
they use, are kept on district machines and are not part of this repo.

## Related pages

- [Respiratory activity levels by county](https://southwest-district-health.github.io/respiratory-activity-levels/)

## Contact

Maintained by Lekshmi Rita-Venugopal, MD, MPH, Epidemiologist Program Manager,
Southwest District Health.

This page is for general awareness and is not a medical diagnosis.

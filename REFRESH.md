# Weekly refresh

How to update the forecast page each Monday, September through May.

Project folder: `Documents\Respiratory forecast` (in R, `~/Respiratory forecast`).

---

## 1. Update the activity levels page first

The forecast page takes its "current activity" sentence from the activity
levels page's log, so rebuild that page first for the same week. If you skip
this, the forecast page still builds but falls back to the hand-typed values
in `CURRENT_ACTIVITY`, and the console tells you so.

## 2. Run the forecast pipeline

In this order, as usual:

1. `wasterwater_API.R`
2. `merge_ww_essence.R`
3. `Resp_forecast_ensemble_final.R`

The last one saves `~/Data/Essence/precomputed_forecasts_new.rds`, which is
all the page needs.

## 3. Build the page

```r
source("build_forecast.R", echo = TRUE)
```

The console prints, for each illness, the last week of data, last week's
visits, the four-week average, and the trend, followed by the activity
sentence that will appear on the page. Read them before publishing.

The script warns you if:

- the forecast file is more than six days old, meaning the pipeline wasn't
  rerun
- the activity levels log is from a different week than the forecast data

## 4. Check the page

Open `index.html` in a browser and click through Flu, RSV, and COVID-19.
Check the date in the top right corner and that the chart, table, and summary
all agree with the console output.

## 5. Publish

In this repo, click **Add file**, then **Upload files**, and drag in the new
`index.html`. Commit with a message like `Week ending 09/19/26`.

GitHub updates the live page within a couple of minutes. Check it in a
private browser window, which shows exactly what visitors see instead of a
copy your browser saved earlier.

Only upload `index.html`. The data files, the `.rds` file, and the log stay on
your machine.

---

## Good to know

**Forecast log.** Every run adds that week's published forecast to
`published_forecast_log.csv` in the project folder: the expected visits, the
range, and the trend for each of the next four weeks. Rerunning a week
replaces that week's rows. Over a season, this is what lets you check how
accurate the forecasts were.

**Model settings** live in `Resp_forecast_ensemble_final.R`, not here. This
script only turns the results into a page.

**Off-season.** Between June and August nothing needs to happen. The page
keeps showing the last forecast, with its date in the top corner.

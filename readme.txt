Kp forecast
===========

The forecast uses the previous 24 hours of Kp and solar wind to predict the
3-hour Kp bins from 3 hours through 24 hours ahead. The issued model is one
24-hour VARX. In-situ measurements and the AWSoM Earth model share those
coefficients; only the solar-wind block changes. The two pages are
self-contained English HTML files and open in a browser with no other files.

On the scored issues so far (about 2026-02-25 through 2026-10-15, n = 1856),
24-hour endpoint RMSE is:

  Persistence            1.579
  In-situ measurements   1.218
  AWSoM Earth model      1.264
  Boosted                1.212

Boosted changes the forecast only when the last observed Kp is already about
5 to 6.5, so the curves otherwise match in-situ. The charts mark only those
changed points.


Model
-----

A forecast is issued at time tau, the end of the 3-hour Kp bin that just
finished (00, 03, 06, 09, 12, 15, 18, 21 UTC).

Each lead is its own regression and predicts that Kp bin directly. Predicted
values are never fed back in. The leads are 3, 6, 9, 12, 15, 18, 21, and 24
hours, the next eight Kp bins. The fit is unpenalized least squares (thin
SVD). Predictors are standardized on the training rows. Coefficients are
cached in Rcode/Main/output/Main_VARX_Realtime/VARX_Kp_coef.RData and refit
only when Kp, OMNI, or the training cutoff changes.

Predictors contain only data that had already ended by tau:

  * The Date on a Kp row is the start of that bin. The bin is known only at
    its end.
  * The Date on a solar-wind hour is the start of that hour. The hour is
    known only at its end.
  * Solar wind inside the forecast window is not a predictor. In-situ
    measurements of that window do not exist at issue time, and the AWSoM
    file currently ends near the present hour. It does not contain the next
    24 hours.

The state is the previous 24 hours: 8 completed Kp bins, and 24 completed
hours of each of B, Bz, Vs, Pd, and Pt.

A solar-rotation block is added. For each issue and each solar-wind input,
lags from 25 to 29 days are tried in 12-hour steps. The lag is the one whose
speed, over a 5-day window, correlates best with the most recent 5 days.
The model uses Kp and mean speed at that lag. It also finds windows in the
same 25-to-29-day range whose Kp shape matches the last 48 hours, and uses
the Kp rise or fall over the 24 hours after those windows, weighted by the
match. Those changes had already ended, so they are not future data.

Persistence repeats the last observed Kp at every lead.
Skill = 1 - RMSE / RMSE of persistence.


How each piece is built
-----------------------

1. Source data

   Kp is the GFZ 3-hour index, https://kp.gfz.de/en/data
   Raw files are in Data/Origin/Kp/. The processed file is
   Data/Processed/Kp.csv. Training Kp runs through about January 2026.
   The realtime script uses one Kp source after that file: NOAA planetary Kp,
   https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json
   Only bins that have already ended are appended, and they are rounded to
   the nearest third. Bins already in Kp.csv are left as they are. The GFZ
   nowcast is not mixed into this live history. Training labels are unchanged.

   OMNI is https://omniweb.gsfc.nasa.gov/ow_min.html
   The download is 1-minute data. Preprocessing interpolates it and then
   takes hourly means, which reduces noise and fills short gaps. The hourly
   file is Data/Processed/OMNI_hourly.csv, with columns Date, B, Bz, Vs, Pd,
   Pt. B is |B| in nT, Bz is nT, Vs is km/s, Pd is proton density in cm^-3,
   and Pt is proton temperature in K. The OMNI file used here ends at
   2026-09-17 23:00. VARX training labels stop earlier, before the first
   RTSW model sample.

   RTSW is https://solarwind.engin.umich.edu/rtsw/
   The full series is https://solarwind.engin.umich.edu/rtsw/api/series/?t=earth
   Download_RTSW_2026.R writes Data/Origin/RTSW/Earth_obs.csv and
   Earth_model.csv. Columns are renamed to the OMNI names: b to B, bz to Bz,
   v to Vs, n to Pd, T to Pt.

2. Preprocessing

   Rcode/Preprocessing.R reads Data/Origin and writes the merged tables in
   Data/Aggregated/, plus Data/Processed/Kp.csv and OMNI_hourly.csv.
   The realtime script does not rerun this step. It starts from those two
   processed files and adds the nowcast and RTSW.

3. Offline evaluation

   Rcode/Main/Main_VARX_RawBaseline.R is the offline VARX evaluation on the
   RTSW period, including the January fit of OMNI = a + b * AWSoM.
   Its outputs are in Rcode/Main/output/Main_VARX_RawBaseline/.
   The realtime forecast no longer reads that frozen map. The script is kept
   separate from the realtime pages.

4. Realtime forecast

   Rcode/Main/Main_VARX_Realtime.R issues one forecast:

   * Download the RTSW Earth in-situ series and the AWSoM run.
   * Turn non-finite values and fill sentinels into missing values, interpolate
     onto a 1-minute grid, and average each clock hour. A missing hour that
     sits between two hours already known at the issue is filled by linear
     interpolation. If the series has stopped, the last hour is held for at
     most 3 hours, so one finished Kp bin can still be used.
   * Extend Kp after the processed file from NOAA planetary Kp only. Only
     ended bins are used, rounded to the nearest third.
   * Load or fit the 24-hour VARX above.
   * Walk back from the newest Kp bin until both in-situ and AWSoM have a
     complete 24-hour driver.

   The two inputs are:

   * In-situ measurements: the hourly mean of Earth_obs, with no adjustment.
     A gap of a few hours is interpolated or held, as above. A longer gap
     drops that issue.
   * AWSoM Earth model: the hourly mean of Earth_model, then
     OMNI = a + b * AWSoM. For each issue, a and b are fit on the latest 30
     days of paired hours that ended before this driver starts. Hours that
     are being scored are not in the fit. The slopes of |B| and density are
     near zero, so the map pulls those series toward the OMNI mean. It does
     not restore storm peaks.

   Rcode/Main/Write_VARX_Pages.R writes the English HTML. The realtime page
   shows this one 24-hour VARX, including all eight leads. The 6-hour
   endpoint is lead 2 of that model, not a second model.

   The next 27 days on that page is the same model, stepped forward one day
   at a time, in 3-hour bins. The first day is the 24-hour forecast. Later
   days read the solar wind that followed the closest aligned OMNI cycle,
   and the Kp just predicted becomes the next day's recent Kp. In-situ and
   AWSoM are the same two solar-wind inputs. Boosted is the same blend as
   the 24-hour scores: AWSoM is mixed in only when the latest observed Kp
   is high enough that earlier storms of that lead favored it. The chart
   also draws Kp already observed since the issue, and the Kp that followed
   the closest cycle. The match uses the last 27 days of in-situ solar wind.
   An hour with no in-situ sample uses that same OMNI hour when the hour is
   already in the OMNI file. A cycle is eligible only after its own next 27
   days of Kp have ended.

5. Baseline

   The same script scores issues whose target Kp has already been observed.
   The first 30 days after the RTSW model starts are used only to fit the
   AWSoM map and are not scored. The headline score is the bin that ends
   exactly 6 or 24 hours after issue time. Leads are not pooled. The
   by-lead table is on the same page. Below that, the 27-day outlook uses
   the same criteria on each 3-hour bin: RMSE, MAE, correlation, and skill.
   The reference is the Kp that followed the closest aligned OMNI cycle, not
   persistence. Skill = 1 - RMSE / RMSE of that cycle. The endpoint is the
   bin that ends exactly 27 days after issue time. A lead is scored as soon
   as that bin has ended. The scored forecast is the same 24-hour model
   stepped forward one day at a time, in in-situ, AWSoM, and boosted form.
   Later days read the solar wind that followed the closest aligned OMNI
   cycle. The chart shows Boosted. Buttons add the other two lines.
   Scoring uses Kp
   through the latest ended bin, the same extension as the realtime page.
   Day d is scored as soon as that day has ended, so the charts reach the
   current Kp. The hourly run rewrites this page with the realtime page.

6. Boosted

   Boosted starts from the in-situ VARX. When the last observed Kp is above
   a threshold, it blends in the AWSoM forecast. The threshold and the weight
   are chosen from earlier issues of the same lead whose targets have already
   ended, and only if that blend had a lower error there. On quiet issues the
   two forecasts are the same. At 24 hours the blend changes a few dozen
   issues and lowers RMSE from 1.218 to 1.212. Every forecast chart shows
   Boosted. Buttons on the chart add the in-situ and AWSoM lines. If the
   current Kp is below the threshold, boosted stays on the in-situ forecast.


Why the forecast rises more slowly than observed Kp
---------------------------------------------------

On the 24-hour endpoint chart, forecast peaks fall about 18 hours after
observed peaks. Bz, speed, and Kp available at issue time line up with Kp
at that time, not with Kp a day later. Current Bz correlates about 0.19 with
Kp 24 hours ahead and about 0.54 with current Kp. The 3-hour equation already
puts its largest solar-wind weights on the newest hour of Bz and speed. At
24 hours the Bz weights are near zero.

Pulling the forecast toward the latest Kp, toward the previous rotation's
Kp change, or by a recent mean bias all raised the 24-hour RMSE, and the
peaks did not move earlier. The issued VARX does not apply that stretch.
When observed Kp later reaches 5 or more, the forecast is still usually near
2, because the solar wind of that storm was not yet in the inputs.


Possible later changes
----------------------

  * If RTSW later publishes AWSoM solar wind for the 24 hours after issue
    time, that block can be added as a predictor. That is the route that
    could show a storm before Kp rises. The file now ends near the current
    hour, so it cannot be used.
  * Do not amplify the forecast so that the curve matches the observed
    spikes. The part that is not predictable raises RMSE when it is inflated.
  * Boosted separates from in-situ only at high Kp. Beating in-situ on
    ordinary issues needs a leading signal of the kind above, not another
    tree model on the same previous 24 hours of solar wind. Residual
    boosting on those inputs was worse than in-situ.
  * A neural VARX has not been trained. If it is, compare it with this
    linear model on the 24-hour endpoint, and look at Kp >= 6 separately.
  * The frozen January values of a and b remain only in the offline
    evaluation. Realtime refits them on the latest 30 days. Updating the
    map did not improve Kp skill.


Run it on a schedule and publish the pages
------------------------------------------

Kp arrives every three hours. The script walks back to the latest issue for
which both drivers are complete, so cron does not have to hit the bin
boundary. Every three hours is enough.

These commands use the full path, so they work from any directory, including
the home directory:

  chmod +x /Users/yuyian/pphd-files/Kp-Forecasting/Rcode/Main/run_realtime.sh
  /Users/yuyian/pphd-files/Kp-Forecasting/Rcode/Main/run_realtime.sh

The script runs Main_VARX_Realtime.R, then copies the pages and tables to
publish/. The copied files are

  publish/VARX_Kp_realtime.html          current forecast
  publish/VARX_Kp_baseline.html          scored comparison
  publish/VARX_Kp_realtime_latest.csv
  publish/VARX_Kp_realtime_27d.csv       next 27 days, every 3 hours
  publish/VARX_Kp_baseline_metrics.csv
  publish/VARX_Kp_baseline_endpoint.csv
  publish/VARX_Kp_baseline_27d_metrics.csv
  publish/VARX_Kp_baseline_27d.csv       scored 27-day outlook

The two HTML files link to each other by relative path and need no other
files, so they have to stay in the same directory. The working outputs,
including coefficients, drivers, and the run log, stay in
Rcode/Main/output/Main_VARX_Realtime/. The directory other people should
open is publish/.

GitHub Actions runs the same script every three hours, at minutes 27 and 52
after each Kp bin end (00, 03, 06, 09, 12, 15, 18, 21 UTC). The second start
covers a delayed or dropped first start. The laptop can be off.
The workflow is .github/workflows/forecast.yml. It publishes publish/ with
GitHub Pages. The page is

  https://yuyian.github.io/Kp-Forecasting/

and that address redirects to VARX_Kp_realtime.html. A run can also be
started by hand from the Actions tab. The raw one-minute OMNI files are not
in the git repository.

On this Mac the hourly job is still a LaunchAgent. It runs at minute 12 of
every hour while the user is logged in. If the Mac was asleep, the missed
hour runs once after wake. The log is publish/cron.log. A lock stops a
second run from starting while the first is still downloading.

  cp Rcode/Main/local.kp-forecast.plist ~/Library/LaunchAgents/
  launchctl bootout gui/$(id -u) ~/Library/LaunchAgents/local.kp-forecast.plist
  launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/local.kp-forecast.plist

Google Drive is not mounted on this Mac. The shared copy is the University
of Michigan Dropbox folder that already holds this project:

  ~/Library/CloudStorage/Dropbox-UniversityofMichigan/Yu Yian/Kp-Forecasting/publish

Dropbox syncs that folder. Share the publish folder from the Dropbox app
with the people who should open the pages. Only the pages and tables are
copied there. The R code in that Dropbox tree is not updated by this job.
Set KP_SHARE_DIR to use a different shared folder, including a Google Drive
folder after Drive for desktop is installed and signed in. KP_PUBLISH_DIR
still receives an extra copy when it is set.

publish/VARX_Kp_realtime.html is the current forecast, and it links to the
baseline page.

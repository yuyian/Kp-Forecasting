Kp, ap, and F10.7 forecast
==========================

The issued forecast is three quantities, not three outputs of one model.
The two pages are self-contained English HTML files and open in a browser
with no other files.

  https://yuyian.github.io/Kp-Forecasting/

That address opens the current forecast. The same directory has the
baseline page, which scores issues whose target has already been observed.


What is issued
--------------

A forecast is issued at the end of a finished 3-hour Kp bin (00, 03, 06,
09, 12, 15, 18, 21 UTC). Only data that had already ended by that time are
used. The Date on a Kp row is the start of that bin, so the bin is known
at its end. The Date on a solar-wind hour is the start of that hour, so
the hour is known at its end. Solar wind inside the forecast window is
not a predictor.

Kp, next 24 hours. Eight vanilla XGBoost models, one for each 3-hour lead
from 3 hours through 24 hours. Each model predicts that bin directly.
Predicted Kp is not fed back in. The inputs are the last eight Kp bins,
24 hours of near-Earth solar wind (B, Bz, speed, proton density, and
proton temperature), and a solar-rotation match from 25 to 29 days earlier.
The saved models are in
Rcode/Main/output/surtai_vanilla_xgboost_allhorizons/models/.
Training ends on 2026-01-25.

Kp, from 27 hours through 27 days. The in-situ regression, stepped forward
one day at a time. The first eight bins of that series are the XGBoost
values above. Later days read the solar wind that followed the closest
matching OMNI cycle, and the Kp just predicted becomes the next day's
recent Kp. AWSoM is the same regression with the AWSoM solar-wind block.
It is a comparison, not the issued line. On the baseline page that
comparison is labeled Issued AWSoM, next to Issued in-situ.

ap. The Bartels planetary amplitude of the issued Kp. Between scale points
it is interpolated in log(ap). It is not a separate model, so it follows
the XGBoost for the first 24 hours and the in-situ regression after that.

F10.7. A daily ridge regression on the 10.7 cm flux itself. It does not
use the Kp XGBoost. The Kp model reads geomagnetic and solar-wind inputs
and predicts Kp. The flux is a daily solar measurement, so the ridge uses
the flux's own history and the 27-day solar rotation. On 2019-2022 that
ridge beat repeating the latest measured flux, a 27-day recurrence, and a
blend of the two. Each run downloads the CelesTrak file SW-All.txt and
fits through the last measured day. Penticton measures at 20:00 UT, so at
a 12:00 UTC issue the flux for that UT day is still a forecast.


Pages
-----

The realtime page has one row of buttons: Kp, ap, F10.7, and solar wind
over the past day. The solar-wind button shows the five measured and
AWSoM series for the 24 hours before the issue. It is not a forecast.
Below that, Next 24 hours and Next 27 days switch the chart. One quantity
and one range are visible at a time. History of the selected quantity
stays on the page when the range changes.

The baseline page scores the same issued models on 2026 issues whose
target has already ended. Kp, ap, and F10.7 each have their own table.
For Kp and ap after 24 hours, the rows are Issued in-situ, Issued AWSoM,
and the best matched solar cycle. Skill for those rows is against that
cycle. F10.7 is compared with persistence, which repeats the last
measured flux.


Data
----

Kp training is the GFZ 3-hour index, https://kp.gfz.de/en/data
The processed file is Data/Processed/Kp.csv. Bins after that file come
from NOAA planetary Kp,
https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json
rounded to the nearest third. A bin NOAA no longer has is filled from
the GFZ nowcast. Training labels are unchanged.

OMNI is https://omniweb.gsfc.nasa.gov/ow_min.html
The hourly file is Data/Processed/OMNI_hourly.csv. B is |B| in nT, Bz is
nT, Vs is km/s, Pd is proton density in cm^-3, and Pt is proton
temperature in K. Flow pressure is not Pd.

RTSW is https://solarwind.engin.umich.edu/rtsw/
The realtime script downloads the Earth in-situ series and the AWSoM run.
In-situ measurements are the hourly mean of Earth_obs, with no adjustment.
AWSoM is the hourly mean of Earth_model, then OMNI = a + b * AWSoM. For
each issue, a and b are fit on the latest 30 days of paired hours that
ended before the driver starts.

F10.7 is https://celestrak.org/SpaceData/SW-All.txt
The observed section is the training record. A copy is kept at
Data/Origin/F107/SW-All.txt and replaced on each run.


Run it and publish the pages
----------------------------

Kp arrives every three hours. The script walks back to the latest issue
for which both solar-wind drivers are complete, so the timer does not
have to hit the bin boundary. GitHub Actions runs the same script every
hour at minute 17 UTC. The workflow is .github/workflows/forecast.yml.
It publishes the pages with GitHub Pages.

  chmod +x Rcode/Main/run_realtime.sh
  Rcode/Main/run_realtime.sh

The script runs Main_VARX_Realtime.R, then copies the pages and tables
to publish/. The two HTML files link to each other and have to stay in
the same directory. Working outputs stay in
Rcode/Main/output/Main_VARX_Realtime/.

On this Mac a LaunchAgent runs at minute 12 of every hour while the user
is logged in, and that start also launches the GitHub run. If the Mac
was asleep, the missed hour runs once after wake. The log is
publish/cron.log.

  cp Rcode/Main/local.kp-forecast.plist ~/Library/LaunchAgents/
  launchctl bootout gui/$(id -u) ~/Library/LaunchAgents/local.kp-forecast.plist
  launchctl bootstrap gui/$(id -u) ~/Library/LaunchAgents/local.kp-forecast.plist

A second copy goes to the shared folder when that disk is mounted.
Set KP_SHARE_DIR to choose the folder. KP_PUBLISH_DIR receives an extra
copy when it is set.

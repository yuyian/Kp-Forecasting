OMNI-lead experiment
=====================

Code:  Rcode/Main/OmniLead/Main_VARX_OmniLead.R
Scores: Rcode/Main/OmniLead/results/
Run from the repository root:

  /usr/local/bin/Rscript Rcode/Main/OmniLead/Main_VARX_OmniLead.R

This does not change the realtime forecast.

What a perfect OMNI forecast would be worth
--------------------------------------------
The Kp model trained on the true next 24 hours of OMNI, and tested on those
same true hours, has 24 h RMSE 0.664 (correlation 0.859) on the 1627 issues
where OMNI still exists. The published in-situ VARX on those issues is 1.212.
That gap is the value of knowing tomorrow's solar wind. It is not an
operational forecast.

How tomorrow's OMNI is predicted
--------------------------------
Three forecasts are built. The one used for Kp is the stack.

  Linear. One least-squares model maps the previous 24 hours of B, Bz, Vs,
  Pd and Pt to each feature at each lead from 1 to 24 hours. Fit on OMNI
  hours whose target ends before 2016-01-01 (184033 hours, rank 121/121).

  Solar rotation. The lag from 25 to 29 days whose last 5 days of speed best
  match the origin is selected, and the solar wind that followed that window
  is copied in as the next 24 hours. Negative correlations are discarded.

  Stack. On 2016-2018 (8768 origins), each feature and lead is a regression
  on the linear forecast, the rotation path and persistence of the last hour.
  A depth-3 boosted model for Bz and speed was fit on the same years. On
  2019 it was worse at lead 24 (Bz 1.965 vs 1.960, speed 54.2 vs 52.8), so
  the stack is kept and the boosted model is not used.

Lead-24 stack weights (intercept, linear, rotation, persistence):

  B    0.44    0.77    0.14    -0.04
  Bz  -0.06    0.56    0.03    -0.01
  Vs   -68     1.18    0.20    -0.22
  Pd  -0.59    1.29    0.02    -0.16
  Pt  -20263   1.10    0.17    -0.05

The rotation term is small. At a 24 h horizon the previous rotation is a
worse speed forecast than repeating the last hour, because speed changes
more slowly than the error in a 27-day lag. Rotation would matter at leads
of several days. It does not matter inside one day.

On calendar 2024, 24 hours ahead, paired RMSE:

  feature   persistence   linear   rotation   stack
  B            4.90        3.81      5.30      3.91
  Bz           6.12        4.29      5.87      4.28
  Vs          65.8        62.1     107.2      65.0
  Pd           7.27        5.47      7.27      5.54
  Pt         132600      101500    125900    102400

The stack does not beat the linear forecast. Speed at 24 h is already close
to persistence, and Bz skill stays near +0.30.

Why the forecast is not pasted on as extra columns
--------------------------------------------------
Each linear predicted hour is a linear function of the past 24 hours. A Kp
model that already contains those hours has the same fit after that forecast
is appended (rank stays 147). The rotation path is not a linear function of
the last 24 hours, but by itself it is a poor 24 h forecast, and the stack
barely uses it.

The Kp model is therefore trained on the true next 24 hours of OMNI
(29130 issues, 2016 through the RTSW start, rank 267). At issue time those
columns are filled with a forecast. A second model adds nonlinear summaries
(minimum Bz, fraction of hours with Bz below 0, maximum speed, maximum |B|)
of the stack forecast. That summary model is trained only from 2020, so it
does not see the years that chose the stack weights (17703 issues, rank 155).

How AWSoM is used inside the forecast window
---------------------------------------------
Where the archive contains the 24 hours after an issue, those hours are
mapped with OMNI = a + b * AWSoM, fit on the first 30 days of overlap
(721 hours, through 2026-02-24 16:14 UTC) and then held fixed. They fill the
same columns as the statistical forecast.

On the overlap, lead-24 error against OMNI is:

  source                         Bz RMSE    Vs RMSE
  repeat the last hour            4.19       77.6
  linear forecast                 3.09       62.6
  solar-rotation path             4.44      107.4
  rotation stack                  3.09       62.6
  calibrated AWSoM                3.08       94.8

AWSoM matches the linear forecast on Bz and is worse than persistence on
speed. The stored series is an Earth-time nowcast that ends near the current
hour, so this is an upper bound for a feed that issued those hours at the
issue time.

24 h Kp endpoint, same 1856 issues as the published baseline
-------------------------------------------------------------
Kp.csv ends 2026-10-15 12:00 UTC, so 23 later baseline issues are left out.

  model                                      n     RMSE    corr
  Previous boosted                         1856    1.212   0.376
  Previous in-situ VARX                    1856    1.218   0.370
  Perfect prog, rotation stack             1856    1.211   0.379
  Trained on the rotation-stack forecast   1856    1.210   0.378
  Summaries of the rotation stack          1856    1.212   0.374
  Trained on the linear forecast           1856    1.223   0.359
  Retrained, past solar wind only          1856    1.223   0.359
  Perfect prog, linear 24 h forecast       1856    1.227   0.359
  Perfect prog, AWSoM next 24 h            1856    1.306   0.229
  Oracle, true OMNI over the next 24 h     1627    0.664   0.859

The stack ties the published boosted blend. Training on the stack forecast
from 2020, so the forecast is already out of sample, gives 1.210 at 24 h and
a worse error at 3-12 h. Training on the linear forecast stays at rank 147
and matches past-only at 1.223, because that forecast is a linear function of
the past 24 hours. The realtime model is unchanged.

Files written under results/
----------------------------
  kp_24h_comparison.csv        the table above
  kp_by_lead.csv               RMSE at 3, 6, ..., 24 h
  kp_24h_predictions.csv       one row per issue
  omni_forecast_skill_2024.csv persistence, linear, rotation, stack, by lead
  tomorrow_block_vs_omni.csv   each forecast against OMNI on the test issues
  OmniLead_coef.RData          linear model, stack weights, Kp fits

Longer memory for the solar-wind forecast
-----------------------------------------
Main_SolarForecast.R adds what the 24 h model cannot see: daily means of |B|,
Bz, speed and density over the previous six days, and the change the best
25-29 day rotation made over the same lead (not the level). Fit before 2016
on 88980 hours, scored on 2024. Rank of the full model is 218/218. Lead-24
absolute-weight shares are about 0.60-0.68 on the last day, 0.08-0.10 on the
daily means, and 0.24-0.30 on the rotation increment.

2024 RMSE, linear versus daily means versus daily+rotation:

  feature  lead   linear   daily   daily+rotation
  Bz        3 h    3.32     3.33      3.33
  Bz       12 h    3.98     3.98      3.98
  Bz       24 h    4.06     4.07      4.07
  Vs        3 h   22.49    22.47     22.49
  Vs       12 h   45.76    45.66     45.65
  Vs       24 h   61.84    61.76     61.88
  B        24 h    3.88     3.83      3.85

The extra columns do not beat the 24 h linear forecast at 3, 6, 9, 12 or 24 h.
The one visible difference is |B| at 24 h, 3.88 to 3.83. That is too small to
move Kp. Full table: results/solar_forecast_extended_2024.csv.

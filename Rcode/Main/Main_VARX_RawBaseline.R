# =============================================================================
# Main_VARX_RawBaseline.R
# Direct VARX baseline for the 3-hourly Kp index
# =============================================================================
#
# Forecast the processed Kp series from
#   Data/Processed/Kp.csv            (3-hourly Kp; Date is the bin START)
#   Data/Processed/OMNI_hourly.csv   (hourly B, Bz, Vs, Pd, Pt; Date is the hour START)
#
# Estimator: the same unpenalized VARX used by the geomagnetic raw baseline
# (thin-SVD least squares, rank cutoff ols_tol, no ridge, no CV). One model
# per horizon. The forecast is direct: future Kp is never fed back in.
#
# Information set (no look-ahead)
# --------------------------------
# A forecast is issued at time tau, the END of a completed 3-hour Kp bin.
# Predictors at tau may contain only:
#   * Kp bins whose interval has already ended (end time <= tau)
#   * OMNI hours whose interval has already ended (hour start + 1h <= tau)
# The response is Kp strictly after tau. Nothing from the forecast window
# (neither future Kp nor future OMNI) enters X.
#
# Because Kp is stamped at the bin start, the row dated t covers [t, t+3h)
# and is not treated as known until t+3h. An OMNI row dated h covers [h, h+1h)
# and is not treated as known until h+1h.
#
# Horizons
# --------
#   6h  -> the next 2 Kp bins, ending at tau+3h and tau+6h
#   24h -> the next 8 Kp bins, ending at tau+3h, tau+6h, ..., tau+24h
# Each horizon has its own model. Metrics are reported on the endpoint bin,
# the Kp interval that ends exactly H hours after tau. Intermediate leads
# are in the by-lead table.
#
# State
# -----
# The lag-1 state is the preceding `history_hr` hours (default 24):
#   8 completed Kp bins + 24 completed hours of each OMNI feature.
# history_hr must be a multiple of 3. This is the raw-curve analogue of the
# source script's lag-1 window: no FPCA and no averaging of OMNI inside a
# Kp bin.
#
# Split
# -----
# Training labels are Kp bins that end at or before the first RTSW model
# sample (2026-01-25 16:14 UTC in the current file). Everything before
# RTSW starts is training. The test window is the RTSW
# model span (Data/Origin/RTSW/Earth_model.csv), and only issue times whose
# full history_hr driver window also exists in OMNI_hourly. That is the
# stretch on which the same coefficients can be fed either the original
# OMNI values or the RTSW model values.
#
# RTSW is reduced the same way as OMNI: non-finite values and OMNI fill
# sentinels (9999.99 for B, Bz, Vs, Pd; 9999999 for Pt) become missing,
# gaps are linearly interpolated in time, then each clock hour is the mean
# of that hour. Kp history is identical across test passes; only
# B, Bz, Vs, Pd, Pt change.
#
# Calibration
# -----------
# RTSW is not used until it starts, so the map from RTSW to OMNI can only
# be learned on the overlap. The first `cal_days` days of that overlap are
# the calibration window. An affine map (OMNI = a + b * RTSW) and a
# quantile map are fit there, per feature, on hours that contain a real
# RTSW sample. Kp scores for those maps use only later issue times, whose
# 24 h driver window starts after the calibration window. An affine map
# fit on the whole overlap is also reported; it sees the OMNI values of
# the hours it is applied to, so it is a ceiling, not an operational result.
#
# Dependencies: base R only.
# Run from anywhere:  Rscript Rcode/Main/Main_VARX_RawBaseline.R
# =============================================================================

rm(list = ls())

script_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path = normalizePath(sub("^--file=", "", script_arg[1]))
  repo_root = normalizePath(file.path(dirname(script_path), "..", ".."))
} else {
  repo_root = normalizePath("~/pphd-files/Kp-Forecasting")
}
setwd(repo_root)

# ---- parameters --------------------------------------------------------------
history_hr   = 24L                 # length of the lag-1 state, in hours
horizons_hr  = c(6L, 24L)
omni_features = c("B", "Bz", "Vs", "Pd", "Pt")
rtsw_model_file = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_model.csv")
cal_days = 30L
ols_tol      = 1e-10
kp_step_hr   = 3L

stopifnot(history_hr %% kp_step_hr == 0)
stopifnot(all(horizons_hr %% kp_step_hr == 0))
stopifnot(all(horizons_hr > 0))

out_dir = file.path(repo_root, "Rcode", "Main", "output", "Main_VARX_RawBaseline")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- load --------------------------------------------------------------------
Kp_df = read.csv(file.path(repo_root, "Data", "Processed", "Kp.csv"),
                 stringsAsFactors = FALSE)
OMNI_df = read.csv(file.path(repo_root, "Data", "Processed", "OMNI_hourly.csv"),
                   stringsAsFactors = FALSE)

# Midnight stamps are written as "YYYY-MM-DD" with no clock time. A mixed
# vector makes as.POSIXct() guess "%Y-%m-%d" and drop every hour. Pad those
# rows, then parse the whole column with an explicit UTC format.
parse_utc = function(x) {
  x = trimws(as.character(x))
  date_only = grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  x[date_only] = paste0(x[date_only], " 00:00:00")
  out = as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (anyNA(out)) stop(sprintf("Unparsed timestamps, first: %s", x[which(is.na(out))[1]]))
  out
}
Kp_df$Date   = parse_utc(Kp_df$Date)
OMNI_df$Date = parse_utc(OMNI_df$Date)
Kp_df   = Kp_df[order(Kp_df$Date), ]
OMNI_df = OMNI_df[order(OMNI_df$Date), ]
stopifnot(!anyDuplicated(Kp_df$Date), !anyDuplicated(OMNI_df$Date))
stopifnot(all(omni_features %in% names(OMNI_df)))
stopifnot("Data" %in% names(Kp_df))

kp_time   = Kp_df$Date
kp_value  = Kp_df$Data
omni_time = OMNI_df$Date
OMNI_mat  = as.matrix(OMNI_df[, omni_features, drop = FALSE])

kp_dt_hr   = as.numeric(diff(kp_time), units = "hours")
omni_dt_hr = as.numeric(diff(omni_time), units = "hours")
stopifnot(length(kp_dt_hr) == 0 || all(abs(kp_dt_hr - kp_step_hr) < 1e-6))
stopifnot(length(omni_dt_hr) == 0 || all(abs(omni_dt_hr - 1) < 1e-6))
stopifnot(abs(as.numeric(kp_time[1] - omni_time[1], units = "secs")) < 1)

# Integer alignment: Kp index i (1-based) ends at kp_time[i] + 3h.
# The last OMNI hour that has ended by then starts at that instant minus 1h,
# which is OMNI row 3*i when both series share an epoch and have no gaps.
N_kp = length(kp_time)
N_omni = nrow(OMNI_mat)
i_probe = min(8L, N_kp)
stopifnot(abs(as.numeric(omni_time[3 * i_probe] - (kp_time[i_probe] + 2 * 3600),
                         units = "secs")) < 1)

cat(sprintf("Kp:   %d bins, %s  ->  %s\n", N_kp, kp_time[1], kp_time[N_kp]))
cat(sprintf("OMNI: %d hours, %s  ->  %s\n", N_omni, omni_time[1], omni_time[N_omni]))
cat(sprintf("State: previous %dh (%d Kp bins + %d hours x %d OMNI features)\n",
            history_hr, history_hr / kp_step_hr, history_hr, length(omni_features)))

# RTSW model -> hourly, same two steps as Utils.R / Preprocessing.R:
# fill sentinels become NA, linear interpolation in time, then the mean
# of each clock hour. Hours with no native sample are the interpolated ones.
floor_hour = function(t) {
  as.POSIXct(as.numeric(t) - (as.numeric(t) %% 3600), origin = "1970-01-01", tz = "UTC")
}
make_rtsw_hourly = function(path) {
  raw = read.csv(path, stringsAsFactors = FALSE)
  raw$Date = parse_utc(raw$Date)
  raw = raw[order(raw$Date), ]
  stopifnot(all(omni_features %in% names(raw)))
  time = raw$Date
  mat = as.matrix(raw[, omni_features, drop = FALSE])
  # OMNI high-resolution fill is 9999.99 on B, Bz, Vs, and Pd, and 9999999
  # on proton temperature. RTSW temperatures of a few 1e5 K are real.
  sentinel = c(B = 9999, Bz = 9999, Vs = 9999, Pd = 9999, Pt = 9999999)
  for (nm in omni_features) {
    mat[!is.finite(mat[, nm]) | mat[, nm] >= sentinel[nm], nm] = NA
  }
  ok_any = apply(mat, 1, function(r) any(is.finite(r)))
  time = time[ok_any]
  mat = mat[ok_any, , drop = FALSE]

  hour0 = floor_hour(min(time))
  hour1 = floor_hour(max(time))
  hours = seq(hour0, hour1, by = "hour")
  breaks = c(as.numeric(hours), as.numeric(hours[length(hours)]) + 3600)
  native_bin = findInterval(as.numeric(time), breaks)
  n_native = tabulate(native_bin, nbins = length(hours))

  grid = seq(min(time), max(time), by = "min")
  grid_bin = findInterval(as.numeric(grid), breaks)
  out = matrix(NA_real_, length(hours), length(omni_features))
  colnames(out) = omni_features
  for (j in seq_along(omni_features)) {
    y = mat[, j]
    ok = is.finite(y)
    yg = approx(as.numeric(time[ok]), y[ok], xout = as.numeric(grid),
                method = "linear", rule = 1)$y
    keep = grid_bin >= 1L & grid_bin <= length(hours) & is.finite(yg)
    means = tapply(yg[keep], grid_bin[keep], mean)
    out[as.integer(names(means)), j] = as.numeric(means)
  }
  list(time = hours, mat = out, n_native = n_native,
       native_start = min(time), native_end = max(time))
}

RTSW_h = make_rtsw_hourly(rtsw_model_file)
rtsw_ok = apply(RTSW_h$mat, 1, function(r) all(is.finite(r)))
# A training target must finish before the first RTSW sample, so January
# 2026 up to the RTSW start is training and the RTSW span stays test-only.
train_end = RTSW_h$native_start
cat(sprintf("Train labels end at or before RTSW start: %s\n", train_end))
cat(sprintf("RTSW model native: %s -> %s\n", RTSW_h$native_start, RTSW_h$native_end))
cat(sprintf("RTSW hourly: %d hours, %s -> %s, %d hours filled by interpolation\n",
            sum(rtsw_ok), RTSW_h$time[which(rtsw_ok)[1]],
            RTSW_h$time[max(which(rtsw_ok))],
            sum(rtsw_ok & RTSW_h$n_native == 0)))
write.csv(data.frame(Date = format(RTSW_h$time, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
                     RTSW_h$mat, n_native = RTSW_h$n_native,
                     filled = RTSW_h$n_native == 0,
                     stringsAsFactors = FALSE),
          file.path(repo_root, "Data", "Origin", "RTSW", "Earth_model_hourly.csv"),
          row.names = FALSE)

# ---- design helpers ----------------------------------------------------------
n_kp_hist = history_hr / kp_step_hr

# Issue time of the forecast that has just observed Kp bin i.
issue_time = function(i) kp_time[i] + kp_step_hr * 3600

# X rows for last-observed Kp indices `i`: intercept, Kp history, OMNI history.
# Columns run oldest -> newest. Newest Kp bin ends at tau; newest OMNI hour
# ends at tau. The next OMNI hour (index 3*i + 1) starts at tau and is excluded.
build_X = function(i) {
  i = as.integer(i)
  kp_end_age = seq(from = history_hr - kp_step_hr, to = 0, by = -kp_step_hr)
  Kp_hist = sapply(seq_len(n_kp_hist), function(j) kp_value[i - n_kp_hist + j])
  colnames(Kp_hist) = sprintf("Kp_ended_%02dh_before_tau", kp_end_age)

  omni_end = 3L * i
  omni_age = (history_hr - 1L):0L
  omni_blocks = lapply(seq_along(omni_features), function(f) {
    v = OMNI_mat[, f]
    block = sapply(omni_age, function(age) v[omni_end - age])
    colnames(block) = sprintf("%s_ended_%02dh_before_tau", omni_features[f], omni_age)
    block
  })
  cbind(Intercept = 1, Kp_hist, do.call(cbind, omni_blocks))
}

build_Y = function(i, horizon_hr) {
  n_ahead = horizon_hr / kp_step_hr
  leads = seq(kp_step_hr, horizon_hr, by = kp_step_hr)   # bin END offset from tau
  Y = sapply(seq_len(n_ahead), function(k) kp_value[i + k])
  colnames(Y) = sprintf("Kp_ending_tau_plus_%02dh", leads)
  Y
}

standardize_fit = function(M) {
  ctr = colMeans(M)
  scl = apply(M, 2, sd)
  scl[scl == 0 | is.na(scl)] = 1
  list(center = ctr, scale = scl,
       M = sweep(sweep(M, 2, ctr, "-"), 2, scl, "/"))
}

apply_scale = function(M, center, scale) {
  sweep(sweep(M, 2, center, "-"), 2, scale, "/")
}

ols_coef_from_svd = function(X, Y, tol = 1e-10) {
  k = min(nrow(X), ncol(X))
  sv = svd(X, nu = k, nv = k)
  keep = sv$d > tol * max(sv$d)
  Dinv = ifelse(keep, 1 / sv$d, 0)
  UtY = crossprod(sv$u, Y)
  B = t(sv$v %*% (Dinv * UtY))          # K x p
  list(B = B, n_dropped = sum(!keep), rank = sum(keep))
}

rmse_ = function(y, yhat) sqrt(mean((y - yhat)^2, na.rm = TRUE))
mae_  = function(y, yhat) mean(abs(y - yhat), na.rm = TRUE)
cor_  = function(y, yhat) suppressWarnings(cor(as.vector(y), as.vector(yhat), use = "complete.obs"))

metric_row = function(horizon_hr, scope, input, y, yhat, ypers) {
  data.frame(
    horizon_hr = horizon_hr,
    input = input,
    scope = scope,
    n = sum(is.finite(y) & is.finite(yhat)),
    RMSE = rmse_(y, yhat),
    MAE = mae_(y, yhat),
    Correlation = cor_(y, yhat),
    RMSE_persistence = rmse_(y, ypers),
    MAE_persistence = mae_(y, ypers),
    skill_RMSE = 1 - rmse_(y, yhat) / rmse_(y, ypers),
    stringsAsFactors = FALSE
  )
}

# For each issue index, the history_hr driver hours (oldest -> newest) and
# whether every one of them exists in both the RTSW hourly series and OMNI.
# exog follows `hour_values` (default: raw RTSW). Column order matches
# build_X: one feature block, oldest hour to newest.
rtsw_driver = function(i, hour_values = RTSW_h$mat) {
  tau_i = issue_time(i)
  age = history_hr:1L
  hour_mat = outer(as.numeric(tau_i), age * 3600, `-`)
  rtsw_num = as.numeric(RTSW_h$time)
  omni_num = as.numeric(omni_time)
  n = length(i)
  exog = matrix(NA_real_, n, history_hr * length(omni_features))
  all_hit = rep(TRUE, n)
  n_filled = integer(n)
  col = 1L
  for (f in seq_along(omni_features)) {
    for (a in seq_along(age)) {
      hit = match(hour_mat[, a], rtsw_num)
      hit_o = match(hour_mat[, a], omni_num)
      ok = !is.na(hit) & !is.na(hit_o) & rtsw_ok[hit]
      ok[is.na(ok)] = FALSE
      all_hit = all_hit & ok
      exog[, col] = hour_values[hit, f]
      n_filled = n_filled + as.integer(!is.na(hit) & RTSW_h$n_native[hit] == 0)
      col = col + 1L
    }
  }
  list(exog = exog, all_hit = all_hit, n_filled = n_filled)
}

# Hours where both series exist. Calibration uses only real RTSW samples
# (n_native > 0) in the first cal_days after RTSW starts.
pair_idx = match(as.numeric(RTSW_h$time), as.numeric(omni_time))
paired = rtsw_ok & !is.na(pair_idx)
cal_end = train_end + cal_days * 86400
in_cal = paired & RTSW_h$n_native > 0 & RTSW_h$time < cal_end
in_eval = paired & RTSW_h$n_native > 0 & RTSW_h$time >= cal_end

affine_fit = function(x, y) {
  ok = is.finite(x) & is.finite(y)
  x = x[ok]; y = y[ok]
  if (length(x) < 5 || var(x) == 0) return(c(a = mean(y) - mean(x), b = 1))
  b = cov(x, y) / var(x)
  c(a = mean(y) - b * mean(x), b = b)
}
qmap_apply = function(x, x_cal, y_cal) {
  ok = is.finite(x_cal) & is.finite(y_cal)
  xc = x_cal[ok]
  yc = y_cal[ok]
  n = length(xc)
  xs = sort(xc)
  r = findInterval(x, xs)
  p = pmin(pmax(r / n, 0.5 / n), 1 - 0.5 / n)
  as.numeric(quantile(yc, p, names = FALSE, type = 7))
}

shift_rows = data.frame()
ab_cal = data.frame(feature = omni_features, a = NA_real_, b = NA_real_, n_cal = NA_integer_)
ab_all = ab_cal
RTSW_affine = RTSW_h$mat
RTSW_qmap = RTSW_h$mat
RTSW_affine_all = RTSW_h$mat
for (f in seq_along(omni_features)) {
  nm = omni_features[f]
  x = RTSW_h$mat[, f]
  y = rep(NA_real_, length(x))
  y[paired] = OMNI_mat[pair_idx[paired], f]
  fit_c = affine_fit(x[in_cal], y[in_cal])
  fit_a = affine_fit(x[paired & RTSW_h$n_native > 0], y[paired & RTSW_h$n_native > 0])
  ab_cal$a[f] = fit_c["a"]; ab_cal$b[f] = fit_c["b"]; ab_cal$n_cal[f] = sum(in_cal)
  ab_all$a[f] = fit_a["a"]; ab_all$b[f] = fit_a["b"]
  RTSW_affine[, f] = fit_c["a"] + fit_c["b"] * x
  RTSW_affine_all[, f] = fit_a["a"] + fit_a["b"] * x
  RTSW_qmap[, f] = qmap_apply(x, x[in_cal], y[in_cal])
  for (period in c("calibration", "evaluation")) {
    m = if (period == "calibration") in_cal else in_eval
    shift_rows = rbind(shift_rows, data.frame(
      feature = nm, period = period, n = sum(m),
      mean_OMNI = mean(y[m]), mean_RTSW = mean(x[m]),
      bias_RTSW_minus_OMNI = mean(x[m] - y[m]),
      sd_OMNI = sd(y[m]), sd_RTSW = sd(x[m]),
      cor = cor(x[m], y[m]),
      stringsAsFactors = FALSE))
  }
}
cat(sprintf("\nCalibration window: RTSW hours with a native sample before %s\n", cal_end))
cat("Per-feature shift (RTSW minus OMNI) and held-out affine map OMNI = a + b * RTSW:\n")
print(cbind(shift_rows[shift_rows$period == "evaluation", c("feature", "n", "bias_RTSW_minus_OMNI", "sd_OMNI", "sd_RTSW", "cor")],
            a = ab_cal$a, b = ab_cal$b), row.names = FALSE, digits = 3)
write.csv(shift_rows, file.path(out_dir, "RTSW_distribution_shift.csv"), row.names = FALSE)
write.csv(ab_cal, file.path(out_dir, "RTSW_affine_calibration.csv"), row.names = FALSE)

# Eligible issue indices: full Kp history, full OMNI history, and at least
# the longest horizon still inside the Kp record. Per-horizon masks below
# trim the target end and the train/test cut.
n_ahead_max = max(horizons_hr) / kp_step_hr
i_all = seq.int(n_kp_hist, N_kp - n_ahead_max)
i_all = i_all[3L * i_all <= N_omni]
i_all = i_all[(3L * i_all - history_hr + 1L) >= 1L]

# Causality of the index rule itself (every eligible origin, not a sample).
# History bin i ends at tau; target bin i+1 starts at tau; OMNI row 3*i
# starts at tau-1h and OMNI row 3*i+1 starts at tau.
tau_all = issue_time(i_all)
stopifnot(all(abs(as.numeric(kp_time[i_all] + kp_step_hr * 3600 - tau_all, units = "secs")) < 1))
stopifnot(all(kp_time[i_all + 1L] == tau_all))
stopifnot(all(abs(as.numeric(omni_time[3L * i_all] + 3600 - tau_all, units = "secs")) < 1))
# The hour AFTER the state starts at tau. It is not a predictor. The final
# origin can sit on the last OMNI hour, so that next row need not exist.
has_next_omni = (3L * i_all + 1L) <= N_omni
stopifnot(all(omni_time[3L * i_all[has_next_omni] + 1L] == tau_all[has_next_omni]))

cat(sprintf("Eligible origins with a full %dh target still in the Kp record: %d\n",
            max(horizons_hr), length(i_all)))
cat(sprintf("  first issue %s   last issue %s\n", tau_all[1], tau_all[length(tau_all)]))

Results = data.frame()
Results_lead = data.frame()
Results_extreme = data.frame()

for (horizon_hr in horizons_hr) {
  n_ahead = horizon_hr / kp_step_hr
  # Origins that still have this horizon inside the Kp record. i_all was
  # cut for the longest horizon, so every i_all is valid for shorter ones
  # as well; recompute the upper bound explicitly so a single-horizon edit
  # stays correct.
  i_h = seq.int(n_kp_hist, N_kp - n_ahead)
  i_h = i_h[3L * i_h <= N_omni]
  i_h = i_h[(3L * i_h - history_hr + 1L) >= 1L]

  tau = issue_time(i_h)
  target_end = tau + horizon_hr * 3600
  is_train = target_end <= train_end
  driver = rtsw_driver(i_h)
  # Paired test: the whole driver window sits in the RTSW model hourly
  # series and in OMNI, so the two inputs are scored on the same issues.
  is_test = driver$all_hit & !is_train
  stopifnot(!any(is_train & is_test))

  if (!any(is_train)) stop(sprintf("Horizon %dh: no training origins.", horizon_hr))
  if (!any(is_test))  stop(sprintf("Horizon %dh: no test origins in the RTSW model window.", horizon_hr))

  cat(sprintf("\n=== Horizon %dh ===\n", horizon_hr))
  cat(sprintf("Train origins: %d   issue %s -> %s   (targets end by %s)\n",
              sum(is_train), tau[which(is_train)[1]], tau[max(which(is_train))],
              train_end))
  cat(sprintf("Test  origins: %d   issue %s -> %s\n",
              sum(is_test), tau[which(is_test)[1]], tau[max(which(is_test))]))
  cat(sprintf("Test origins whose 24h RTSW input contains a filled hour: %d\n",
              sum(is_test & driver$n_filled > 0)))

  X = build_X(i_h)
  Y = build_Y(i_h, horizon_hr)
  Y_pers = matrix(kp_value[i_h], nrow = length(i_h), ncol = n_ahead)

  stopifnot(nrow(X) == length(i_h), nrow(Y) == length(i_h))
  stopifnot(!anyNA(X[is_train | is_test, ]), !anyNA(Y[is_train | is_test, ]))
  stopifnot(!anyNA(driver$exog[is_test, ]))

  pred_cols = setdiff(colnames(X), "Intercept")
  sc = standardize_fit(X[is_train, pred_cols, drop = FALSE])
  X_train = cbind(Intercept = 1, apply_scale(X[is_train, pred_cols, drop = FALSE],
                                             sc$center, sc$scale))
  X_omni  = cbind(Intercept = 1, apply_scale(X[is_test, pred_cols, drop = FALSE],
                                             sc$center, sc$scale))
  # Same Kp history and same training scale. Only the solar-wind block changes.
  exog_ix = (n_kp_hist + 1L):length(pred_cols)
  rtsw_scaled = apply_scale(driver$exog[is_test, , drop = FALSE],
                            sc$center[exog_ix], sc$scale[exog_ix])
  X_rtsw = X_omni
  X_rtsw[, exog_ix + 1L] = rtsw_scaled

  predict_exog = function(exog_raw) {
    scaled = apply_scale(exog_raw, sc$center[exog_ix], sc$scale[exog_ix])
    Xc = X_omni
    Xc[, exog_ix + 1L] = scaled
    Xc %*% t(fit$B)
  }
  map_affine = function(exog, a, b) {
    out = exog
    for (f in seq_along(omni_features)) {
      ix = ((f - 1L) * history_hr + 1L):(f * history_hr)
      out[, ix] = a[f] + b[f] * exog[, ix]
    }
    out
  }
  map_q = function(exog) {
    out = exog
    for (f in seq_along(omni_features)) {
      ix = ((f - 1L) * history_hr + 1L):(f * history_hr)
      out[, ix] = qmap_apply(as.vector(exog[, ix]),
                             RTSW_h$mat[in_cal, f], OMNI_mat[pair_idx[in_cal], f])
    }
    out
  }

  fit = ols_coef_from_svd(X_train, Y[is_train, , drop = FALSE], ols_tol)
  raw_exog = driver$exog[is_test, , drop = FALSE]
  Yhat_omni = X_omni %*% t(fit$B)
  Yhat_rtsw = X_rtsw %*% t(fit$B)
  Yhat_affine = predict_exog(map_affine(raw_exog, ab_cal$a, ab_cal$b))
  Yhat_qmap = predict_exog(map_q(raw_exog))
  Yhat_affine_all = predict_exog(map_affine(raw_exog, ab_all$a, ab_all$b))
  colnames(Yhat_omni) = colnames(Y)

  cat(sprintf("OLS rank %d / %d predictors (dropped %d singular directions)\n",
              fit$rank, ncol(X_train), fit$n_dropped))

  y_true = Y[is_test, , drop = FALSE]
  y_pers = Y_pers[is_test, , drop = FALSE]
  leads = seq(kp_step_hr, horizon_hr, by = kp_step_hr)
  n_filled_test = driver$n_filled[is_test]
  # Driver window starts at tau - history_hr, so this excludes every hour
  # that was used to fit the held-out maps.
  after_cal = tau[is_test] - history_hr * 3600 >= cal_end
  cat(sprintf("Issues whose driver window is entirely after calibration: %d\n",
              sum(after_cal)))

  score_named = function(scope, y, pers, preds) {
    out = data.frame()
    for (nm in names(preds))
      out = rbind(out, metric_row(horizon_hr, scope, nm, y, preds[[nm]], pers))
    out
  }
  preds_all = list(OMNI = Yhat_omni, RTSW_model = Yhat_rtsw,
                   RTSW_affine = Yhat_affine, RTSW_qmap = Yhat_qmap,
                   RTSW_affine_all = Yhat_affine_all)

  Results = rbind(Results,
                  score_named(sprintf("lead_%dh", horizon_hr),
                              y_true[, n_ahead], y_pers[, n_ahead],
                              lapply(preds_all[c("OMNI", "RTSW_model")], `[`, , n_ahead)))
  if (any(after_cal)) {
    Results = rbind(Results,
                    score_named(sprintf("lead_%dh_after_cal", horizon_hr),
                                y_true[after_cal, n_ahead], y_pers[after_cal, n_ahead],
                                lapply(preds_all, function(M) M[after_cal, n_ahead])))
  }

  for (k in seq_len(n_ahead)) {
    Results_lead = rbind(Results_lead,
                         score_named(sprintf("lead_%dh", leads[k]),
                                     y_true[, k], y_pers[, k],
                                     lapply(preds_all[c("OMNI", "RTSW_model")], `[`, , k)))
    if (any(after_cal)) {
      Results_lead = rbind(Results_lead,
                           score_named(sprintf("lead_%dh_after_cal", leads[k]),
                                       y_true[after_cal, k], y_pers[after_cal, k],
                                       lapply(preds_all, function(M) M[after_cal, k])))
    }
  }

  for (input_name in names(preds_all)) {
    yhat_i = preds_all[[input_name]]
    scopes_ext = list(all_test = rep(TRUE, nrow(y_true)))
    if (any(after_cal)) scopes_ext$after_cal = after_cal
    for (scope_name in names(scopes_ext)) {
      if (scope_name == "all_test" && !input_name %in% c("OMNI", "RTSW_model")) next
      keep = scopes_ext[[scope_name]]
      yt = y_true[keep, , drop = FALSE]
      yh = yhat_i[keep, , drop = FALSE]
      yp = y_pers[keep, , drop = FALSE]
      extreme = yt >= 6
      n_ext = sum(extreme)
      if (n_ext == 0) {
        Results_extreme = rbind(Results_extreme, data.frame(
          horizon_hr = horizon_hr, input = input_name, scope = scope_name,
          n_extreme = 0, pct = 0,
          RMSE = NA_real_, MAE = NA_real_, Correlation = NA_real_,
          RMSE_persistence = NA_real_, stringsAsFactors = FALSE))
      } else {
        Results_extreme = rbind(Results_extreme, data.frame(
          horizon_hr = horizon_hr, input = input_name, scope = scope_name,
          n_extreme = n_ext,
          pct = round(100 * n_ext / length(yt), 2),
          RMSE = rmse_(yt[extreme], yh[extreme]),
          MAE = mae_(yt[extreme], yh[extreme]),
          Correlation = cor_(yt[extreme], yh[extreme]),
          RMSE_persistence = rmse_(yt[extreme], yp[extreme]),
          stringsAsFactors = FALSE))
      }
    }
  }

  tau_test = tau[is_test]
  target_start = sapply(seq_len(n_ahead), function(k) kp_time[i_h[is_test] + k])
  target_start = as.POSIXct(as.vector(target_start), origin = "1970-01-01", tz = "UTC")
  issue_long = rep(tau_test, times = n_ahead)

  stopifnot(all(target_start >= issue_long))
  stopifnot(all(target_start + kp_step_hr * 3600
                <= issue_long + horizon_hr * 3600 + 1))

  fmt_utc = function(x) format(x, "%Y-%m-%d %H:%M:%S", tz = "UTC")
  pred_long = data.frame(
    horizon_hr = horizon_hr,
    issue_time = fmt_utc(issue_long),
    target_bin_start = fmt_utc(target_start),
    lead_end_hr = rep(leads, each = length(tau_test)),
    y_true = as.vector(y_true),
    y_pred_omni = as.vector(Yhat_omni),
    y_pred_rtsw = as.vector(Yhat_rtsw),
    y_pred_rtsw_affine = as.vector(Yhat_affine),
    y_pred_rtsw_qmap = as.vector(Yhat_qmap),
    y_persistence = as.vector(y_pers),
    n_filled_hours = rep(n_filled_test, times = n_ahead),
    stringsAsFactors = FALSE
  )
  pred_long = pred_long[order(pred_long$issue_time, pred_long$lead_end_hr), ]

  write.csv(pred_long,
            file.path(out_dir, sprintf("VARX_Kp_h%02d_predictions.csv", horizon_hr)),
            row.names = FALSE)

  result_obj = list(
    horizon_hr = horizon_hr,
    history_hr = history_hr,
    omni_features = omni_features,
    train_end = train_end,
    issue_time = tau_test,
    y_true = y_true,
    y_pred_omni = Yhat_omni,
    y_pred_rtsw = Yhat_rtsw,
    y_persistence = y_pers,
    n_filled_hours = n_filled_test,
    lead_end_hr = leads,
    n_dropped = fit$n_dropped,
    rank = fit$rank,
    n_predictors = ncol(X_train)
  )
  save(result_obj,
       file = file.path(out_dir, sprintf("VARX_Kp_h%02d.RData", horizon_hr)))
}

omni_rows = Results$input == "OMNI"
gap_rows = data.frame()
for (alt in c("RTSW_model", "RTSW_affine", "RTSW_qmap", "RTSW_affine_all")) {
  alt_rows = Results$input == alt
  for (i in which(omni_rows)) {
    key_h = Results$horizon_hr[i]
    key_s = Results$scope[i]
    j = which(alt_rows & Results$horizon_hr == key_h & Results$scope == key_s)
    if (length(j) != 1) next
    gap_rows = rbind(gap_rows, data.frame(
      horizon_hr = key_h, input = paste0(alt, "_minus_OMNI"), scope = key_s,
      n = Results$n[i],
      RMSE = Results$RMSE[j] - Results$RMSE[i],
      MAE = Results$MAE[j] - Results$MAE[i],
      Correlation = Results$Correlation[j] - Results$Correlation[i],
      RMSE_persistence = Results$RMSE_persistence[i],
      MAE_persistence = Results$MAE_persistence[i],
      skill_RMSE = Results$skill_RMSE[j] - Results$skill_RMSE[i],
      stringsAsFactors = FALSE))
  }
}
Results = rbind(Results, gap_rows)

write.csv(Results, file.path(out_dir, "VARX_Kp_metrics.csv"), row.names = FALSE)
write.csv(Results_lead, file.path(out_dir, "VARX_Kp_metrics_by_lead.csv"), row.names = FALSE)
write.csv(Results_extreme, file.path(out_dir, "VARX_Kp_extreme_Kpge6.csv"), row.names = FALSE)

cat("\n=== Same model, OMNI vs raw RTSW vs calibrated RTSW ===\n")
cat("after_cal scores use only issue times whose 24 h driver starts after the calibration window.\n")
cat("RTSW_affine_all is fit on the whole overlap, including the hours it is applied to.\n")
cat("A *_minus_OMNI row is that input's metric minus OMNI. Positive RMSE means worse than OMNI.\n")
print(Results, row.names = FALSE, digits = 4)
cat("\n=== By lead ===\n")
print(Results_lead, row.names = FALSE, digits = 4)
cat("\n=== Points with true Kp >= 6 ===\n")
print(Results_extreme, row.names = FALSE, digits = 4)

# Each forecast is drawn alone against the observed series. Five traces on
# one axis hide the color; a single saturated line does not.
col_obs = "#111111"
line_methods = list(
  list(key = "y_pred_omni", name = "OMNI", col = "#0057FF", input = "OMNI"),
  list(key = "y_pred_rtsw", name = "RTSW raw", col = "#E10600", input = "RTSW_model"),
  list(key = "y_pred_rtsw_affine", name = "RTSW affine", col = "#087F23", input = "RTSW_affine")
)
metric_of = function(h, input, scope, field) {
  hit = Results$horizon_hr == h & Results$input == input & Results$scope == scope
  Results[[field]][which(hit)[1]]
}

png(file.path(out_dir, "VARX_Kp_RTSW_comparison.png"), width = 2100, height = 1120, res = 140)
par(mfrow = c(2, 3), mar = c(3.1, 3.7, 2.7, 0.7), oma = c(1.8, 0.5, 0.4, 0.3),
    mgp = c(2.05, 0.55, 0), tcl = -0.3, las = 1)
for (horizon_hr in horizons_hr) {
  d = read.csv(file.path(out_dir, sprintf("VARX_Kp_h%02d_predictions.csv", horizon_hr)),
               stringsAsFactors = FALSE)
  d = d[d$lead_end_hr == horizon_hr, ]
  d$issue = as.POSIXct(d$issue_time, tz = "UTC")
  d = d[d$issue >= cal_end + history_hr * 3600, ]
  d = d[order(d$issue), ]
  d$day = as.Date(d$issue, tz = "UTC")
  daily = aggregate(d[, c("y_true", "y_pred_omni", "y_pred_rtsw", "y_pred_rtsw_affine")],
                    by = list(day = d$day), FUN = mean)
  scope = sprintf("lead_%dh_after_cal", horizon_hr)
  for (spec in line_methods) {
    plot(daily$day, daily$y_true, type = "l", col = col_obs, lwd = 2.3,
         ylim = c(0, 9), xlab = "", ylab = "Kp",
         main = sprintf("%dh %s    RMSE %.3f",
                        horizon_hr, spec$name,
                        metric_of(horizon_hr, spec$input, scope, "RMSE")))
    lines(daily$day, daily[[spec$key]], col = spec$col, lwd = 2.3)
    legend("topright", c("Observed", spec$name),
           col = c(col_obs, spec$col), lwd = 3.2, bty = "n", cex = 0.85)
  }
}
mtext("Issue date (UTC). Scores use every 3 h issue after the calibration window; lines are the daily mean of that endpoint forecast.",
      side = 1, outer = TRUE, line = 0.35, cex = 0.72)
dev.off()

# Numerical comparison on the post-calibration issues only.
bar_scopes = list(
  list(h = 6, scope = "lead_6h_after_cal", label = "6h endpoint"),
  list(h = 24, scope = "lead_24h_after_cal", label = "24h endpoint")
)
bar_inputs = c("OMNI", "RTSW_model", "RTSW_affine", "RTSW_qmap")
bar_names = c("Persistence", "OMNI", "RTSW raw", "RTSW affine", "RTSW quantile")
bar_cols = c("#444444", "#0057FF", "#E10600", "#087F23", "#6C2BD9")
rmse_mat = sapply(bar_scopes, function(s) {
  pers = metric_of(s$h, "OMNI", s$scope, "RMSE_persistence")
  vals = sapply(bar_inputs, function(inp) metric_of(s$h, inp, s$scope, "RMSE"))
  c(pers, vals)
})
skill_mat = sapply(bar_scopes, function(s) {
  vals = sapply(bar_inputs, function(inp) metric_of(s$h, inp, s$scope, "skill_RMSE"))
  c(0, vals)
})
colnames(rmse_mat) = vapply(bar_scopes, function(s) s$label, "")
colnames(skill_mat) = colnames(rmse_mat)
rownames(rmse_mat) = bar_names
rownames(skill_mat) = bar_names

png(file.path(out_dir, "VARX_Kp_RTSW_metrics.png"), width = 1680, height = 1280, res = 140)
par(mfrow = c(2, 1), mar = c(4.2, 4.6, 2.6, 1.2), oma = c(0.3, 0, 0.2, 0),
    mgp = c(2.4, 0.6, 0), tcl = -0.3, las = 1)
bp = barplot(rmse_mat, beside = TRUE, col = bar_cols, ylim = c(0, max(rmse_mat) * 1.22),
             ylab = "RMSE (Kp)", border = NA,
             main = "Kp forecast RMSE after the calibration window")
text(bp, rmse_mat, labels = sprintf("%.2f", rmse_mat), pos = 3, cex = 0.62, xpd = NA)
legend("topleft", bar_names, fill = bar_cols, border = NA, bty = "n", cex = 0.82, horiz = TRUE)
mtext("Lower is closer to the observed Kp. Persistence repeats the last observed Kp. n = 1635 issue times.",
      side = 1, line = 2.6, cex = 0.72, adj = 0)

sk_pad = 0.06
bp = barplot(skill_mat, beside = TRUE, col = bar_cols, border = NA,
             ylim = c(min(skill_mat) - sk_pad, max(skill_mat) + sk_pad),
             ylab = "skill  =  1 - RMSE / RMSE_persistence",
             main = "Skill against persistence, same issue times")
abline(h = 0, col = "#222222", lwd = 1)
text(as.vector(bp), as.vector(skill_mat),
     labels = sprintf("%+.2f", as.vector(skill_mat)),
     pos = ifelse(as.vector(skill_mat) >= 0, 3, 1), cex = 0.62, xpd = NA)
mtext("Positive skill beats persistence. The gray bar is persistence itself, so its skill is 0.",
      side = 1, line = 2.6, cex = 0.72, adj = 0)
dev.off()
cat(sprintf("\nWrote metrics, predictions, VARX_Kp_RTSW_comparison.png, and VARX_Kp_RTSW_metrics.png under %s\n", out_dir))

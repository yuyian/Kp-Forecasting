# =============================================================================
# Predict the next 24 hours of OMNI solar-wind features, then forecast Kp
# with those features included. Does not change the realtime model.
#
# How tomorrow's OMNI is predicted
# --------------------------------
# Three forecasts of the next 24 hours are combined.
#
#   Linear. One least-squares model maps the previous 24 hours of B, Bz, Vs,
#   Pd and Pt to each feature at each lead. Fit on OMNI hours whose target
#   ends before 2016-01-01.
#
#   Solar rotation. The lag between 25 and 29 days that best matches the last
#   5 days of speed is selected, and the solar wind that followed that window
#   is copied forward as the next 24 hours. This is the part the linear model
#   misses: speed structures recur with the Sun's rotation, and that path is
#   not a linear function of the last 24 hours.
#
#   Stack. On 2016-2018, each feature and lead is a regression on the linear
#   forecast, the rotation path and persistence of the last hour. A shallow
#   boosted model of Bz and speed is fit on the same years and kept only when
#   it beats the stack on 2019. The Kp summary model then sees this forecast
#   only from 2020 onward.
#
# A linear Kp model that already contains those past 24 hours cannot gain
# from appending the linear forecast: the forecast is an exact function of
# columns already in the design. Two uses are not redundant.
#
#   Perfect prog. The Kp model is trained on the true next 24 hours of OMNI,
#   so its coefficients are the map from tomorrow's solar wind to Kp. At
#   issue time those columns are filled with the statistical forecast.
#
#   Forecast summaries. Minimum Bz, mean speed, the fraction of hours with
#   Bz below 0, and a few other nonlinear summaries of the predicted path
#   are added beside the past hours. They are trained on the statistical
#   forecast, so training and the operational test use the same kind of input.
#
# How AWSoM is used where the archive has it
# ------------------------------------------
# The stored AWSoM series is a model estimate of solar wind at Earth, on the
# same timestamps as OMNI. Inside the forecast window those hours are turned
# onto the OMNI scale with OMNI = a + b * AWSoM, fit once on the first 30
# days of overlap and then held fixed. The calibrated hours fill the same
# "tomorrow" columns as the statistical forecast, and their summaries fill
# the summary columns. The live RTSW file ends near the current hour, so
# this is what the Kp model would do if that next-24-hour AWSoM block were
# available at issue time. It is not the operational feed.
#
# Comparison
# ----------
# Scored on issues in VARX_Kp_baseline_endpoint.csv whose Kp bin is in
# Data/Processed/Kp.csv.
#   previous in-situ / boosted : the current realtime VARX
#   past only                  : same state, retrained here, no future solar wind
#   perfect prog, forecast     : trained on true future OMNI, tested on the forecast
#   perfect prog, AWSoM        : same coefficients, tested on calibrated AWSoM
#   summaries                  : past state plus nonlinear summaries of the forecast
#   oracle                     : trained and tested on true future OMNI
#                                (ceiling; not an operational forecast)
# =============================================================================

script_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path = normalizePath(sub("^--file=", "", script_arg[1]))
  repo_root = normalizePath(file.path(dirname(script_path), "..", "..", ".."))
} else {
  repo_root = normalizePath("~/pphd-files/Kp-Forecasting")
}

history_hr = 24L
omni_features = c("B", "Bz", "Vs", "Pd", "Pt")
ols_tol = 1e-10
kp_step_hr = 3L
cal_days = 30L
match_hr = 5L * 24L
lags_hr = seq(25L * 24L, 29L * 24L, by = 12L)
kp_match_bins = 16L
sentinel = c(B = 9999, Bz = 9999, Vs = 9999, Pd = 9999, Pt = 9999999)
omni_fit_end = as.POSIXct("2016-01-01 00:00:00", tz = "UTC")

out_dir = file.path(repo_root, "Rcode", "Main", "OmniLead", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
kp_path = file.path(repo_root, "Data", "Processed", "Kp.csv")
omni_path = file.path(repo_root, "Data", "Processed", "OMNI_hourly.csv")
obs_path = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_obs.csv")
model_path = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_model.csv")
endpoint_path = file.path(repo_root, "Rcode", "Main", "output", "Main_VARX_Realtime",
                          "VARX_Kp_baseline_endpoint.csv")
stopifnot(file.exists(kp_path), file.exists(omni_path), file.exists(obs_path),
          file.exists(model_path), file.exists(endpoint_path))

fmt_utc = function(x) format(x, "%Y-%m-%d %H:%M:%S", tz = "UTC")
parse_utc = function(x) {
  x = trimws(as.character(x))
  date_only = grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  x[date_only] = paste0(x[date_only], " 00:00:00")
  out = as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (anyNA(out)) stop(sprintf("Unparsed timestamps, first: %s", x[which(is.na(out))[1]]))
  out
}
floor_hour = function(t) {
  as.POSIXct(as.numeric(t) - (as.numeric(t) %% 3600), origin = "1970-01-01", tz = "UTC")
}
finite_rows = function(M) rowSums(is.finite(M)) == ncol(M)
rmse = function(y, yhat) {
  m = is.finite(y) & is.finite(yhat)
  if (!any(m)) return(NA_real_)
  sqrt(mean((y[m] - yhat[m])^2))
}
standardize_fit = function(M) {
  ctr = colMeans(M)
  scl = apply(M, 2, sd)
  scl[scl == 0 | is.na(scl)] = 1
  list(center = ctr, scale = scl)
}
apply_scale = function(M, center, scale) sweep(sweep(M, 2, center, "-"), 2, scale, "/")
ols_coef_from_svd = function(X, Y, tol = 1e-10) {
  k = min(nrow(X), ncol(X))
  sv = svd(X, nu = k, nv = k)
  keep = sv$d > tol * max(sv$d)
  Dinv = ifelse(keep, 1 / sv$d, 0)
  UtY = crossprod(sv$u, Y)
  B = t(sv$v %*% (Dinv * UtY))
  list(B = B, n_dropped = sum(!keep), rank = sum(keep))
}
predict_ols = function(M, fit) {
  cbind(1, apply_scale(M, fit$center, fit$scale)) %*% t(fit$B)
}

# Feature-major block, oldest hour to newest, ending at end_idx.
sw_block = function(mat, end_idx) {
  n = length(end_idx)
  p = ncol(mat)
  out = matrix(NA_real_, n, history_hr * p)
  for (f in seq_len(p)) {
    v = mat[, f]
    for (a in seq_len(history_hr))
      out[, (f - 1L) * history_hr + a] = v[end_idx - history_hr + a]
  }
  out
}

cat("Loading Kp and OMNI\n")
Kp_df = read.csv(kp_path, stringsAsFactors = FALSE)
OMNI_df = read.csv(omni_path, stringsAsFactors = FALSE)
Kp_df$Date = parse_utc(Kp_df$Date)
OMNI_df$Date = parse_utc(OMNI_df$Date)
Kp_df = Kp_df[order(Kp_df$Date), ]
OMNI_df = OMNI_df[order(OMNI_df$Date), ]
kp_time = Kp_df$Date
kp_value = Kp_df$Data
omni_time = OMNI_df$Date
OMNI_mat = as.matrix(OMNI_df[, omni_features, drop = FALSE])
tn = as.numeric(omni_time)
stopifnot(median(diff(tn)) == 3600)
n_kp_hist = history_hr / kp_step_hr
N_omni = nrow(OMNI_mat)
N_kp = length(kp_time)

# ---- OMNI feature forecast: previous 24 h -> next 24 h ----------------------
cat("Fitting the OMNI lead model on hours before 2016-01-01\n")
end_idx = seq.int(history_hr, N_omni - history_hr)
contig = (tn[end_idx] - tn[end_idx - history_hr + 1L] == (history_hr - 1L) * 3600) &
  (tn[end_idx + history_hr] - tn[end_idx] == history_hr * 3600)
end_idx = end_idx[contig]
in_fit = omni_time[end_idx + history_hr] < omni_fit_end
end_fit = end_idx[in_fit]
X_om = sw_block(OMNI_mat, end_fit)
Y_om = sw_block(OMNI_mat, end_fit + history_hr)
keep = finite_rows(X_om) & finite_rows(Y_om)
X_om = X_om[keep, , drop = FALSE]
Y_om = Y_om[keep, , drop = FALSE]
cat(sprintf("  OMNI training hours %d\n", nrow(X_om)))
sc_om = standardize_fit(X_om)
fit_om = ols_coef_from_svd(cbind(1, apply_scale(X_om, sc_om$center, sc_om$scale)), Y_om, ols_tol)
omni_model = list(B = fit_om$B, center = sc_om$center, scale = sc_om$scale,
                  rank = fit_om$rank, n_train = nrow(X_om))
cat(sprintf("  rank %d / %d\n", fit_om$rank, ncol(X_om) + 1L))
rm(X_om, Y_om)
predict_tomorrow = function(hist) predict_ols(hist, omni_model)

# Recurrence reads this series. OMNI is the prefix; in-situ hours after the
# last OMNI hour are appended once RTSW is loaded.
sw_time = omni_time
sw_num = tn
sw_mat = OMNI_mat
row_cor = function(A, B) {
  A = A - rowMeans(A); B = B - rowMeans(B)
  num = rowSums(A * B)
  den = sqrt(rowSums(A * A) * rowSums(B * B))
  out = num / den
  out[!is.finite(out)] = 0
  out
}
contiguous_end = function(end_idx, n_back) {
  ok = end_idx >= n_back & end_idx <= length(sw_num)
  out = rep(FALSE, length(end_idx))
  if (!any(ok)) return(out)
  e = end_idx[ok]
  out[ok] = sw_num[e] - sw_num[e - n_back + 1L] == (n_back - 1L) * 3600
  out
}
window_from = function(v, end_idx, n_back) {
  M = matrix(NA_real_, length(end_idx), n_back)
  for (a in seq_len(n_back)) M[, a] = v[end_idx - n_back + a]
  M
}
repeat_last = function(hist) {
  hist[, rep(seq(history_hr, ncol(hist), by = history_hr), each = history_hr), drop = FALSE]
}
# Next 24 h copied from the previous rotation whose last 5 days of speed
# match the origin. Lags with negative correlation are not used.
recur_path = function(end_idx) {
  n = length(end_idx)
  vs = sw_mat[, match("Vs", omni_features)]
  cors = matrix(NA_real_, n, length(lags_hr))
  for (j in seq_along(lags_hr)) {
    L = lags_hr[j]
    past_end = end_idx - L
    fut_end = end_idx - L + history_hr
    ok = contiguous_end(end_idx, match_hr) & contiguous_end(past_end, match_hr) &
      contiguous_end(fut_end, history_hr)
    if (!any(ok)) next
    cj = row_cor(window_from(vs, end_idx[ok], match_hr), window_from(vs, past_end[ok], match_hr))
    cj[cj <= 0] = NA_real_
    cors[ok, j] = cj
  }
  finite_cor = is.finite(cors)
  cors[!finite_cor] = -Inf
  best = max.col(cors, ties.method = "first")
  has = rowSums(finite_cor) > 0
  out = matrix(NA_real_, n, history_hr * length(omni_features))
  for (j in seq_along(lags_hr)) {
    use = has & best == j
    if (!any(use)) next
    out[use, ] = sw_block(sw_mat, end_idx[use] - lags_hr[j] + history_hr)
  }
  out
}
xgb_side = function(hist) {
  bz_cols = history_hr + seq_len(history_hr)
  vs_cols = 2L * history_hr + seq_len(history_hr)
  cbind(
    hist[, bz_cols[history_hr]],
    apply(hist[, bz_cols[(history_hr - 5L):history_hr], drop = FALSE], 1, min),
    hist[, vs_cols[history_hr]],
    rowMeans(hist[, vs_cols, drop = FALSE]),
    hist[, history_hr]
  )
}
fit_stack_col = function(y, lin, rec, per) {
  ok = is.finite(y) & is.finite(lin) & is.finite(rec) & is.finite(per)
  beta = c(0, 1, 0, 0)
  if (sum(ok) < 100) return(beta)
  beta = tryCatch(qr.coef(qr(cbind(1, lin[ok], rec[ok], per[ok])), y[ok]),
                  error = function(e) rep(NA_real_, 4))
  if (any(!is.finite(beta))) beta = c(0, 1, 0, 0)
  beta
}
apply_stack = function(lin, rec, per) {
  n = nrow(lin)
  m = ncol(lin)
  out = matrix(stack_coef[, 1], n, m, byrow = TRUE) +
    lin * matrix(stack_coef[, 2], n, m, byrow = TRUE) +
    rec * matrix(stack_coef[, 3], n, m, byrow = TRUE) +
    per * matrix(stack_coef[, 4], n, m, byrow = TRUE)
  miss = !is.finite(out)
  out[miss] = lin[miss]
  out
}

cat("Fitting the rotation stack on 2016-2018\n")
take_origins = function(t0, t1, step) {
  ix = end_idx[omni_time[end_idx] >= t0 & omni_time[end_idx] < t1]
  if (step > 1L && length(ix) > step) ix = ix[seq(1L, length(ix), by = step)]
  ix
}
build_trio = function(ix) {
  hist = sw_block(sw_mat, ix)
  fut = sw_block(sw_mat, ix + history_hr)
  list(hist = hist, fut = fut, lin = predict_tomorrow(hist),
       rec = recur_path(ix), per = repeat_last(hist))
}
fit_ix = take_origins(as.POSIXct("2016-01-01", tz = "UTC"),
                      as.POSIXct("2019-01-01", tz = "UTC"), 3L)
trio = build_trio(fit_ix)
stack_coef = t(sapply(seq_len(ncol(trio$fut)), function(j) {
  fit_stack_col(trio$fut[, j], trio$lin[, j], trio$rec[, j], trio$per[, j])
}))
cat(sprintf("  stack hours %d\n", nrow(trio$fut)))
cat("  lead-24 weights (intercept, linear, rotation, persistence)\n")
for (f in seq_along(omni_features)) {
  b = stack_coef[f * history_hr, ]
  cat(sprintf("  %-2s  %8.3f  %6.3f  %6.3f  %6.3f\n", omni_features[f], b[1], b[2], b[3], b[4]))
}

use_xgb = rep(FALSE, length(omni_features))
xgb_models = vector("list", ncol(trio$fut))
xgb_ok = requireNamespace("xgboost", quietly = TRUE)
if (xgb_ok) {
  cat("Fitting boosted residuals for Bz and speed\n")
  side = xgb_side(trio$hist)
  xgb_features = match(c("Bz", "Vs"), omni_features)
  for (f in xgb_features) {
    cols = (f - 1L) * history_hr + seq_len(history_hr)
    for (h in seq_len(history_hr)) {
      j = cols[h]
      X = cbind(trio$lin[, j], trio$rec[, j], trio$per[, j], side)
      okj = finite_rows(X) & is.finite(trio$fut[, j])
      if (sum(okj) < 500) next
      xgb_models[[j]] = xgboost::xgb.train(
        params = list(max_depth = 3, eta = 0.08, min_child_weight = 30,
                      subsample = 0.8, objective = "reg:squarederror"),
        data = xgboost::xgb.DMatrix(X[okj, , drop = FALSE], label = trio$fut[okj, j]),
        nrounds = 40, verbose = 0, nthread = 2
      )
    }
  }
  sel_ix = take_origins(as.POSIXct("2019-01-01", tz = "UTC"),
                        as.POSIXct("2020-01-01", tz = "UTC"), 3L)
  sel = build_trio(sel_ix)
  sel_stack = apply_stack(sel$lin, sel$rec, sel$per)
  sel_xgb = sel_stack
  sel_side = xgb_side(sel$hist)
  for (f in xgb_features) {
    cols = (f - 1L) * history_hr + seq_len(history_hr)
    for (h in seq_len(history_hr)) {
      j = cols[h]
      if (is.null(xgb_models[[j]])) next
      X = cbind(sel$lin[, j], sel$rec[, j], sel$per[, j], sel_side)
      okj = finite_rows(X)
      sel_xgb[okj, j] = predict(xgb_models[[j]], X[okj, , drop = FALSE])
    }
    j24 = f * history_hr
    r_stack = rmse(sel$fut[, j24], sel_stack[, j24])
    m = is.finite(sel$fut[, j24]) & is.finite(sel_xgb[, j24])
    r_xgb = rmse(sel$fut[m, j24], sel_xgb[m, j24])
    use_xgb[f] = is.finite(r_xgb) && r_xgb < r_stack
    cat(sprintf("  2019 lead-24 %-2s  stack %.3f  boosted %.3f  %s\n",
                omni_features[f], r_stack, r_xgb,
                if (use_xgb[f]) "keep boosted" else "keep stack"))
  }
  rm(sel, sel_stack, sel_xgb, sel_side, side)
} else {
  cat("  xgboost is not installed; the rotation stack is the forecast\n")
}
rm(trio)
forecast_name = if (any(use_xgb)) "rotation stack + boosted Bz/Vs" else "rotation stack"

forecast_chosen = function(hist, end_idx) {
  lin = predict_tomorrow(hist)
  rec = recur_path(end_idx)
  per = repeat_last(hist)
  out = apply_stack(lin, rec, per)
  if (any(use_xgb)) {
    side = xgb_side(hist)
    for (f in which(use_xgb)) {
      cols = (f - 1L) * history_hr + seq_len(history_hr)
      for (h in seq_len(history_hr)) {
        j = cols[h]
        if (is.null(xgb_models[[j]])) next
        X = cbind(lin[, j], rec[, j], per[, j], side)
        okj = finite_rows(X)
        if (!any(okj)) next
        out[okj, j] = predict(xgb_models[[j]], X[okj, , drop = FALSE])
      }
    }
  }
  miss = !is.finite(out)
  out[miss] = lin[miss]
  out
}

cat("OMNI forecast skill on 2024\n")
chk_ix = take_origins(as.POSIXct("2024-01-01", tz = "UTC"),
                      as.POSIXct("2025-01-01", tz = "UTC"), 3L)
chk = build_trio(chk_ix)
chk_stack = apply_stack(chk$lin, chk$rec, chk$per)
chk_best = forecast_chosen(chk$hist, chk_ix)
omni_skill = do.call(rbind, lapply(seq_along(omni_features), function(f) {
  cols = (f - 1L) * history_hr + seq_len(history_hr)
  do.call(rbind, lapply(seq_along(cols), function(h) {
    j = cols[h]
    m = is.finite(chk$fut[, j]) & is.finite(chk$per[, j]) & is.finite(chk$lin[, j]) &
      is.finite(chk$rec[, j]) & is.finite(chk_stack[, j]) & is.finite(chk_best[, j])
    data.frame(
      feature = omni_features[f], lead_hr = h, n = sum(m),
      RMSE_persistence = rmse(chk$fut[m, j], chk$per[m, j]),
      RMSE_linear = rmse(chk$fut[m, j], chk$lin[m, j]),
      RMSE_rotation = rmse(chk$fut[m, j], chk$rec[m, j]),
      RMSE_stack = rmse(chk$fut[m, j], chk_stack[m, j]),
      RMSE_forecast = rmse(chk$fut[m, j], chk_best[m, j]),
      stringsAsFactors = FALSE
    )
  }))
}))
omni_skill$skill = 1 - omni_skill$RMSE_forecast / omni_skill$RMSE_persistence
write.csv(omni_skill, file.path(out_dir, "omni_forecast_skill_2024.csv"), row.names = FALSE)
for (f in omni_features) {
  sub = omni_skill[omni_skill$feature == f & omni_skill$lead_hr == 24L, ]
  cat(sprintf("  %-2s  24h  pers %.2f  linear %.2f  rotation %.2f  stack %.2f  forecast %.2f  skill %+.3f\n",
              f, sub$RMSE_persistence, sub$RMSE_linear, sub$RMSE_rotation,
              sub$RMSE_stack, sub$RMSE_forecast, sub$skill))
}
rm(chk, chk_stack, chk_best)

# ---- RTSW hourly means (archive already on disk; no download) ---------------
make_rtsw_hourly = function(path) {
  raw = read.csv(path, stringsAsFactors = FALSE)
  raw$Date = parse_utc(raw$Date)
  raw = raw[order(raw$Date), ]
  time = raw$Date
  mat = as.matrix(raw[, omni_features, drop = FALSE])
  for (nm in omni_features)
    mat[!is.finite(mat[, nm]) | mat[, nm] >= sentinel[nm], nm] = NA
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
  list(time = hours, mat = out, n_native = n_native, native_start = min(time))
}
cat("Reducing RTSW obs and model to hourly means\n")
obs_h = make_rtsw_hourly(obs_path)
model_h = make_rtsw_hourly(model_path)
extra_h = obs_h$time > max(sw_time)
if (any(extra_h)) {
  sw_time = c(sw_time, obs_h$time[extra_h])
  sw_num = c(sw_num, as.numeric(obs_h$time[extra_h]))
  sw_mat = rbind(sw_mat, obs_h$mat[extra_h, , drop = FALSE])
}
cat(sprintf("  recurrence series runs through %s\n", fmt_utc(max(sw_time))))
train_end = model_h$native_start
cal_end = train_end + cal_days * 86400
cat(sprintf("  RTSW model starts %s ; calibration ends %s\n", fmt_utc(train_end), fmt_utc(cal_end)))

# Affine map from the first 30 days only. Later hours, including hours inside
# a forecast window, are not used to fit a or b.
pair_omni = match(as.numeric(model_h$time), tn)
model_pair_ok = model_h$n_native > 0 & apply(model_h$mat, 1, function(r) all(is.finite(r)))
omni_pair_ok = !is.na(pair_omni)
omni_pair_ok[omni_pair_ok] = apply(OMNI_mat[pair_omni[omni_pair_ok], , drop = FALSE], 1,
                                   function(r) all(is.finite(r)))
elig = which(model_pair_ok & omni_pair_ok & model_h$time < cal_end)
affine_fit = function(x, y) {
  ok = is.finite(x) & is.finite(y)
  x = x[ok]; y = y[ok]
  if (length(x) < 5 || var(x) == 0) return(c(a = mean(y) - mean(x), b = 1))
  b = cov(x, y) / var(x)
  c(a = mean(y) - b * mean(x), b = b)
}
aw_a = numeric(length(omni_features))
aw_b = aw_a
for (f in seq_along(omni_features)) {
  ab = affine_fit(model_h$mat[elig, f], OMNI_mat[pair_omni[elig], f])
  aw_a[f] = unname(ab["a"])
  aw_b[f] = unname(ab["b"])
}
cat(sprintf("  AWSoM map fit on %d hours through %s\n", length(elig), fmt_utc(cal_end)))

# ---- recurrence, same definition as the realtime VARX -----------------------
row_cor = function(A, B) {
  A = A - rowMeans(A); B = B - rowMeans(B)
  num = rowSums(A * B)
  den = sqrt(rowSums(A * A) * rowSums(B * B))
  out = num / den
  out[!is.finite(out)] = 0
  out
}
pull_hours = function(hour_num, grid, values, latest = NULL) {
  ix = match(as.vector(hour_num), as.numeric(grid))
  bad = is.na(ix)
  if (!is.null(latest)) bad = bad | as.vector(hour_num) > as.numeric(latest)
  ix[bad] = 1L
  vals = values[ix]
  dim(vals) = dim(hour_num)
  vals[bad] = NA_real_
  vals
}
hours_ending = function(t_end, n_back) outer(as.numeric(t_end), (n_back:1L) * 3600, "-")
kp_matched_growth = function(i_kp, kp_src) {
  n = length(i_kp)
  W = kp_match_bins
  n_ahead = 24L / kp_step_hr
  lag_bins = lags_hr / kp_step_hr
  pull_kp = function(ix) {
    good = ix >= 1L & ix <= length(kp_src)
    out = rep(NA_real_, length(ix))
    out[good] = kp_src[ix[good]]
    out
  }
  as_rows = function(M) { if (is.null(dim(M))) M = matrix(M, nrow = n); M }
  recent = as_rows(sapply(seq_len(W), function(j) pull_kp(i_kp - W + j)))
  cors = matrix(0, n, length(lag_bins))
  deltas = array(NA_real_, c(n, length(lag_bins), n_ahead))
  for (j in seq_along(lag_bins)) {
    L = lag_bins[j]
    past = as_rows(sapply(seq_len(W), function(t) pull_kp(i_kp - L - W + t)))
    okj = finite_rows(recent) & finite_rows(past)
    if (any(okj)) {
      cj = row_cor(recent[okj, , drop = FALSE], past[okj, , drop = FALSE])
      cors[okj, j] = pmax(cj, 0)
    }
    for (k in seq_len(n_ahead))
      deltas[, j, k] = pull_kp(i_kp - L + k) - pull_kp(i_kp - L)
    cors[rowSums(is.finite(deltas[, j, , drop = FALSE])) < n_ahead, j] = 0
  }
  sw = rowSums(cors)
  grow = matrix(NA_real_, n, n_ahead)
  for (k in seq_len(n_ahead))
    grow[, k] = rowSums(cors * deltas[, , k], na.rm = TRUE) / sw
  grow[sw <= 0 | !is.finite(sw)] = NA_real_
  grow
}
recur_batch = function(tau, i_kp, kp_src, kp_hist_mat, exog,
                       speed_time, speed_vs, vs_time, vs_values, speed_latest = NULL) {
  n = length(tau)
  recent = pull_hours(hours_ending(tau, match_hr), speed_time, speed_vs, speed_latest)
  cors = matrix(NA_real_, n, length(lags_hr))
  for (j in seq_along(lags_hr)) {
    past = pull_hours(hours_ending(tau - lags_hr[j] * 3600, match_hr),
                      speed_time, speed_vs, speed_latest)
    okj = finite_rows(recent) & finite_rows(past)
    if (any(okj))
      cors[okj, j] = row_cor(recent[okj, , drop = FALSE], past[okj, , drop = FALSE])
  }
  finite_cor = is.finite(cors)
  cors[!finite_cor] = -Inf
  best_j = max.col(cors, ties.method = "first")
  has = rowSums(finite_cor) > 0
  best_lag = lags_hr[best_j]
  rcor = cors[cbind(seq_len(n), best_j)]
  rcor[!has] = NA_real_
  lag_bins = best_lag / kp_step_hr
  rkp = matrix(NA_real_, n, n_kp_hist)
  for (j in seq_len(n_kp_hist)) {
    ix = i_kp - lag_bins - n_kp_hist + j
    good = ix >= 1L & ix <= length(kp_src)
    rkp[good, j] = kp_src[ix[good]]
  }
  raw_vs = pull_hours(hours_ending(tau - best_lag * 3600, history_hr), vs_time, vs_values, speed_latest)
  rvs = rowMeans(raw_vs)
  grow = kp_matched_growth(i_kp, kp_src)
  ok = has & finite_rows(rkp) & is.finite(rvs) & is.finite(rcor) & finite_rows(grow)
  list(ok = ok, extra = cbind(rkp, rvs, rcor, grow))
}

fit_kp = function(X, Y) {
  sc = standardize_fit(X)
  fit = ols_coef_from_svd(cbind(1, apply_scale(X, sc$center, sc$scale)), Y, ols_tol)
  list(B = fit$B, center = sc$center, scale = sc$scale, rank = fit$rank, n_train = nrow(X))
}

# Nonlinear summaries of a 24 h solar-wind block. A linear function of the
# past hours cannot reproduce a minimum or a threshold count.
sw_summary = function(block) {
  sl = function(f) block[, (f - 1L) * history_hr + seq_len(history_hr), drop = FALSE]
  B = sl(1L); Bz = sl(2L); Vs = sl(3L); Pd = sl(4L)
  cbind(
    Bz_min = apply(Bz, 1, min),
    Bz_mean = rowMeans(Bz),
    Bz_neg = rowMeans(Bz < 0),
    Vs_mean = rowMeans(Vs),
    Vs_max = apply(Vs, 1, max),
    Pd_mean = rowMeans(Pd),
    B_mean = rowMeans(B),
    B_max = apply(B, 1, max)
  )
}

# Training issues: Kp target ends before the first RTSW sample, and the OMNI
# lead model is not asked to predict hours it was trained on.
cat("Building Kp training rows from 2016 to the RTSW start\n")
issue_of = function(i) kp_time[i] + kp_step_hr * 3600
n_ahead = 24L / kp_step_hr
i_train = seq.int(n_kp_hist + kp_match_bins + max(lags_hr) / kp_step_hr, N_kp - n_ahead)
tau_train = issue_of(i_train)
keep_t = tau_train >= omni_fit_end + 30 * 86400 & tau_train + 24 * 3600 <= train_end
i_train = i_train[keep_t]
tau_train = tau_train[keep_t]
last_hour = match(as.numeric(tau_train - 3600), tn)
contig_tr = !is.na(last_hour) &
  (tn[last_hour] - tn[last_hour - history_hr + 1L] == (history_hr - 1L) * 3600)
i_train = i_train[contig_tr]
tau_train = tau_train[contig_tr]
last_hour = last_hour[contig_tr]
past_tr = sw_block(OMNI_mat, last_hour)
fut_true = sw_block(OMNI_mat, last_hour + history_hr)
fut_hat = forecast_chosen(past_tr, last_hour)
kp_hist = sapply(seq_len(n_kp_hist), function(j) kp_value[i_train - n_kp_hist + j])
if (is.null(dim(kp_hist))) kp_hist = matrix(kp_hist, ncol = n_kp_hist)
Y = sapply(seq_len(n_ahead), function(k) kp_value[i_train + k])
rb = recur_batch(tau_train, i_train, kp_value, kp_hist, past_tr,
                 omni_time, OMNI_mat[, match("Vs", omni_features)],
                 omni_time, OMNI_mat[, match("Vs", omni_features)])
ok = rb$ok & finite_rows(past_tr) & finite_rows(fut_true) &
  finite_rows(kp_hist) & finite_rows(Y)
ok_sum = ok & finite_rows(fut_hat) & tau_train >= as.POSIXct("2020-01-01", tz = "UTC")
cat(sprintf("  Kp training issues %d ; summary issues from 2020: %d\n", sum(ok), sum(ok_sum)))
sum_hat = sw_summary(fut_hat[ok_sum, , drop = FALSE])
fut_lin_tr = predict_tomorrow(past_tr)
ok_lintr = ok & finite_rows(fut_lin_tr)
fit_past = fit_kp(cbind(kp_hist[ok, ], past_tr[ok, ], rb$extra[ok, ]), Y[ok, ])
fit_sum = fit_kp(cbind(kp_hist[ok_sum, ], past_tr[ok_sum, ], sum_hat, rb$extra[ok_sum, ]), Y[ok_sum, ])
fit_prog = fit_kp(cbind(kp_hist[ok, ], past_tr[ok, ], fut_true[ok, ], rb$extra[ok, ]), Y[ok, ])
# Same future columns as perfect-prog, but filled with the forecast in training
# as well as at issue time. Linear forecasts are a function of the past 24 h,
# so that fit should not beat past-only. The stack forecast adds the rotation
# path. Its weights were fit on 2016-2018, so this Kp fit starts in 2020.
fit_on_linear = fit_kp(cbind(kp_hist[ok_lintr, ], past_tr[ok_lintr, ], fut_lin_tr[ok_lintr, ],
                             rb$extra[ok_lintr, ]), Y[ok_lintr, ])
fit_on_forecast = fit_kp(cbind(kp_hist[ok_sum, ], past_tr[ok_sum, ], fut_hat[ok_sum, ],
                               rb$extra[ok_sum, ]), Y[ok_sum, ])
cat(sprintf("  ranks  past %d   summaries %d   perfect-prog %d\n",
            fit_past$rank, fit_sum$rank, fit_prog$rank))
cat(sprintf("  ranks  trained on linear forecast %d   trained on stack forecast %d\n",
            fit_on_linear$rank, fit_on_forecast$rank))
save(omni_model, stack_coef, use_xgb, fit_past, fit_sum, fit_prog,
     fit_on_linear, fit_on_forecast, aw_a, aw_b,
     file = file.path(out_dir, "OmniLead_coef.RData"))
rm(past_tr, fut_true, fut_hat, kp_hist, Y, rb, sum_hat)

# ---- score the same issues as the published baseline ------------------------
cat("Scoring the published baseline issues\n")
ep = read.csv(endpoint_path, stringsAsFactors = FALSE)
ep$issue = parse_utc(ep$issue_time)
bin_start = ep$issue - kp_step_hr * 3600
i_ep = match(as.numeric(bin_start), as.numeric(kp_time))
in_kp = !is.na(i_ep) & i_ep >= n_kp_hist
cat(sprintf("  %d of %d issues have their Kp bin in Kp.csv\n", sum(in_kp), length(in_kp)))
ep = ep[in_kp, ]
i_ep = i_ep[in_kp]
tau = ep$issue
last_obs_hour = match(as.numeric(tau - 3600), as.numeric(obs_h$time))
past_te = matrix(NA_real_, length(tau), history_hr * length(omni_features))
has_past = !is.na(last_obs_hour) & last_obs_hour >= history_hr
past_te[has_past, ] = sw_block(obs_h$mat, last_obs_hour[has_past])
# Drop a driver if any in-situ hour is missing or after the last real sample.
obs_real = obs_h$n_native > 0
last_real = max(obs_h$time[obs_real])
hour_ok = rep(TRUE, length(tau))
for (i in which(has_past)) {
  hrs = as.numeric(tau[i]) - (history_hr:1L) * 3600
  ix = match(hrs, as.numeric(obs_h$time))
  hour_ok[i] = !anyNA(ix) && all(obs_h$time[ix] <= last_real) && all(is.finite(obs_h$mat[ix, ]))
}
kp_hist_te = sapply(seq_len(n_kp_hist), function(j) kp_value[i_ep - n_kp_hist + j])
rb_te = recur_batch(tau, i_ep, kp_value, kp_hist_te, past_te,
                    obs_h$time, obs_h$mat[, match("Vs", omni_features)],
                    obs_h$time, obs_h$mat[, match("Vs", omni_features)],
                    speed_latest = last_real)
base_ok = hour_ok & rb_te$ok & finite_rows(past_te) & finite_rows(kp_hist_te)

fut_lin = matrix(NA_real_, length(tau), history_hr * length(omni_features))
fut_lin[base_ok, ] = predict_tomorrow(past_te[base_ok, , drop = FALSE])
end_sw = match(as.numeric(tau - 3600), sw_num)
has_sw = base_ok & !is.na(end_sw)
fut_stat = fut_lin
fut_rec = matrix(NA_real_, length(tau), history_hr * length(omni_features))
if (any(has_sw)) {
  fut_rec[has_sw, ] = recur_path(end_sw[has_sw])
  fut_stat[has_sw, ] = forecast_chosen(past_te[has_sw, , drop = FALSE], end_sw[has_sw])
}

# Archived AWSoM over the next 24 hours, mapped onto OMNI.
fut_aw = matrix(NA_real_, length(tau), history_hr * length(omni_features))
aw_ok = rep(FALSE, length(tau))
for (i in which(base_ok)) {
  hrs = as.numeric(tau[i]) + (0:(history_hr - 1L)) * 3600
  ix = match(hrs, as.numeric(model_h$time))
  if (anyNA(ix) || !all(is.finite(model_h$mat[ix, ]))) next
  block = matrix(NA_real_, 1L, history_hr * length(omni_features))
  for (f in seq_along(omni_features)) {
    cols = (f - 1L) * history_hr + seq_len(history_hr)
    block[1, cols] = aw_a[f] + aw_b[f] * model_h$mat[ix, f]
  }
  fut_aw[i, ] = block
  aw_ok[i] = TRUE
}

fut_om = matrix(NA_real_, length(tau), history_hr * length(omni_features))
om_ok = rep(FALSE, length(tau))
for (i in which(base_ok)) {
  hrs = as.numeric(tau[i]) + (0:(history_hr - 1L)) * 3600
  ix = match(hrs, tn)
  if (anyNA(ix) || !all(is.finite(OMNI_mat[ix, ]))) next
  if (tn[ix[history_hr]] - tn[ix[1L]] != (history_hr - 1L) * 3600) next
  fut_om[i, ] = sw_block(OMNI_mat, ix[history_hr])
  om_ok[i] = TRUE
}

y24 = sapply(seq_len(n_ahead), function(k) {
  ix = i_ep + k
  out = rep(NA_real_, length(ix))
  good = ix >= 1L & ix <= N_kp
  out[good] = kp_value[ix[good]]
  out
})
pred_prog = function(fit, fut, use) {
  out = matrix(NA_real_, length(tau), n_ahead)
  if (!any(use)) return(out)
  X = cbind(kp_hist_te[use, , drop = FALSE], past_te[use, , drop = FALSE],
            fut[use, , drop = FALSE], rb_te$extra[use, , drop = FALSE])
  out[use, ] = predict_ols(X, fit)
  out
}
pred_past_only = function(fit, use) {
  out = matrix(NA_real_, length(tau), n_ahead)
  X = cbind(kp_hist_te[use, , drop = FALSE], past_te[use, , drop = FALSE],
            rb_te$extra[use, , drop = FALSE])
  out[use, ] = predict_ols(X, fit)
  out
}
pred_from_summary = function(fit, fut, use) {
  out = matrix(NA_real_, length(tau), n_ahead)
  if (!any(use)) return(out)
  X = cbind(kp_hist_te[use, , drop = FALSE], past_te[use, , drop = FALSE],
            sw_summary(fut[use, , drop = FALSE]), rb_te$extra[use, , drop = FALSE])
  out[use, ] = predict_ols(X, fit)
  out
}

P_past = pred_past_only(fit_past, base_ok)
use_lin = base_ok & finite_rows(fut_lin)
use_stat = base_ok & finite_rows(fut_stat)
P_prog_lin = pred_prog(fit_prog, fut_lin, use_lin)
P_prog_stat = pred_prog(fit_prog, fut_stat, use_stat)
P_on_linear = pred_prog(fit_on_linear, fut_lin, use_lin)
P_on_forecast = pred_prog(fit_on_forecast, fut_stat, use_stat)
P_prog_aw = pred_prog(fit_prog, fut_aw, aw_ok)
P_oracle = pred_prog(fit_prog, fut_om, om_ok)
P_sum_stat = pred_from_summary(fit_sum, fut_stat, use_stat)
P_sum_aw = pred_from_summary(fit_sum, fut_aw, aw_ok)
both = aw_ok & use_stat
fut_blend = fut_stat
fut_blend[both, ] = 0.5 * fut_stat[both, ] + 0.5 * fut_aw[both, ]
P_prog_blend = pred_prog(fit_prog, fut_blend, both)

metric_row = function(name, y, yhat, pers, subset) {
  m = subset & is.finite(yhat) & is.finite(y) & is.finite(pers)
  data.frame(
    model = name, n = sum(m),
    RMSE = if (any(m)) rmse(y[m], yhat[m]) else NA_real_,
    MAE = if (any(m)) mean(abs(y[m] - yhat[m])) else NA_real_,
    Correlation = if (sum(m) > 5) suppressWarnings(cor(y[m], yhat[m])) else NA_real_,
    skill_vs_persistence = if (any(m)) 1 - rmse(y[m], yhat[m]) / rmse(y[m], pers[m]) else NA_real_,
    stringsAsFactors = FALSE
  )
}

# Paired 24 h comparison against the published truth and the published forecasts.
y_pub = ep$y_true_24h
pers = ep$y_persistence
prev_in = ep$y_Earth_obs_24h
prev_aw = ep$y_model_affine_24h
prev_bo = ep$y_boost_24h
common = base_ok & is.finite(y_pub) & is.finite(prev_in)
rows = list(
  metric_row("Previous in-situ VARX", y_pub, prev_in, pers, common),
  metric_row("Previous AWSoM VARX", y_pub, prev_aw, pers, common),
  metric_row("Previous boosted", y_pub, prev_bo, pers, common),
  metric_row("Retrained, past solar wind only", y_pub, P_past[, 8], pers, common),
  metric_row("Perfect prog, linear 24 h forecast", y_pub, P_prog_lin[, 8], pers, common),
  metric_row(sprintf("Perfect prog, %s", forecast_name), y_pub, P_prog_stat[, 8], pers, common),
  metric_row("Trained on the linear forecast", y_pub, P_on_linear[, 8], pers, common),
  metric_row(sprintf("Trained on the %s", forecast_name), y_pub, P_on_forecast[, 8], pers, common),
  metric_row(sprintf("Summaries of the %s", forecast_name), y_pub, P_sum_stat[, 8], pers, common),
  metric_row("Perfect prog, AWSoM next 24 h", y_pub, P_prog_aw[, 8], pers, common & aw_ok),
  metric_row("Perfect prog, blend of forecast and AWSoM", y_pub, P_prog_blend[, 8], pers, common & both),
  metric_row("Summaries of AWSoM next 24 h", y_pub, P_sum_aw[, 8], pers, common & aw_ok),
  metric_row("Previous in-situ, OMNI-overlap issues", y_pub, prev_in, pers, common & om_ok),
  metric_row("Oracle, true OMNI over the next 24 h", y_pub, P_oracle[, 8], pers, common & om_ok)
)
metrics = do.call(rbind, rows)
write.csv(metrics, file.path(out_dir, "kp_24h_comparison.csv"), row.names = FALSE)

rmse_on = function(y, yhat, m) {
  m = m & is.finite(y) & is.finite(yhat)
  if (!any(m)) return(c(n = 0, RMSE = NA_real_))
  c(n = sum(m), RMSE = rmse(y[m], yhat[m]))
}
lead_rows = list()
for (k in seq_len(n_ahead)) {
  yk = y24[, k]
  a = rmse_on(yk, P_past[, k], base_ok)
  b = rmse_on(yk, P_prog_stat[, k], use_stat)
  c = rmse_on(yk, P_sum_stat[, k], use_stat)
  d = rmse_on(yk, P_prog_aw[, k], aw_ok)
  e = rmse_on(yk, P_on_forecast[, k], use_stat)
  lead_rows[[length(lead_rows) + 1L]] = data.frame(
    lead_hr = k * kp_step_hr,
    past_only = a["RMSE"], predicted_OMNI = b["RMSE"],
    trained_on_forecast = e["RMSE"],
    summaries = c["RMSE"], AWSoM_horizon = d["RMSE"],
    n_past = a["n"], n_awsom = d["n"],
    stringsAsFactors = FALSE
  )
}
by_lead = do.call(rbind, lead_rows)
write.csv(by_lead, file.path(out_dir, "kp_by_lead.csv"), row.names = FALSE)

detail = data.frame(
  issue_time = fmt_utc(tau),
  y_true_24h = y_pub,
  y_previous_insitu = prev_in,
  y_previous_boosted = prev_bo,
  y_past_only = P_past[, 8],
  y_prog_linear = P_prog_lin[, 8],
  y_prog_predicted = P_prog_stat[, 8],
  y_prog_awsom = P_prog_aw[, 8],
  y_prog_blend = P_prog_blend[, 8],
  y_summary_predicted = P_sum_stat[, 8],
  y_summary_awsom = P_sum_aw[, 8],
  y_oracle = P_oracle[, 8],
  stringsAsFactors = FALSE
)
write.csv(detail, file.path(out_dir, "kp_24h_predictions.csv"), row.names = FALSE)

both_truth = om_ok & aw_ok & finite_rows(fut_stat) & finite_rows(fut_om)
pers_block = past_te
for (f in seq_along(omni_features)) {
  cols = (f - 1L) * history_hr + seq_len(history_hr)
  pers_block[, cols] = past_te[, cols[history_hr]]
}
block_skill = function(name, block) {
  rows = lapply(seq_along(omni_features), function(f) {
    cols = (f - 1L) * history_hr + seq_len(history_hr)
    j24 = cols[history_hr]
    err = fut_om[, j24] - block[, j24]
    err_all = fut_om[, cols] - block[, cols]
    data.frame(
      source = name, feature = omni_features[f],
      n = sum(both_truth & is.finite(err)),
      RMSE_lead24 = sqrt(mean(err[both_truth]^2, na.rm = TRUE)),
      RMSE_1to24 = sqrt(mean(err_all[both_truth, ]^2, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
sw_vs_omni = rbind(
  block_skill("Persistence of the last hour", pers_block),
  block_skill("Linear 24 h forecast", fut_lin),
  block_skill("Solar-rotation path", fut_rec),
  block_skill(forecast_name, fut_stat),
  block_skill("Calibrated AWSoM", fut_aw)
)
write.csv(sw_vs_omni, file.path(out_dir, "tomorrow_block_vs_omni.csv"), row.names = FALSE)

cat("\n24 h endpoint, paired with the published baseline\n")
for (i in seq_len(nrow(metrics))) {
  cat(sprintf("  %-54s  n %4d  RMSE %.3f  corr %.3f  skill %+.3f\n",
              metrics$model[i], metrics$n[i], metrics$RMSE[i],
              metrics$Correlation[i], metrics$skill_vs_persistence[i]))
}
cat("\nNext-24h solar wind against OMNI, on issues where OMNI and AWSoM both exist\n")
for (i in seq_len(nrow(sw_vs_omni))) {
  if (sw_vs_omni$feature[i] %in% c("Bz", "Vs"))
    cat(sprintf("  %-32s %-2s  lead-24 RMSE %8.2f   all-leads %8.2f\n",
                sw_vs_omni$source[i], sw_vs_omni$feature[i],
                sw_vs_omni$RMSE_lead24[i], sw_vs_omni$RMSE_1to24[i]))
}
cat("\nRMSE by lead (past / perfect-prog / trained on forecast / summaries / AWSoM)\n")
for (i in seq_len(nrow(by_lead))) {
  cat(sprintf("  %2dh   %.3f   %.3f   %.3f   %.3f   %.3f\n",
              by_lead$lead_hr[i], by_lead$past_only[i], by_lead$predicted_OMNI[i],
              by_lead$trained_on_forecast[i], by_lead$summaries[i], by_lead$AWSoM_horizon[i]))
}
cat(sprintf("\nWrote %s\n", out_dir))

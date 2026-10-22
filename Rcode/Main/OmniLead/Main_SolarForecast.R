# =============================================================================
# A solar-wind forecast that can beat the 24 h linear model.
#
# The 24 h model already has full rank on the last 24 hours, so another
# function of those same hours cannot help. Two blocks sit outside that window:
#
#   Daily means of |B|, Bz, speed and density over the six days before the
#   last 24 hours. Speed and sector structure live on that scale.
#
#   Rotation increment. Among lags of 25-29 days, pick the one whose last five
#   days of speed match the origin (negative correlations are dropped). The
#   predictor at lead h is the change that rotation made over h hours,
#   value(t+h-L) - value(t-L), for |B|, Bz and speed. Copying the level from
#   the previous rotation loses to persistence inside one day; the change is
#   the piece that can still be useful.
#
# Both models are least squares, fit on OMNI hours whose target ends before
# 2016-01-01, and scored on calendar 2024. The linear 24 h model is refit on
# the same hours so the only difference at test time is the extra columns.
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
daily_features = c("B", "Bz", "Vs", "Pd")
rot_features = c("B", "Bz", "Vs")
ols_tol = 1e-10
match_hr = 5L * 24L
lags_hr = seq(25L * 24L, 29L * 24L, by = 12L)
daily_days = 2:7
report_leads = c(3L, 6L, 9L, 12L, 24L)
omni_fit_end = as.POSIXct("2016-01-01 00:00:00", tz = "UTC")

out_dir = file.path(repo_root, "Rcode", "Main", "OmniLead", "results")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
omni_path = file.path(repo_root, "Data", "Processed", "OMNI_hourly.csv")
stopifnot(file.exists(omni_path))

parse_utc = function(x) {
  x = trimws(as.character(x))
  date_only = grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  x[date_only] = paste0(x[date_only], " 00:00:00")
  out = as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (anyNA(out)) stop(sprintf("Unparsed timestamps, first: %s", x[which(is.na(out))[1]]))
  out
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
fit_ols = function(X, Y) {
  sc = standardize_fit(X)
  fit = ols_coef_from_svd(cbind(1, apply_scale(X, sc$center, sc$scale)), Y, ols_tol)
  list(B = fit$B, center = sc$center, scale = sc$scale, rank = fit$rank, n_train = nrow(X))
}

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
row_cor = function(A, B) {
  A = A - rowMeans(A); B = B - rowMeans(B)
  num = rowSums(A * B)
  den = sqrt(rowSums(A * A) * rowSums(B * B))
  out = num / den
  out[!is.finite(out)] = 0
  out
}
contiguous_end = function(end_idx, n_back, n_time) {
  ok = end_idx >= n_back & end_idx <= n_time
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

cat("Loading OMNI\n")
OMNI_df = read.csv(omni_path, stringsAsFactors = FALSE)
OMNI_df$Date = parse_utc(OMNI_df$Date)
OMNI_df = OMNI_df[order(OMNI_df$Date), ]
sw_time = OMNI_df$Date
sw_num = as.numeric(sw_time)
sw_mat = as.matrix(OMNI_df[, omni_features, drop = FALSE])
stopifnot(median(diff(sw_num)) == 3600)
N = nrow(sw_mat)
feat_col = match(omni_features, omni_features)
daily_col = match(daily_features, omni_features)
rot_col = match(rot_features, omni_features)

# Origins with a full past week and a full next day on the clock.
week_hr = 7L * 24L
need = max(lags_hr) + match_hr
end_idx = seq.int(need, N - history_hr)
ok_span = contiguous_end(end_idx, week_hr, N) &
  contiguous_end(end_idx + history_hr, history_hr, N)
end_idx = end_idx[ok_span]
cat(sprintf("  origins with a full week and a full next day: %d\n", length(end_idx)))

choose_lag = function(end_idx) {
  n = length(end_idx)
  vs = sw_mat[, match("Vs", omni_features)]
  cors = matrix(NA_real_, n, length(lags_hr))
  for (j in seq_along(lags_hr)) {
    L = lags_hr[j]
    past_end = end_idx - L
    fut_end = end_idx - L + history_hr
    ok = contiguous_end(end_idx, match_hr, N) &
      contiguous_end(past_end, match_hr, N) &
      contiguous_end(fut_end, history_hr, N)
    if (!any(ok)) next
    cj = row_cor(window_from(vs, end_idx[ok], match_hr),
                 window_from(vs, past_end[ok], match_hr))
    cj[cj <= 0] = NA_real_
    cors[ok, j] = cj
  }
  finite_cor = is.finite(cors)
  cors[!finite_cor] = -Inf
  best = max.col(cors, ties.method = "first")
  has = rowSums(finite_cor) > 0
  lag = rep(NA_real_, n)
  rcor = rep(NA_real_, n)
  lag[has] = lags_hr[best[has]]
  rcor[has] = cors[cbind(which(has), best[has])]
  list(lag = lag, cor = rcor, ok = has)
}

daily_means = function(end_idx) {
  out = matrix(NA_real_, length(end_idx), length(daily_days) * length(daily_col))
  col = 1L
  for (f in daily_col) {
    v = sw_mat[, f]
    for (d in daily_days) {
      dest = end_idx - (d - 1L) * 24L
      acc = rep(0, length(end_idx))
      for (a in 0:23) acc = acc + v[dest - a]
      out[, col] = acc / 24
      col = col + 1L
    }
  }
  out
}

# Change over the next h hours on the chosen previous rotation.
rotation_increment = function(end_idx, lag) {
  n = length(end_idx)
  out = matrix(NA_real_, n, history_hr * length(rot_col))
  has = is.finite(lag)
  if (!any(has)) return(out)
  origin = end_idx[has] - lag[has]
  for (i in seq_along(rot_col)) {
    v = sw_mat[, rot_col[i]]
    base = v[origin]
    for (h in seq_len(history_hr)) {
      out[has, (i - 1L) * history_hr + h] = v[origin + h] - base
    }
  }
  out
}

cat("Choosing the solar-rotation lag\n")
lag_info = choose_lag(end_idx)
cat(sprintf("  origins with a positive speed match: %d\n", sum(lag_info$ok)))

in_fit = sw_time[end_idx + history_hr] < omni_fit_end & lag_info$ok
in_chk = sw_time[end_idx] >= as.POSIXct("2024-01-01", tz = "UTC") &
  sw_time[end_idx] < as.POSIXct("2025-01-01", tz = "UTC") & lag_info$ok
fit_ix = end_idx[in_fit]
chk_ix = end_idx[in_chk]
fit_lag = lag_info$lag[in_fit]
chk_lag = lag_info$lag[in_chk]
fit_cor = lag_info$cor[in_fit]
chk_cor = lag_info$cor[in_chk]
if (length(fit_ix) > 80000L) {
  keep = seq(1L, length(fit_ix), by = 2L)
  fit_ix = fit_ix[keep]
  fit_lag = fit_lag[keep]
  fit_cor = fit_cor[keep]
}
if (length(chk_ix) > 4000L) {
  keep = seq(1L, length(chk_ix), by = 3L)
  chk_ix = chk_ix[keep]
  chk_lag = chk_lag[keep]
  chk_cor = chk_cor[keep]
}

build_X = function(ix, lag, rcor, with_rotation) {
  X = cbind(sw_block(sw_mat, ix), daily_means(ix))
  if (with_rotation) X = cbind(X, rotation_increment(ix, lag), rcor)
  X
}

cat("Building training matrices\n")
Y_fit = sw_block(sw_mat, fit_ix + history_hr)
X_lin = sw_block(sw_mat, fit_ix)
X_day = build_X(fit_ix, fit_lag, fit_cor, FALSE)
X_ext = build_X(fit_ix, fit_lag, fit_cor, TRUE)
keep = finite_rows(Y_fit) & finite_rows(X_lin) & finite_rows(X_day) & finite_rows(X_ext)
Y_fit = Y_fit[keep, , drop = FALSE]
X_lin = X_lin[keep, , drop = FALSE]
X_day = X_day[keep, , drop = FALSE]
X_ext = X_ext[keep, , drop = FALSE]
cat(sprintf("  training hours %d\n", nrow(Y_fit)))

cat("Fitting the 24 h model, the daily-mean model, and the rotation-increment model\n")
fit_lin = fit_ols(X_lin, Y_fit)
fit_day = fit_ols(X_day, Y_fit)
fit_ext = fit_ols(X_ext, Y_fit)
cat(sprintf("  rank  24h %d / %d    daily %d / %d    rotation %d / %d\n",
            fit_lin$rank, ncol(X_lin) + 1L,
            fit_day$rank, ncol(X_day) + 1L,
            fit_ext$rank, ncol(X_ext) + 1L))
rm(Y_fit, X_lin, X_day, X_ext)

# Where the lead-24 speed and Bz equations put their weight.
coef_share = function(fit, target_index, n_hourly, n_daily, n_rot) {
  b = fit$B[target_index, -1L]
  b = abs(b)
  s = sum(b)
  if (!is.finite(s) || s == 0) return(c(hourly = NA, daily = NA, rotation = NA))
  c(hourly = sum(b[seq_len(n_hourly)]) / s,
    daily = sum(b[n_hourly + seq_len(n_daily)]) / s,
    rotation = sum(b[n_hourly + n_daily + seq_len(n_rot + 1L)]) / s)
}
n_hourly = history_hr * length(omni_features)
n_daily = length(daily_days) * length(daily_col)
n_rot = history_hr * length(rot_col)
for (f in c("Bz", "Vs")) {
  j = match(f, omni_features) * history_hr
  sh = coef_share(fit_ext, j, n_hourly, n_daily, n_rot)
  cat(sprintf("  lead-24 %-2s weight share  hourly %.2f  daily %.2f  rotation %.2f\n",
              f, sh["hourly"], sh["daily"], sh["rotation"]))
}

cat("Scoring calendar 2024\n")
Y = sw_block(sw_mat, chk_ix + history_hr)
P_lin = predict_ols(sw_block(sw_mat, chk_ix), fit_lin)
P_day = predict_ols(build_X(chk_ix, chk_lag, chk_cor, FALSE), fit_day)
P_ext = predict_ols(build_X(chk_ix, chk_lag, chk_cor, TRUE), fit_ext)
last_rep = sw_block(sw_mat, chk_ix)
for (f in seq_along(omni_features)) {
  cols = (f - 1L) * history_hr + seq_len(history_hr)
  last_rep[, cols] = last_rep[, cols[history_hr]]
}

skill = do.call(rbind, lapply(seq_along(omni_features), function(f) {
  cols = (f - 1L) * history_hr + seq_len(history_hr)
  do.call(rbind, lapply(seq_along(cols), function(h) {
    j = cols[h]
    m = is.finite(Y[, j]) & is.finite(last_rep[, j]) & is.finite(P_lin[, j]) &
      is.finite(P_day[, j]) & is.finite(P_ext[, j])
    data.frame(
      feature = omni_features[f], lead_hr = h, n = sum(m),
      RMSE_persistence = rmse(Y[m, j], last_rep[m, j]),
      RMSE_linear24 = rmse(Y[m, j], P_lin[m, j]),
      RMSE_daily = rmse(Y[m, j], P_day[m, j]),
      RMSE_rotation = rmse(Y[m, j], P_ext[m, j]),
      stringsAsFactors = FALSE
    )
  }))
}))
skill$skill_daily = 1 - skill$RMSE_daily / skill$RMSE_linear24
skill$skill_rotation = 1 - skill$RMSE_rotation / skill$RMSE_linear24
write.csv(skill, file.path(out_dir, "solar_forecast_extended_2024.csv"), row.names = FALSE)
save(fit_lin, fit_day, fit_ext, file = file.path(out_dir, "SolarForecast_coef.RData"))

cat("\n2024 RMSE at 3, 6, 9, 12 and 24 h\n")
cat(sprintf("  %-4s %5s %8s %8s %8s %8s %8s\n",
            "feat", "lead", "persist", "linear", "daily", "rot+day", "vs lin"))
for (f in omni_features) {
  for (h in report_leads) {
    sub = skill[skill$feature == f & skill$lead_hr == h, ]
    cat(sprintf("  %-4s %5d %8.2f %8.2f %8.2f %8.2f %+8.3f\n",
                f, h, sub$RMSE_persistence, sub$RMSE_linear24,
                sub$RMSE_daily, sub$RMSE_rotation, sub$skill_rotation))
  }
}
cat(sprintf("\nWrote %s\n", file.path(out_dir, "solar_forecast_extended_2024.csv")))

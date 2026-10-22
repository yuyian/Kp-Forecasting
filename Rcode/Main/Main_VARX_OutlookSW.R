# Compare the 27-day Kp outlook with solar wind that has already ended.
#
# The issued outlook is not changed. At issue time the next 27 days of solar
# wind have not reached Earth. The RTSW AWSoM file is a nowcast at Earth
# time, so hours after the issue are not used. Two blocks are added, both
# known at issue time:
#   * means of B, Bz, Vs, Pd, Pt, and southward Bz over the previous 24 hours
#   * the previous rotation's daily mean Vs, Bz, and southward Bz
#     (mean of max(-Bz, 0)), on the same 27 days as the Kp analogue
# A day is used when at least 12 of its 24 hours have a real sample.
# Requiring every hour drops almost every RTSW issue, because both series
# have occasional gaps.
# Coefficients are fit on OMNI, on issues whose day-27 target ends before
# the RTSW model starts. The same coefficients are then filled with RTSW
# in-situ, with AWSoM calibrated to OMNI, or with OMNI itself.

script_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path = normalizePath(sub("^--file=", "", script_arg[1]))
  repo_root = normalizePath(file.path(dirname(script_path), "..", ".."))
} else {
  repo_root = normalizePath("~/pphd-files/Kp-Forecasting")
}

omni_features = c("B", "Bz", "Vs", "Pd", "Pt")
sentinel = c(B = 9999, Bz = 9999, Vs = 9999, Pd = 9999, Pt = 9999999)
ols_tol = 1e-10
kp_step_hr = 3L
history_hr = 24L
cal_days = 30L
outlook_days = 27L
outlook_lags_hr = seq(25L * 24L, 29L * 24L, by = kp_step_hr)
outlook_match_bins = 5L * (24L / kp_step_hr)
bins_per_day = 24L / kp_step_hr
n_kp_hist = history_hr / kp_step_hr
min_hours = 12L
recent_names = c(omni_features, "Bs")
rot_names = c("Vs", "Bz", "Bs")

kp_path = file.path(repo_root, "Data", "Processed", "Kp.csv")
omni_path = file.path(repo_root, "Data", "Processed", "OMNI_hourly.csv")
obs_path = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_obs.csv")
model_path = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_model.csv")
out_dir = file.path(repo_root, "Rcode", "Main", "output", "Main_VARX_Realtime")
coef_path = file.path(out_dir, "VARX_Kp_outlook27_coef.RData")

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
finite_rows = function(M) {
  if (!length(M)) return(logical())
  if (is.null(dim(M))) M = matrix(M, nrow = 1)
  rowSums(is.finite(M)) == ncol(M)
}
row_cor = function(A, B) {
  A = A - rowMeans(A)
  B = B - rowMeans(B)
  num = rowSums(A * B)
  den = sqrt(rowSums(A * A) * rowSums(B * B))
  out = num / den
  out[!is.finite(out)] = 0
  out
}
standardize_fit = function(M) {
  ctr = colMeans(M)
  scl = apply(M, 2, sd)
  scl[scl == 0 | is.na(scl)] = 1
  list(center = ctr, scale = scl)
}
apply_scale = function(M, center, scale) {
  if (is.null(dim(M))) M = matrix(M, nrow = 1)
  sweep(sweep(M, 2, center, "-"), 2, scale, "/")
}
ols_coef_from_svd = function(X, Y, tol = 1e-10) {
  k = min(nrow(X), ncol(X))
  sv = svd(X, nu = k, nv = k)
  keep = sv$d > tol * max(sv$d)
  Dinv = ifelse(keep, 1 / sv$d, 0)
  UtY = crossprod(sv$u, Y)
  B = t(sv$v %*% (Dinv * UtY))
  list(B = B, n_dropped = sum(!keep), rank = sum(keep))
}
predict_y = function(X, center, scale, B) {
  cbind(1, apply_scale(X, center, scale)) %*% t(B)
}

message("Loading Kp, OMNI, and the issued 27-day coefficients")
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
for (nm in omni_features)
  OMNI_mat[!is.finite(OMNI_mat[, nm]) | OMNI_mat[, nm] >= sentinel[nm], nm] = NA
env = new.env()
load(coef_path, envir = env)
fit_kp = env$fit_outlook
train_end = as.POSIXct(env$outlook_sig$train_end, origin = "1970-01-01", tz = "UTC")
stopifnot(length(fit_kp$center) == n_kp_hist + outlook_days + 1L)
stopifnot(nrow(fit_kp$B) == outlook_days, ncol(fit_kp$B) == length(fit_kp$center) + 1L)

# Same extension as the live script: GFZ if it answers, otherwise NOAA.
# Existing file bins are not replaced by NOAA, and the file is not written.
extend_kp = function(kp_time, kp_value) {
  added = 0L
  source = "file"
  tmp = tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  gfz_ok = identical(system2("curl", c(
    "--ipv4", "--max-time", "20", "-fsSL", "-A", "Mozilla/5.0",
    "https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_nowcast.txt", "-o", tmp)), 0L)
  if (gfz_ok) {
    nw_lines = readLines(tmp, warn = FALSE)
    nw_lines = nw_lines[!grepl("^\\s*#", nw_lines) & nzchar(trimws(nw_lines))]
    nw = read.table(text = nw_lines, stringsAsFactors = FALSE)
    if (ncol(nw) >= 28) {
      kp8 = as.matrix(nw[, 8:15, drop = FALSE])
      hour_head = seq(0, 21, by = kp_step_hr)
      for (r in seq_len(nrow(nw))) {
        if (nw[r, 28] != 0) next
        for (k in seq_along(hour_head)) {
          val = kp8[r, k]
          if (!is.finite(val) || val < 0) next
          start = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00",
                                     nw[r, 1], nw[r, 2], nw[r, 3], hour_head[k]),
                             format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
          if (start + kp_step_hr * 3600 > Sys.time()) next
          j = match(as.numeric(start), as.numeric(kp_time))
          if (!is.na(j)) {
            kp_value[j] = val
          } else if (length(kp_time) &&
                     abs(as.numeric(start - kp_time[length(kp_time)], units = "secs") -
                         kp_step_hr * 3600) < 1) {
            kp_time = c(kp_time, start)
            kp_value = c(kp_value, val)
            added = added + 1L
          }
        }
      }
      source = "GFZ"
    }
  }
  if (source == "file") {
    noaa_tmp = tempfile(fileext = ".json")
    on.exit(unlink(noaa_tmp), add = TRUE)
    noaa_ok = identical(system2("curl", c(
      "--max-time", "25", "-fsSL", "-A", "Mozilla/5.0",
      "https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json",
      "-o", noaa_tmp)), 0L)
    if (noaa_ok) {
      noaa = tryCatch(jsonlite::fromJSON(noaa_tmp), error = function(e) NULL)
      if (!is.null(noaa) && !is.null(noaa$time_tag)) {
        noaa_start = parse_utc(gsub("T", " ", noaa$time_tag, fixed = TRUE))
        noaa_kp = round(as.numeric(noaa$Kp) * 3) / 3
        ord = order(noaa_start)
        noaa_start = noaa_start[ord]
        noaa_kp = noaa_kp[ord]
        for (k in seq_along(noaa_start)) {
          start = noaa_start[k]
          val = noaa_kp[k]
          if (!is.finite(val) || val < 0 || val > 9) next
          if (start + kp_step_hr * 3600 > Sys.time()) next
          if (!is.na(match(as.numeric(start), as.numeric(kp_time)))) next
          if (length(kp_time) &&
              abs(as.numeric(start - kp_time[length(kp_time)], units = "secs") -
                  kp_step_hr * 3600) < 1) {
            kp_time = c(kp_time, start)
            kp_value = c(kp_value, val)
            added = added + 1L
          }
        }
        if (added > 0L) source = "NOAA"
      }
    }
  }
  message(sprintf("Kp extension %s, appended %d bins, last %s",
                  source, added, fmt_utc(max(kp_time))))
  list(time = kp_time, value = kp_value, source = source)
}
live = extend_kp(kp_time, kp_value)
kp_time = live$time
kp_value = live$value
N_kp = length(kp_time)
N_file = nrow(Kp_df)

make_rtsw_hourly = function(path) {
  raw = read.csv(path, stringsAsFactors = FALSE)
  raw$Date = parse_utc(raw$Date)
  raw = raw[order(raw$Date), ]
  stopifnot(all(omni_features %in% names(raw)))
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
  real = n_native > 0 & apply(out, 1, function(r) all(is.finite(r)))
  list(time = hours, mat = out, ok_hour = real,
       native_start = min(time), native_end = max(hours[real]))
}
message("Reducing RTSW obs and model to hourly means")
obs_h = make_rtsw_hourly(obs_path)
model_h = make_rtsw_hourly(model_path)
message(sprintf("In-situ hours with a real sample end at %s", fmt_utc(obs_h$native_end)))
message(sprintf("AWSoM hours with a real sample end at %s", fmt_utc(model_h$native_end)))

hour_pos = function(origin, n_hours, hour_num) {
  ix = (as.numeric(hour_num) - origin) / 3600 + 1
  r = round(ix)
  bad = !is.finite(ix) | abs(ix - r) > 1e-3 | r < 1 | r > n_hours
  out = as.integer(r)
  out[bad] = NA_integer_
  out
}
require_regular = function(time, label) {
  d = diff(as.numeric(time))
  if (any(abs(d - 3600) >= 1))
    stop(label, " is not a continuous hourly grid.")
  invisible(TRUE)
}
require_regular(omni_time, "OMNI")
require_regular(obs_h$time, "In-situ hourly")
require_regular(model_h$time, "AWSoM hourly")

day_means = function(kp, i, end_off, n_days) {
  n = length(i)
  out = matrix(NA_real_, n, n_days)
  if (!n) return(out)
  end_off = rep(end_off, length.out = n)
  for (d in seq_len(n_days)) {
    acc = rep(0, n)
    last = i + end_off + d * bins_per_day
    good = last >= bins_per_day & last <= length(kp)
    acc[!good] = NA_real_
    for (a in 0:(bins_per_day - 1L)) acc[good] = acc[good] + kp[last[good] - a]
    out[, d] = acc / bins_per_day
  }
  out
}
choose_kp_lag = function(kp, i) {
  n = length(i)
  W = outlook_match_bins
  recent = matrix(NA_real_, n, W)
  for (a in seq_len(W)) recent[, a] = kp[i - W + a]
  cors = matrix(NA_real_, n, length(outlook_lags_hr))
  for (j in seq_along(outlook_lags_hr)) {
    L = outlook_lags_hr[j] / kp_step_hr
    past = matrix(NA_real_, n, W)
    for (a in seq_len(W)) past[, a] = kp[i - L - W + a]
    okj = finite_rows(recent) & finite_rows(past)
    if (any(okj)) cors[okj, j] = row_cor(recent[okj, , drop = FALSE], past[okj, , drop = FALSE])
  }
  cors[, outlook_lags_hr < outlook_days * 24] = NA_real_
  finite_cor = is.finite(cors)
  cors[!finite_cor] = -Inf
  best = max.col(cors, ties.method = "first")
  has = rowSums(finite_cor) > 0
  lag = rep(NA_real_, n)
  rcor = rep(NA_real_, n)
  lag[has] = outlook_lags_hr[best[has]]
  rcor[has] = pmax(cors[cbind(which(has), best[has])], -1)
  list(lag = lag, cor = rcor, ok = has)
}
kp_block = function(kp, i) {
  info = choose_kp_lag(kp, i)
  recent = matrix(NA_real_, length(i), n_kp_hist)
  for (a in seq_len(n_kp_hist)) recent[, a] = kp[i - n_kp_hist + a]
  recur = matrix(NA_real_, length(i), outlook_days)
  if (any(info$ok))
    recur[info$ok, ] = day_means(kp, i[info$ok], -info$lag[info$ok] / kp_step_hr, outlook_days)
  Y = day_means(kp, i, 0, outlook_days)
  pers = day_means(kp, i, -bins_per_day, 1L)[, 1]
  ok = info$ok & finite_rows(recent) & finite_rows(recur) & is.finite(info$cor)
  list(X = cbind(recent, recur, info$cor), Y = Y, ok = ok,
       lag = info$lag, recur = recur, pers = pers)
}

# Means over the hours that actually have a sample. A window is kept when
# at least min_hours of the 24 are present. Recent columns are B, Bz, Vs,
# Pd, Pt, then southward Bz. Each rotation day is mean Vs, mean Bz, mean
# southward Bz, for the same clock hours as that day's Kp analogue.
hour_mask = function(base, n, n_hours, ok_hour, h) {
  ix = base + h
  known = !is.na(ix)
  good = known
  good[known] = ok_hour[ix[known]]
  good[is.na(good)] = FALSE
  good
}
solar_blocks = function(tau, lag_hr, mat, ok_hour, origin) {
  n = length(tau)
  n_hours = nrow(mat)
  vs_i = match("Vs", colnames(mat))
  bz_i = match("Bz", colnames(mat))
  recent = matrix(NA_real_, n, length(recent_names))
  base = hour_pos(origin, n_hours, as.numeric(tau) - history_hr * 3600)
  sum_r = matrix(0, n, length(omni_features))
  sum_bs = rep(0, n)
  n_r = rep(0L, n)
  for (h in 0:23) {
    good = hour_mask(base, n, n_hours, ok_hour, h)
    use = which(good)
    if (!length(use)) next
    n_r[use] = n_r[use] + 1L
    for (f in seq_along(omni_features))
      sum_r[use, f] = sum_r[use, f] + mat[ix_of(base, h)[use], f]
    sum_bs[use] = sum_bs[use] + pmax(-mat[ix_of(base, h)[use], bz_i], 0)
  }
  recent_ok = n_r >= min_hours
  recent[recent_ok, seq_along(omni_features)] = sum_r[recent_ok, , drop = FALSE] / n_r[recent_ok]
  recent[recent_ok, length(recent_names)] = sum_bs[recent_ok] / n_r[recent_ok]
  rot = matrix(NA_real_, n, outlook_days * length(rot_names))
  rot_ok = is.finite(lag_hr)
  if (n > 20000L) message(sprintf("Solar-wind block for %d issues", n))
  for (d in seq_len(outlook_days)) {
    hour0 = as.numeric(tau) - lag_hr * 3600 + (d - 1) * 86400
    base_d = hour_pos(origin, n_hours, hour0)
    sum_vs = rep(0, n)
    sum_bz = rep(0, n)
    sum_bs_d = rep(0, n)
    n_d = rep(0L, n)
    for (h in 0:23) {
      good = hour_mask(base_d, n, n_hours, ok_hour, h)
      use = which(good)
      if (!length(use)) next
      ix = ix_of(base_d, h)[use]
      n_d[use] = n_d[use] + 1L
      bz = mat[ix, bz_i]
      sum_bz[use] = sum_bz[use] + bz
      sum_vs[use] = sum_vs[use] + mat[ix, vs_i]
      sum_bs_d[use] = sum_bs_d[use] + pmax(-bz, 0)
    }
    cols = (d - 1L) * length(rot_names) + seq_along(rot_names)
    keep = rot_ok & n_d >= min_hours
    rot[keep, cols[1]] = sum_vs[keep] / n_d[keep]
    rot[keep, cols[2]] = sum_bz[keep] / n_d[keep]
    rot[keep, cols[3]] = sum_bs_d[keep] / n_d[keep]
    rot_ok = keep
  }
  rot[!rot_ok, ] = NA_real_
  list(recent = recent, rot = rot, ok = recent_ok & rot_ok)
}
ix_of = function(base, h) base + h
# Hourly AWSoM Bz is calibrated before taking max(-Bz, 0), because that
# rectification does not commute with OMNI = a + b * AWSoM.
daily_southward = function(tau, lag_hr, bz, ok_hour, origin, a, b_slope) {
  n = length(tau)
  n_hours = length(bz)
  out = matrix(NA_real_, n, outlook_days)
  ok = is.finite(lag_hr) & is.finite(a) & is.finite(b_slope)
  for (d in seq_len(outlook_days)) {
    hour0 = as.numeric(tau) - lag_hr * 3600 + (d - 1L) * 86400
    base_d = hour_pos(origin, n_hours, hour0)
    sum_bs = rep(0, n)
    n_d = rep(0L, n)
    for (h in 0:23) {
      good = hour_mask(base_d, n, n_hours, ok_hour, h)
      use = which(good)
      if (!length(use)) next
      n_d[use] = n_d[use] + 1L
      cal = a[use] + b_slope[use] * bz[ix_of(base_d, h)[use]]
      sum_bs[use] = sum_bs[use] + pmax(-cal, 0)
    }
    keep = ok & n_d >= min_hours
    out[keep, d] = sum_bs[keep] / n_d[keep]
    ok = keep
  }
  out[!ok, ] = NA_real_
  list(Bs = out, ok = ok)
}
recent_southward = function(tau, bz, ok_hour, origin, a, b_slope) {
  n = length(tau)
  n_hours = length(bz)
  base = hour_pos(origin, n_hours, as.numeric(tau) - history_hr * 3600)
  sum_bs = rep(0, n)
  n_r = rep(0L, n)
  for (h in 0:23) {
    good = hour_mask(base, n, n_hours, ok_hour, h)
    use = which(good)
    if (!length(use)) next
    n_r[use] = n_r[use] + 1L
    cal = a[use] + b_slope[use] * bz[ix_of(base, h)[use]]
    sum_bs[use] = sum_bs[use] + pmax(-cal, 0)
  }
  out = rep(NA_real_, n)
  keep = is.finite(a) & is.finite(b_slope) & n_r >= min_hours
  out[keep] = sum_bs[keep] / n_r[keep]
  list(Bs = out, ok = keep)
}

affine_fit = function(x, y) {
  ok = is.finite(x) & is.finite(y)
  x = x[ok]
  y = y[ok]
  if (length(x) < 5 || var(x) == 0) return(c(a = mean(y) - mean(x), b = 1))
  b = cov(x, y) / var(x)
  c(a = mean(y) - b * mean(x), b = b)
}
pair_row = match(as.numeric(model_h$time), as.numeric(omni_time))
model_pair = model_h$ok_hour & !is.na(pair_row)
omni_pair = model_pair
omni_pair[model_pair] = apply(OMNI_mat[pair_row[model_pair], , drop = FALSE], 1,
                              function(r) all(is.finite(r)))
elig = which(omni_pair)
if (length(elig) < history_hr * 7L) stop("Fewer than 7 days of AWSoM-OMNI pairs.")
pair_time = model_h$time[elig]
# For an issue at tau the recent driver starts at tau - 24 h, so the map
# uses pairs whose hour ended strictly before that.
map_ab = function(tau) {
  n = length(tau)
  last = findInterval(as.numeric(tau) - (history_hr + 1) * 3600, as.numeric(pair_time))
  cache = new.env(parent = emptyenv())
  out = matrix(NA_real_, n, 2L * length(omni_features))
  for (i in seq_len(n)) {
    L = last[i]
    if (L < 1L) next
    key = as.character(L)
    if (!exists(key, envir = cache, inherits = FALSE)) {
      newest = pair_time[L]
      ix = elig[seq_len(L)]
      ix = ix[model_h$time[ix] > newest - cal_days * 86400]
      ab = rep(NA_real_, 2L * length(omni_features))
      if (length(ix) >= 1L) {
        for (f in seq_along(omni_features)) {
          fit = affine_fit(model_h$mat[ix, f], OMNI_mat[pair_row[ix], f])
          ab[f] = unname(fit["a"])
          ab[length(omni_features) + f] = unname(fit["b"])
        }
      }
      assign(key, ab, envir = cache)
    }
    out[i, ] = get(key, envir = cache, inherits = FALSE)
  }
  list(a = out[, seq_along(omni_features), drop = FALSE],
       b = out[, length(omni_features) + seq_along(omni_features), drop = FALSE],
       ok = apply(out, 1, function(r) all(is.finite(r))))
}
calibrate_blocks = function(blocks, a, b) {
  recent = blocks$recent
  rot = blocks$rot
  for (f in seq_along(omni_features))
    recent[, f] = a[, f] + b[, f] * recent[, f]
  vs_f = match("Vs", omni_features)
  bz_f = match("Bz", omni_features)
  for (d in seq_len(outlook_days)) {
    cols = (d - 1L) * length(rot_names) + seq_along(rot_names)
    rot[, cols[1]] = a[, vs_f] + b[, vs_f] * rot[, cols[1]]
    rot[, cols[2]] = a[, bz_f] + b[, bz_f] * rot[, cols[2]]
    # Column 3 is southward Bz. The mean of max(-Bz, 0) is not linear, so
    # it is recomputed from hourly calibrated Bz after this returns.
  }
  blocks$recent = recent
  blocks$rot = rot
  blocks
}

message(sprintf("Training cutoff %s", fmt_utc(train_end)))
i_min = max(outlook_lags_hr) / kp_step_hr + outlook_match_bins
i_train = seq.int(i_min, N_file - outlook_days * bins_per_day)
tau_train = kp_time[i_train] + kp_step_hr * 3600
i_train = i_train[tau_train + outlook_days * 86400 <= train_end]
tau_train = kp_time[i_train] + kp_step_hr * 3600
message(sprintf("Building the training design, %d issues", length(i_train)))
tr_kp = kp_block(kp_value, i_train)
tr_omni = solar_blocks(tau_train, tr_kp$lag, OMNI_mat,
                       apply(OMNI_mat, 1, function(r) all(is.finite(r))),
                       as.numeric(omni_time[1]))
tr_keep = tr_kp$ok & tr_omni$ok & finite_rows(tr_kp$Y)
message(sprintf("Training rows with complete OMNI solar wind: %d", sum(tr_keep)))
X_kp = tr_kp$X[tr_keep, , drop = FALSE]
X_full = cbind(X_kp, tr_omni$recent[tr_keep, , drop = FALSE], tr_omni$rot[tr_keep, , drop = FALSE])
X_rot = cbind(X_kp, tr_omni$rot[tr_keep, , drop = FALSE])
Y_tr = tr_kp$Y[tr_keep, , drop = FALSE]
stopifnot(all(is.finite(X_full)), all(is.finite(X_rot)), all(is.finite(Y_tr)))
fit_one = function(X, label) {
  message(sprintf("Fitting %s, %d columns", label, ncol(X) + 1L))
  sc = standardize_fit(X)
  fit = ols_coef_from_svd(cbind(1, apply_scale(X, sc$center, sc$scale)), Y_tr, ols_tol)
  message(sprintf("  rank %d / %d", fit$rank, ncol(X) + 1L))
  list(center = sc$center, scale = sc$scale, B = fit$B, rank = fit$rank, label = label)
}
fit_full = fit_one(X_full, "Kp plus recent solar wind and the previous rotation")
fit_rot = fit_one(X_rot, "Kp plus the previous rotation's solar wind")
rm(X_full, X_rot, X_kp, Y_tr)
gc()

message("Building the held-out issues")
i_te = seq.int(i_min, N_kp - bins_per_day)
tau_te = kp_time[i_te] + kp_step_hr * 3600
keep_te = tau_te > train_end
i_te = i_te[keep_te]
tau_te = tau_te[keep_te]
te_kp = kp_block(kp_value, i_te)
te_obs = solar_blocks(tau_te, te_kp$lag, obs_h$mat, obs_h$ok_hour, as.numeric(obs_h$time[1]))
te_mod = solar_blocks(tau_te, te_kp$lag, model_h$mat, model_h$ok_hour, as.numeric(model_h$time[1]))
ab = map_ab(tau_te)
te_mod$ok = te_mod$ok & ab$ok
te_mod = calibrate_blocks(te_mod, ab$a, ab$b)
bz_f = match("Bz", omni_features)
bs = daily_southward(tau_te, te_kp$lag, model_h$mat[, bz_f], model_h$ok_hour,
                     as.numeric(model_h$time[1]), ab$a[, bz_f], ab$b[, bz_f])
for (d in seq_len(outlook_days))
  te_mod$rot[, (d - 1L) * length(rot_names) + 3L] = bs$Bs[, d]
rb = recent_southward(tau_te, model_h$mat[, bz_f], model_h$ok_hour,
                      as.numeric(model_h$time[1]), ab$a[, bz_f], ab$b[, bz_f])
te_mod$recent[, length(recent_names)] = rb$Bs
te_mod$ok = te_mod$ok & bs$ok & rb$ok
te_omni = solar_blocks(tau_te, te_kp$lag, OMNI_mat,
                       apply(OMNI_mat, 1, function(r) all(is.finite(r))),
                       as.numeric(omni_time[1]))

score_sample = function(sample, row_ok, models) {
  rows = list()
  r = 0L
  for (d in seq_len(outlook_days)) {
    use = row_ok & is.finite(te_kp$Y[, d]) & is.finite(te_kp$recur[, d]) & is.finite(te_kp$pers)
    if (sum(use) < 2L) next
    y = te_kp$Y[use, d]
    y_rec = te_kp$recur[use, d]
    y_pers = te_kp$pers[use]
    Xk = te_kp$X[use, , drop = FALSE]
    y_kp = predict_y(Xk, fit_kp$center, fit_kp$scale, fit_kp$B)[, d]
    built = list(
      kp = y_kp,
      insitu = predict_y(cbind(Xk, te_obs$recent[use, , drop = FALSE], te_obs$rot[use, , drop = FALSE]),
                         fit_full$center, fit_full$scale, fit_full$B)[, d],
      awsom = predict_y(cbind(Xk, te_mod$recent[use, , drop = FALSE], te_mod$rot[use, , drop = FALSE]),
                        fit_full$center, fit_full$scale, fit_full$B)[, d],
      insitu_rotation = predict_y(cbind(Xk, te_obs$rot[use, , drop = FALSE]),
                                  fit_rot$center, fit_rot$scale, fit_rot$B)[, d],
      awsom_rotation = predict_y(cbind(Xk, te_mod$rot[use, , drop = FALSE]),
                                 fit_rot$center, fit_rot$scale, fit_rot$B)[, d],
      omni = predict_y(cbind(Xk, te_omni$recent[use, , drop = FALSE], te_omni$rot[use, , drop = FALSE]),
                       fit_full$center, fit_full$scale, fit_full$B)[, d],
      omni_rotation = predict_y(cbind(Xk, te_omni$rot[use, , drop = FALSE]),
                                fit_rot$center, fit_rot$scale, fit_rot$B)[, d]
    )
    rmse = function(a, b) sqrt(mean((a - b)^2))
    rk = rmse(y, y_kp)
    for (nm in models) {
      yh = built[[nm]]
      rm = rmse(y, yh)
      r = r + 1L
      rows[[r]] = data.frame(
        sample = sample, model = nm, day = d, n = length(y),
        RMSE = rm, MAE = mean(abs(y - yh)),
        Correlation = suppressWarnings(cor(y, yh)),
        RMSE_kp = rk,
        skill_vs_kp = if (rk == 0) NA_real_ else 1 - rm / rk,
        RMSE_recur = rmse(y, y_rec),
        RMSE_pers = rmse(y, y_pers),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, rows)
}

paired_ok = te_kp$ok & te_obs$ok & te_mod$ok
omni_ok = te_kp$ok & te_omni$ok
message(sprintf("Paired RTSW issues %d, %s through %s",
                sum(paired_ok), fmt_utc(min(tau_te[paired_ok])), fmt_utc(max(tau_te[paired_ok]))))
message(sprintf("OMNI held-out issues %d, %s through %s",
                sum(omni_ok), fmt_utc(min(tau_te[omni_ok])), fmt_utc(max(tau_te[omni_ok]))))
tab = rbind(
  score_sample("paired_rtsw", paired_ok, c("kp", "insitu", "awsom", "insitu_rotation", "awsom_rotation")),
  score_sample("omni_heldout", omni_ok, c("kp", "omni", "omni_rotation"))
)
out_path = file.path(out_dir, "outlook_sw_comparison.csv")
write.csv(tab, out_path, row.names = FALSE)
message(sprintf("Wrote %s", out_path))

show_days = c(1L, 7L, 14L, 27L)
for (s in unique(tab$sample)) {
  message(sprintf("--- %s ---", s))
  sub = tab[tab$sample == s & tab$day %in% show_days, ]
  for (d in show_days) {
    sd = sub[sub$day == d, ]
    if (!nrow(sd)) next
    message(sprintf("day %d  n %d  recurrence %.3f  previous day %.3f",
                    d, sd$n[1], sd$RMSE_recur[1], sd$RMSE_pers[1]))
    for (i in seq_len(nrow(sd))) {
      message(sprintf("  %-18s RMSE %.3f  MAE %.3f  corr %.3f  vs Kp-only %+.3f",
                      sd$model[i], sd$RMSE[i], sd$MAE[i], sd$Correlation[i], sd$skill_vs_kp[i]))
    }
  }
}

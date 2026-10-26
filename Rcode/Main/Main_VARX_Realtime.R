# =============================================================================
# Main_VARX_Realtime.R
# One-shot realtime Kp forecast from the latest RTSW obs and model
# =============================================================================
#
# Run every three hours:
#   Rscript Rcode/Main/Main_VARX_Realtime.R
#
# What it grabs
# -------------
# RTSW Earth series (obs and model) from
#   https://solarwind.engin.umich.edu/rtsw/api/series/?t=earth
# Bins after that file come from NOAA planetary Kp, rounded to the nearest
# third. A bin NOAA no longer has is filled from the GFZ nowcast. A bin
# neither source has is left missing. Training Kp is unchanged.
#
# The issued forecast is the better model at each range. The next 24 hours are
# the original vanilla XGBoost, loaded by Issue_Best.R. From 27 hours through
# 27 days the issued value is the in-situ VARX below, rolled forward. F10.7 is
# a separate daily ridge regression. In-situ and AWSoM still share one VARX;
# only the solar-wind block changes. Besides the previous 24 hours, it uses
# Kp and mean speed from one earlier window. The lag is chosen, per issue
# and per input, from 25 to 29 days at 12 h steps, by how well speed over
# the last 5 days matches that earlier window. It also uses how Kp changed
# over the 24 hours after earlier windows whose own Kp shape matches the
# last 48 hours. Those changes are 25 to 29 days old, so they have already
# ended. Only hours that have already ended are used (hour start <= tau - 1h):
#   * Earth_obs:     Data/Origin/RTSW/Earth_obs.csv, hourly mean, no adjustment.
#                    Shown as "In-situ measurements".
#   * model_affine:  Data/Origin/RTSW/Earth_model.csv, hourly mean, then
#                    OMNI = a + b * AWSoM. Shown as "AWSoM Earth model".
#                    a and b are refit on the latest 30 days of AWSoM–OMNI
#                    pairs that ended before this issue's 24 h driver.
# Solar wind inside the forecast window is not a predictor.
#
# Training labels end before the first RTSW model sample. Coefficients are
# cached and refit only when Kp, OMNI, or that training cutoff changes.
# Pages are written by Write_VARX_Pages.R.
#
# VARX_Kp_baseline.html scores this model where the 24 h target Kp has
# already been observed. The first 30 days after RTSW starts are the affine
# calibration window and are not scored. Headline criteria are the endpoint
# bin only. Intermediate leads are on the same page. The same page scores
# the 27-day outlook on the live Kp series, through the latest ended bin.
# Day d is scored as soon as that day has ended, so day 1 reaches the
# current Kp and day 27 reaches the Kp from about four weeks earlier.
#
# A second product on the same page is the next 27 days. Each day is the
# mean of its eight Kp bins. The predictors are the last 24 hours of Kp and
# the daily means from one earlier window of 25 to 29 days, every 3 hours,
# chosen by the correlation of the last 5 days of Kp. A lag shorter than 27
# days would put the last analogue day after the issue, so it is not chosen.
# That regression is then mixed with the Kp that followed the ten closest
# 27-day OMNI cycles from 1995 through the training cutoff. The query is the
# last 27 days of in-situ solar wind. The weight uses only earlier issues
# whose target day has ended, and it stays between the two forecasts.
# Coefficients of the regression are cached separately from the 24 h VARX.
# =============================================================================

script_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (Sys.getlocale("LC_CTYPE") == "C") {
  loc = tryCatch(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), error = function(e) "")
  if (!nzchar(loc)) try(Sys.setlocale("LC_CTYPE", "UTF-8"), silent = TRUE)
}
if (length(script_arg)) {
  script_path = normalizePath(sub("^--file=", "", script_arg[1]))
  repo_root = normalizePath(file.path(dirname(script_path), "..", ".."))
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
sentinel = c(B = 9999, Bz = 9999, Vs = 9999, Pd = 9999, Pt = 9999999)

kp_path = file.path(repo_root, "Data", "Processed", "Kp.csv")
omni_path = file.path(repo_root, "Data", "Processed", "OMNI_hourly.csv")
out_dir = file.path(repo_root, "Rcode", "Main", "output", "Main_VARX_Realtime")
coef_path = file.path(out_dir, "VARX_Kp_coef.RData")
# Size and checksum, not modification time. A checkout on GitHub changes
# mtime, and that must not throw away a cache of the same file.
file_stamp = function(path) {
  info = file.info(path)
  list(size = unname(info$size), md5 = unname(tools::md5sum(path)))
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(repo_root, "Data", "Origin", "RTSW"), recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(kp_path), file.exists(omni_path))

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

curl_to = function(url, dest, extra = character()) {
  status = system2("curl", c("--max-time", "90", "-fsSL", extra, url, "-o", dest))
  if (!identical(status, 0L)) {
    message("Verified download failed; retrying without certificate check.")
    status = system2("curl", c("--max-time", "90", "-kfsSL", extra, url, "-o", dest))
  }
  if (!identical(status, 0L)) stop("curl download failed for ", url, ", status ", status)
  invisible(dest)
}

# ---- grab RTSW ---------------------------------------------------------------
# KP_SCORE_OUTLOOK=1 rebuilds only the scored 27-day section from Kp on disk.
score_outlook_only = identical(Sys.getenv("KP_SCORE_OUTLOOK"), "1")
obs_path = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_obs.csv")
model_path = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_model.csv")
if (!score_outlook_only) {
  cat("Downloading RTSW Earth obs and model\n")
  status = system2(Sys.which("Rscript"),
                   file.path(repo_root, "Rcode", "Main", "Download_RTSW_2026.R"))
  if (!identical(status, 0L)) stop("RTSW download failed, status ", status)
  stopifnot(file.exists(obs_path), file.exists(model_path))
}

# ---- load archives -----------------------------------------------------------
Kp_df = read.csv(kp_path, stringsAsFactors = FALSE)
OMNI_df = read.csv(omni_path, stringsAsFactors = FALSE)
Kp_df$Date = parse_utc(Kp_df$Date)
OMNI_df$Date = parse_utc(OMNI_df$Date)
Kp_df = Kp_df[order(Kp_df$Date), ]
OMNI_df = OMNI_df[order(OMNI_df$Date), ]
stopifnot(!anyDuplicated(Kp_df$Date), !anyDuplicated(OMNI_df$Date))
stopifnot(all(omni_features %in% names(OMNI_df)), "Data" %in% names(Kp_df))

kp_time = Kp_df$Date
kp_value = Kp_df$Data
omni_time = OMNI_df$Date
OMNI_mat = as.matrix(OMNI_df[, omni_features, drop = FALSE])
N_kp = length(kp_time)
N_omni = nrow(OMNI_mat)
n_kp_hist = history_hr / kp_step_hr
vs_col = match("Vs", omni_features)

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
  list(time = hours, mat = out, n_native = n_native,
       native_start = min(time), native_end = max(time))
}
# Hours inside the driver that are missing. A gap with a known hour on both
# sides is linear interpolation. Those hours have to have ended by the issue,
# so a sample after tau is not used. If the series has already stopped, the
# last hour is held for at most one Kp bin. That is what lets a finished Kp
# bin be used when AWSoM or in-situ is a few hours behind.
causal_solar_wind = function(hour_num, series_time, series_mat, hold_hr) {
  # A vector is one issue, oldest hour first. A matrix is issues by hours.
  if (is.null(dim(hour_num))) {
    hour_mat = matrix(as.numeric(hour_num), nrow = 1)
  } else {
    hour_mat = matrix(as.numeric(hour_num), nrow(hour_num), ncol(hour_num))
  }
  know_row = apply(hour_mat, 1, max)
  know = rep(know_row, times = ncol(hour_mat))
  h = as.vector(hour_mat)
  st = as.numeric(series_time)
  out = matrix(NA_real_, length(h), ncol(series_mat))
  colnames(out) = colnames(series_mat)
  hold_sec = hold_hr * 3600
  for (j in seq_len(ncol(series_mat))) {
    ok = is.finite(series_mat[, j])
    if (!any(ok)) next
    sx = st[ok]
    sy = series_mat[ok, j]
    ord = order(sx)
    sx = sx[ord]
    sy = sy[ord]
    ix = match(h, sx)
    hit = !is.na(ix)
    out[hit, j] = sy[ix[hit]]
    miss = which(!hit)
    if (!length(miss)) next
    prev = findInterval(h[miss], sx)
    has_prev = prev >= 1L
    has_next = prev < length(sx)
    prev_t = rep(NA_real_, length(miss))
    next_t = prev_t
    prev_t[has_prev] = sx[prev[has_prev]]
    next_i = prev + 1L
    next_t[has_next] = sx[next_i[has_next]]
    can_interp = has_prev & has_next & next_t <= know[miss]
    can_interp[is.na(can_interp)] = FALSE
    can_hold = has_prev & !can_interp & (h[miss] - prev_t) <= hold_sec &
      (h[miss] - prev_t) > 0 & h[miss] <= know[miss]
    can_hold[is.na(can_hold)] = FALSE
    ii = which(can_interp)
    if (length(ii)) {
      w = (h[miss][ii] - prev_t[ii]) / (next_t[ii] - prev_t[ii])
      out[miss[ii], j] = (1 - w) * sy[prev[ii]] + w * sy[next_i[ii]]
    }
    ii = which(can_hold)
    if (length(ii)) out[miss[ii], j] = sy[prev[ii]]
  }
  out
}

if (!score_outlook_only) {
cat("Reducing RTSW obs and model to hourly means\n")
obs_h = make_rtsw_hourly(obs_path)
model_h = make_rtsw_hourly(model_path)
# The lag match reads the same hourly series. Hold the last complete hour for
# one Kp bin so a finished bin is not dropped when the feed is a few hours short.
# A hole earlier in the series is filled inside the lag match, same rule as
# the 24 h driver: interpolate when a later sample is already known, otherwise
# leave it missing.
hold_series_tail = function(h, hold_hr) {
  ok = apply(h$mat, 1, function(r) all(is.finite(r)))
  if (!any(ok) || hold_hr < 1L) return(h)
  last = max(which(ok))
  add = seq_len(hold_hr)
  h$time = c(h$time, h$time[last] + add * 3600)
  h$mat = rbind(h$mat, h$mat[rep(last, hold_hr), , drop = FALSE])
  h$n_native = c(h$n_native, integer(hold_hr))
  h
}
obs_h = hold_series_tail(obs_h, kp_step_hr)
model_h = hold_series_tail(model_h, kp_step_hr)
obs_real = obs_h$n_native > 0 & apply(obs_h$mat, 1, function(r) all(is.finite(r)))
if (!any(obs_real)) stop("RTSW obs has no hour with a real sample.")
last_obs = max(obs_h$time[obs_real])
cat(sprintf("Obs hours with a real sample end at %s\n", fmt_utc(last_obs)))
cat(sprintf("Model hourly ends at %s\n", fmt_utc(max(model_h$time))))

# Updated AWSoM calibration. Pairs are hours with a real AWSoM sample and a
# finite OMNI hour. For an issue at tau the driver starts at tau - 24h, so
# the map uses only pairs that ended strictly before that, and only the
# latest cal_days of those pairs. The same hour can therefore get a
# different map at a later issue.
affine_fit = function(x, y) {
  ok = is.finite(x) & is.finite(y)
  x = x[ok]
  y = y[ok]
  if (length(x) < 5 || var(x) == 0) return(c(a = mean(y) - mean(x), b = 1))
  b = cov(x, y) / var(x)
  c(a = mean(y) - b * mean(x), b = b)
}
pair_omni = match(as.numeric(model_h$time), as.numeric(omni_time))
model_pair_ok = model_h$n_native > 0 & apply(model_h$mat, 1, function(r) all(is.finite(r)))
omni_pair_ok = !is.na(pair_omni)
omni_pair_ok[omni_pair_ok] = apply(OMNI_mat[pair_omni[omni_pair_ok], , drop = FALSE], 1,
                                   function(r) all(is.finite(r)))
elig_ix = which(model_pair_ok & omni_pair_ok)
if (length(elig_ix) < 24L * 7L) stop("Fewer than 7 days of AWSoM–OMNI pairs.")
affine_cache = new.env(parent = emptyenv())
affine_updated = function(cutoff) {
  ended = model_h$time[elig_ix] <= cutoff - 3600
  if (!any(ended)) stop("No AWSoM–OMNI pair ended before ", fmt_utc(cutoff))
  last = max(which(ended))
  key = as.character(last)
  if (exists(key, envir = affine_cache, inherits = FALSE))
    return(get(key, envir = affine_cache, inherits = FALSE))
  ix = elig_ix[seq_len(last)]
  latest = model_h$time[ix[length(ix)]]
  ix = ix[model_h$time[ix] > latest - cal_days * 86400]
  a = numeric(length(omni_features))
  b = a
  for (f in seq_along(omni_features)) {
    fit = affine_fit(model_h$mat[ix, f], OMNI_mat[pair_omni[ix], f])
    a[f] = unname(fit["a"])
    b[f] = unname(fit["b"])
  }
  out = list(a = a, b = b, n = length(ix),
             start = min(model_h$time[ix]), end = max(model_h$time[ix]))
  assign(key, out, envir = affine_cache)
  out
}
cat(sprintf("AWSoM–OMNI pairs available through %s\n",
            fmt_utc(max(model_h$time[elig_ix]))))

# ---- extend Kp after the file, without touching the training series -------
# Training keeps Data/Processed/Kp.csv. Bins in the file are not overwritten.
# A missing bin is filled from NOAA, then from the GFZ nowcast. A bin that
# neither source has, or that has not ended, stays missing.
read_noaa_kp = function() {
  noaa_tmp = tempfile(fileext = ".json")
  on.exit(unlink(noaa_tmp), add = TRUE)
  ok = tryCatch({
    curl_to("https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json",
            noaa_tmp, c("-A", "Mozilla/5.0"))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(NULL)
  noaa = tryCatch(jsonlite::fromJSON(noaa_tmp), error = function(e) NULL)
  if (is.null(noaa) || is.null(noaa$time_tag)) return(NULL)
  start = parse_utc(gsub("T", " ", noaa$time_tag, fixed = TRUE))
  kp = round(as.numeric(noaa$Kp) * 3) / 3
  ord = order(start)
  start = start[ord]
  kp = kp[ord]
  ended = which(is.finite(kp) & kp >= 0 & kp <= 9 & start + kp_step_hr * 3600 <= Sys.time())
  if (!length(ended)) return(NULL)
  data.frame(start = start[ended], kp = kp[ended], stringsAsFactors = FALSE)
}
read_gfz_kp = function() {
  tmp = tempfile(fileext = ".txt")
  on.exit(unlink(tmp), add = TRUE)
  ok = tryCatch({
    curl_to("https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_nowcast.txt",
            tmp, c("--ipv4", "-A", "Mozilla/5.0"))
    TRUE
  }, error = function(e) FALSE)
  if (!ok) return(NULL)
  lines = readLines(tmp, warn = FALSE)
  lines = lines[!grepl("^\\s*#", lines) & nzchar(trimws(lines))]
  if (!length(lines)) return(NULL)
  nw = tryCatch(read.table(text = lines, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(nw) || ncol(nw) < 15L) return(NULL)
  hour_head = seq(0, 21, by = kp_step_hr)
  start = rep(as.POSIXct(NA), nrow(nw) * length(hour_head))
  kp = rep(NA_real_, length(start))
  n = 0L
  for (r in seq_len(nrow(nw))) {
    for (k in seq_along(hour_head)) {
      val = nw[r, 7L + k]
      if (!is.finite(val) || val < 0 || val > 9) next
      t0 = as.POSIXct(sprintf("%04d-%02d-%02d %02d:00:00",
                              nw[r, 1], nw[r, 2], nw[r, 3], hour_head[k]),
                      format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
      if (is.na(t0) || t0 + kp_step_hr * 3600 > Sys.time()) next
      n = n + 1L
      start[n] = t0
      kp[n] = round(val * 3) / 3
    }
  }
  if (!n) return(NULL)
  data.frame(start = start[seq_len(n)], kp = kp[seq_len(n)], stringsAsFactors = FALSE)
}
refresh_noaa_tail = function() {
  noaa = read_noaa_kp()
  gfz = read_gfz_kp()
  if (is.null(noaa) && is.null(gfz))
    return(list(ok = FALSE, updated = 0L, added = 0L, from_gfz = 0L, open = 0L))
  last_set = c(
    if (!is.null(noaa) && nrow(noaa)) as.numeric(max(noaa$start)),
    if (!is.null(gfz) && nrow(gfz)) as.numeric(max(gfz$start))
  )
  last_set = last_set[is.finite(last_set) & last_set > as.numeric(kp_file_end)]
  if (!length(last_set))
    return(list(ok = TRUE, updated = 0L, added = 0L, from_gfz = 0L, open = 0L))
  updated = 0L
  added = 0L
  from_gfz = 0L
  open = 0L
  step = kp_step_hr * 3600
  t = kp_file_end + step
  last = max(last_set)
  while (as.numeric(t) <= last + 1) {
    jn = if (is.null(noaa)) NA_integer_ else match(as.numeric(t), as.numeric(noaa$start))
    jg = if (is.null(gfz)) NA_integer_ else match(as.numeric(t), as.numeric(gfz$start))
    val = NA_real_
    used_gfz = FALSE
    if (!is.na(jn)) {
      val = noaa$kp[jn]
    } else if (!is.na(jg)) {
      val = gfz$kp[jg]
      used_gfz = TRUE
    }
    k = match(as.numeric(t), as.numeric(kp_time_live))
    if (!is.na(k)) {
      if (is.finite(val) && !isTRUE(all.equal(kp_value_live[k], val))) {
        kp_value_live[k] <<- val
        updated = updated + 1L
      }
    } else {
      kp_time_live <<- c(kp_time_live, t)
      kp_value_live <<- c(kp_value_live, val)
      added = added + 1L
    }
    if (used_gfz && is.finite(val)) from_gfz = from_gfz + 1L
    if (!is.finite(val)) open = open + 1L
    t = t + step
  }
  list(ok = TRUE, updated = updated, added = added, from_gfz = from_gfz, open = open)
}
kp_time_live = kp_time
kp_value_live = kp_value
kp_file_end = max(kp_time)
kp_extend_source = NULL
nowcast_note = "NOAA planetary Kp was not applied. Kp history is the processed file only."
noaa_tail = refresh_noaa_tail()
if (isTRUE(noaa_tail$ok)) {
  last_i = length(kp_time_live)
  kp_extend_source = "NOAA"
  fill_note = ""
  if (noaa_tail$from_gfz > 0L)
    fill_note = paste0(fill_note, sprintf(
      " %d bins missing from NOAA were filled from the GFZ nowcast.", noaa_tail$from_gfz))
  if (noaa_tail$open > 0L)
    fill_note = paste0(fill_note, sprintf(
      " %d bins are still missing and are left open.", noaa_tail$open))
  nowcast_note = sprintf(
    "Kp after %s is NOAA planetary Kp, rounded to the nearest third.%s Latest ended bin %s to %s, Kp %.3f. Training Kp is unchanged.",
    fmt_utc(kp_file_end), fill_note,
    fmt_utc(kp_time_live[last_i]),
    fmt_utc(kp_time_live[last_i] + kp_step_hr * 3600),
    kp_value_live[last_i])
} else {
  message("NOAA planetary Kp was not applied.")
}
cat(nowcast_note, "\n")
} else {
  cached_outlook = new.env()
  load(file.path(out_dir, "VARX_Kp_outlook27_coef.RData"), envir = cached_outlook)
  train_end = as.POSIXct(cached_outlook$outlook_sig$train_end, origin = "1970-01-01", tz = "UTC")
  kp_time_live = kp_time
  kp_value_live = kp_value
  cat("Scoring the 27-day outlook only. RTSW and the GFZ nowcast were not downloaded.\n")
}

# ---- one 24 h VARX: previous 24 h, plus a searched 25–29 day lag ----------
# Extra columns: 8 Kp bins ending at tau - L, mean Vs over the 24 h ending
# at tau - L, and the match correlation. L is the lag in lags_hr with the
# highest speed correlation. Coupling of Bz and speed is not used.
if (!score_outlook_only) train_end = model_h$native_start
issue_time = function(i) kp_time[i] + kp_step_hr * 3600

standardize_fit = function(M) {
  ctr = colMeans(M)
  scl = apply(M, 2, sd)
  scl[scl == 0 | is.na(scl)] = 1
  list(center = ctr, scale = scl)
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
  B = t(sv$v %*% (Dinv * UtY))
  list(B = B, n_dropped = sum(!keep), rank = sum(keep))
}

# ---- recurrence state --------------------------------------------------------
finite_rows = function(M) rowSums(is.finite(M)) == ncol(M)
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
row_cor = function(A, B) {
  A = A - rowMeans(A)
  B = B - rowMeans(B)
  num = rowSums(A * B)
  den = sqrt(rowSums(A * A) * rowSums(B * B))
  out = num / den
  out[!is.finite(out)] = 0
  out
}
hours_ending = function(t_end, n_back) {
  outer(as.numeric(t_end), (n_back:1L) * 3600, "-")
}
# Growth of Kp after earlier windows, weighted by how well that window's Kp
# shape matches the last kp_match_bins. This is the previous rotation's
# rise or fall, not the Kp level at the lag.
kp_match_bins = 16L
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
  as_rows = function(M) {
    if (is.null(dim(M))) M = matrix(M, nrow = n)
    M
  }
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
                       speed_time, speed_vs, vs_time, vs_values,
                       scale_vs = NULL, speed_latest = NULL) {
  n = length(tau)
  stopifnot(nrow(kp_hist_mat) == n, nrow(exog) == n)
  stopifnot(ncol(exog) == history_hr * length(omni_features))
  # A short hole in speed is filled the same way as the 24 h driver. Hours
  # after speed_latest stay missing, so a sample that has not ended is not used.
  pull_speed = function(hour_num, grid, values, latest) {
    filled = causal_solar_wind(hour_num, grid, matrix(values, ncol = 1), kp_step_hr)
    if (is.null(dim(hour_num))) {
      out = matrix(filled, nrow = 1)
    } else {
      out = matrix(filled, nrow = nrow(hour_num), ncol = ncol(hour_num))
    }
    if (!is.null(latest)) {
      bad = matrix(as.vector(hour_num) > as.numeric(latest),
                   nrow = nrow(out), ncol = ncol(out))
      out[bad] = NA_real_
    }
    out
  }
  recent = pull_speed(hours_ending(tau, match_hr), speed_time, speed_vs, speed_latest)
  cors = matrix(NA_real_, n, length(lags_hr))
  for (j in seq_along(lags_hr)) {
    past = pull_speed(hours_ending(tau - lags_hr[j] * 3600, match_hr),
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
  raw_vs = pull_speed(hours_ending(tau - best_lag * 3600, history_hr),
                      vs_time, vs_values, speed_latest)
  rvs = rowMeans(raw_vs)
  if (!is.null(scale_vs)) rvs = scale_vs[, 1] + scale_vs[, 2] * rvs
  grow = kp_matched_growth(i_kp, kp_src)
  ok = has & finite_rows(rkp) & is.finite(rvs) & is.finite(rcor) & finite_rows(grow)
  pred = cbind(kp_hist_mat, exog, rkp, rvs, rcor, grow)
  list(ok = ok, pred = pred, lag = best_lag, grow = grow)
}
recur_predict = function(pred) {
  scaled = sweep(sweep(pred, 2, fit_recur$center, "-"), 2, fit_recur$scale, "/")
  cbind(1, scaled) %*% t(fit_recur$B)
}

recur_sig = list(
  kp = file_stamp(kp_path),
  omni = file_stamp(omni_path),
  train_end = as.numeric(train_end),
  history_hr = history_hr,
  horizon_hr = 24L,
  omni_features = omni_features,
  ols_tol = ols_tol,
  match_hr = match_hr,
  lags_hr = lags_hr,
  kp_match_bins = kp_match_bins,
  model = "period_kpgrowth"
)
fit_recur = NULL
if (file.exists(coef_path)) {
  env = new.env()
  load(coef_path, envir = env)
  if (isTRUE(all.equal(env$recur_sig, recur_sig))) fit_recur = env$fit_recur
}
if (is.null(fit_recur)) {
  cat("Fitting the 24 h VARX with a searched 25-29 day lag\n")
  n_ahead = 24L / kp_step_hr
  look_hr = max(lags_hr) + match_hr
  i_h = seq.int(n_kp_hist, N_kp - n_ahead)
  i_h = i_h[3L * i_h <= N_omni]
  i_h = i_h[(3L * i_h - look_hr) >= 1L]
  i_h = i_h[(i_h - max(lags_hr) / kp_step_hr - kp_match_bins) >= 1L]
  i_train = i_h[issue_time(i_h) + 24 * 3600 <= train_end]
  end_idx = 3L * i_train
  stopifnot(max(abs(as.numeric(omni_time[end_idx] - (issue_time(i_train) - 3600)))) < 1)
  hour_block = function(v, end_at, n_back) {
    sapply((n_back - 1L):0L, function(age) v[end_at - age])
  }
  kp_hist_tr = sapply(seq_len(n_kp_hist), function(j) kp_value[i_train - n_kp_hist + j])
  feat = lapply(seq_along(omni_features), function(f) hour_block(OMNI_mat[, f], end_idx, history_hr))
  recent = hour_block(OMNI_mat[, vs_col], end_idx, match_hr)
  cors = sapply(lags_hr, function(L) {
    row_cor(recent, hour_block(OMNI_mat[, vs_col], end_idx - L, match_hr))
  })
  best_j = max.col(cors, ties.method = "first")
  best_lag = lags_hr[best_j]
  recur_cor = cors[cbind(seq_along(best_j), best_j)]
  lag_bins = best_lag / kp_step_hr
  recur_kp = sapply(seq_len(n_kp_hist), function(j) kp_value[i_train - lag_bins - n_kp_hist + j])
  recur_vs = rowMeans(sapply(seq_len(history_hr), function(j) {
    OMNI_mat[end_idx - best_lag - history_hr + j, vs_col]
  }))
  grow_tr = kp_matched_growth(i_train, kp_value)
  aug_pred = cbind(kp_hist_tr, do.call(cbind, feat), recur_kp, recur_vs, recur_cor, grow_tr)
  keep_pred = rowSums(is.finite(aug_pred)) == ncol(aug_pred)
  if (!all(keep_pred)) {
    cat(sprintf("  dropping %d training rows without a Kp-shape match\n", sum(!keep_pred)))
    i_train = i_train[keep_pred]
    aug_pred = aug_pred[keep_pred, , drop = FALSE]
    best_lag = best_lag[keep_pred]
  }
  stopifnot(all(is.finite(aug_pred)))
  check = recur_batch(
    issue_time(i_train[1]), i_train[1], kp_value,
    matrix(aug_pred[1, seq_len(n_kp_hist)], nrow = 1),
    matrix(aug_pred[1, n_kp_hist + seq_len(history_hr * length(omni_features))], nrow = 1),
    omni_time, OMNI_mat[, vs_col], omni_time, OMNI_mat[, vs_col])
  stopifnot(check$ok)
  if (max(abs(check$pred - aug_pred[1, ])) > 1e-6)
    stop("Recurrence features from timestamps do not match the training row.")
  Y = sapply(seq_len(n_ahead), function(k) kp_value[i_train + k])
  sc = standardize_fit(aug_pred)
  fit = ols_coef_from_svd(cbind(1, apply_scale(aug_pred, sc$center, sc$scale)), Y, ols_tol)
  fit_recur = list(
    B = fit$B, center = sc$center, scale = sc$scale,
    rank = fit$rank, n_dropped = fit$n_dropped, n_train = length(i_train))
  cat(sprintf("  24h  train %d  rank %d / %d  median lag %.1f days\n",
              length(i_train), fit$rank, ncol(aug_pred) + 1L, median(best_lag) / 24))
  save(fit_recur, recur_sig, file = coef_path)
} else {
  cat(sprintf("Loaded cached coefficients (%s)\n", coef_path))
}

# ---- 1 to 27 day outlook: daily mean Kp, previous rotation already observed -
outlook_days = 27L
outlook_lags_hr = seq(25L * 24L, 29L * 24L, by = kp_step_hr)
outlook_match_bins = 5L * (24L / kp_step_hr)
bins_per_day = 24L / kp_step_hr
outlook_coef_path = file.path(out_dir, "VARX_Kp_outlook27_coef.RData")
outlook_sig = list(
  kp = file_stamp(kp_path),
  train_end = as.numeric(train_end),
  n_days = outlook_days,
  lags_hr = outlook_lags_hr,
  match_bins = outlook_match_bins,
  model = "kp_daily_recur_27"
)
day_means = function(kp, i, end_off, n_days) {
  n = length(i)
  out = matrix(NA_real_, n, n_days)
  if (!n) return(out)
  for (d in seq_len(n_days)) {
    acc = rep(0, n)
    last = i + end_off + d * bins_per_day
    good = last >= bins_per_day & last <= length(kp)
    acc[!good] = NA_real_
    for (a in 0:(bins_per_day - 1L)) {
      ix = last - a
      acc[good] = acc[good] + kp[ix[good]]
    }
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
  # Day 27 of the analogue ends at tau - L + 27 days. L must be at least 27
  # days, otherwise that day is still in the future.
  cors[, outlook_lags_hr < outlook_days * 24] = NA_real_
  finite_cor = is.finite(cors)
  cors[!finite_cor] = -Inf
  best = max.col(cors, ties.method = "first")
  has = rowSums(finite_cor) > 0
  lag = rep(NA_real_, n)
  rcor = rep(NA_real_, n)
  lag[has] = outlook_lags_hr[best[has]]
  rcor[has] = cors[cbind(which(has), best[has])]
  # max.col wrote -Inf back into cors; restore the correlation that was chosen
  rcor[has] = pmax(rcor[has], -1)
  list(lag = lag, cor = rcor, ok = has)
}
outlook_design = function(kp, i, with_future) {
  info = choose_kp_lag(kp, i)
  recent = matrix(NA_real_, length(i), n_kp_hist)
  for (a in seq_len(n_kp_hist)) recent[, a] = kp[i - n_kp_hist + a]
  recur = matrix(NA_real_, length(i), outlook_days)
  fut = NULL
  use = info$ok
  if (any(use)) {
    recur[use, ] = day_means(kp, i[use], -info$lag[use] / kp_step_hr, outlook_days)
    if (with_future) fut = day_means(kp, i, 0, outlook_days)
  } else if (with_future) {
    fut = matrix(NA_real_, length(i), outlook_days)
  }
  ok = info$ok & finite_rows(recent) & finite_rows(recur) & is.finite(info$cor)
  if (with_future) ok = ok & finite_rows(fut)
  list(X = cbind(recent, recur, info$cor), Y = fut, ok = ok,
       lag = info$lag, cor = info$cor, recur = recur)
}
cat_outlook_skill = function(skill) {
  if (is.null(skill) || !nrow(skill)) return(invisible())
  show = skill[skill$day %in% c(1L, 3L, 7L, 14L, 21L, 27L), ]
  cat("27-day outlook, issues after the training cutoff\n")
  for (i in seq_len(nrow(show))) {
    cat(sprintf("  day %2d  n %4d  RMSE %.3f  recurrence %.3f  persistence %.3f\n",
                show$day[i], show$n[i], show$RMSE[i], show$RMSE_recur[i], show$RMSE_pers[i]))
  }
}
fit_outlook = NULL
outlook_skill = NULL
if (file.exists(outlook_coef_path)) {
  env = new.env()
  load(outlook_coef_path, envir = env)
  if (isTRUE(all.equal(env$outlook_sig, outlook_sig))) {
    fit_outlook = env$fit_outlook
    outlook_skill = env$outlook_skill
  }
}
if (is.null(fit_outlook)) {
  cat("Fitting the 27-day Kp outlook\n")
  i_min = max(outlook_lags_hr) / kp_step_hr + outlook_match_bins
  i_out = seq.int(i_min, N_kp - outlook_days * bins_per_day)
  i_out = i_out[issue_time(i_out) + outlook_days * 24 * 3600 <= train_end]
  des = outlook_design(kp_value, i_out, TRUE)
  keep = des$ok
  cat(sprintf("  training issues %d\n", sum(keep)))
  X = des$X[keep, , drop = FALSE]
  Y = des$Y[keep, , drop = FALSE]
  sc = standardize_fit(X)
  fit = ols_coef_from_svd(cbind(1, apply_scale(X, sc$center, sc$scale)), Y, ols_tol)
  fit_outlook = list(
    B = fit$B, center = sc$center, scale = sc$scale,
    rank = fit$rank, n_train = nrow(X))
  cat(sprintf("  rank %d / %d\n", fit$rank, ncol(X) + 1L))
  i_te = seq.int(i_min, N_kp - outlook_days * bins_per_day)
  i_te = i_te[issue_time(i_te) > train_end]
  outlook_skill = NULL
  if (length(i_te)) {
    te = outlook_design(kp_value, i_te, TRUE)
    m = te$ok
    yhat = matrix(NA_real_, length(i_te), outlook_days)
    yhat[m, ] = cbind(1, apply_scale(te$X[m, , drop = FALSE], sc$center, sc$scale)) %*% t(fit$B)
    pers = day_means(kp_value, i_te, -bins_per_day, 1L)[, 1]
    outlook_skill = do.call(rbind, lapply(seq_len(outlook_days), function(d) {
      md = m & is.finite(te$Y[, d]) & is.finite(yhat[, d]) & is.finite(te$recur[, d]) & is.finite(pers)
      data.frame(
        day = d, n = sum(md),
        RMSE = sqrt(mean((te$Y[md, d] - yhat[md, d])^2)),
        RMSE_recur = sqrt(mean((te$Y[md, d] - te$recur[md, d])^2)),
        RMSE_pers = sqrt(mean((te$Y[md, d] - pers[md])^2)),
        stringsAsFactors = FALSE
      )
    }))
  }
  save(fit_outlook, outlook_skill, outlook_sig, file = outlook_coef_path)
} else {
  cat(sprintf("Loaded cached 27-day outlook (%s)\n", outlook_coef_path))
}
cat_outlook_skill(outlook_skill)

# Historical OMNI cycles, 1995 through the training cutoff. A cycle may be
# used only after the 27 days of Kp that followed it have ended. The query
# is the last 27 days of in-situ measurements, and a day counts when at
# least 12 of its hours have a real sample.
archive_cycle_days = 27L
archive_min_hours = 12L
archive_top_k = 10L
archive_min_past = 80L
archive_env = new.env(parent = emptyenv())

archive_hour_index = function(origin, hour_num, n_hours) {
  ix = (as.numeric(hour_num) - origin) / 3600 + 1
  r = round(ix)
  bad = !is.finite(ix) | abs(ix - r) > 1e-3 | r < 1 | r > n_hours
  out = as.integer(r)
  out[bad] = NA_integer_
  out
}

archive_cycle_means = function(mat, last_idx, strict) {
  n = length(last_idx)
  days = vector("list", archive_cycle_days)
  keep = rep(TRUE, n)
  for (j in 0:(archive_cycle_days - 1L)) {
    acc = matrix(0, n, length(omni_features))
    n_ok = rep(0L, n)
    for (h in 0:23) {
      ix = last_idx - j * 24L - h
      known = !is.na(ix) & ix >= 1L & ix <= nrow(mat)
      good = rep(FALSE, n)
      if (any(known)) {
        rows = mat[ix[known], , drop = FALSE]
        row_ok = rowSums(is.finite(rows)) == ncol(rows)
        pos = which(known)
        good[pos[row_ok]] = TRUE
        if (any(row_ok))
          acc[pos[row_ok], ] = acc[pos[row_ok], ] + rows[row_ok, , drop = FALSE]
      }
      n_ok = n_ok + good
    }
    enough = if (strict) n_ok == 24L else n_ok >= archive_min_hours
    keep = keep & enough
    days[[archive_cycle_days - j]] = acc / pmax(n_ok, 1L)
  }
  list(days = days, keep = keep)
}

archive_pack = function(days, keep) {
  packed = vector("list", length(omni_features))
  for (f in seq_along(omni_features)) {
    M = matrix(NA_real_, length(keep), archive_cycle_days)
    for (j in seq_len(archive_cycle_days)) M[, j] = days[[j]][, f]
    packed[[f]] = M[keep, , drop = FALSE]
  }
  packed
}

archive_cors = function(feat, query) {
  n = nrow(feat[[1]])
  cors = rep(0, n)
  used = 0L
  for (f in seq_along(omni_features)) {
    q = query[, f]
    if (sum(is.finite(q)) < archive_cycle_days || sd(q) == 0) next
    M = feat[[f]]
    q = q - mean(q)
    qn = sqrt(sum(q * q))
    L = M - rowMeans(M)
    ln = sqrt(rowSums(L * L))
    c = as.numeric(L %*% q) / (ln * qn)
    c[!is.finite(c)] = 0
    cors = cors + c
    used = used + 1L
  }
  if (!used) return(rep(NA_real_, n))
  cors / used
}

archive_kp_after = function(end_time) {
  i0 = match(as.numeric(end_time), as.numeric(kp_time))
  n = length(end_time)
  out = matrix(NA_real_, n, outlook_days)
  known = !is.na(i0)
  for (d in seq_len(outlook_days)) {
    acc = rep(0, n)
    n_ok = rep(0L, n)
    for (a in 0:7) {
      ix = i0 + (d - 1L) * 8L + a
      good = known & !is.na(ix) & ix >= 1L & ix <= length(kp_value)
      acc[good] = acc[good] + kp_value[ix[good]]
      n_ok = n_ok + good
    }
    out[n_ok == 8L, d] = acc[n_ok == 8L] / 8
  }
  out
}

archive_library = function() {
  if (exists("library", envir = archive_env, inherits = FALSE))
    return(get("library", envir = archive_env, inherits = FALSE))
  cache_path = file.path(out_dir, "VARX_Kp_cycle_library.RData")
  sig = list(
    omni = file_stamp(omni_path),
    train_end = as.numeric(train_end),
    days = archive_cycle_days)
  if (file.exists(cache_path)) {
    env = new.env()
    load(cache_path, envir = env)
    if (isTRUE(all.equal(env$sig, sig))) {
      assign("library", env$library, envir = archive_env)
      return(env$library)
    }
  }
  cat("Building OMNI cycle library, 1995 through the training cutoff\n")
  omni_use = OMNI_mat
  for (nm in omni_features)
    omni_use[!is.finite(omni_use[, nm]) | omni_use[, nm] >= sentinel[nm], nm] = NA
  origin = as.numeric(omni_time[1])
  library = vector("list", 8L)
  names(library) = sprintf("%02d", seq(0, 21, by = kp_step_hr))
  for (ph in seq(0, 21, by = kp_step_hr)) {
    start = as.POSIXct(sprintf("1995-02-01 %02d:00:00", ph), tz = "UTC")
    ends = seq(start, train_end, by = "24 hours")
    ends = ends[ends <= train_end]
    last_idx = archive_hour_index(origin, as.numeric(ends) - 3600, nrow(omni_use))
    built = archive_cycle_means(omni_use, last_idx, TRUE)
    outcome = archive_kp_after(ends)
    keep = built$keep & rowSums(is.finite(outcome)) == outlook_days
    library[[sprintf("%02d", ph)]] = list(
      end = ends[keep],
      feat = archive_pack(built$days, keep),
      outcome = outcome[keep, , drop = FALSE])
    cat(sprintf("  phase %02d UTC  %d cycles\n", ph, sum(keep)))
  }
  save(library, sig, file = cache_path)
  assign("library", library, envir = archive_env)
  library
}

archive_query = function(tau) {
  if (!exists("obs_match", envir = archive_env, inherits = FALSE)) {
    mat = obs_h$mat
    drop = obs_h$n_native == 0 | obs_h$time > last_obs
    mat[drop, ] = NA
    assign("obs_match", mat, envir = archive_env)
    assign("obs_origin", as.numeric(obs_h$time[1]), envir = archive_env)
    omni_use = OMNI_mat
    for (nm in omni_features)
      omni_use[!is.finite(omni_use[, nm]) | omni_use[, nm] >= sentinel[nm], nm] = NA
    assign("omni_match", omni_use, envir = archive_env)
    assign("omni_origin", as.numeric(omni_time[1]), envir = archive_env)
  }
  obs_mat = get("obs_match", envir = archive_env, inherits = FALSE)
  obs_origin = get("obs_origin", envir = archive_env, inherits = FALSE)
  omni_use = get("omni_match", envir = archive_env, inherits = FALSE)
  omni_origin = get("omni_origin", envir = archive_env, inherits = FALSE)
  # Oldest hour first. A missing in-situ hour uses OMNI when that hour is
  # already in the file, so a gap in RTSW does not drop the whole match.
  hours = as.numeric(tau) - (archive_cycle_days * 24L):1L * 3600
  block = matrix(NA_real_, length(hours), length(omni_features))
  ix_obs = archive_hour_index(obs_origin, hours, nrow(obs_mat))
  got = !is.na(ix_obs)
  if (any(got)) block[got, ] = obs_mat[ix_obs[got], , drop = FALSE]
  miss = rowSums(is.finite(block)) < ncol(block)
  ix_om = archive_hour_index(omni_origin, hours, nrow(omni_use))
  fill = miss & !is.na(ix_om)
  if (any(fill)) block[fill, ] = omni_use[ix_om[fill], , drop = FALSE]
  built = archive_cycle_means(block, nrow(block), FALSE)
  if (!all(built$keep)) return(NULL)
  query = matrix(NA_real_, archive_cycle_days, length(omni_features))
  for (j in seq_len(archive_cycle_days)) query[j, ] = built$days[[j]][1, ]
  query
}

archive_top10 = function(tau) {
  library = archive_library()
  n_lead = outlook_days * 8L
  out = matrix(NA_real_, length(tau), outlook_days)
  best_bin = matrix(NA_real_, length(tau), n_lead)
  best_end = rep(NA_real_, length(tau))
  hour = as.integer(format(tau, "%H", tz = "UTC"))
  for (i in seq_along(tau)) {
    query = archive_query(tau[i])
    if (is.null(query)) next
    lib = library[[sprintf("%02d", hour[i])]]
    eligible = lib$end + outlook_days * 86400 <= tau[i]
    if (!any(eligible)) next
    cors = archive_cors(lapply(lib$feat, function(M) M[eligible, , drop = FALSE]), query)
    if (!any(is.finite(cors))) next
    take = min(archive_top_k, sum(is.finite(cors)))
    ord = order(cors, decreasing = TRUE)[seq_len(take)]
    w = pmax(cors[ord], 0)
    if (sum(w) == 0) w = rep(1, length(ord))
    y = lib$outcome[eligible, , drop = FALSE]
    out[i, ] = colSums(w * y[ord, , drop = FALSE]) / sum(w)
    s = lib$end[eligible][ord[1]]
    i0 = match(as.numeric(s), as.numeric(kp_time))
    if (!is.na(i0) && i0 + n_lead - 1L <= length(kp_value))
      best_bin[i, ] = kp_value[i0:(i0 + n_lead - 1L)]
    best_end[i] = as.numeric(s)
    if (i %% 500L == 0L) cat(sprintf("  matched %d / %d issues\n", i, length(tau)))
  }
  list(top = out, best_bin = best_bin, best_end = best_end)
}

# y = w * current + (1 - w) * archive. w is fit on earlier issues only.
archive_mix = function(Y, current, top) {
  yhat = current
  w_next = rep(1, ncol(current))
  for (d in seq_len(ncol(current))) {
    y = Y[, d]
    c = current[, d]
    t = top[, d]
    num = 0
    den = 0
    n_ok = 0L
    for (i in seq_along(c)) {
      if (n_ok >= archive_min_past && den >= 1e-8 &&
          is.finite(c[i]) && is.finite(t[i])) {
        ww = max(0, min(1, num / den))
        yhat[i, d] = ww * c[i] + (1 - ww) * t[i]
      }
      if (is.finite(y[i]) && is.finite(c[i]) && is.finite(t[i])) {
        gap = c[i] - t[i]
        num = num + (y[i] - t[i]) * gap
        den = den + gap * gap
        n_ok = n_ok + 1L
      }
    }
    if (n_ok >= archive_min_past && den >= 1e-8)
      w_next[d] = max(0, min(1, num / den))
  }
  list(yhat = yhat, w_next = w_next)
}

score_outlook_heldout = function() {
  # Predictors come from Kp already known at issue time. Day d is scored
  # once its eight bins have ended, including bins appended after Kp.csv.
  # Training still uses the file only.
  i_min = max(outlook_lags_hr) / kp_step_hr + outlook_match_bins
  N_live = length(kp_time_live)
  i_last = N_live - bins_per_day
  if (i_last < i_min) stop("Not enough Kp to score the 27-day outlook.")
  i_te = seq.int(i_min, i_last)
  tau = kp_time_live[i_te] + kp_step_hr * 3600
  keep_issue = tau > train_end
  i_te = i_te[keep_issue]
  tau = tau[keep_issue]
  te = outlook_design(kp_value_live, i_te, FALSE)
  m = te$ok
  i_te = i_te[m]
  tau = tau[m]
  yhat = cbind(1, apply_scale(te$X[m, , drop = FALSE], fit_outlook$center, fit_outlook$scale)) %*% t(fit_outlook$B)
  Y = day_means(kp_value_live, i_te, 0, outlook_days)
  y_kp = yhat
  y_archive = matrix(NA_real_, nrow(yhat), outlook_days)
  y_match_bin = matrix(NA_real_, nrow(yhat), outlook_days * 8L)
  best_end = rep(NA_real_, nrow(yhat))
  w_next = rep(1, outlook_days)
  if (exists("obs_h")) {
    cat("Matching the last 27 days of in-situ solar wind to OMNI cycles\n")
    hit = archive_top10(tau)
    y_archive = hit$top
    y_match_bin = hit$best_bin
    best_end = hit$best_end
    mixed = archive_mix(Y, y_kp, y_archive)
    yhat = mixed$yhat
    w_next = mixed$w_next
    cat(sprintf("Archive weight on the current regression, day 1 %.2f, day 7 %.2f, day 27 %.2f\n",
                w_next[1], w_next[7], w_next[27]))
  }
  recur = te$recur[m, , drop = FALSE]
  # Each lead is one 3-hour Kp bin. The reference is the Kp that followed
  # the closest aligned OMNI cycle, not persistence. Inside a day the eight
  # forecast bins follow the previous rotation and average to the daily value.
  n_iss = length(i_te)
  n_lead = outlook_days * 8L
  leads_hr = seq(kp_step_hr, outlook_days * 24L, by = kp_step_hr)
  lag_bins = as.integer(round(te$lag[m] / kp_step_hr))
  pers = kp_value_live[i_te]
  y_true_bin = matrix(NA_real_, n_iss, n_lead)
  y_recur_bin = y_true_bin
  y_pred_bin = y_true_bin
  for (d in seq_len(outlook_days)) {
    cols = (d - 1L) * 8L + seq_len(8L)
    for (k in seq_len(8L)) {
      lead = cols[k]
      true_ix = i_te + lead
      rot_ix = i_te - lag_bins + lead
      ok_true = true_ix >= 1L & true_ix <= N_live
      ok_rot = rot_ix >= 1L & rot_ix <= N_live
      y_true_bin[ok_true, lead] = kp_value_live[true_ix[ok_true]]
      y_recur_bin[ok_rot, lead] = kp_value_live[rot_ix[ok_rot]]
    }
    mu = rowMeans(y_recur_bin[, cols, drop = FALSE])
    for (k in seq_len(8L))
      y_pred_bin[, cols[k]] = yhat[, d] + (y_recur_bin[, cols[k]] - mu)
  }
  bin_row = function(lead_hr, input, y, yhat_v, ref) {
    md = is.finite(y) & is.finite(yhat_v) & is.finite(ref)
    n = sum(md)
    rm = if (n) sqrt(mean((y[md] - yhat_v[md])^2)) else NA_real_
    rr = if (n) sqrt(mean((y[md] - ref[md])^2)) else NA_real_
    data.frame(
      horizon_hr = outlook_days * 24L,
      lead_end_hr = lead_hr,
      input = input,
      n = n,
      RMSE = rm,
      MAE = if (n) mean(abs(y[md] - yhat_v[md])) else NA_real_,
      Correlation = if (n > 1L) suppressWarnings(cor(y[md], yhat_v[md])) else NA_real_,
      RMSE_reference = rr,
      skill_RMSE = if (!n || !is.finite(rr) || rr == 0) NA_real_ else 1 - rm / rr,
      stringsAsFactors = FALSE
    )
  }
  metrics = do.call(rbind, lapply(seq_len(n_lead), function(j) {
    ref = y_match_bin[, j]
    rbind(
      bin_row(leads_hr[j], "Best matched cycle", y_true_bin[, j], ref, ref),
      bin_row(leads_hr[j], "Forecast", y_true_bin[, j], y_pred_bin[, j], ref)
    )
  }))
  list(tau = tau, y_true = Y, y_pred = yhat, y_kp = y_kp, y_archive = y_archive,
       w_next = w_next, y_recur = recur, pers = pers, metrics = metrics,
       y_true_bin = y_true_bin, y_pred_bin = y_pred_bin, y_recur_bin = y_recur_bin,
       y_match_bin = y_match_bin, best_end = best_end, leads_hr = leads_hr)
}
cat("Scoring the 27-day outlook where those days have already been observed\n")
outlook_eval = score_outlook_heldout()
outlook_metrics_path = file.path(out_dir, "VARX_Kp_baseline_27d_metrics.csv")
outlook_series_path = file.path(out_dir, "VARX_Kp_baseline_27d.csv")
write.csv(outlook_eval$metrics, outlook_metrics_path, row.names = FALSE)
outlook_series = do.call(rbind, lapply(seq_along(outlook_eval$leads_hr), function(j) {
  keep = is.finite(outlook_eval$y_true_bin[, j])
  if (!any(keep)) return(NULL)
  lead = outlook_eval$leads_hr[j]
  data.frame(
    issue_time = fmt_utc(outlook_eval$tau[keep]),
    lead_end_hr = lead,
    target_end = fmt_utc(outlook_eval$tau[keep] + lead * 3600),
    y_true = outlook_eval$y_true_bin[keep, j],
    y_pred = outlook_eval$y_pred_bin[keep, j],
    y_recur = outlook_eval$y_recur_bin[keep, j],
    y_matched = outlook_eval$y_match_bin[keep, j],
    stringsAsFactors = FALSE
  )
}))
write.csv(outlook_series, outlook_series_path, row.names = FALSE)
cat(sprintf("Wrote %s\n", outlook_metrics_path))
cat(sprintf("Wrote %s\n", outlook_series_path))
for (lead in c(24L, 168L, 648L)) {
  sub = outlook_eval$metrics[outlook_eval$metrics$lead_end_hr == lead, ]
  fc = sub[sub$input == "Forecast", ]
  ref = sub[sub$input == "Best matched cycle", ]
  scored = which(is.finite(outlook_eval$y_true_bin[, lead / kp_step_hr]))
  last_target = if (length(scored))
    fmt_utc(outlook_eval$tau[max(scored)] + lead * 3600) else "none"
  cat(sprintf("  lead %3dh  n %d  through %s  forecast %.3f  best cycle %.3f  skill %+.3f\n",
              lead, fc$n, last_target, fc$RMSE, ref$RMSE, fc$skill_RMSE))
}
if (score_outlook_only) {
  source(file.path(repo_root, "Rcode", "Main", "Write_VARX_Pages.R"))
  section = paste0("<!-- outlook-27-start -->", outlook_baseline_section(outlook_eval), "<!-- outlook-27-end -->")
  baseline_page = file.path(out_dir, "VARX_Kp_baseline.html")
  page = readChar(baseline_page, file.info(baseline_page)$size, useBytes = TRUE)
  Encoding(page) = "UTF-8"
  if (grepl("<!-- outlook-27-start -->", page, fixed = TRUE)) {
    page = sub("(?s)<!-- outlook-27-start -->.*<!-- outlook-27-end -->", section, page, perl = TRUE)
  } else {
    page = sub("</main>", paste0(section, "</main>"), page, fixed = TRUE)
  }
  write_html(baseline_page, page)
  publish_dir = Sys.getenv("KP_PUBLISH_DIR", file.path(repo_root, "publish"))
  dir.create(publish_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(baseline_page, file.path(publish_dir, "VARX_Kp_baseline.html"), overwrite = TRUE)
  file.copy(outlook_metrics_path, file.path(publish_dir, basename(outlook_metrics_path)), overwrite = TRUE)
  file.copy(outlook_series_path, file.path(publish_dir, basename(outlook_series_path)), overwrite = TRUE)
  cat(sprintf("Updated the 27-day section of %s\n", baseline_page))
  quit(save = "no", status = 0)
}

# ---- latest issue whose 24h window exists in both obs and model -------------
# One row per source. Exog is feature-major, oldest hour to newest.
build_driver = function(i_live, input) {
  tau = kp_time_live[i_live] + kp_step_hr * 3600
  hours = as.POSIXct(as.numeric(tau) - (history_hr:1L) * 3600,
                     origin = "1970-01-01", tz = "UTC")
  if (input == "Earth_obs") {
    used = causal_solar_wind(hours, obs_h$time, obs_h$mat, kp_step_hr)
    src_n = obs_h$n_native
    src_t = obs_h$time
  } else {
    raw = causal_solar_wind(hours, model_h$time, model_h$mat, kp_step_hr)
    if (!all(is.finite(raw))) return(NULL)
    ab = affine_updated(tau - history_hr * 3600)
    used = sweep(raw, 2, ab$b, "*")
    used = sweep(used, 2, ab$a, "+")
    src_n = model_h$n_native
    src_t = model_h$time
  }
  if (!all(is.finite(used))) return(NULL)
  ix = match(as.numeric(hours), as.numeric(src_t))
  n_native = integer(length(hours))
  known = !is.na(ix)
  n_native[known] = src_n[ix[known]]
  list(input = input, tau = tau, hours = hours, used = used, n_native = n_native)
}

inputs = c("Earth_obs", "model_affine")
i_live = length(kp_time_live)
drivers_by_input = NULL
while (i_live >= n_kp_hist) {
  hist_ix = (i_live - n_kp_hist + 1L):i_live
  hist_ok = all(is.finite(kp_value_live[hist_ix])) &&
    all(abs(diff(as.numeric(kp_time_live[hist_ix])) - kp_step_hr * 3600) < 1)
  built = if (!hist_ok) NULL else lapply(inputs, function(input) build_driver(i_live, input))
  if (!is.null(built) && !any(vapply(built, is.null, logical(1)))) {
    drivers_try = setNames(built, inputs)
    kp_row = matrix(kp_value_live[hist_ix], nrow = 1)
    tau_try = drivers_try[[1]]$tau
    lag_ok = TRUE
    for (input in inputs) {
      driver = drivers_try[[input]]
      exog = matrix(as.vector(driver$used), nrow = 1)
      if (input == "Earth_obs") {
        rb = recur_batch(tau_try, i_live, kp_value_live, kp_row, exog,
                         obs_h$time, obs_h$mat[, vs_col],
                         obs_h$time, obs_h$mat[, vs_col],
                         speed_latest = last_obs)
      } else {
        ab_now = affine_updated(tau_try - history_hr * 3600)
        rb = recur_batch(tau_try, i_live, kp_value_live, kp_row, exog,
                         model_h$time, model_h$mat[, vs_col],
                         model_h$time, model_h$mat[, vs_col],
                         scale_vs = matrix(c(ab_now$a[vs_col], ab_now$b[vs_col]), nrow = 1))
      }
      if (!isTRUE(rb$ok)) {
        cat(sprintf("Lag features are still incomplete for %s at %s. Trying the previous issue.\n",
                    input, fmt_utc(tau_try)))
        lag_ok = FALSE
        break
      }
    }
    if (lag_ok) {
      drivers_by_input = drivers_try
      break
    }
  }
  i_live = i_live - 1L
}
if (is.null(drivers_by_input))
  stop("No issue time has a complete 24h window in both RTSW obs and the RTSW model.")

kp_hist = kp_value_live[(i_live - n_kp_hist + 1L):i_live]
generated_at = fmt_utc(Sys.time())
tau_live = drivers_by_input[[1]]$tau
forecast_rows = list()
driver_rows = list()
recur_lag = list()
xgb_X = NULL
leads = seq(kp_step_hr, 24L, by = kp_step_hr)
for (input in inputs) {
  driver = drivers_by_input[[input]]
  exog = matrix(as.vector(driver$used), nrow = 1)
  kp_row = matrix(kp_hist, nrow = 1)
  if (input == "Earth_obs") {
    rb = recur_batch(tau_live, i_live, kp_value_live, kp_row, exog,
                     obs_h$time, obs_h$mat[, vs_col],
                     obs_h$time, obs_h$mat[, vs_col],
                     speed_latest = last_obs)
  } else {
    ab_now = affine_updated(tau_live - history_hr * 3600)
    rb = recur_batch(tau_live, i_live, kp_value_live, kp_row, exog,
                     model_h$time, model_h$mat[, vs_col],
                     model_h$time, model_h$mat[, vs_col],
                     scale_vs = matrix(c(ab_now$a[vs_col], ab_now$b[vs_col]), nrow = 1))
  }
  if (!isTRUE(rb$ok))
    stop("Lag features are incomplete for ", input, " at ", fmt_utc(tau_live))
  if (input == "Earth_obs") xgb_X = rb$pred
  yhat = as.numeric(recur_predict(rb$pred))
  one = function(x) unname(rep(as.vector(x), length(leads)))
  forecast_rows[[length(forecast_rows) + 1L]] = data.frame(
    generated_at = one(generated_at),
    issue_time = one(fmt_utc(tau_live)),
    input = one(input),
    horizon_hr = one(24L),
    lead_end_hr = leads,
    target_bin_start = fmt_utc(tau_live + (leads - kp_step_hr) * 3600),
    target_bin_end = fmt_utc(tau_live + leads * 3600),
    y_pred = yhat,
    y_persistence = one(kp_value_live[i_live]),
    n_native_hours = one(sum(driver$n_native > 0)),
    n_filled_hours = one(sum(driver$n_native == 0)),
    stringsAsFactors = FALSE
  )
  recur_lag[[input]] = rb$lag
  driver_rows[[input]] = data.frame(
    issue_time = fmt_utc(tau_live),
    input = input,
    hour_start = fmt_utc(driver$hours),
    n_native = driver$n_native,
    driver$used,
    stringsAsFactors = FALSE
  )
}
forecast = do.call(rbind, forecast_rows)
drivers = do.call(rbind, driver_rows)

# The next 27 days uses the 24 h VARX, stepped one day at a time. The first
# day is the forecast already issued above. Later days read OMNI from the
# closest aligned cycle, and the Kp just predicted becomes the next history.
roll_outlook = function(y0, past_time, past_mat, transform) {
  n_lead = outlook_days * 8L
  y = rep(NA_real_, n_lead)
  y[seq_len(8L)] = y0
  if (!is.finite(live_cycle$best_end) || i_live <= 260L) return(y)
  omni_use = get("omni_match", envir = archive_env, inherits = FALSE)
  omni_origin = get("omni_origin", envir = archive_env, inherits = FALSE)
  kp_back = 260L
  hnum = seq(as.numeric(tau_live) - 40 * 86400,
             as.numeric(tau_live) + outlook_days * 86400, by = 3600)
  sw = matrix(NA_real_, length(hnum), length(omni_features))
  past = hnum < as.numeric(tau_live)
  ix = match(hnum[past], as.numeric(past_time))
  ok = !is.na(ix)
  if (any(ok)) {
    rows = past_mat[ix[ok], , drop = FALSE]
    if (!is.null(transform))
      rows = sweep(sweep(rows, 2, transform$b, "*"), 2, transform$a, "+")
    sw[which(past)[ok], ] = rows
  }
  fut = !past
  ix_om = archive_hour_index(
    omni_origin, live_cycle$best_end + hnum[fut] - as.numeric(tau_live), nrow(omni_use))
  ok_om = !is.na(ix_om)
  if (any(ok_om)) sw[which(fut)[ok_om], ] = omni_use[ix_om[ok_om], , drop = FALSE]
  # A short hole in the working wind table is filled from the hours on either
  # side. The hours after the issue are the chosen OMNI cycle, already known.
  filled = causal_solar_wind(hnum, hnum, sw, kp_step_hr)
  sw = matrix(filled, nrow = length(hnum), ncol = ncol(sw))
  colnames(sw) = omni_features
  kp_work = c(kp_value_live[(i_live - kp_back + 1L):i_live], y)
  for (day in 1:(outlook_days - 1L)) {
    i_now = kp_back + day * 8L
    tau_d = tau_live + day * 24 * 3600
    hour_start = as.numeric(tau_d) - (history_hr:1L) * 3600
    ixh = match(hour_start, hnum)
    if (anyNA(ixh)) {
      cat(sprintf("27-day roll stopped on day %d: a driver hour is outside the wind table.\n", day + 1L))
      break
    }
    block = sw[ixh, , drop = FALSE]
    if (!all(is.finite(block))) {
      cat(sprintf("27-day roll stopped on day %d: %d driver values are missing.\n",
                  day + 1L, sum(!is.finite(block))))
      break
    }
    rb = recur_batch(tau_d, i_now, kp_work,
                     matrix(kp_work[(i_now - n_kp_hist + 1L):i_now], nrow = 1),
                     matrix(as.vector(block), nrow = 1),
                     hnum, sw[, vs_col], hnum, sw[, vs_col])
    if (!isTRUE(rb$ok)) {
      # Speed match is usable, but no earlier Kp window has the same shape,
      # so there is no rise or fall to add. Continue with that growth at zero.
      if (!all(is.finite(rb$pred[, seq_len(138), drop = FALSE]))) {
        cat(sprintf("27-day roll stopped on day %d: the solar-rotation match is incomplete.\n", day + 1L))
        break
      }
      rb$pred[, 139:146] = 0
    }
    yhat = as.numeric(recur_predict(rb$pred))
    if (!all(is.finite(yhat))) break
    slots = day * 8L + seq_len(8L)
    kp_work[i_now + seq_len(8L)] = yhat
    y[slots] = yhat
  }
  y
}
ab_live = affine_updated(tau_live - history_hr * 3600)
live_cycle = archive_top10(tau_live)
f_obs = forecast[forecast$input == "Earth_obs", ]
f_mod = forecast[forecast$input == "model_affine", ]
f_obs = f_obs[order(f_obs$lead_end_hr), ]
f_mod = f_mod[order(f_mod$lead_end_hr), ]
y_insitu = roll_outlook(f_obs$y_pred, obs_h$time, obs_h$mat, NULL)
y_awsom = roll_outlook(f_mod$y_pred, model_h$time, model_h$mat,
                       list(a = ab_live$a, b = ab_live$b))
n_bin = outlook_days * 8L
y_matched = as.numeric(live_cycle$best_bin[1, ])
if (length(y_matched) != n_bin) y_matched = rep(NA_real_, n_bin)
outlook = data.frame(
  day = rep(seq_len(outlook_days), each = 8L),
  day_end = fmt_utc(tau_live + seq_len(n_bin) * kp_step_hr * 3600),
  y_insitu = y_insitu,
  y_awsom = y_awsom,
  y_boost = y_insitu,
  y_matched = y_matched,
  stringsAsFactors = FALSE
)
cat(sprintf("27-day roll  in-situ %d/216  AWSoM %d/216  matched cycle %s\n",
            sum(is.finite(y_insitu)), sum(is.finite(y_awsom)),
            if (is.finite(live_cycle$best_end[1]))
              fmt_utc(as.POSIXct(live_cycle$best_end[1], origin = "1970-01-01", tz = "UTC"))
            else "none"))

source(file.path(repo_root, "Rcode", "Main", "Write_VARX_Pages.R"))
source(file.path(repo_root, "Rcode", "Main", "Issue_Best.R"))

score_baseline = function() {
  cal_end = train_end + cal_days * 86400
  n24 = 24L / kp_step_hr
  N_live = length(kp_time_live)
  i_cand = seq.int(n_kp_hist, N_live - n24)
  bad_step = abs(diff(as.numeric(kp_time_live)) - kp_step_hr * 3600) >= 1
  cs = c(0L, cumsum(bad_step))
  i_cand = i_cand[(cs[i_cand + n24] - cs[i_cand]) == 0]
  tau = kp_time_live[i_cand] + kp_step_hr * 3600
  after_cal = as.numeric(tau) - history_hr * 3600 >= as.numeric(cal_end)
  i_cand = i_cand[after_cal]
  tau = tau[after_cal]
  if (!length(i_cand))
    stop("No issue time is after the calibration window with a 24 h target still on file.")

  ages = history_hr:1L
  hour_num = outer(as.numeric(tau), ages * 3600, "-")
  raw_obs = causal_solar_wind(hour_num, obs_h$time, obs_h$mat, kp_step_hr)
  raw_mod = causal_solar_wind(hour_num, model_h$time, model_h$mat, kp_step_hr)
  obs_ok = matrix(rowSums(is.finite(raw_obs)) == ncol(raw_obs), nrow(hour_num), ncol(hour_num))
  mod_ok = matrix(rowSums(is.finite(raw_mod)) == ncol(raw_mod), nrow(hour_num), ncol(hour_num))
  # In-situ hours after the last real sample are only the short hold. Do not
  # score an issue whose driver reaches past that last real hour.
  obs_time = matrix(hour_num, nrow(hour_num), ncol(hour_num))
  live_tail = rowSums(obs_time > as.numeric(last_obs)) > 0
  keep = rowSums(obs_ok) == ncol(hour_num) & rowSums(mod_ok) == ncol(hour_num) & !live_tail
  i_keep = i_cand[keep]
  tau_keep = tau[keep]
  raw_obs = matrix(raw_obs, nrow(hour_num) * ncol(hour_num), ncol(raw_obs))
  raw_mod = matrix(raw_mod, nrow(hour_num) * ncol(hour_num), ncol(raw_mod))
  # Rebuild as issue x hour after the keep filter. causal_solar_wind stored
  # hours in column-major order, one column of hour_num at a time.
  take = function(raw) {
    cols = vector("list", length(omni_features))
    for (f in seq_along(omni_features)) {
      m = matrix(raw[, f], nrow(hour_num), ncol(hour_num))
      cols[[f]] = m[keep, , drop = FALSE]
    }
    cols
  }
  obs_cols = take(raw_obs)
  mod_cols = take(raw_mod)
  if (!length(i_keep))
    stop("No issue time has a complete 24 h window in both Earth_obs and Earth_model.")

  kp_hist_mat = vapply(seq_len(n_kp_hist), function(j) {
    kp_value_live[i_keep - n_kp_hist + j]
  }, numeric(length(i_keep)))
  y_true = vapply(seq_len(n24), function(k) kp_value_live[i_keep + k], numeric(length(i_keep)))
  finite_row = rowSums(is.finite(kp_hist_mat)) == n_kp_hist & rowSums(is.finite(y_true)) == n24
  i_keep = i_keep[finite_row]
  tau_keep = tau_keep[finite_row]
  obs_cols = lapply(obs_cols, function(m) m[finite_row, , drop = FALSE])
  mod_cols = lapply(mod_cols, function(m) m[finite_row, , drop = FALSE])
  kp_hist_mat = kp_hist_mat[finite_row, , drop = FALSE]
  y_true = y_true[finite_row, , drop = FALSE]
  if (!length(i_keep)) stop("Scored issues have non-finite Kp history or targets.")

  exog_from_cols = function(cols) {
    ex = matrix(NA_real_, nrow(cols[[1]]), ncol(cols[[1]]) * length(cols))
    col0 = 1L
    for (f in seq_along(cols)) {
      ex[, col0:(col0 + ncol(cols[[f]]) - 1L)] = cols[[f]]
      col0 = col0 + ncol(cols[[f]])
    }
    ex
  }
  ex_obs = exog_from_cols(obs_cols)
  ex_mod = matrix(NA_real_, length(i_keep), history_hr * length(omni_features))
  scale_vs = matrix(NA_real_, length(i_keep), 2)
  for (i in seq_len(length(i_keep))) {
    ab = affine_updated(tau_keep[i] - history_hr * 3600)
    scale_vs[i, ] = c(ab$a[vs_col], ab$b[vs_col])
    col0 = 1L
    for (f in seq_along(omni_features)) {
      ex_mod[i, col0:(col0 + history_hr - 1L)] = ab$a[f] + ab$b[f] * mod_cols[[f]][i, ]
      col0 = col0 + history_hr
    }
  }
  rb_obs = recur_batch(tau_keep, i_keep, kp_value_live, kp_hist_mat, ex_obs,
                       obs_h$time, obs_h$mat[, vs_col],
                       obs_h$time, obs_h$mat[, vs_col],
                       speed_latest = last_obs)
  rb_mod = recur_batch(tau_keep, i_keep, kp_value_live, kp_hist_mat, ex_mod,
                       model_h$time, model_h$mat[, vs_col],
                       model_h$time, model_h$mat[, vs_col],
                       scale_vs = scale_vs)
  scored = rb_obs$ok & rb_mod$ok
  if (!any(scored)) stop("No issue has a complete 25-29 day lag in both inputs.")
  y24o = matrix(NA_real_, length(i_keep), n24)
  y24m = matrix(NA_real_, length(i_keep), n24)
  y24o[scored, ] = recur_predict(rb_obs$pred[scored, , drop = FALSE])
  y24m[scored, ] = recur_predict(rb_mod$pred[scored, , drop = FALSE])

  # When the last observed Kp is already high, the in-situ VARX does not
  # fall back as fast as the storm does, and the AWSoM VARX is closer.
  # For each lead, earlier ended issues choose the Kp level and the blend.
  # The blend is used only if it lowered RMSE on those issues.
  pers = kp_value_live[i_keep]
  blend_lead = function(y, y_o, y_m, pers, tau_num, lead_hr) {
    out = y_o
    ok = is.finite(y_o) & is.finite(y_m) & is.finite(pers)
    thr_use = Inf
    w_use = 1
    next_fit = -Inf
    n_switch = 0L
    for (i in which(ok)) {
      if (tau_num[i] >= next_fit) {
        tr = which(ok & is.finite(y) & (tau_num + lead_hr * 3600) <= tau_num[i])
        thr_use = Inf
        w_use = 1
        if (length(tr) >= 80L) {
          best_r = sqrt(mean((y[tr] - y_o[tr])^2))
          for (thr in c(4, 4.5, 5, 5.5, 6, 6.5)) {
            if (sum(pers[tr] >= thr) < 8L) next
            hi = pers[tr] >= thr
            for (w in c(0, 0.25, 0.5, 0.75)) {
              pr = y_o[tr]
              pr[hi] = w * y_o[tr][hi] + (1 - w) * y_m[tr][hi]
              r = sqrt(mean((y[tr] - pr)^2))
              if (r < best_r - 1e-4) {
                best_r = r
                thr_use = thr
                w_use = w
              }
            }
          }
        }
        next_fit = tau_num[i] + 5 * 24 * 3600
      }
      if (pers[i] >= thr_use && w_use < 1) {
        out[i] = w_use * y_o[i] + (1 - w_use) * y_m[i]
        n_switch = n_switch + 1L
      }
    }
    cat(sprintf("    blended %d issues   Kp level %s   in-situ weight %.2f\n",
                n_switch,
                if (is.finite(thr_use)) sprintf("%.1f", thr_use) else "none",
                w_use))
    out
  }
  # Rule for one new issue, from every earlier issue whose target has ended.
  boost_rule_at = function(y, y_o, y_m, pers, tau_num, lead_hr, tau_new) {
    tr = which(is.finite(y) & is.finite(y_o) & is.finite(y_m) & is.finite(pers) &
               (tau_num + lead_hr * 3600) <= tau_new)
    thr_use = Inf
    w_use = 1
    if (length(tr) >= 80L) {
      best_r = sqrt(mean((y[tr] - y_o[tr])^2))
      for (thr in c(4, 4.5, 5, 5.5, 6, 6.5)) {
        if (sum(pers[tr] >= thr) < 8L) next
        hi = pers[tr] >= thr
        for (w in c(0, 0.25, 0.5, 0.75)) {
          pr = y_o[tr]
          pr[hi] = w * y_o[tr][hi] + (1 - w) * y_m[tr][hi]
          r = sqrt(mean((y[tr] - pr)^2))
          if (r < best_r - 1e-4) {
            best_r = r
            thr_use = thr
            w_use = w
          }
        }
      }
    }
    c(thr = thr_use, w = w_use)
  }
  yboost = matrix(NA_real_, length(i_keep), n24)
  tau_num = as.numeric(tau_keep)
  leads = seq(kp_step_hr, 24L, by = kp_step_hr)
  cat("Blending in AWSoM when earlier storms say it beats in-situ\n")
  for (k in seq_along(leads)) {
    yboost[, k] = blend_lead(y_true[, k], y24o[, k], y24m[, k], pers, tau_num, leads[k])
    cat(sprintf("  lead %2dh  finite %d\n", leads[k], sum(is.finite(yboost[, k]))))
  }

  metric_one = function(horizon_hr, lead_hr, input, y, yhat, ypers) {
    rm = sqrt(mean((y - yhat)^2))
    rp = sqrt(mean((y - ypers)^2))
    data.frame(
      horizon_hr = horizon_hr,
      lead_end_hr = lead_hr,
      input = input,
      n = length(y),
      RMSE = rm,
      MAE = mean(abs(y - yhat)),
      Correlation = suppressWarnings(cor(y, yhat)),
      RMSE_persistence = rp,
      skill_RMSE = if (rp == 0) NA_real_ else 1 - rm / rp,
      stringsAsFactors = FALSE
    )
  }
  rows = list()
  for (hz in 24L) {
    n_ahead = hz / kp_step_hr
    Yo = y24o
    Ym = y24m
    leads = seq(kp_step_hr, hz, by = kp_step_hr)
    for (k in seq_along(leads)) {
      rows[[length(rows) + 1L]] = metric_one(hz, leads[k], "Persistence", y_true[scored, k], pers[scored], pers[scored])
      rows[[length(rows) + 1L]] = metric_one(hz, leads[k], lab_meas, y_true[scored, k], Yo[scored, k], pers[scored])
      rows[[length(rows) + 1L]] = metric_one(hz, leads[k], lab_awso, y_true[scored, k], Ym[scored, k], pers[scored])
      m = is.finite(yboost[, k])
      if (any(m)) {
        rows[[length(rows) + 1L]] = metric_one(hz, leads[k], "Boosted", y_true[m, k], yboost[m, k], pers[m])
      }
    }
  }
  metrics = do.call(rbind, rows)
  # The 27-day chart uses the same blend. The decision reads only the latest
  # observed Kp, and each 3-hour slot uses the rule of that lead.
  pers_now = unname(kp_value_live[i_live])
  y_boost = outlook$y_insitu
  cat("27-day boost, same rule as the 24 h forecast\n")
  for (k in seq_along(leads)) {
    rule = boost_rule_at(y_true[, k], y24o[, k], y24m[, k], pers, tau_num, leads[k],
                         as.numeric(tau_live))
    cat(sprintf("    lead %2dh  Kp level %s  in-situ weight %.2f\n",
                leads[k],
                if (is.finite(rule["thr"])) sprintf("%.1f", rule["thr"]) else "none",
                unname(rule["w"])))
    if (pers_now >= unname(rule["thr"]) && unname(rule["w"]) < 1) {
      idx = seq(k, length(y_boost), by = 8L)
      y_boost[idx] = unname(rule["w"]) * outlook$y_insitu[idx] +
        (1 - unname(rule["w"])) * outlook$y_awsom[idx]
    }
  }
  outlook$y_boost <<- y_boost
  list(
    cal_end = cal_end, tau = tau_keep, pers = pers, y_true = y_true,
    y24o = y24o, y24m = y24m, yboost = yboost, scored = scored,
    lag_obs = rb_obs$lag, lag_mod = rb_mod$lag, metrics = metrics
  )
}


kp_hist_time = kp_time_live[(i_live - n_kp_hist + 1L):i_live]
ab_live = affine_updated(tau_live - history_hr * 3600)
# Kp that arrived after the forecast was issued. Prefer the GFZ series the
# model uses. If that file still ends at the issue, draw NOAA planetary Kp
# so the chart can show the days that have already happened.
kp_chart_time = NULL
kp_chart_value = NULL
kp_chart_source = NULL
day0 = as.POSIXct(format(tau_live, "%Y-%m-%d"), tz = "UTC")
gfz_end = kp_time_live + kp_step_hr * 3600
gfz_after = kp_time_live >= day0 & gfz_end <= Sys.time() & gfz_end > tau_live
if (any(gfz_after)) {
  kp_chart_time = gfz_end[gfz_after]
  kp_chart_value = kp_value_live[gfz_after]
  kp_chart_source = if (identical(kp_extend_source, "NOAA")) "NOAA" else "GFZ"
} else {
  noaa_tmp = tempfile(fileext = ".json")
  noaa_ok = tryCatch({
    curl_to("https://services.swpc.noaa.gov/products/noaa-planetary-k-index.json",
            noaa_tmp, c("-A", "Mozilla/5.0"))
    TRUE
  }, error = function(e) FALSE)
  if (noaa_ok) {
    noaa = tryCatch(jsonlite::fromJSON(noaa_tmp), error = function(e) NULL)
    if (!is.null(noaa) && !is.null(noaa$time_tag)) {
      noaa_start = parse_utc(gsub("T", " ", noaa$time_tag, fixed = TRUE))
      noaa_end = noaa_start + kp_step_hr * 3600
      noaa_kp = as.numeric(noaa$Kp)
      keep = noaa_start >= day0 & noaa_end <= Sys.time() & is.finite(noaa_kp)
      if (any(keep)) {
        kp_chart_time = noaa_end[keep]
        kp_chart_value = noaa_kp[keep]
        kp_chart_source = "NOAA"
      }
    }
  }
  unlink(noaa_tmp)
}
page_path = file.path(out_dir, "VARX_Kp_realtime.html")
latest_path = file.path(out_dir, "VARX_Kp_realtime_latest.csv")
driver_path = file.path(out_dir, "VARX_Kp_realtime_drivers.csv")
log_path = file.path(out_dir, "VARX_Kp_realtime_log.csv")
forecast_live = forecast
write.csv(forecast, latest_path, row.names = FALSE)
write.csv(drivers, driver_path, row.names = FALSE)
if (file.exists(log_path)) {
  prev = read.csv(log_path, stringsAsFactors = FALSE)
  if (identical(names(prev), names(forecast))) {
    forecast = rbind(prev, forecast)
  }
}
write.csv(forecast, log_path, row.names = FALSE)

cat(sprintf("\nIssue time (UTC): %s\n", fmt_utc(tau_live)))
cat(sprintf("Last observed Kp, bin %s to %s: %.3f\n",
            fmt_utc(kp_time_live[i_live]), fmt_utc(tau_live), kp_value_live[i_live]))
cat(sprintf("Earth_obs: the whole 24h driver is the hourly mean of %s, with no adjustment.\n", obs_path))
cat(sprintf("model_affine: hourly AWSoM after OMNI = a + b * AWSoM, fit on %d hours %s through %s.\n",
            ab_live$n, fmt_utc(ab_live$start), fmt_utc(ab_live$end)))
cat(sprintf("AWSoM map fit on %d hours %s through %s.\n",
            ab_live$n, fmt_utc(ab_live$start), fmt_utc(ab_live$end)))
for (f in seq_along(omni_features)) {
  cat(sprintf("  %-2s  a %12.4f  b %9.4f\n", omni_features[f], ab_live$a[f], ab_live$b[f]))
}
cat("Kp history is the same for both inputs and is not adjusted.\n\n")
sub = forecast[forecast$generated_at == generated_at, ]
obs_sub = sub[sub$input == "Earth_obs", ]
mod_sub = sub[sub$input == "model_affine", ]
obs_sub = obs_sub[order(obs_sub$lead_end_hr), ]
mod_sub = mod_sub[order(mod_sub$lead_end_hr), ]
cat("24-hour forecast\n")
for (r in seq_len(nrow(obs_sub))) {
  cat(sprintf("  Kp ending %s   in-situ %.3f   AWSoM %.3f   persistence %.3f\n",
              obs_sub$target_bin_end[r], obs_sub$y_pred[r], mod_sub$y_pred[r],
              obs_sub$y_persistence[r]))
}
cat(sprintf("  lag  in-situ %.1f days   AWSoM %.1f days\n\n",
            recur_lag[["Earth_obs"]] / 24, recur_lag[["model_affine"]] / 24))
cat(sprintf("Wrote %s\n", latest_path))
cat(sprintf("Wrote %s\n", driver_path))
cat(sprintf("Appended %s\n", log_path))

cat("\nScoring baseline where the 24 h target Kp has already been observed\n")
base = score_baseline()
f107_path = file.path(repo_root, "Data", "Origin", "F107", "SW-All.txt")
dir.create(dirname(f107_path), recursive = TRUE, showWarnings = FALSE)
f107_tmp = tempfile(fileext = ".txt")
f107_ok = tryCatch({
  curl_to("https://celestrak.org/SpaceData/SW-All.txt", f107_tmp)
  TRUE
}, error = function(e) FALSE)
if (isTRUE(f107_ok)) {
  f107_head = readLines(f107_tmp, warn = FALSE, n = 40)
  if (any(grepl("^BEGIN OBSERVED", f107_head))) {
    file.copy(f107_tmp, f107_path, overwrite = TRUE)
    cat("Updated F10.7 from CelesTrak SW-All.txt\n")
  } else {
    cat("CelesTrak download was not SW-All.txt. Keeping the file already saved.\n")
  }
} else {
  cat("F10.7 download failed. Keeping the file already saved.\n")
}
unlink(f107_tmp)
issued = issue_best_forecast(repo_root, out_dir, xgb_X, outlook, tau_live)
outlook = issued$outlook
outlook_path = file.path(out_dir, "VARX_Kp_realtime_27d.csv")
write.csv(outlook, outlook_path, row.names = FALSE)
cat(sprintf("Wrote %s\n", outlook_path))
f107_raw = load_f107_observed(file.path(repo_root, "Data", "Origin", "F107", "SW-All.txt"))
f107_measured = data.frame(
  date = f107_raw$date, f107_obs = f107_raw$obs, f107_adj = f107_raw$adj,
  stringsAsFactors = FALSE)
f107_measured = f107_measured[f107_raw$q %in% c(0L, 1L), ]
write_realtime_page(page_path, forecast_live, drivers, kp_hist_time, kp_hist, tau_live, generated_at, ab_live,
                    recur_lag[["Earth_obs"]], recur_lag[["model_affine"]],
                    outlook = outlook, f107 = issued$f107, f107_measured = f107_measured,
                    kp_obs_time = kp_chart_time, kp_obs_value = kp_chart_value,
                    kp_obs_source = kp_chart_source,
                    kp_history_note = if (identical(kp_extend_source, "NOAA")) nowcast_note else NULL)
cat(sprintf("Wrote %s\n", page_path))
# Day 1 of the 27-day outlook is the mean of the next eight Kp bins. The
# 24 h VARX already predicts those bins from solar wind. Compare that mean
# with the Kp-only day-1 outlook on the issues both can score.
day1_ix = match(as.numeric(base$tau), as.numeric(outlook_eval$tau))
day1_kp_ok = rep(FALSE, length(day1_ix))
day1_known = !is.na(day1_ix)
day1_kp_ok[day1_known] = is.finite(outlook_eval$y_pred[day1_ix[day1_known], 1])
day1_ok = base$scored & day1_known & day1_kp_ok &
  rowSums(is.finite(base$y_true)) == ncol(base$y_true) &
  rowSums(is.finite(base$y24o)) == ncol(base$y24o) &
  rowSums(is.finite(base$y24m)) == ncol(base$y24m) &
  rowSums(is.finite(base$yboost)) == ncol(base$yboost)
day1_obs = rowMeans(base$y_true[day1_ok, , drop = FALSE])
day1_kp = outlook_eval$y_kp[day1_ix[day1_ok], 1]
day1_mix = outlook_eval$y_pred[day1_ix[day1_ok], 1]
day1_rmse = function(yhat) sqrt(mean((day1_obs - yhat)^2))
cat(sprintf("Day-1 daily mean on %d shared issues\n", sum(day1_ok)))
cat(sprintf("  27-day Kp only          RMSE %.3f\n", day1_rmse(day1_kp)))
cat(sprintf("  27-day with archive     RMSE %.3f\n", day1_rmse(day1_mix)))
cat(sprintf("  mean of in-situ 24 h    RMSE %.3f\n", day1_rmse(rowMeans(base$y24o[day1_ok, , drop = FALSE]))))
cat(sprintf("  mean of AWSoM 24 h      RMSE %.3f\n", day1_rmse(rowMeans(base$y24m[day1_ok, , drop = FALSE]))))
cat(sprintf("  mean of boosted 24 h    RMSE %.3f\n", day1_rmse(rowMeans(base$yboost[day1_ok, , drop = FALSE]))))
# The realtime 27-day forecast, scored on every past issue. In-situ and AWSoM
# differ in the solar wind already observed. Later days both read the OMNI
# that followed the closest aligned cycle. Boosted uses the 24-hour rule.
score_outlook_rolls = function() {
  cat("Scoring the 27-day roll for in-situ, AWSoM, and boosted\n")
  omni_use = get("omni_match", envir = archive_env, inherits = FALSE)
  omni_origin = get("omni_origin", envir = archive_env, inherits = FALSE)
  n_lead = outlook_days * 8L
  n_iss = length(outlook_eval$tau)
  y_insitu = matrix(NA_real_, n_iss, n_lead)
  y_awsom = matrix(NA_real_, n_iss, n_lead)
  kp_back = 260L
  roll_one = function(i, kind) {
    if (!is.finite(outlook_eval$best_end[i])) return(rep(NA_real_, n_lead))
    tau = outlook_eval$tau[i]
    i_now0 = match(as.numeric(tau) - kp_step_hr * 3600, as.numeric(kp_time_live))
    if (is.na(i_now0) || i_now0 <= kp_back) return(rep(NA_real_, n_lead))
    hnum = seq(as.numeric(tau) - 40 * 86400, as.numeric(tau) + outlook_days * 86400, by = 3600)
    sw = matrix(NA_real_, length(hnum), length(omni_features))
    past = hnum < as.numeric(tau)
    if (kind == "obs") {
      ix = match(hnum[past], as.numeric(obs_h$time))
      ok = !is.na(ix)
      if (any(ok)) sw[which(past)[ok], ] = obs_h$mat[ix[ok], , drop = FALSE]
    } else {
      ab = tryCatch(affine_updated(tau - history_hr * 3600), error = function(e) NULL)
      if (is.null(ab)) return(rep(NA_real_, n_lead))
      ix = match(hnum[past], as.numeric(model_h$time))
      ok = !is.na(ix)
      if (any(ok)) {
        rows = model_h$mat[ix[ok], , drop = FALSE]
        sw[which(past)[ok], ] = sweep(sweep(rows, 2, ab$b, "*"), 2, ab$a, "+")
      }
    }
    fut = !past
    ix_om = archive_hour_index(
      omni_origin, outlook_eval$best_end[i] + hnum[fut] - as.numeric(tau), nrow(omni_use))
    ok_om = !is.na(ix_om)
    if (any(ok_om)) sw[which(fut)[ok_om], ] = omni_use[ix_om[ok_om], , drop = FALSE]
    y = rep(NA_real_, n_lead)
    kp_work = c(kp_value_live[(i_now0 - kp_back + 1L):i_now0], y)
    for (day in 0:(outlook_days - 1L)) {
      i_now = kp_back + day * 8L
      tau_d = tau + day * 24 * 3600
      hour_start = as.numeric(tau_d) - (history_hr:1L) * 3600
      ixh = match(hour_start, hnum)
      if (anyNA(ixh)) break
      block = sw[ixh, , drop = FALSE]
      if (!all(is.finite(block))) break
      rb = recur_batch(tau_d, i_now, kp_work,
                       matrix(kp_work[(i_now - n_kp_hist + 1L):i_now], nrow = 1),
                       matrix(as.vector(block), nrow = 1),
                       hnum, sw[, vs_col], hnum, sw[, vs_col])
      if (!isTRUE(rb$ok)) break
      yhat = as.numeric(recur_predict(rb$pred))
      if (!all(is.finite(yhat))) break
      slots = day * 8L + seq_len(8L)
      kp_work[i_now + seq_len(8L)] = yhat
      y[slots] = yhat
    }
    y
  }
  for (i in seq_len(n_iss)) {
    y_insitu[i, ] = roll_one(i, "obs")
    y_awsom[i, ] = roll_one(i, "model")
    if (i %% 200L == 0L) cat(sprintf("  rolled %d / %d\n", i, n_iss))
  }
  rule_for = function(y, y_o, y_m, pers, tau_num, lead_hr) {
    ok = is.finite(y_o) & is.finite(y_m) & is.finite(pers)
    thr_use = Inf
    w_use = 1
    next_fit = -Inf
    thr = rep(Inf, length(tau_num))
    w = rep(1, length(tau_num))
    for (i in which(ok)) {
      if (tau_num[i] >= next_fit) {
        tr = which(ok & is.finite(y) & (tau_num + lead_hr * 3600) <= tau_num[i])
        thr_use = Inf
        w_use = 1
        if (length(tr) >= 80L) {
          best_r = sqrt(mean((y[tr] - y_o[tr])^2))
          for (thr_c in c(4, 4.5, 5, 5.5, 6, 6.5)) {
            if (sum(pers[tr] >= thr_c) < 8L) next
            hi = pers[tr] >= thr_c
            for (ww in c(0, 0.25, 0.5, 0.75)) {
              pr = y_o[tr]
              pr[hi] = ww * y_o[tr][hi] + (1 - ww) * y_m[tr][hi]
              r = sqrt(mean((y[tr] - pr)^2))
              if (r < best_r - 1e-4) {
                best_r = r
                thr_use = thr_c
                w_use = ww
              }
            }
          }
        }
        next_fit = tau_num[i] + 5 * 24 * 3600
      }
      thr[i] = thr_use
      w[i] = w_use
    }
    list(thr = thr, w = w)
  }
  y_boost = y_insitu
  leads8 = seq(kp_step_hr, 24L, by = kp_step_hr)
  tau_num = as.numeric(base$tau)
  issue_num = as.numeric(outlook_eval$tau)
  for (k in seq_along(leads8)) {
    rule = rule_for(base$y_true[, k], base$y24o[, k], base$y24m[, k], base$pers, tau_num, leads8[k])
    pos = findInterval(issue_num, tau_num)
    idx = seq(k, n_lead, by = 8L)
    for (i in which(pos > 0L)) {
      if (outlook_eval$pers[i] >= rule$thr[pos[i]] && rule$w[pos[i]] < 1) {
        y_boost[i, idx] = rule$w[pos[i]] * y_insitu[i, idx] +
          (1 - rule$w[pos[i]]) * y_awsom[i, idx]
      }
    }
  }
  both = match(issue_num, tau_num)
  known = !is.na(both) & is.finite(y_insitu[, 8]) & is.finite(base$y24o[both, 8])
  if (any(known)) {
    cat(sprintf("  step 0 vs 24 h endpoint  in-situ n %d median abs %.4f\n",
                sum(known), median(abs(y_insitu[known, 8] - base$y24o[both[known], 8]))))
    known_m = known & is.finite(y_awsom[, 8]) & is.finite(base$y24m[both, 8])
    if (any(known_m))
      cat(sprintf("  step 0 vs 24 h endpoint  AWSoM   n %d median abs %.4f\n",
                  sum(known_m), median(abs(y_awsom[known_m, 8] - base$y24m[both[known_m], 8]))))
  }
  bin_row = function(lead_hr, input, y, yhat_v, ref) {
    md = is.finite(y) & is.finite(yhat_v) & is.finite(ref)
    n = sum(md)
    rm = if (n) sqrt(mean((y[md] - yhat_v[md])^2)) else NA_real_
    rr = if (n) sqrt(mean((y[md] - ref[md])^2)) else NA_real_
    data.frame(
      horizon_hr = outlook_days * 24L, lead_end_hr = lead_hr, input = input, n = n,
      RMSE = rm,
      MAE = if (n) mean(abs(y[md] - yhat_v[md])) else NA_real_,
      Correlation = if (n > 1L) suppressWarnings(cor(y[md], yhat_v[md])) else NA_real_,
      RMSE_reference = rr,
      skill_RMSE = if (!n || !is.finite(rr) || rr == 0) NA_real_ else 1 - rm / rr,
      stringsAsFactors = FALSE
    )
  }
  metrics = do.call(rbind, lapply(seq_len(n_lead), function(j) {
    ref = outlook_eval$y_match_bin[, j]
    rbind(
      bin_row(outlook_eval$leads_hr[j], "Best matched cycle", outlook_eval$y_true_bin[, j], ref, ref),
      bin_row(outlook_eval$leads_hr[j], "Boosted", outlook_eval$y_true_bin[, j], y_boost[, j], ref),
      bin_row(outlook_eval$leads_hr[j], "In-situ", outlook_eval$y_true_bin[, j], y_insitu[, j], ref),
      bin_row(outlook_eval$leads_hr[j], "AWSoM", outlook_eval$y_true_bin[, j], y_awsom[, j], ref)
    )
  }))
  outlook_eval$y_insitu_bin <<- y_insitu
  outlook_eval$y_awsom_bin <<- y_awsom
  outlook_eval$y_boost_bin <<- y_boost
  outlook_eval$y_pred_bin <<- y_boost
  outlook_eval$metrics <<- metrics
  for (lead in c(24L, 168L, 648L)) {
    sub = metrics[metrics$lead_end_hr == lead, ]
    ref = sub[sub$input == "Best matched cycle", ]
    cat(sprintf("  lead %3dh  n %d  best cycle %.3f  boosted %.3f  in-situ %.3f  AWSoM %.3f\n",
                lead, ref$n, ref$RMSE, sub$RMSE[sub$input == "Boosted"],
                sub$RMSE[sub$input == "In-situ"], sub$RMSE[sub$input == "AWSoM"]))
  }
  invisible(metrics)
}
score_outlook_rolls()
outlook_series = do.call(rbind, lapply(seq_along(outlook_eval$leads_hr), function(j) {
  keep = is.finite(outlook_eval$y_true_bin[, j])
  if (!any(keep)) return(NULL)
  lead = outlook_eval$leads_hr[j]
  data.frame(
    issue_time = fmt_utc(outlook_eval$tau[keep]),
    lead_end_hr = lead,
    target_end = fmt_utc(outlook_eval$tau[keep] + lead * 3600),
    y_true = outlook_eval$y_true_bin[keep, j],
    y_boost = outlook_eval$y_boost_bin[keep, j],
    y_insitu = outlook_eval$y_insitu_bin[keep, j],
    y_awsom = outlook_eval$y_awsom_bin[keep, j],
    y_matched = outlook_eval$y_match_bin[keep, j],
    stringsAsFactors = FALSE
  )
}))
write.csv(outlook_eval$metrics, outlook_metrics_path, row.names = FALSE)
write.csv(outlook_series, outlook_series_path, row.names = FALSE)
baseline_page = file.path(out_dir, "VARX_Kp_baseline.html")
baseline_metrics_path = file.path(out_dir, "VARX_Kp_baseline_metrics.csv")
baseline_endpoint_path = file.path(out_dir, "VARX_Kp_baseline_endpoint.csv")
xgb_pred_path = file.path(repo_root, "Rcode", "Main", "output",
                           "surtai_vanilla_xgboost_allhorizons", "external_predictions.csv")
issued_section = issued_quantity_baseline(
  xgb = if (file.exists(xgb_pred_path)) read.csv(xgb_pred_path, stringsAsFactors = FALSE) else NULL,
  roll = outlook_series,
  flux = score_f107_recent(repo_root)
)
write_baseline_page(baseline_page, base, generated_at, outlook_eval = outlook_eval,
                    issued_section = issued_section)
write.csv(base$metrics, baseline_metrics_path, row.names = FALSE)
write.csv(data.frame(
  issue_time = fmt_utc(base$tau),
  y_persistence = base$pers,
  y_true_24h = base$y_true[, 8],
  y_Earth_obs_24h = base$y24o[, 8],
  y_model_affine_24h = base$y24m[, 8],
  y_boost_24h = base$yboost[, 8],
  lag_obs_hr = base$lag_obs,
  lag_model_hr = base$lag_mod,
  stringsAsFactors = FALSE
), baseline_endpoint_path, row.names = FALSE)
ep = base$metrics[base$metrics$lead_end_hr == base$metrics$horizon_hr, ]
cat(sprintf("Baseline issues %d, %s through %s\n",
            length(base$tau), fmt_utc(min(base$tau)), fmt_utc(max(base$tau))))
for (i in seq_len(nrow(ep))) {
  cat(sprintf("  %2dh endpoint  %-32s  n %d  RMSE %.3f  MAE %.3f  corr %.3f  skill %+.3f\n",
              ep$horizon_hr[i], ep$input[i], ep$n[i], ep$RMSE[i], ep$MAE[i],
              ep$Correlation[i], ep$skill_RMSE[i]))
}
same = is.finite(base$yboost[, 8])
if (any(same)) {
  fair = function(tag, yhat) {
    rm = sqrt(mean((base$y_true[same, 8] - yhat[same])^2))
    cat(sprintf("  24h endpoint on boosted issues  %-20s  n %d  RMSE %.3f\n",
                tag, sum(same), rm))
  }
  fair("Persistence", base$pers)
  fair("In-situ", base$y24o[, 8])
  fair("AWSoM", base$y24m[, 8])
  fair("Boosted", base$yboost[, 8])
  slope = function(yh) coef(lm(yh[same] ~ base$y_true[same, 8]))[2]
  cat(sprintf("  slope of forecast on observed Kp   in-situ %.3f   AWSoM %.3f   boosted %.3f\n",
              slope(base$y24o[, 8]), slope(base$y24m[, 8]), slope(base$yboost[, 8])))
  cat(sprintf("  sd of 3 h changes                  observed %.3f   in-situ %.3f   boosted %.3f\n",
              sd(diff(base$y_true[same, 8])), sd(diff(base$y24o[same, 8])), sd(diff(base$yboost[same, 8]))))
}
updated_path = file.path(out_dir, "RTSW_affine_updated.csv")
write.csv(data.frame(
  feature = omni_features, a = ab_live$a, b = ab_live$b, n_cal = ab_live$n,
  window_start = fmt_utc(ab_live$start), window_end = fmt_utc(ab_live$end),
  stringsAsFactors = FALSE
), updated_path, row.names = FALSE)
cat(sprintf("Wrote %s\n", baseline_page))
cat(sprintf("Wrote %s\n", baseline_metrics_path))
cat(sprintf("Wrote %s\n", baseline_endpoint_path))
cat(sprintf("Wrote %s\n", updated_path))

if (identical(Sys.info()[["sysname"]], "Darwin")) {
  system2("open", c(page_path, baseline_page))
}

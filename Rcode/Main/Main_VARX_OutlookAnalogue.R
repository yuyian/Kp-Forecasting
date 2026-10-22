# 27-day Kp by analogue, not by linear regression.
#
# The library is OMNI from 1995 through the RTSW training cutoff in January
# 2026. Each candidate is one 27-day cycle of B, Bz, Vs, Pd and Pt. At an
# issue, the last 27 days of solar wind are matched to that library. The
# forecast is the Kp that followed the closest cycles. A cycle is eligible
# only when its following 27 days of Kp had already ended. The issued
# outlook is not changed.

script_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path = normalizePath(sub("^--file=", "", script_arg[1]))
  repo_root = normalizePath(file.path(dirname(script_path), "..", ".."))
} else {
  repo_root = normalizePath("~/pphd-files/Kp-Forecasting")
}

omni_features = c("B", "Bz", "Vs", "Pd", "Pt")
sentinel = c(B = 9999, Bz = 9999, Vs = 9999, Pd = 9999, Pt = 9999999)
kp_step_hr = 3L
outlook_days = 27L
cycle_days = 27L
min_hours = 12L
top_k = 10L
show_days = c(1L, 7L, 14L, 27L)

kp_path = file.path(repo_root, "Data", "Processed", "Kp.csv")
omni_path = file.path(repo_root, "Data", "Processed", "OMNI_hourly.csv")
obs_path = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_obs.csv")
model_path = file.path(repo_root, "Data", "Origin", "RTSW", "Earth_model.csv")
out_dir = file.path(repo_root, "Rcode", "Main", "output", "Main_VARX_Realtime")
coef_path = file.path(out_dir, "VARX_Kp_outlook27_coef.RData")
base_path = file.path(out_dir, "VARX_Kp_baseline_27d.csv")

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
hour_index = function(origin, hour_num, n_hours) {
  ix = (as.numeric(hour_num) - origin) / 3600 + 1
  r = round(ix)
  bad = !is.finite(ix) | abs(ix - r) > 1e-3 | r < 1 | r > n_hours
  out = as.integer(r)
  out[bad] = NA_integer_
  out
}

message("Loading Kp, OMNI, and the scored 27-day outlook")
env = new.env()
load(coef_path, envir = env)
train_end = as.POSIXct(env$outlook_sig$train_end, origin = "1970-01-01", tz = "UTC")
Kp = read.csv(kp_path, stringsAsFactors = FALSE)
Kp$Date = parse_utc(Kp$Date)
Kp = Kp[order(Kp$Date), ]
kp_time = Kp$Date
kp_value = Kp$Data
if (any(abs(diff(as.numeric(kp_time)) - kp_step_hr * 3600) >= 1))
  stop("Kp is not a continuous 3-hour grid.")
OMNI = read.csv(omni_path, stringsAsFactors = FALSE)
OMNI$Date = parse_utc(OMNI$Date)
OMNI = OMNI[order(OMNI$Date), ]
omni_time = OMNI$Date
if (any(abs(diff(as.numeric(omni_time)) - 3600) >= 1))
  stop("OMNI is not a continuous hourly grid.")
OMNI_mat = as.matrix(OMNI[, omni_features, drop = FALSE])
for (nm in omni_features)
  OMNI_mat[!is.finite(OMNI_mat[, nm]) | OMNI_mat[, nm] >= sentinel[nm], nm] = NA
base = read.csv(base_path, stringsAsFactors = FALSE)
base$issue_time = parse_utc(base$issue_time)
message(sprintf("Library ends at the training cutoff %s", fmt_utc(train_end)))

# One 27-day cycle ending at s, oldest day first. last_idx is the hour
# starting at s - 1 h. strict = TRUE keeps a day only when all 24 hours exist.
cycle_means = function(mat, last_idx, strict) {
  n = length(last_idx)
  days = vector("list", cycle_days)
  keep = rep(TRUE, n)
  for (j in 0:(cycle_days - 1L)) {
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
    enough = if (strict) n_ok == 24L else n_ok >= min_hours
    keep = keep & enough
    days[[cycle_days - j]] = acc / pmax(n_ok, 1L)
  }
  list(days = days, keep = keep)
}

feature_cors = function(days, query) {
  n = nrow(days[[1]])
  cors = rep(0, n)
  used = 0L
  for (f in seq_along(omni_features)) {
    q = query[, f]
    if (sum(is.finite(q)) < cycle_days || sd(q) == 0) next
    M = days[[f]]
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

# Library matrices are stored per feature, n x 27, so the correlation does
# not rebuild them for every issue.
pack_days = function(days, keep) {
  packed = vector("list", length(omni_features))
  for (f in seq_along(omni_features)) {
    M = matrix(NA_real_, length(keep), cycle_days)
    for (j in seq_len(cycle_days)) M[, j] = days[[j]][, f]
    packed[[f]] = M[keep, , drop = FALSE]
  }
  packed
}

kp_after = function(end_time) {
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

message("Building the 1995 to January 2026 OMNI cycle library")
phases = seq(0, 21, by = kp_step_hr)
library = vector("list", length(phases))
names(library) = sprintf("%02d", phases)
omni_origin = as.numeric(omni_time[1])
for (ph in phases) {
  start = as.POSIXct(sprintf("1995-02-01 %02d:00:00", ph), tz = "UTC")
  ends = seq(start, train_end, by = "24 hours")
  ends = ends[ends + outlook_days * 86400 <= train_end + outlook_days * 86400]
  ends = ends[ends <= train_end]
  last_idx = hour_index(omni_origin, as.numeric(ends) - 3600, nrow(OMNI_mat))
  built = cycle_means(OMNI_mat, last_idx, TRUE)
  outcome = kp_after(ends)
  keep = built$keep & rowSums(is.finite(outcome)) == outlook_days
  key = sprintf("%02d", ph)
  library[[key]] = list(
    end = ends[keep],
    feat = pack_days(built$days, keep),
    outcome = outcome[keep, , drop = FALSE]
  )
  message(sprintf("  phase %02d UTC  %d cycles", ph, sum(keep)))
}

native_hourly = function(path) {
  raw = read.csv(path, stringsAsFactors = FALSE)
  raw$Date = parse_utc(raw$Date)
  raw = raw[order(raw$Date), ]
  mat = as.matrix(raw[, omni_features, drop = FALSE])
  for (nm in omni_features)
    mat[!is.finite(mat[, nm]) | mat[, nm] >= sentinel[nm], nm] = NA
  hour = floor_hour(raw$Date)
  hours = seq(min(hour), max(hour), by = "hour")
  out = matrix(NA_real_, length(hours), length(omni_features))
  colnames(out) = omni_features
  bin = match(as.numeric(hour), as.numeric(hours))
  for (f in seq_along(omni_features)) {
    ok = is.finite(mat[, f]) & !is.na(bin)
    means = tapply(mat[ok, f], bin[ok], mean)
    out[as.integer(names(means)), f] = as.numeric(means)
  }
  list(time = hours, mat = out, origin = as.numeric(hours[1]))
}

analogue = function(query, feat, outcome, end, tau) {
  eligible = end + outlook_days * 86400 <= tau
  if (!any(eligible)) return(NULL)
  cors = feature_cors(lapply(feat, function(M) M[eligible, , drop = FALSE]), query)
  if (!any(is.finite(cors))) return(NULL)
  take = min(top_k, sum(is.finite(cors)))
  ord = order(cors, decreasing = TRUE)[seq_len(take)]
  w = pmax(cors[ord], 0)
  if (sum(w) == 0) w = rep(1, length(ord))
  y = outcome[eligible, , drop = FALSE]
  list(
    best = y[ord[1], ],
    top = colSums(w * y[ord, , drop = FALSE]) / sum(w),
    cor = cors[ord[1]],
    end = end[eligible][ord[1]]
  )
}

issue_query = function(hourly, tau) {
  last_idx = hour_index(hourly$origin, as.numeric(tau) - 3600, nrow(hourly$mat))
  if (is.na(last_idx)) return(NULL)
  built = cycle_means(hourly$mat, last_idx, FALSE)
  if (!built$keep) return(NULL)
  query = matrix(NA_real_, cycle_days, length(omni_features))
  for (j in seq_len(cycle_days)) query[j, ] = built$days[[j]][1, ]
  query
}

message("Reading in-situ and AWSoM")
obs_h = native_hourly(obs_path)
model_h = native_hourly(model_path)

issues = sort(unique(base$issue_time))
phase_of = function(tau) as.integer(format(tau, "%H", tz = "UTC"))
message(sprintf("Matching %d issues", length(issues)))
best_omni = matrix(NA_real_, length(issues), outlook_days)
top_omni = best_omni
best_obs = best_omni
top_obs = best_omni
best_aw = best_omni
top_aw = best_omni
cor_obs = rep(NA_real_, length(issues))
end_obs = rep(NA_real_, length(issues))
pair_row = hour_index(omni_origin, model_h$time, nrow(OMNI_mat))
model_ok = rowSums(is.finite(model_h$mat)) == length(omni_features)
omni_at_model = rep(FALSE, nrow(model_h$mat))
pair_known = !is.na(pair_row)
omni_at_model[pair_known] = rowSums(is.finite(OMNI_mat[pair_row[pair_known], , drop = FALSE])) == length(omni_features)
elig = which(model_ok & omni_at_model)
pair_time = model_h$time[elig]
affine_cache = new.env(parent = emptyenv())
affine_at = function(cutoff) {
  key = format(cutoff, "%Y-%m-%d", tz = "UTC")
  if (exists(key, envir = affine_cache, inherits = FALSE))
    return(get(key, envir = affine_cache, inherits = FALSE))
  last = findInterval(as.numeric(cutoff) - 3600, as.numeric(pair_time))
  ab = rep(NA_real_, 2L * length(omni_features))
  if (last >= 24L * 7L) {
    newest = pair_time[last]
    ix = elig[seq_len(last)]
    ix = ix[model_h$time[ix] > newest - 30 * 86400]
    if (length(ix) >= 24L) {
      for (f in seq_along(omni_features)) {
        x = model_h$mat[ix, f]
        y = OMNI_mat[pair_row[ix], f]
        ok = is.finite(x) & is.finite(y)
        x = x[ok]
        y = y[ok]
        b = if (length(x) < 5 || var(x) == 0) 1 else cov(x, y) / var(x)
        ab[f] = mean(y) - b * mean(x)
        ab[length(omni_features) + f] = b
      }
    }
  }
  assign(key, ab, envir = affine_cache)
  ab
}
for (i in seq_along(issues)) {
  tau = issues[i]
  key = sprintf("%02d", phase_of(tau))
  lib = library[[key]]
  q_omni_idx = hour_index(omni_origin, as.numeric(tau) - 3600, nrow(OMNI_mat))
  if (!is.na(q_omni_idx)) {
    built = cycle_means(OMNI_mat, q_omni_idx, TRUE)
    if (built$keep && as.numeric(tau) - cycle_days * 86400 <= as.numeric(max(omni_time))) {
      query = matrix(NA_real_, cycle_days, length(omni_features))
      for (j in seq_len(cycle_days)) query[j, ] = built$days[[j]][1, ]
      hit = analogue(query, lib$feat, lib$outcome, lib$end, tau)
      if (!is.null(hit)) {
        best_omni[i, ] = hit$best
        top_omni[i, ] = hit$top
      }
    }
  }
  query = issue_query(obs_h, tau)
  if (!is.null(query)) {
    hit = analogue(query, lib$feat, lib$outcome, lib$end, tau)
    if (!is.null(hit)) {
      best_obs[i, ] = hit$best
      top_obs[i, ] = hit$top
      cor_obs[i] = hit$cor
      end_obs[i] = as.numeric(hit$end)
    }
  }
  query = issue_query(model_h, tau)
  if (!is.null(query)) {
    ab = affine_at(tau - cycle_days * 86400)
    if (all(is.finite(ab))) {
      for (f in seq_along(omni_features))
        query[, f] = ab[f] + ab[length(omni_features) + f] * query[, f]
      hit = analogue(query, lib$feat, lib$outcome, lib$end, tau)
      if (!is.null(hit)) {
        best_aw[i, ] = hit$best
        top_aw[i, ] = hit$top
      }
    }
  }
  if (i %% 500L == 0L) message(sprintf("  %d / %d", i, length(issues)))
}

score_rows = function(sample, model, y, yhat, y_kp, y_rec) {
  ok = is.finite(y) & is.finite(yhat) & is.finite(y_kp)
  if (sum(ok) < 2L) return(NULL)
  rmse = function(a, b) sqrt(mean((a[ok] - b[ok])^2))
  data.frame(
    sample = sample, model = model, day = NA_integer_, n = sum(ok),
    RMSE = rmse(y, yhat),
    RMSE_kp = rmse(y, y_kp),
    RMSE_recur = rmse(y, y_rec),
    skill_vs_kp = 1 - rmse(y, yhat) / rmse(y, y_kp),
    stringsAsFactors = FALSE
  )
}

rows = list()
for (d in show_days) {
  sub = base[base$day == d, ]
  ix = match(as.numeric(sub$issue_time), as.numeric(issues))
  add = function(sample, model, yhat) {
    r = score_rows(sample, model, sub$y_true, yhat, sub$y_pred, sub$y_recur)
    if (!is.null(r)) {
      r$day = d
      rows[[length(rows) + 1L]] <<- r
    }
  }
  add("omni_query", "best_cycle", best_omni[ix, d])
  add("omni_query", "top10_cycles", top_omni[ix, d])
  add("insitu_query", "best_cycle", best_obs[ix, d])
  add("insitu_query", "top10_cycles", top_obs[ix, d])
  add("awsom_query", "best_cycle", best_aw[ix, d])
  add("awsom_query", "top10_cycles", top_aw[ix, d])
  same = is.finite(best_obs[ix, d])
  add("same_issues", "current_outlook", ifelse(same, sub$y_pred, NA))
  add("same_issues", "previous_rotation", ifelse(same, sub$y_recur, NA))
  add("same_issues", "insitu_best_cycle", best_obs[ix, d])
  add("same_issues", "insitu_top10", top_obs[ix, d])
  add("same_issues", "awsom_best_cycle", ifelse(same, best_aw[ix, d], NA))
  add("same_issues", "awsom_top10", ifelse(same, top_aw[ix, d], NA))
}
# Current outlook plus the top 10. The weight uses only issues whose target
# day has already ended, and it stays between the two forecasts.
convex_blend = function(y, current, top, min_past) {
  n = length(y)
  yhat = rep(NA_real_, n)
  w = rep(NA_real_, n)
  for (i in seq_len(n)) {
    if (!is.finite(current[i]) || !is.finite(top[i])) next
    if (i == 1L) {
      yhat[i] = current[i]
      w[i] = 1
      next
    }
    past = which(is.finite(y[seq_len(i - 1L)]) &
                   is.finite(current[seq_len(i - 1L)]) &
                   is.finite(top[seq_len(i - 1L)]))
    if (length(past) < min_past) {
      yhat[i] = current[i]
      w[i] = 1
      next
    }
    gap = current[past] - top[past]
    den = sum(gap * gap)
    ww = if (den < 1e-8) 1 else sum((y[past] - top[past]) * gap) / den
    ww = max(0, min(1, ww))
    yhat[i] = ww * current[i] + (1 - ww) * top[i]
    w[i] = ww
  }
  list(yhat = yhat, w = w)
}

oracle_weight = function(y, current, top) {
  ok = is.finite(y) & is.finite(current) & is.finite(top)
  gap = current[ok] - top[ok]
  den = sum(gap * gap)
  if (den < 1e-8) return(1)
  max(0, min(1, sum((y[ok] - top[ok]) * gap) / den))
}

blend_rows = list()
for (d in seq_len(outlook_days)) {
  sub = base[base$day == d, ]
  sub = sub[order(sub$issue_time), ]
  ix = match(as.numeric(sub$issue_time), as.numeric(issues))
  current = sub$y_pred
  top = top_obs[ix, d]
  known = is.finite(sub$y_true) & is.finite(current) & is.finite(top)
  if (sum(known) < 2L) next
  ww = oracle_weight(sub$y_true, current, top)
  mixed = ww * current + (1 - ww) * top
  walk = convex_blend(sub$y_true, current, top, 80L)
  ready = known & is.finite(walk$w) & walk$w < 1 - 1e-8
  # The first issues have no past weight yet, so both methods equal the
  # current forecast there. Compare only issues that actually mix the two.
  if (sum(ready) < 2L) ready = known
  rmse = function(a) sqrt(mean((sub$y_true[ready] - a[ready])^2))
  blend_rows[[length(blend_rows) + 1L]] = data.frame(
    day = d, n = sum(ready),
    RMSE_current = rmse(current),
    RMSE_top10 = rmse(top),
    RMSE_half = rmse(0.5 * current + 0.5 * top),
    RMSE_walk = rmse(walk$yhat),
    RMSE_oracle = rmse(mixed),
    oracle_w_current = ww,
    mean_walk_w = mean(walk$w[ready]),
    stringsAsFactors = FALSE
  )
}
blends = do.call(rbind, blend_rows)
blend_path = file.path(out_dir, "outlook_analogue_blend.csv")
write.csv(blends, blend_path, row.names = FALSE)
message(sprintf("Wrote %s", blend_path))
for (d in show_days) {
  r = blends[blends$day == d, ]
  message(sprintf(
    "blend day %d  n %d  current %.3f  top10 %.3f  half %.3f  walk %.3f  oracle %.3f (w %.2f)",
    d, r$n, r$RMSE_current, r$RMSE_top10, r$RMSE_half, r$RMSE_walk,
    r$RMSE_oracle, r$oracle_w_current))
}
message(sprintf(
  "all 27 days, mean RMSE  current %.3f  top10 %.3f  half %.3f  walk %.3f  oracle %.3f",
  mean(blends$RMSE_current), mean(blends$RMSE_top10), mean(blends$RMSE_half),
  mean(blends$RMSE_walk), mean(blends$RMSE_oracle)))

tab = do.call(rbind, rows)
out_path = file.path(out_dir, "outlook_analogue_comparison.csv")
write.csv(tab, out_path, row.names = FALSE)
message(sprintf("Wrote %s", out_path))
matched_end = as.POSIXct(end_obs, origin = "1970-01-01", tz = "UTC")
message(sprintf("In-situ matches %d issues, median correlation %.3f, match ends %s to %s",
                sum(is.finite(cor_obs)),
                median(cor_obs, na.rm = TRUE),
                fmt_utc(min(matched_end, na.rm = TRUE)),
                fmt_utc(max(matched_end, na.rm = TRUE))))
for (s in unique(tab$sample)) {
  message(sprintf("--- %s ---", s))
  sub = tab[tab$sample == s, ]
  for (d in show_days) {
    sd = sub[sub$day == d, ]
    if (!nrow(sd)) next
    message(sprintf("day %d  n %d", d, sd$n[1]))
    for (i in seq_len(nrow(sd))) {
      message(sprintf("  %-20s RMSE %.3f  vs current %+.3f  rotation %.3f",
                      sd$model[i], sd$RMSE[i], sd$skill_vs_kp[i], sd$RMSE_recur[i]))
    }
  }
}

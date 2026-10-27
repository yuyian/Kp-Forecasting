# Issued forecast for the webpage and the Swift driver files.
# Next 24 hours: the original vanilla XGBoost (one model per 3-hour lead).
# After that, through 27 days: the in-situ VARX rolled on the closest OMNI cycle.
# F10.7: one ridge regression per lead day, refit through the last measured flux.

kp_to_ap = function(kp) {
  knots = c(0, 1/3, 2/3, 1, 4/3, 5/3, 2, 7/3, 8/3, 3, 10/3, 11/3,
            4, 13/3, 14/3, 5, 16/3, 17/3, 6, 19/3, 20/3, 7, 22/3, 23/3,
            8, 25/3, 26/3, 9)
  apv = c(0, 2, 3, 4, 5, 6, 7, 9, 12, 15, 18, 22, 27, 32, 39, 48, 56, 67,
          80, 94, 111, 132, 154, 179, 207, 236, 300, 400)
  sapply(kp, function(v) {
    if (!is.finite(v) || v <= knots[1]) return(apv[1])
    if (v >= knots[length(knots)]) return(apv[length(apv)])
    i = which(knots <= v + 1e-9)
    i = i[length(i)]
    if (i >= length(knots)) return(apv[length(apv)])
    k0 = knots[i]; k1 = knots[i + 1]; a0 = apv[i]; a1 = apv[i + 1]
    if (a0 == 0) return(a0 + (v - k0) / (k1 - k0) * (a1 - a0))
    w = (v - k0) / (k1 - k0)
    exp(log(a0) + w * (log(a1) - log(a0)))
  })
}

predict_student_xgb = function(repo_root, X) {
  if (!requireNamespace("xgboost", quietly = TRUE))
    stop("The issued 24-hour forecast needs the xgboost package.")
  model_dir = file.path(repo_root, "Rcode", "Main", "output",
                        "surtai_vanilla_xgboost_allhorizons", "models")
  X = as.matrix(X)
  if (nrow(X) != 1L) stop("The issued XGBoost forecast expects one feature row.")
  leads = seq(3L, 24L, by = 3L)
  pred = numeric(length(leads))
  for (k in seq_along(leads)) {
    stem = sprintf("vanilla_xgboost_%02dh", leads[k])
    meta_path = file.path(model_dir, paste0(stem, "_metadata.rds"))
    model_path = file.path(model_dir, paste0(stem, ".ubj"))
    if (!file.exists(meta_path) || !file.exists(model_path))
      stop("Missing saved XGBoost model for the ", leads[k], " h lead in ", model_dir)
    meta = readRDS(meta_path)
    if (length(meta$scaler$center) != ncol(X))
      stop("XGBoost feature count does not match the saved scaler.")
    got = unname(tools::md5sum(model_path))
    if (!identical(got, unname(meta$booster_md5)))
      stop("Saved XGBoost file for ", leads[k], " h does not match its checksum.")
    z = sweep(sweep(X, 2L, meta$scaler$center, "-"), 2L, meta$scaler$scale, "/")
    if (any(!is.finite(z))) stop("The XGBoost inputs are not finite.")
    booster = xgboost::xgb.load(model_path)
    one = as.numeric(predict(booster, xgboost::xgb.DMatrix(z)))
    if (length(one) != 1L || !is.finite(one))
      stop("XGBoost did not return a finite forecast at ", leads[k], " h.")
    pred[k] = one
  }
  pred
}

load_f107_observed = function(path) {
  lines = readLines(path, warn = FALSE)
  a = match("BEGIN OBSERVED", lines)
  b = match("END OBSERVED", lines)
  if (is.na(a) || is.na(b) || b <= a + 1L) stop("SW-All.txt has no OBSERVED section.")
  body = lines[(a + 1L):(b - 1L)]
  body = body[nzchar(body) & grepl("^[0-9]{4}", body)]
  p = strsplit(body, "\\s+")
  date = as.Date(vapply(p, function(z) sprintf("%04d-%02d-%02d", as.integer(z[1]), as.integer(z[2]), as.integer(z[3])), ""))
  if (any(diff(date) != 1)) stop("Observed F10.7 is not one row per day.")
  list(
    date = date,
    adj = vapply(p, function(z) as.numeric(z[27]), 0),
    q = vapply(p, function(z) as.integer(z[28]), 0L),
    obs = vapply(p, function(z) as.numeric(z[31]), 0)
  )
}

forecast_f107 = function(repo_root, issue_time) {
  src = file.path(repo_root, "Data", "Origin", "F107", "SW-All.txt")
  if (!file.exists(src)) stop("Missing historical F10.7 file: ", src)
  raw = load_f107_observed(src)
  date = raw$date
  obs = raw$obs
  adj = raw$adj
  real = raw$q %in% c(0L, 1L)
  n = length(obs)
  origin = date[n]
  issue_day = as.Date(format(issue_time, "%Y-%m-%d", tz = "UTC"))
  if (origin != issue_day - 1)
    stop("The last measured F10.7 day is ", format(origin),
         ", which is not the day before the Kp issue ", format(issue_day), ".")
  cache_path = file.path(repo_root, "Rcode", "Main", "output", "Main_VARX_Realtime", "f107_fit.rds")
  cache = if (file.exists(cache_path)) tryCatch(readRDS(cache_path), error = function(e) NULL) else NULL
  if (!is.null(cache) && identical(cache$origin, origin) && identical(cache$n, n)) {
    beta = cache$beta
    mu = cache$mu
    sd = cache$sd
  } else {
    ok = real
    for (h in 1:28) {
      target_ok = c(real[seq.int(h + 1L, n)], rep(FALSE, h))
      ok = ok & target_ok
    }
    ok[seq_len(81L)] = FALSE
    ok[seq.int(n - 27L, n)] = FALSE
    usable = which(ok)
    issue_date = date[usable]
    train = usable[issue_date < as.Date("2019-01-01")]
    inner_fit = train[date[train] < as.Date("2015-01-01")]
    inner = train[date[train] >= as.Date("2015-01-01")]
    final_issues = usable[date[usable + 28L] <= origin]
    cs = c(0, cumsum(obs))
    trail_of = function(idx) (cs[idx + 1L] - cs[idx - 80L]) / 81
    design = function(idx, h) {
      nrot = if (h > 27L) 2L else 1L
      rec = idx + h - 27L * nrot
      doy = as.integer(format(date[idx + h], "%j"))
      ang = 2 * pi * doy / 365.25
      cbind(obs[idx], obs[idx - 1L], obs[idx - 2L], obs[idx - 6L], obs[rec],
            obs[idx] - obs[idx - 1L], obs[idx] - obs[idx - 27L], trail_of(idx),
            sin(ang), cos(ang), 1)
    }
    scale_fit = function(X) {
      mu = colMeans(X[, -ncol(X), drop = FALSE])
      xc = sweep(X[, -ncol(X), drop = FALSE], 2, mu, "-")
      s = sqrt(colSums(xc * xc) / nrow(X))
      s[!is.finite(s) | s == 0] = 1
      list(mu = mu, sd = s)
    }
    scale_apply = function(X, sc) {
      Z = X
      Z[, -ncol(X)] = sweep(sweep(X[, -ncol(X), drop = FALSE], 2, sc$mu, "-"), 2, sc$sd, "/")
      Z
    }
    fit_ridge = function(X, y, alpha) {
      p = ncol(X)
      A = crossprod(X) + alpha * diag(p)
      A[p, p] = A[p, p] - alpha
      as.numeric(solve(A, crossprod(X, y)))
    }
    alphas = c(0.1, 1, 10, 100, 1000)
    beta = vector("list", 28L)
    mu = vector("list", 28L)
    sd = vector("list", 28L)
    for (h in 1:28) {
      Xf = design(inner_fit, h)
      sc = scale_fit(Xf)
      yf = obs[inner_fit + h]
      Xi = scale_apply(design(inner, h), sc)
      yi = obs[inner + h]
      best_a = alphas[1]
      best_r = Inf
      for (a in alphas) {
        b = fit_ridge(scale_apply(Xf, sc), yf, a)
        hat = as.numeric(Xi %*% b)
        r = sqrt(mean((yi - hat)^2))
        if (is.finite(r) && r < best_r) {
          best_a = a
          best_r = r
        }
      }
      Xtr = design(final_issues, h)
      sc = scale_fit(Xtr)
      beta[[h]] = fit_ridge(scale_apply(Xtr, sc), obs[final_issues + h], best_a)
      mu[[h]] = sc$mu
      sd[[h]] = sc$sd
    }
    dir.create(dirname(cache_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(list(origin = origin, n = n, beta = beta, mu = mu, sd = sd), cache_path)
  }
  yday = as.integer(format(date, "%j"))
  ratio = tapply(adj[real & obs > 0] / obs[real & obs > 0], yday[real & obs > 0], median)
  i0 = n
  out_date = origin + seq_len(28L)
  hat = numeric(28L)
  for (h in 1:28) {
    nrot = if (h > 27L) 2L else 1L
    rec = i0 + h - 27L * nrot
    doy = as.integer(format(origin + h, "%j"))
    ang = 2 * pi * doy / 365.25
    x = c(obs[i0], obs[i0 - 1L], obs[i0 - 2L], obs[i0 - 6L], obs[rec],
          obs[i0] - obs[i0 - 1L], obs[i0] - obs[i0 - 27L],
          mean(obs[(i0 - 80L):i0]), sin(ang), cos(ang))
    z = c((x - mu[[h]]) / sd[[h]], 1)
    hat[h] = max(60, sum(z * beta[[h]]))
  }
  data.frame(
    date = out_date,
    f107_obs = hat,
    f107_adj = hat * as.numeric(ratio[as.character(as.integer(format(out_date, "%j")))]),
    stringsAsFactors = FALSE
  )
}

# Ridge scores for targets from `from` through the last measured day.
# Alpha is chosen on issues before 2019, the same choice as the issued forecast.
# Coefficients are fit only on issues whose 28-day target ended before `from`.
score_f107_recent = function(repo_root, from = as.Date("2026-01-01")) {
  src = file.path(repo_root, "Data", "Origin", "F107", "SW-All.txt")
  raw = load_f107_observed(src)
  date = raw$date
  obs = raw$obs
  real = raw$q %in% c(0L, 1L)
  n = length(obs)
  origin = date[n]
  ok = real
  for (h in 1:28) {
    target_ok = c(real[seq.int(h + 1L, n)], rep(FALSE, h))
    ok = ok & target_ok
  }
  ok[seq_len(81L)] = FALSE
  ok[seq.int(n - 27L, n)] = FALSE
  usable = which(ok)
  cs = c(0, cumsum(obs))
  trail_of = function(idx) (cs[idx + 1L] - cs[idx - 80L]) / 81
  design = function(idx, h) {
    nrot = if (h > 27L) 2L else 1L
    rec = idx + h - 27L * nrot
    doy = as.integer(format(date[idx + h], "%j"))
    ang = 2 * pi * doy / 365.25
    cbind(obs[idx], obs[idx - 1L], obs[idx - 2L], obs[idx - 6L], obs[rec],
          obs[idx] - obs[idx - 1L], obs[idx] - obs[idx - 27L], trail_of(idx),
          sin(ang), cos(ang), 1)
  }
  scale_fit = function(X) {
    mu = colMeans(X[, -ncol(X), drop = FALSE])
    xc = sweep(X[, -ncol(X), drop = FALSE], 2, mu, "-")
    s = sqrt(colSums(xc * xc) / nrow(X))
    s[!is.finite(s) | s == 0] = 1
    list(mu = mu, sd = s)
  }
  scale_apply = function(X, sc) {
    Z = X
    Z[, -ncol(X)] = sweep(sweep(X[, -ncol(X), drop = FALSE], 2, sc$mu, "-"), 2, sc$sd, "/")
    Z
  }
  fit_ridge = function(X, y, alpha) {
    p = ncol(X)
    A = crossprod(X) + alpha * diag(p)
    A[p, p] = A[p, p] - alpha
    as.numeric(solve(A, crossprod(X, y)))
  }
  train = usable[date[usable] < as.Date("2019-01-01")]
  inner_fit = train[date[train] < as.Date("2015-01-01")]
  inner = train[date[train] >= as.Date("2015-01-01")]
  fit_issues = usable[date[usable + 28L] < from]
  alphas = c(0.1, 1, 10, 100, 1000)
  rows = list()
  for (h in 1:28) {
    Xf = design(inner_fit, h)
    sc0 = scale_fit(Xf)
    yf = obs[inner_fit + h]
    Xi = scale_apply(design(inner, h), sc0)
    yi = obs[inner + h]
    best_a = alphas[1]
    best_r = Inf
    for (a in alphas) {
      b = fit_ridge(scale_apply(Xf, sc0), yf, a)
      hat = as.numeric(Xi %*% b)
      r = sqrt(mean((yi - hat)^2))
      if (is.finite(r) && r < best_r) {
        best_a = a
        best_r = r
      }
    }
    sc = scale_fit(design(fit_issues, h))
    beta = fit_ridge(scale_apply(design(fit_issues, h), sc), obs[fit_issues + h], best_a)
    idx = which(real & seq_along(date) > 81L & (date + h) >= from & (date + h) <= origin)
    idx = idx[real[idx + h]]
    x = design(idx, h)
    z = x
    z[, -ncol(x)] = sweep(sweep(x[, -ncol(x), drop = FALSE], 2, sc$mu, "-"), 2, sc$sd, "/")
    hat = pmax(60, as.numeric(z %*% beta))
    y = obs[idx + h]
    pers = obs[idx]
    rm = sqrt(mean((y - hat)^2))
    rp = sqrt(mean((y - pers)^2))
    rows[[length(rows) + 1L]] = data.frame(
      split = sprintf("%s to %s", format(from), format(origin)),
      lead_day = h, model = "ridge", n = length(y),
      RMSE = rm, skill_vs_persistence = 1 - rm / rp,
      stringsAsFactors = FALSE)
    rows[[length(rows) + 1L]] = data.frame(
      split = sprintf("%s to %s", format(from), format(origin)),
      lead_day = h, model = "persistence", n = length(y),
      RMSE = rp, skill_vs_persistence = 0,
      stringsAsFactors = FALSE)
  }
  tab = do.call(rbind, rows)
  rownames(tab) = NULL
  tab
}

write_submission_annex = function(repo_root, sub_dir, issue_time, xgb_X = NULL,
                                    solar_wind = NULL, kp_time = NULL, kp_value = NULL) {
  model_dir = file.path(repo_root, "Rcode", "Main", "output",
                        "surtai_vanilla_xgboost_allhorizons", "models")
  leads = seq(3L, 24L, by = 3L)
  metas = lapply(leads, function(h) {
    readRDS(file.path(model_dir, sprintf("vanilla_xgboost_%02dh_metadata.rds", h)))
  })
  feat = metas[[1]]$feature_names
  if (any(vapply(metas, function(m) !identical(m$feature_names, feat), logical(1))))
    stop("The eight XGBoost models do not share one feature list.")
  xgb_rows = lapply(metas, function(m) {
    p = m$selected_params
    data.frame(
      lead_end_hr = m$lead_end_hr,
      objective = p$objective,
      booster = p$booster,
      tree_method = p$tree_method,
      max_depth = p$max_depth,
      eta = p$eta,
      min_child_weight = p$min_child_weight,
      subsample = p$subsample,
      colsample_bytree = p$colsample_bytree,
      rounds = m$selected_rounds,
      training_cutoff_utc = format(m$historical_cutoff, "%Y-%m-%d %H:%M", tz = "UTC"),
      training_rows = m$refit_rows,
      stringsAsFactors = FALSE
    )
  })
  write.csv(do.call(rbind, xgb_rows),
            file.path(sub_dir, "parameters_xgboost.csv"),
            row.names = FALSE, quote = FALSE)
  write.csv(data.frame(position = seq_along(feat), feature = feat),
            file.path(sub_dir, "xgboost_features.csv"),
            row.names = FALSE, quote = FALSE)
  if (!is.null(xgb_X)) {
    x = as.numeric(xgb_X)
    if (length(x) != length(feat)) stop("The issued XGBoost row does not have 146 inputs.")
    write.csv(data.frame(position = seq_along(feat), feature = feat, value = x),
              file.path(sub_dir, "xgboost_input_at_issue.csv"),
              row.names = FALSE, quote = FALSE)
  }
  coef_path = file.path(repo_root, "Rcode", "Main", "output",
                        "Main_VARX_Realtime", "VARX_Kp_coef.RData")
  env = new.env()
  load(coef_path, envir = env)
  B = env$fit_recur$B
  terms = c("intercept", feat)
  if (ncol(B) != length(terms)) stop("The in-situ regression width does not match the feature list.")
  varx = do.call(rbind, lapply(seq_len(nrow(B)), function(i) {
    data.frame(lead_end_hr = leads[i], term = terms, coefficient = B[i, ],
               stringsAsFactors = FALSE)
  }))
  write.csv(varx, file.path(sub_dir, "parameters_varx_insitu.csv"),
            row.names = FALSE, quote = FALSE)
  fit = readRDS(file.path(repo_root, "Rcode", "Main", "output",
                          "Main_VARX_Realtime", "f107_fit.rds"))
  ridge_names = c("flux", "flux_lag_1d", "flux_lag_2d", "flux_lag_6d",
                  "recurrence", "delta_1d", "delta_27d", "trail_81d",
                  "sin_target_doy", "cos_target_doy", "intercept")
  ridge = do.call(rbind, lapply(seq_along(fit$beta), function(h) {
    b = fit$beta[[h]]
    center = c(fit$mu[[h]], NA_real_)
    scale = c(fit$sd[[h]], NA_real_)
    data.frame(lead_day = h, term = ridge_names, center = center, scale = scale,
               coefficient = b, stringsAsFactors = FALSE)
  }))
  write.csv(ridge, file.path(sub_dir, "parameters_f107_ridge.csv"),
            row.names = FALSE, quote = FALSE)
  if (!is.null(kp_time) && !is.null(kp_value)) {
    write.csv(data.frame(
      bin_start_utc = format(kp_time, "%Y-%m-%d %H:%M", tz = "UTC"),
      bin_end_utc = format(kp_time + 3 * 3600, "%Y-%m-%d %H:%M", tz = "UTC"),
      kp = kp_value, stringsAsFactors = FALSE
    ), file.path(sub_dir, "kp_input_at_issue.csv"), row.names = FALSE, quote = FALSE)
  }
  if (!is.null(solar_wind) && nrow(solar_wind)) {
    wind = solar_wind
    if ("input" %in% names(wind)) wind = wind[wind$input == "Earth_obs", ]
    keep = intersect(c("hour_start", "n_native", "B", "Bz", "Vs", "Pd", "Pt"), names(wind))
    write.csv(wind[, keep, drop = FALSE],
              file.path(sub_dir, "solar_wind_input_at_issue.csv"),
              row.names = FALSE, quote = FALSE)
  }
  raw = load_f107_observed(file.path(repo_root, "Data", "Origin", "F107", "SW-All.txt"))
  origin = as.Date(format(issue_time, "%Y-%m-%d", tz = "UTC")) - 1
  ix = which(raw$date == origin)
  if (length(ix) == 1L && ix >= 81L) {
    sl = (ix - 80L):ix
    write.csv(data.frame(
      date = format(raw$date[sl]),
      f107_obs = raw$obs[sl],
      f107_adj = raw$adj[sl],
      stringsAsFactors = FALSE
    ), file.path(sub_dir, "f107_input_at_issue.csv"), row.names = FALSE, quote = FALSE)
  }
  issue_s = format(issue_time, "%Y-%m-%d %H:%M", tz = "UTC")
  writeLines(c(
    "University of Michigan submission for the NASA CCMC Swift re-entry challenge.",
    paste("Issued", issue_s, "UTC."),
    "",
    "This folder is an F10.7 and Kp forecast. ap is the Bartels amplitude of that Kp.",
    "These series replace the predictive F10.7 and Kp/Ap supplied by CCMC.",
    "They are not those provided drivers.",
    "This folder has no thermosphere density.",
    "This folder has no reentry time and no trajectory.",
    "No ballistic coefficient and no initial state vector were used.",
    "",
    "Files",
    "  UM_Kp_ap_forecast_*.csv       3-hour Kp, ap, and F10.7",
    "  UM_Kp_ap_daily_*.csv          daily ap, Kp, and F10.7",
    "  f107_forecast_*.csv           daily F10.7 from the ridge regression",
    "  kp_input_at_issue.csv         eight Kp bins known at issue time",
    "  solar_wind_input_at_issue.csv 24 hours of near-Earth solar wind known at issue time",
    "  f107_input_at_issue.csv       last 81 measured flux days used by the ridge",
    "  xgboost_input_at_issue.csv    the 146 numbers read by the 24-hour models, when saved",
    "  xgboost_features.csv          names of those 146 inputs",
    "  parameters_xgboost.csv        depth, learning rate, and rounds for each lead",
    "  parameters_varx_insitu.csv    in-situ regression used from 27 hours through 27 days",
    "  parameters_f107_ridge.csv     ridge coefficients for each lead day",
    "",
    "The 24-hour Kp models read the eight Kp bins, 24 hours of measured |B|, Bz,",
    "speed, proton density, and proton temperature, and a solar-rotation match",
    "from 25 to 29 days earlier. AWSoM solar wind is not an input to the issued Kp.",
    "From 27 hours through 27 days, Kp is the in-situ regression stepped forward",
    "on the closest OMNI cycle. F10.7 is a separate ridge regression on the flux",
    "history. Its penalty was chosen on 2015-2018 from 0.1, 1, 10, 100, and 1000,",
    "and the stored coefficients are the ones applied. Forecast flux is kept at",
    "or above 60 sfu. f107_adj multiplies that forecast by the historical median",
    "ratio of adjusted to observed flux for the target day of the year.",
    "",
    "What may be made public",
    "The solar and geomagnetic numbers in this folder may be made public.",
    "That includes the forecasts, the inputs at issue time, and the parameters.",
    "",
    "What is not submitted and should not be made public",
    "Nothing outside this folder was uploaded. Please do not publish serialized",
    "model files, the training archive, validation tables, or any contact address.",
    "No thermosphere density and no reentry prediction are included."
  ), file.path(sub_dir, "SUBMISSION.txt"))
  invisible(sub_dir)
}

write_issued_package = function(repo_root, out_dir, issue_time, outlook, f107,
                                xgb_X = NULL, solar_wind = NULL,
                                kp_time = NULL, kp_value = NULL) {
  if (!"y_issued" %in% names(outlook)) stop("The outlook has no issued column.")
  end = as.POSIXct(outlook$day_end, tz = "UTC")
  start = end - 3 * 3600
  lead = as.integer(round(as.numeric(difftime(end, issue_time, units = "hours"))))
  kp = round(outlook$y_issued, 3)
  ap = round(kp_to_ap(kp), 2)
  flux_obs = f107$f107_obs[match(as.Date(format(start, "%Y-%m-%d", tz = "UTC")), f107$date)]
  flux_adj = f107$f107_adj[match(as.Date(format(start, "%Y-%m-%d", tz = "UTC")), f107$date)]
  if (any(!is.finite(flux_obs))) stop("An issued Kp bin has no F10.7 forecast.")
  issue_s = format(issue_time, "%Y-%m-%d %H:%M", tz = "UTC")
  three = data.frame(
    issue_time_utc = issue_s,
    forecast_time_utc = format(end, "%Y-%m-%d %H:%M", tz = "UTC"),
    lead_h = lead,
    kp = sprintf("%.3f", kp),
    ap = sprintf("%.2f", ap),
    f107_obs = sprintf("%.1f", flux_obs),
    f107_adj = sprintf("%.1f", flux_adj),
    stringsAsFactors = FALSE
  )
  day = format(start, "%Y-%m-%d", tz = "UTC")
  daily = do.call(rbind, lapply(split(seq_along(day), day), function(ix) {
    data.frame(
      date = day[ix[1]],
      n_bins = length(ix),
      ap_mean = sprintf("%.2f", mean(ap[ix])),
      ap_max = sprintf("%.2f", max(ap[ix])),
      kp_mean = sprintf("%.3f", mean(kp[ix])),
      kp_max = sprintf("%.3f", max(kp[ix])),
      f107_obs = three$f107_obs[ix[1]],
      f107_adj = three$f107_adj[ix[1]],
      stringsAsFactors = FALSE
    )
  }))
  daily = daily[order(daily$date), ]
  rownames(daily) = NULL
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(three, file.path(out_dir, "issued_forecast_3h.csv"), row.names = FALSE, quote = FALSE)
  write.csv(daily, file.path(out_dir, "issued_forecast_daily.csv"), row.names = FALSE, quote = FALSE)
  if (!all(c("y_issued", "y_boost", "y_insitu", "y_awsom", "y_matched", "day", "day_end") %in% names(outlook)))
    stop("The 27-day outlook is missing issued or boosted values.")
  end27 = as.POSIXct(outlook$day_end, tz = "UTC")
  start27 = end27 - 3 * 3600
  lead27 = as.integer(round(as.numeric(difftime(end27, issue_time, units = "hours"))))
  issued_kp = round(outlook$y_issued, 3)
  boosted_kp = round(outlook$y_boost, 3)
  insitu_kp = round(outlook$y_insitu, 3)
  awsom_kp = round(outlook$y_awsom, 3)
  matched_kp = round(outlook$y_matched, 3)
  flux27 = match(as.Date(format(start27, "%Y-%m-%d", tz = "UTC")), f107$date)
  if (any(is.na(flux27))) stop("A 27-day bin has no F10.7 forecast.")
  twentyseven = data.frame(
    forecast_time_utc = format(end27, "%Y-%m-%d %H:%M", tz = "UTC"),
    day = outlook$day,
    lead_h = lead27,
    issued_kp = sprintf("%.3f", issued_kp),
    boosted_kp = sprintf("%.3f", boosted_kp),
    insitu_kp = sprintf("%.3f", insitu_kp),
    awsom_kp = sprintf("%.3f", awsom_kp),
    matched_kp = sprintf("%.3f", matched_kp),
    issued_ap = sprintf("%.2f", round(kp_to_ap(issued_kp), 2)),
    boosted_ap = sprintf("%.2f", round(kp_to_ap(boosted_kp), 2)),
    insitu_ap = sprintf("%.2f", round(kp_to_ap(insitu_kp), 2)),
    awsom_ap = sprintf("%.2f", round(kp_to_ap(awsom_kp), 2)),
    f107_obs = sprintf("%.1f", f107$f107_obs[flux27]),
    f107_adj = sprintf("%.1f", f107$f107_adj[flux27]),
    stringsAsFactors = FALSE
  )
  write.csv(twentyseven, file.path(out_dir, "issued_forecast_27d.csv"), row.names = FALSE, quote = FALSE)
  stamp = format(issue_time, "%Y%m%dT%H%M", tz = "UTC")
  day_dir = file.path(repo_root, "swift_challenge", format(issue_time, "%Y%m%d", tz = "UTC"))
  dir.create(day_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(three, file.path(day_dir, paste0("UM_Kp_ap_forecast_", stamp, ".csv")), row.names = FALSE, quote = FALSE)
  write.csv(daily, file.path(day_dir, paste0("UM_Kp_ap_daily_", stamp, ".csv")), row.names = FALSE, quote = FALSE)
  write.csv(data.frame(
    date = format(f107$date), lead_day = seq_len(nrow(f107)),
    f107_obs = sprintf("%.1f", f107$f107_obs),
    f107_adj = sprintf("%.1f", f107$f107_adj),
    model = "ridge", stringsAsFactors = FALSE
  ), file.path(day_dir, paste0("f107_forecast_", format(issue_time, "%Y%m%d", tz = "UTC"), ".csv")),
  row.names = FALSE, quote = FALSE)
  hi = daily[order(-as.numeric(daily$ap_mean)), ]
  readme = c(
    "University of Michigan Kp and F10.7 forecast for the NASA CCMC Swift re-entry challenge.",
    paste("Issued", issue_s, "UTC."),
    "",
    "The issued Kp uses the better model at each range. The first 24 hours are the",
    "original vanilla XGBoost: one model per 3-hour lead, trained on OMNI and Kp",
    "through 2026-01-25. Each model sees the last eight Kp bins, 24 hours of measured",
    "B, Bz, Vs, proton density and proton temperature, and the solar-rotation match",
    "from that training recipe. The eight leads are forecast directly. From the bin",
    "that ends 27 hours after issue through 27 days, Kp is the in-situ regression",
    "stepped forward one day at a time. Later days use the solar wind that followed",
    "the closest matching OMNI cycle from 1995 through January 2026, and the predicted",
    "Kp is fed forward. AWSoM is not part of the issued Kp.",
    "",
    "ap is the 3-hour planetary amplitude for that kp, from the Bartels scale used by",
    "GFZ and NOAA. Between scale points, ap is interpolated in log(ap). The daily ap",
    "is the mean of the 3-hour ap values whose bin starts on that UT day.",
    paste0("The first day, ", daily$date[1], ", has ", daily$n_bins[1],
           if (daily$n_bins[1] == 1) " bin" else " bins",
           " because the forecast starts at ", issue_s, " UTC. The last day, ",
           daily$date[nrow(daily)], ", has ", daily$n_bins[nrow(daily)],
           if (daily$n_bins[nrow(daily)] == 1) " bin." else " bins."),
    "",
    "forecast_time_utc is the end of the 3-hour period that Kp describes.",
    "lead_h is the number of hours from the issue time to forecast_time_utc.",
    "The first period ends 3 hours after the issue, so its lead_h is 3.",
    "",
    "F10.7 is the observed 10.7 cm flux. Penticton measures at 20:00 UT, so the flux",
    "on the issue date was not yet known at 12:00 UTC. Each lead day has its own ridge",
    "regression, fit on the measured flux from 1957 through the day before issue.",
    "On 2019-2022 the mean RMSE over leads 1 to 27 was 13.20 sfu, against 16.61 for",
    "repeating the latest flux. f107_adj multiplies the forecast by the median historical",
    "ratio of adjusted to observed flux for that day of the year. On the three-hour",
    "file both flux columns repeat the value of the UT day on which the bin starts.",
    "",
    paste0("The highest daily Ap values are ",
           paste(sprintf("%s (mean Ap %s, maximum Ap %s, maximum Kp %s)",
                         hi$date[1:2], hi$ap_mean[1:2], hi$ap_max[1:2], hi$kp_max[1:2]),
                 collapse = " and "),
           ".")
  )
  writeLines(readme, file.path(day_dir, "README.txt"))
  sub_dir = file.path(repo_root, "swift_challenge", "submission")
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  old = list.files(sub_dir, full.names = TRUE)
  if (length(old)) unlink(old)
  file.copy(file.path(day_dir, paste0("UM_Kp_ap_forecast_", stamp, ".csv")), sub_dir)
  file.copy(file.path(day_dir, paste0("UM_Kp_ap_daily_", stamp, ".csv")), sub_dir)
  file.copy(file.path(day_dir, paste0("f107_forecast_", format(issue_time, "%Y%m%d", tz = "UTC"), ".csv")), sub_dir)
  file.copy(file.path(day_dir, "README.txt"), sub_dir)
  write_submission_annex(repo_root, sub_dir, issue_time, xgb_X, solar_wind, kp_time, kp_value)
  list(three = three, daily = daily, dir = day_dir, submission = sub_dir)
}

issue_best_forecast = function(repo_root, out_dir, xgb_X, outlook, issue_time,
                               solar_wind = NULL, kp_time = NULL, kp_value = NULL) {
  y_xgb = predict_student_xgb(repo_root, xgb_X)
  if (length(y_xgb) != 8L) stop("The 24-hour XGBoost forecast must have eight values.")
  outlook$y_issued = outlook$y_insitu
  outlook$y_issued[seq_len(8L)] = y_xgb
  if (any(!is.finite(outlook$y_issued))) stop("The issued 27-day Kp forecast is incomplete.")
  f107 = forecast_f107(repo_root, issue_time)
  wrote = write_issued_package(repo_root, out_dir, issue_time, outlook, f107,
                                xgb_X, solar_wind, kp_time, kp_value)
  cat("Issued 24 h XGBoost:")
  cat(sprintf(" %.3f", y_xgb), "\n", sep = "")
  cat(sprintf("Issued package: %s\n", wrote$dir))
  list(outlook = outlook, f107 = f107, y_xgb = y_xgb)
}

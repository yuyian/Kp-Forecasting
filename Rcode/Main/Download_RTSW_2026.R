# =============================================================================
# Download_RTSW_2026.R
# Earth solar-wind series from the Michigan CLEAR realtime AWSoM page
# https://solarwind.engin.umich.edu/rtsw/
# =============================================================================
#
# In-situ observations stay on the JSON feed:
#   https://solarwind.engin.umich.edu/rtsw/api/series/?t=earth
# AWSoM at Earth is the satellite file, not the model block of that JSON:
#   https://solarwind.engin.umich.edu/rtsw/data/sat_earth_realtime.sat
#
# Nothing is filtered by year. The script keeps every obs array that lines
# up with the time vector, and every AWSoM sample in the satellite file.
# cal_until marks obs rows that are still the provisional realtime feed.
#
# Columns that exist in Data/Processed/OMNI_hourly.csv are written in that
# order: Date, B, Bz, Vs, Pd, Pt.
#   obs JSON: b -> B, bz -> Bz, v -> Vs, n -> Pd, T -> Pt, time_utc -> Date
#   AWSoM sat: n -> Pd, sqrt(ux^2+uy^2+uz^2) -> Vs, t -> Pt.
#   bx, by, bz in the sat file are Gauss; they are multiplied by 1e5 to get nT.
#   B is the magnitude of that vector, and Bz is the converted bz.
#
# Dependencies: jsonlite.
# Run: Rscript Rcode/Main/Download_RTSW_2026.R
# =============================================================================

script_arg = grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg)) {
  script_path = normalizePath(sub("^--file=", "", script_arg[1]))
  repo_root = normalizePath(file.path(dirname(script_path), "..", ".."))
} else {
  repo_root = normalizePath("~/pphd-files/Kp-Forecasting")
}

suppressPackageStartupMessages(library(jsonlite))

rtsw_url = "https://solarwind.engin.umich.edu/rtsw/api/series/?t=earth"
awsom_sat_url = "https://solarwind.engin.umich.edu/rtsw/data/sat_earth_realtime.sat"
out_dir = file.path(repo_root, "Data", "Origin", "RTSW")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

curl_download = function(url, dest) {
  status = system2("curl", c("--max-time", "120", "-fsSL", url, "-o", dest))
  if (!identical(status, 0L)) {
    message("Verified download failed; retrying without certificate check.")
    status = system2("curl", c("--max-time", "120", "-kfsSL", url, "-o", dest))
  }
  if (!identical(status, 0L)) stop("curl download failed for ", url, ", status ", status)
  invisible(dest)
}

read_rtsw = function(url) {
  tmp = tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  curl_download(url, tmp)
  fromJSON(tmp)
}

# AWSoM Earth trajectory. Magnetic components are Gauss in this file.
read_awsom_sat = function(url) {
  tmp = tempfile(fileext = ".sat")
  on.exit(unlink(tmp), add = TRUE)
  curl_download(url, tmp)
  raw = utils::read.table(tmp, skip = 1, header = TRUE, check.names = FALSE,
                          stringsAsFactors = FALSE)
  need = c("year", "mo", "dy", "hr", "mn", "sc", "n", "ux", "uy", "uz", "bx", "by", "bz", "t")
  missing = setdiff(need, names(raw))
  if (length(missing))
    stop("AWSoM sat file is missing columns: ", paste(missing, collapse = ", "))
  if (!nrow(raw)) stop("AWSoM sat file has no samples.")
  bx = raw$bx * 1e5
  by = raw$by * 1e5
  bz = raw$bz * 1e5
  df = data.frame(
    Date = sprintf("%04d-%02d-%02d %02d:%02d:%02d",
                   raw$year, raw$mo, raw$dy, raw$hr, raw$mn, raw$sc),
    B = round(sqrt(bx * bx + by * by + bz * bz), 2),
    Bz = round(bz, 2),
    Vs = round(sqrt(raw$ux^2 + raw$uy^2 + raw$uz^2), 1),
    Pd = round(raw$n, 2),
    Pt = round(raw$t, 0),
    stringsAsFactors = FALSE
  )
  if (any(!nzchar(df$Date))) stop("AWSoM sat file has a row with no time.")
  # A few samples share a clock second and differ only in milliseconds.
  # The hourly mean keeps the later one.
  df = df[order(df$Date), ]
  n_drop = sum(duplicated(df$Date))
  if (n_drop)
    df = df[!duplicated(df$Date, fromLast = TRUE), ]
  rownames(df) = NULL
  attr(df, "n_drop") = n_drop
  df
}

# Array fields aligned with t. Scalars (bz_frame, cal_until) are not series.
series_fields = function(block) {
  n = length(block$t)
  names(block)[vapply(block, function(x) {
    (is.atomic(x) || is.list(x)) && length(x) == n
  }, logical(1))]
}

as_numeric_col = function(x) {
  if (is.list(x)) x = unlist(lapply(x, function(z) if (is.null(z)) NA else z))
  as.numeric(x)
}

# API name -> column title used in Data/Processed/OMNI_hourly.csv
omni_name = c(time_utc = "Date", b = "B", bz = "Bz", v = "Vs", n = "Pd", T = "Pt")
omni_order = c("Date", "B", "Bz", "Vs", "Pd", "Pt")

apply_omni_names = function(df) {
  hit = intersect(names(df), names(omni_name))
  names(df)[match(hit, names(df))] = omni_name[hit]
  front = intersect(omni_order, names(df))
  df[, c(front, setdiff(names(df), front)), drop = FALSE]
}

series_table = function(block) {
  fields = series_fields(block)
  if (!"t" %in% fields) stop("Series has no time vector t.")
  fields = unique(c("t", setdiff(fields, "t")))
  df = as.data.frame(lapply(setNames(fields, fields), function(nm) as_numeric_col(block[[nm]])),
                     stringsAsFactors = FALSE)
  df$time_utc = format(
    as.POSIXct(df$t, origin = "1970-01-01", tz = "UTC"),
    "%Y-%m-%d %H:%M:%S", tz = "UTC")
  df$t = NULL
  df[, c("time_utc", setdiff(names(df), "time_utc")), drop = FALSE]
}

doc = read_rtsw(rtsw_url)
stopifnot(is.list(doc$obs))

model = read_awsom_sat(awsom_sat_url)

obs = series_table(doc$obs)
if (!is.null(doc$obs$cal_until) && length(doc$obs$cal_until) == 1L) {
  cal_until = as.POSIXct(doc$obs$cal_until, origin = "1970-01-01", tz = "UTC")
  obs$provisional = as.POSIXct(obs$time_utc, tz = "UTC") > cal_until
} else {
  cal_until = NULL
}
obs = apply_omni_names(obs)

target = if (!is.null(doc$target)) doc$target else "earth"
model_path = file.path(out_dir, sprintf("%s_model.csv", target))
obs_path   = file.path(out_dir, sprintf("%s_obs.csv", target))
write.csv(model, model_path, row.names = FALSE)
write.csv(obs, obs_path, row.names = FALSE)

span = function(df) sprintf("%s -> %s", df$Date[1], df$Date[nrow(df)])
cat(sprintf("target: %s\n", target))
cat(sprintf("model  %d rows  %s\n  columns: %s\n  %s\n",
            nrow(model), span(model),
            paste(names(model), collapse = ", "),
            model_path))
if (isTRUE(attr(model, "n_drop") > 0))
  cat(sprintf("model dropped %d repeated clock seconds\n", attr(model, "n_drop")))
cat(sprintf("obs    %d rows  %s\n  columns: %s\n  %s\n",
            nrow(obs), span(obs),
            paste(names(obs), collapse = ", "),
            obs_path))
if (!is.null(cal_until)) {
  cat(sprintf("obs provisional after %s: %d rows\n",
              format(cal_until, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
              sum(obs$provisional)))
}
cat("Units: B |B| nT, Bz nT, Vs km/s, Pd cm^-3, Pt Kelvin.\n")
cat(sprintf("AWSoM source: %s\n", awsom_sat_url))

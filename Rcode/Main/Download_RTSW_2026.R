# =============================================================================
# Download_RTSW_2026.R
# Earth solar-wind series from the Michigan CLEAR realtime AWSoM page
# https://solarwind.engin.umich.edu/rtsw/
# =============================================================================
#
# The page button "all" is only a plot window. The full series is one JSON:
#   https://solarwind.engin.umich.edu/rtsw/api/series/?t=earth
#
# Nothing is filtered by year. The script keeps every array in model and obs
# that lines up with the time vector, and writes the whole span the server
# returns. A one-element field (bz_frame, cal_until) is metadata, not a
# series: bz_frame is stored as a column when present, and cal_until marks
# obs rows that are still the provisional realtime feed.
#
# Columns that exist in Data/Processed/OMNI_hourly.csv are renamed to those
# titles, in that order: Date, B, Bz, Vs, Pd, Pt.
#   b -> B, bz -> Bz, v -> Vs, n -> Pd, T -> Pt, time_utc -> Date
# A series field with no OMNI name (db) stays under the server's name and
# is written after the OMNI columns. bz_frame and provisional are metadata.
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
out_dir = file.path(repo_root, "Data", "Origin", "RTSW")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

read_rtsw = function(url) {
  tmp = tempfile(fileext = ".json")
  on.exit(unlink(tmp), add = TRUE)
  status = system2("curl", c("--max-time", "90", "-fsSL", url, "-o", tmp))
  if (!identical(status, 0L)) {
    message("Verified download failed; retrying without certificate check.")
    status = system2("curl", c("--max-time", "90", "-kfsSL", url, "-o", tmp))
  }
  if (!identical(status, 0L)) stop("curl download failed, status ", status)
  fromJSON(tmp)
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
stopifnot(is.list(doc$model), is.list(doc$obs))

model = apply_omni_names(series_table(doc$model))
if (!is.null(doc$model$bz_frame) && length(doc$model$bz_frame) == 1L)
  model$bz_frame = as.character(doc$model$bz_frame)

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
cat(sprintf("obs    %d rows  %s\n  columns: %s\n  %s\n",
            nrow(obs), span(obs),
            paste(names(obs), collapse = ", "),
            obs_path))
if (!is.null(cal_until)) {
  cat(sprintf("obs provisional after %s: %d rows\n",
              format(cal_until, "%Y-%m-%d %H:%M:%S", tz = "UTC"),
              sum(obs$provisional)))
}
cat("Units: B |B| nT, Bz nT, Vs km/s, Pd cm^-3, Pt Kelvin; db (model only) nT.\n")

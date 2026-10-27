# =============================================================================
# Write_VARX_Pages.R
# HTML for the realtime forecast and the baseline scores.
# The pages show one model: the 24-hour VARX with a searched 25-29 day lag.
# In-situ measurements and the AWSoM Earth model are the two solar-wind inputs.
# Persistence is the reference, not a second forecast model.
# =============================================================================

fmt_utc = function(x) format(x, "%Y-%m-%d %H:%M:%S", tz = "UTC")

parse_utc = function(x) {
  x = trimws(as.character(x))
  date_only = grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  x[date_only] = paste0(x[date_only], " 00:00:00")
  out = as.POSIXct(x, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  if (anyNA(out)) stop(sprintf("Unparsed timestamps, first: %s", x[which(is.na(out))[1]]))
  out
}

html_escape = function(x) {
  x = as.character(x)
  x = gsub("&", "&amp;", x, fixed = TRUE)
  x = gsub("<", "&lt;", x, fixed = TRUE)
  x = gsub(">", "&gt;", x, fixed = TRUE)
  x
}
fmt_num = function(x, digits) sprintf(paste0("%.", digits, "f"), x)

# Names on the RTSW page: "Earth model against in-situ measurements",
# and the coloured line is "the realtime AWSoM run".
lab_meas = "In-situ measurements"
lab_awso = "AWSoM Earth model"
lab_boost = "Boosted"

legend_html = function(items) {
  paste0('<div class="legend">', paste(vapply(seq_along(items), function(i) {
    it = items[[i]]
    fade = if (is.null(it$opacity)) "" else sprintf(";opacity:%.2f", it$opacity)
    sw = if (isTRUE(it$dash)) {
      sprintf('<i class="sw dash" style="border-top-color:%s%s"></i>', it$color, fade)
    } else {
      sprintf('<i class="sw" style="background:%s%s"></i>', it$color, fade)
    }
    paste0("<span>", sw, html_escape(it$label), "</span>")
  }, ""), collapse = ""), "</div>")
}

svg_x = function(x, xmin, xmax, x0, x1) x0 + (x - xmin) / (xmax - xmin) * (x1 - x0)
svg_y = function(y, ymin, ymax, y0, y1) y1 - (y - ymin) / (ymax - ymin) * (y1 - y0)
svg_pts = function(x, y, xmin, xmax, ymin, ymax, x0, x1, y0, y1) {
  paste(sprintf("%.1f,%.1f",
                svg_x(x, xmin, xmax, x0, x1),
                svg_y(y, ymin, ymax, y0, y1)), collapse = " ")
}
svg_runs = function(x, y, xmin, xmax, ymin, ymax, x0, x1, y0, y1, color, width,
                    dash = FALSE, cls = "", opacity = 1) {
  ok = is.finite(x) & is.finite(y)
  if (!any(ok)) return("")
  ix = which(ok)
  gap = rep(FALSE, length(ix))
  if (length(ix) > 1L) {
    dx = diff(x[ix])
    step = dx[is.finite(dx) & dx > 0]
    # A missing stretch is wider than the usual spacing. Do not draw across it.
    gap[-1L] = diff(ix) != 1L
    if (length(step)) gap[-1L] = gap[-1L] | dx > min(step) * 1.5
  }
  run_id = cumsum(c(1L, gap[-1L]))
  dash_attr = if (dash) ' stroke-dasharray="6 3"' else ""
  class_attr = if (nzchar(cls)) sprintf(' class="%s"', cls) else ""
  hide = if (grepl("opt", cls, fixed = TRUE)) ' style="display:none"' else ""
  fade = if (is.finite(opacity) && opacity < 1) sprintf(' stroke-opacity="%.2f"', opacity) else ""
  paste(vapply(split(ix, run_id), function(seg) {
    if (length(seg) < 2L) {
      if (length(seg) != 1L) return("")
      return(sprintf(
        '<circle%s cx="%.1f" cy="%.1f" r="3.5" fill="%s"%s%s/>',
        class_attr,
        svg_x(x[seg], xmin, xmax, x0, x1),
        svg_y(y[seg], ymin, ymax, y0, y1),
        color, fade, hide))
    }
    sprintf('<polyline%s fill="none" stroke="%s" stroke-width="%s"%s%s%s points="%s"/>',
            class_attr, color, width, dash_attr, fade, hide,
            svg_pts(x[seg], y[seg], xmin, xmax, ymin, ymax, x0, x1, y0, y1))
  }, ""), collapse = "")
}

toggle_btns = function(boost = FALSE) paste0(
  '<div class="togs">',
  '<button type="button" class="tog" data-series="insitu" aria-pressed="false">',
  '<i class="sw" style="background:#1d4ed8"></i>In-situ</button>',
  '<button type="button" class="tog" data-series="awsom" aria-pressed="false">',
  '<i class="sw" style="background:#c2410c"></i>AWSoM</button>',
  if (boost) paste0(
    '<button type="button" class="tog" data-series="boost" aria-pressed="false">',
    '<i class="sw" style="background:#0e7490"></i>Boosted</button>') else "",
  '</div>'
)

forecast_svg_chart = function(series, xmin, xmax, ymin, ymax, yticks,
                              xtick_x, xtick_lab, ylab, xlab, aria) {
  W = 760; H = 300; x0 = 52; x1 = 730; y0 = 24; y1 = 248
  if (!is.finite(ymax) || !is.finite(ymin) || ymax <= ymin) ymax = ymin + 1
  grid = paste(vapply(yticks, function(g) {
    yy = svg_y(g, ymin, ymax, y0, y1)
    sprintf('<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="#e7e5e4"/>', x0, yy, x1, yy)
  }, ""), collapse = "")
  ylab_el = paste(vapply(yticks, function(g) {
    yy = svg_y(g, ymin, ymax, y0, y1)
    sprintf('<text x="44" y="%.1f" text-anchor="end" font-size="12" fill="#57534e">%s</text>',
            yy + 4, html_escape(g))
  }, ""), collapse = "")
  xlab_el = paste(vapply(seq_along(xtick_x), function(i) {
    xx = svg_x(xtick_x[i], xmin, xmax, x0, x1)
    sprintf('<text x="%.1f" y="268" text-anchor="middle" font-size="12" fill="#57534e">%s</text>',
            xx, html_escape(xtick_lab[i]))
  }, ""), collapse = "")
  lines = paste(vapply(series, function(s) {
    svg_runs(s$x, s$y, xmin, xmax, ymin, ymax, x0, x1, y0, y1, s$color,
             if (is.null(s$width)) "2.2" else s$width,
             dash = isTRUE(s$dash),
             cls = if (is.null(s$cls)) "" else s$cls,
             opacity = if (is.null(s$opacity)) 1 else s$opacity)
  }, ""), collapse = "")
  sprintf(paste0(
    '<svg viewBox="0 0 %d %d" role="img" aria-label="%s">',
    grid, ylab_el, xlab_el,
    '<text x="14" y="136" text-anchor="middle" font-size="12" fill="#57534e" transform="rotate(-90 14 136)">%s</text>',
    '%s',
    '<text x="%g" y="292" font-size="12" fill="#57534e">%s</text>',
    '</svg>'),
    W, H, html_escape(aria), html_escape(ylab), lines, x0, html_escape(xlab))
}

issued_models_note = function() paste0(
  '<p>Kp, ap, and F10.7 are drawn as three charts, and they are not three outputs of one XGBoost. ',
  'The next 24 hours of Kp are vanilla XGBoost, trained on recent Kp and near-Earth solar wind. ',
  'From 27 hours through 27 days, Kp is the in-situ regression stepped forward on the closest OMNI cycle. ',
  'ap is the Bartels amplitude of that same Kp, so it follows the XGBoost and then the in-situ regression. ',
  'F10.7 does not use that XGBoost. The XGBoost reads geomagnetic and solar-wind inputs and predicts Kp. ',
  'The 10.7 cm flux is a daily solar measurement, so its forecast is a ridge regression on the flux&rsquo;s own history and the 27-day solar rotation. ',
  'On 2019&ndash;2022 that ridge beat repeating the latest measured flux.</p>'
)

range_switch = function() paste0(
  '<div class="togs range-switch">',
  '<button type="button" class="tog range-btn" data-range="24" aria-pressed="true">Next 24 hours</button>',
  '<button type="button" class="tog range-btn" data-range="27" aria-pressed="false">Next 27 days</button>',
  '</div>'
)

view_panel = function(metric, range = NULL, html, hidden = FALSE) paste0(
  '<div class="metric-panel" data-metric="', metric, '"',
  if (is.null(range) || !nzchar(range)) "" else paste0(' data-range="', range, '"'),
  if (hidden) " hidden" else "", ">", html, "</div>"
)

metric_switch = function() paste0(
  '<div class="togs metric-switch">',
  '<button type="button" class="tog metric-btn" data-metric="kp" aria-pressed="true">Kp</button>',
  '<button type="button" class="tog metric-btn" data-metric="ap" aria-pressed="false">ap</button>',
  '<button type="button" class="tog metric-btn" data-metric="f107" aria-pressed="false">F10.7</button>',
  '</div>',
  '<p>One quantity at a time. Kp is the geomagnetic index, ap is its Bartels amplitude, and F10.7 is the 10.7 cm solar flux.</p>',
  issued_models_note()
)

chart_script = function() paste0(
  '<script>',
  'function applyPanels(){',
  'var metricBtn=document.querySelector("button.metric-btn[aria-pressed=\\"true\\"]");',
  'var rangeBtn=document.querySelector("button.range-btn[aria-pressed=\\"true\\"]");',
  'var metric=metricBtn?metricBtn.getAttribute("data-metric"):null;',
  'var range=rangeBtn?rangeBtn.getAttribute("data-range"):null;',
  'var panels=document.querySelectorAll(".metric-panel");',
  'for(var j=0;j<panels.length;j++){',
  'var p=panels[j];',
  'var okM=!!metric&&p.getAttribute("data-metric")===metric;',
  'var pr=p.getAttribute("data-range");',
  'var okR=!range||!pr||pr===range;',
  'p.hidden=!(okM&&okR);}}',
  'document.addEventListener("click",function(ev){',
  'var m=ev.target.closest("button.metric-btn, button.range-btn");',
  'if(m){',
  'var sel=m.classList.contains("range-btn")?"button.range-btn":"button.metric-btn";',
  'var buttons=document.querySelectorAll(sel);',
  'for(var i=0;i<buttons.length;i++)buttons[i].setAttribute("aria-pressed",buttons[i]===m?"true":"false");',
  'applyPanels();',
  'if(m.classList.contains("metric-btn")){',
  'var winds=document.querySelectorAll("button.wind-btn");',
  'for(var k=0;k<winds.length;k++)winds[k].setAttribute("aria-pressed","false");',
  'var wind=document.querySelector(".wind-block");',
  'if(wind)wind.hidden=true;}',
  'return;}',
  'var w=ev.target.closest("button.wind-btn");',
  'if(w){',
  'var on=w.getAttribute("aria-pressed")!=="true";',
  'w.setAttribute("aria-pressed",on?"true":"false");',
  'var block=document.querySelector(".wind-block");',
  'if(block)block.hidden=!on;',
  'if(on){',
  'var metrics=document.querySelectorAll("button.metric-btn");',
  'for(var i=0;i<metrics.length;i++)metrics[i].setAttribute("aria-pressed","false");',
  'applyPanels();}',
  'return;}',
  'var b=ev.target.closest("button.tog");if(!b)return;',
  'var on=b.getAttribute("aria-pressed")!=="true";',
  'b.setAttribute("aria-pressed",on?"true":"false");',
  'var box=b.closest("figure");if(!box)return;',
  'var lines=box.querySelectorAll("."+b.getAttribute("data-series"));',
  'for(var i=0;i<lines.length;i++)lines[i].style.display=on?"":"none";',
  '});',
  '</script>'
)

page_css = function(extra = "") {
  paste0(
    'body{margin:0;background:#fafaf9;color:#1c1917;font:16px/1.45 system-ui,sans-serif;}',
    'main{width:100%;max-width:980px;margin:0 auto;padding:28px 20px 64px;box-sizing:border-box;}',
    'nav{display:flex;gap:16px;margin:0 0 18px;font-size:14px;}',
    'nav a, .downloads a{color:#1d4ed8;text-decoration:none;} nav a:hover, .downloads a:hover{text-decoration:underline;}',
    'nav .here{color:#1c1917;font-weight:650;}',
    'h1{font-size:28px;font-weight:650;margin:0 0 6px;}',
    'h2{font-size:18px;margin:32px 0 8px;}',
    'p,figcaption{color:#57534e;}',
    '.meta{display:flex;flex-wrap:wrap;gap:12px;margin:18px 0 8px;}',
    '.meta div{background:#fff;border:1px solid #e7e5e4;border-radius:10px;padding:12px 14px;min-width:180px;}',
    '.meta b{display:block;font-size:20px;color:#1c1917;}',
    '.legend{display:flex;gap:12px 16px;flex-wrap:wrap;font-size:13px;margin:2px 0 6px;color:#44403c;}',
    '.sw{display:inline-block;width:18px;height:3px;margin-right:6px;vertical-align:middle;}',
    '.sw.dash{height:0;border-top:2px dashed;background:transparent;}',
    'table{border-collapse:collapse;width:100%;background:#fff;font-size:14px;font-variant-numeric:tabular-nums;}',
    'th,td{border-bottom:1px solid #e7e5e4;padding:7px 8px;text-align:right;white-space:nowrap;}',
    'th:first-child,td:first-child{text-align:left;}',
    'th{color:#57534e;font-weight:600;}',
    '.wrap{overflow-x:auto;max-width:100%;border:1px solid #e7e5e4;border-radius:10px;}',
    '.wrap.scroll{max-height:260px;overflow:auto;}',
    '.wrap.scroll table{border-collapse:separate;border-spacing:0;}',
    '.wrap.scroll thead th{position:sticky;top:0;background:#fff;z-index:1;box-shadow:inset 0 -1px 0 #e7e5e4;}',
    'svg{width:100%;max-width:100%;height:auto;display:block;}',
    '.togs{display:flex;flex-wrap:wrap;gap:8px;margin:4px 0 8px;}',
    'button.tog{border:1px solid #e7e5e4;background:#fff;border-radius:999px;padding:3px 10px;font:13px system-ui,sans-serif;color:#44403c;cursor:pointer;}',
    'button.tog[aria-pressed="true"]{background:#1c1917;color:#fff;border-color:#1c1917;}',
    extra)
}

write_html = function(path, page) {
  if (Sys.getlocale("LC_CTYPE") == "C") {
    try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
  }
  raw = charToRaw(enc2utf8(page))
  con = file(path, open = "wb")
  writeBin(raw, con)
  close(con)
  invisible(path)
}

outlook_baseline_section = function(ev) {
  col_fc = "#047857"
  col_obs = "#1d4ed8"
  col_mod = "#c2410c"
  col_rec = "#7c3aed"
  col_true = "#a8a29e"
  met = ev$metrics
  order_in = c("Best matched cycle", "Boosted", "In-situ", "AWSoM")
  met = met[order(met$lead_end_hr, match(met$input, order_in)), ]
  fmt_cor = function(x) if (!is.finite(x)) "" else fmt_num(x, 3)
  skill_td = function(x) {
    if (!is.finite(x)) return("<td></td>")
    color = if (x > 0) "#087F23" else if (x < 0) "#b91c1c" else "#57534e"
    sprintf('<td style="color:%s">%+.3f</td>', color, x)
  }
  metric_rows = function(frame) {
    paste(vapply(seq_len(nrow(frame)), function(i) {
      sprintf("<tr><td>%dh</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td>%s</tr>",
              frame$lead_end_hr[i], html_escape(frame$input[i]), frame$n[i],
              fmt_num(frame$RMSE[i], 3), fmt_num(frame$MAE[i], 3), fmt_cor(frame$Correlation[i]),
              skill_td(frame$skill_RMSE[i]))
    }, ""), collapse = "\n")
  }
  endpoint = met[met$lead_end_hr == max(met$lead_end_hr), ]

  series_svg = function(x, y_true_v, y_hat_v, y_rec_v, y_obs_v, y_mod_v, tick_by) {
    xmin = min(x); xmax = max(x)
    if (xmin == xmax) xmax = xmin + 1
    ymax_s = max(9, ceiling(max(c(y_true_v, y_hat_v, y_rec_v, y_obs_v, y_mod_v), na.rm = TRUE)))
    ymin_s = 0
    yticks_s = if (ymax_s <= 9) c(0, 3, 6, 9) else pretty(c(0, ymax_s), n = 4)
    sw = 760; sh = 270; sx0 = 44; sx1 = 730; sy0 = 16; sy1 = 214
    grid_s = paste(vapply(yticks_s, function(g) {
      yy = svg_y(g, ymin_s, ymax_s, sy0, sy1)
      sprintf('<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="#e7e5e4"/>', sx0, yy, sx1, yy)
    }, ""), collapse = "")
    ylab_s = paste(vapply(yticks_s, function(g) {
      sprintf('<text x="36" y="%.1f" text-anchor="end" font-size="12" fill="#57534e">%g</text>',
              svg_y(g, ymin_s, ymax_s, sy0, sy1) + 4, g)
    }, ""), collapse = "")
    t_min = as.POSIXct(xmin, origin = "1970-01-01", tz = "UTC")
    t_max = as.POSIXct(xmax, origin = "1970-01-01", tz = "UTC")
    if (tick_by == "month") {
      ticks = seq(as.POSIXct(format(t_min, "%Y-%m-01"), tz = "UTC"), t_max, by = "month")
      lab = format(ticks, "%b", tz = "UTC")
    } else {
      ticks = seq(as.POSIXct(format(t_min, "%Y-%m-%d"), tz = "UTC"), t_max, by = "2 days")
      lab = format(ticks, "%m-%d", tz = "UTC")
    }
    keep_tick = ticks >= t_min & ticks <= t_max
    ticks = ticks[keep_tick]
    lab = lab[keep_tick]
    xlab = if (!length(ticks)) "" else paste(vapply(seq_along(ticks), function(i) {
      sprintf('<text x="%.1f" y="236" text-anchor="middle" font-size="11" fill="#57534e">%s</text>',
              svg_x(as.numeric(ticks[i]), xmin, xmax, sx0, sx1), lab[i])
    }, ""), collapse = "")
    sprintf(paste0(
      '<svg viewBox="0 0 %d %d" role="img" aria-label="27 day outlook score">',
      '<text x="14" y="110" text-anchor="middle" font-size="12" fill="#57534e" transform="rotate(-90 14 110)">Kp</text>',
      grid_s, ylab_s, xlab, "%s%s%s%s%s</svg>"),
      sw, sh,
      svg_runs(x, y_true_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_true, "1.2"),
      svg_runs(x, y_hat_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_fc, "2.2"),
      svg_runs(x, y_rec_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_rec, "1.8", opacity = 0.45),
      svg_runs(x, y_obs_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_obs, "1.8", cls = "opt insitu", opacity = 0.55),
      svg_runs(x, y_mod_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_mod, "1.8", cls = "opt awsom", opacity = 0.55))
  }
  legend = paste0(legend_html(list(
    list(color = col_true, label = "Observed Kp"),
    list(color = col_fc, label = "Boosted"),
    list(color = col_rec, label = "Best matched cycle", opacity = 0.7)
  )), toggle_btns())
  chart_leads = c(24L, 168L, 648L)
  chart_title = c("24" = "24 hours ahead", "168" = "7 days ahead", "648" = "27 days ahead")
  chart_note = c(
    "24" = "Each point is the 3-hour Kp bin ending 24 hours after issue. Gray is Kp observed later. Green is Boosted. Violet is the Kp that followed the closest aligned OMNI cycle. The buttons add In-situ and AWSoM.",
    "168" = "The 3-hour period that ends 7 days after the forecast was issued.",
    "648" = "The 3-hour period that ends 27 days after the forecast was issued. This is the endpoint of the 27-day outlook."
  )
  charts = paste(vapply(chart_leads, function(lead) {
    j = match(lead, ev$leads_hr)
    use = is.finite(ev$y_true_bin[, j]) & is.finite(ev$y_boost_bin[, j])
    x = as.numeric(ev$tau[use]) + lead * 3600
    y_true_v = ev$y_true_bin[use, j]
    y_hat_v = ev$y_boost_bin[use, j]
    y_rec_v = ev$y_match_bin[use, j]
    y_obs_v = ev$y_insitu_bin[use, j]
    y_mod_v = ev$y_awsom_bin[use, j]
    recent = x >= (max(x) - 14 * 86400)
    title = chart_title[[as.character(lead)]]
    full = paste0("<figure><figcaption>", title, ", whole test period</figcaption><p>",
                  html_escape(chart_note[[as.character(lead)]]), "</p>", legend,
                  series_svg(x, y_true_v, y_hat_v, y_rec_v, y_obs_v, y_mod_v, "month"),
                  "</figure>")
    zoom = paste0("<figure><figcaption>", title, ", latest 14 days of scored forecasts</figcaption><p>",
                  "The date is the end of that 3-hour period. The same chart, kept to the latest 14 of those periods.</p>",
                  legend,
                  series_svg(x[recent], y_true_v[recent], y_hat_v[recent],
                             y_rec_v[recent], y_obs_v[recent], y_mod_v[recent], "day"),
                  "</figure>")
    paste0(full, zoom)
  }, ""), collapse = "")
  paste0(
    '<h2>Next 27 days, scored after the fact</h2>',
    '<p>Each lead is one 3-hour Kp bin, scored once that period has ended. ',
    'The 3-hour lead reaches the latest Kp. The 27-day lead reaches Kp from about four weeks earlier. ',
    'The reference is the Kp that followed the closest aligned OMNI cycle from 1995 through January 2026. ',
    'Skill = 1 - RMSE / RMSE of that cycle.</p>',
    '<p>The forecast is the 24-hour model stepped forward one day at a time. The first day is that 24-hour forecast. ',
    'Later days use the solar wind that followed the closest cycle, and the predicted Kp becomes the next day&rsquo;s recent Kp. ',
    'Green is Boosted. Violet is the closest cycle. The buttons add In-situ and AWSoM.</p>',
    '<h2>Endpoint criteria</h2>',
    '<div class="wrap"><table><thead><tr>',
    '<th>Lead</th><th>Input</th><th>n</th><th>RMSE</th><th>MAE</th><th>Correlation</th><th>Skill</th>',
    '</tr></thead><tbody>', metric_rows(endpoint), '</tbody></table></div>',
    '<h2>Criteria by lead</h2>',
    '<p>Leads of the 27-day outlook, from 3 hours through 27 days. Later leads have fewer forecasts because those periods have not all ended yet.</p>',
    '<div class="wrap scroll"><table><thead><tr>',
    '<th>Lead</th><th>Input</th><th>n</th><th>RMSE</th><th>MAE</th><th>Correlation</th><th>Skill</th>',
    '</tr></thead><tbody>', metric_rows(met), '</tbody></table></div>',
    charts
  )
}

write_realtime_page = function(path, forecast, drivers, kp_bin_start, kp_bin_value,
                               issue_time, generated_at, ab_live, lag_obs_hr, lag_mod_hr,
                               kp_step_hr = 3L, outlook = NULL,
                               kp_obs_time = NULL, kp_obs_value = NULL, kp_obs_source = NULL,
                               kp_history_note = NULL, f107 = NULL, f107_measured = NULL) {
  col_obs = "#1d4ed8"
  col_mod = "#c2410c"
  col_boost = "#047857"
  col_pers = "#6b7280"
  f24o = forecast[forecast$input == "Earth_obs" & forecast$horizon_hr == 24, ]
  f24m = forecast[forecast$input == "model_affine" & forecast$horizon_hr == 24, ]
  f24o = f24o[order(f24o$lead_end_hr), ]
  f24m = f24m[order(f24m$lead_end_hr), ]
  if (!nrow(f24o) || !nrow(f24m)) stop("Forecast rows are missing.")
  pers = f24o$y_persistence[1]
  has_issued = !is.null(outlook) && "y_issued" %in% names(outlook) &&
    nrow(outlook) >= length(f24o$y_pred) &&
    all(is.finite(outlook$y_issued[seq_along(f24o$y_pred)]))
  y_boost24 = if (!is.null(outlook) && nrow(outlook) >= length(f24o$y_pred))
    outlook$y_boost[seq_along(f24o$y_pred)] else f24o$y_pred
  y_line24 = if (has_issued) outlook$y_issued[seq_along(f24o$y_pred)] else y_boost24
  issued_label = if (has_issued) "Issued" else "Boosted"

  W = 760; H = 300; x0 = 52; x1 = 730; y0 = 24; y1 = 248
  leads = f24o$lead_end_hr
  xmin = 0; xmax = 27; ymin = 0
  ymax = max(c(f24o$y_pred, f24m$y_pred, y_line24, pers), na.rm = TRUE)
  ymax = if (ymax <= 4) 4 else 9
  yticks = if (ymax <= 4) 0:4 else c(0, 3, 6, 9)
  grid = paste(vapply(yticks, function(g) {
    yy = svg_y(g, ymin, ymax, y0, y1)
    sprintf('<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="#e7e5e4"/>', x0, yy, x1, yy)
  }, ""), collapse = "")
  ylab = paste(vapply(yticks, function(g) {
    yy = svg_y(g, ymin, ymax, y0, y1)
    sprintf('<text x="44" y="%.1f" text-anchor="end" font-size="12" fill="#57534e">%g</text>', yy + 4, g)
  }, ""), collapse = "")
  xlab = paste(vapply(leads, function(h) {
    xx = svg_x(h, xmin, xmax, x0, x1)
    sprintf('<text x="%.1f" y="268" text-anchor="middle" font-size="12" fill="#57534e">%dh</text>', xx, h)
  }, ""), collapse = "")
  pers_y = svg_y(pers, ymin, ymax, y0, y1)
  forecast_svg = sprintf(paste0(
    '<svg viewBox="0 0 %d %d" role="img" aria-label="Kp forecast">',
    grid, ylab, xlab,
    '<text x="14" y="136" text-anchor="middle" font-size="12" fill="#57534e" transform="rotate(-90 14 136)">Kp</text>',
    '<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="%s" stroke-dasharray="5 4" stroke-width="1.5"/>',
    '%s%s%s',
    '<text x="%g" y="292" font-size="12" fill="#57534e">Hours after the forecast was issued</text>',
    '</svg>'),
    W, H,
    x0, pers_y, x1, pers_y, col_pers,
    svg_runs(leads, y_line24, xmin, xmax, ymin, ymax, x0, x1, y0, y1, col_boost, "2.4"),
    svg_runs(leads, f24o$y_pred, xmin, xmax, ymin, ymax, x0, x1, y0, y1, col_obs, "2.2", cls = "opt insitu", opacity = 0.55),
    svg_runs(leads, f24m$y_pred, xmin, xmax, ymin, ymax, x0, x1, y0, y1, col_mod, "2.2", cls = "opt awsom", opacity = 0.55),
    x0)

  ap_of = function(kp) round(kp_to_ap(round(kp, 3)), 2)
  ap24 = ap_of(y_line24)
  ap24_obs = ap_of(f24o$y_pred)
  ap24_mod = ap_of(f24m$y_pred)
  ap_pers = ap_of(pers)
  obs_kp_by_end = if (length(kp_obs_time)) {
    stats::setNames(kp_obs_value, sprintf("%.0f", as.numeric(kp_obs_time)))
  } else numeric()
  got24 = vapply(seq_along(leads), function(i) {
    sprintf("%.0f", as.numeric(parse_utc(f24o$target_bin_end[i]))) %in% names(obs_kp_by_end)
  }, logical(1))
  meas24_kp = rep(NA_real_, length(leads))
  meas24_ap = rep(NA_real_, length(leads))
  if (any(got24)) {
    keys = sprintf("%.0f", as.numeric(parse_utc(f24o$target_bin_end)))
    meas24_kp[got24] = vapply(keys[got24], function(key) obs_kp_by_end[[key]], numeric(1))
    meas24_ap[got24] = ap_of(meas24_kp[got24])
  }
  forecast_rows_html = paste(vapply(seq_along(leads), function(i) {
    sprintf("<tr><td>%s</td><td>%dh</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
            html_escape(f24o$target_bin_end[i]), leads[i],
            fmt_num(y_line24[i], 3),
            if (got24[i]) fmt_num(meas24_kp[i], 3) else "",
            fmt_num(f24o$y_pred[i], 3),
            fmt_num(f24m$y_pred[i], 3), fmt_num(pers, 3))
  }, ""), collapse = "\n")
  ap24_rows_html = paste(vapply(seq_along(leads), function(i) {
    sprintf("<tr><td>%s</td><td>%dh</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
            html_escape(f24o$target_bin_end[i]), leads[i],
            fmt_num(ap24[i], 2),
            if (got24[i]) fmt_num(meas24_ap[i], 2) else "",
            fmt_num(ap24_obs[i], 2),
            fmt_num(ap24_mod[i], 2), fmt_num(ap_pers, 2))
  }, ""), collapse = "\n")
  ap24_vals = c(ap24, ap24_obs, ap24_mod, ap_pers, meas24_ap)
  ap24_hi = max(ap24_vals, na.rm = TRUE)
  ap24_ticks = pretty(c(0, ap24_hi), n = 4)
  ap24_ticks = ap24_ticks[ap24_ticks >= 0]
  ap24_hi = max(ap24_ticks)
  ap24_series = list(
    list(x = c(xmin, xmax), y = c(ap_pers, ap_pers), color = col_pers, width = "1.5", dash = TRUE),
    list(x = leads, y = ap24, color = col_boost, width = "2.4"),
    list(x = leads, y = ap24_obs, color = col_obs, width = "2.2", cls = "opt insitu", opacity = 0.55),
    list(x = leads, y = ap24_mod, color = col_mod, width = "2.2", cls = "opt awsom", opacity = 0.55)
  )
  if (any(got24)) {
    ap24_series = c(list(list(
      x = leads[got24], y = meas24_ap[got24], color = "#44403c", width = "2.2"
    )), ap24_series)
  }
  ap24_svg = forecast_svg_chart(
    ap24_series, xmin, xmax, 0, ap24_hi, ap24_ticks,
    leads, paste0(leads, "h"), "ap", "Hours after the forecast was issued",
    "24 hour ap forecast")

  kp27_html = "<p>The 27-day outlook is not available for this issue.</p>"
  ap27_html = kp27_html
  axis_xmin = NULL
  if (!is.null(outlook) && nrow(outlook)) {
    day_end = parse_utc(outlook$day_end)
    x_day = as.numeric(day_end)
    has_obs = length(kp_obs_time) > 0 && any(is.finite(kp_obs_value))
    W2 = 760; H2 = 300; ox0 = 52; ox1 = 730; oy0 = 24; oy1 = 248
    tick_start = as.POSIXct(format(issue_time, "%Y-%m-%d"), tz = "UTC")
    oxmin = as.numeric(tick_start) - 18 * 3600
    oxmax = max(x_day)
    if (has_obs) {
      in_win = is.finite(as.numeric(kp_obs_time)) & is.finite(kp_obs_value) &
        as.numeric(kp_obs_time) >= as.numeric(tick_start) & as.numeric(kp_obs_time) <= oxmax
      kp_obs_time = kp_obs_time[in_win]
      kp_obs_value = kp_obs_value[in_win]
      has_obs = length(kp_obs_time) > 0
    }
    col_boost = "#047857"
    col_match = "#7c3aed"
    col_seen = "#44403c"
    oymin = 0
    y_line27 = if ("y_issued" %in% names(outlook) && all(is.finite(outlook$y_issued)))
      outlook$y_issued else outlook$y_boost
    y_for_scale = c(outlook$y_insitu, outlook$y_awsom, y_line27, outlook$y_boost, outlook$y_matched)
    if (has_obs) y_for_scale = c(y_for_scale, kp_obs_value)
    oymax = max(y_for_scale, na.rm = TRUE)
    oymax = if (!is.finite(oymax) || oymax <= 4) 4 else 9
    oyticks = if (oymax <= 4) 0:4 else c(0, 3, 6, 9)
    ogrid = paste(vapply(oyticks, function(g) {
      yy = svg_y(g, oymin, oymax, oy0, oy1)
      sprintf('<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="#e7e5e4"/>', ox0, yy, ox1, yy)
    }, ""), collapse = "")
    oylab = paste(vapply(oyticks, function(g) {
      yy = svg_y(g, oymin, oymax, oy0, oy1)
      sprintf('<text x="44" y="%.1f" text-anchor="end" font-size="12" fill="#57534e">%g</text>', yy + 4, g)
    }, ""), collapse = "")
    tick_end = as.POSIXct(oxmax, origin = "1970-01-01", tz = "UTC")
    oticks = seq(tick_start, tick_end, by = "7 days")
    oticks = oticks[as.numeric(oticks) >= oxmin & as.numeric(oticks) <= oxmax]
    if (!length(oticks) || as.numeric(tail(oticks, 1)) < oxmax - 2 * 86400) {
      oticks = c(oticks, as.POSIXct(format(tick_end, "%Y-%m-%d"), tz = "UTC"))
    }
    oxlab = paste(vapply(oticks, function(t) {
      xx = svg_x(as.numeric(t), oxmin, oxmax, ox0, ox1)
      sprintf('<text x="%.1f" y="268" text-anchor="middle" font-size="12" fill="#57534e">%s</text>',
              xx, trimws(format(t, "%b %e", tz = "UTC")))
    }, ""), collapse = "")
    line_of = function(y, color, width, cls = "", opacity = 1) svg_runs(
      x_day, y, oxmin, oxmax, oymin, oymax, ox0, ox1, oy0, oy1, color, width,
      cls = cls, opacity = opacity)
    obs_pts = if (!has_obs) "" else svg_runs(
      as.numeric(kp_obs_time), kp_obs_value, oxmin, oxmax, oymin, oymax,
      ox0, ox1, oy0, oy1, col_seen, "2.2")
    outlook_svg = sprintf(paste0(
      '<svg viewBox="0 0 %d %d" role="img" aria-label="27 day Kp outlook">',
      ogrid, oylab, oxlab,
      '<text x="14" y="136" text-anchor="middle" font-size="12" fill="#57534e" transform="rotate(-90 14 136)">Kp</text>',
      '%s%s%s%s%s%s',
      '<text x="%g" y="292" font-size="12" fill="#57534e">Date (UTC)</text>',
      '</svg>'),
      W2, H2,
      obs_pts,
      line_of(y_line27, col_boost, "2.4"),
      line_of(outlook$y_matched, col_match, "1.8", opacity = 0.45),
      line_of(outlook$y_insitu, col_obs, "2.0", "opt insitu", 0.55),
      line_of(outlook$y_awsom, col_mod, "2.0", "opt awsom", 0.55),
      line_of(outlook$y_boost, "#0e7490", "2.0", "opt boost", 0.8),
      ox0)
    kp_r = round(y_line27, 3)
    ap_r = ap_of(y_line27)
    meas_t = as.numeric(kp_bin_start + kp_step_hr * 3600)
    meas_kp = kp_bin_value
    if (has_obs) {
      meas_t = c(meas_t, as.numeric(kp_obs_time))
      meas_kp = c(meas_kp, kp_obs_value)
    }
    meas_ok = is.finite(meas_t) & is.finite(meas_kp)
    meas_t = meas_t[meas_ok]
    meas_kp = meas_kp[meas_ok]
    if (length(meas_t)) {
      ix = vapply(split(seq_along(meas_t), sprintf("%.0f", meas_t)), function(i) i[length(i)], 1L)
      meas_t = meas_t[ix]
      meas_kp = meas_kp[ix]
      ordm = order(meas_t)
      meas_t = meas_t[ordm]
      meas_kp = meas_kp[ordm]
    }
    meas_ap = if (length(meas_t)) ap_of(meas_kp) else numeric()
    on_chart = if (length(meas_t)) meas_t >= oxmin & meas_t <= oxmax else logical()
    ap_scale = c(ap_r, ap_of(outlook$y_insitu), ap_of(outlook$y_awsom),
                 ap_of(outlook$y_boost), meas_ap)
    ap_ymax = max(ap_scale, na.rm = TRUE)
    ap_ticks = pretty(c(0, ap_ymax), n = 4)
    ap_ticks = ap_ticks[ap_ticks >= 0]
    ap_ymax = max(ap_ticks)
    ap_series = list()
    if (any(on_chart)) {
      ap_series = list(list(
        x = meas_t[on_chart], y = meas_ap[on_chart], color = col_seen, width = "2.2"))
    }
    ap_series = c(ap_series, list(
      list(x = x_day, y = ap_r, color = col_boost, width = "2.4")
    ), list(
      list(x = x_day, y = ap_of(outlook$y_insitu), color = col_obs, width = "2.0",
           cls = "opt insitu", opacity = 0.55),
      list(x = x_day, y = ap_of(outlook$y_awsom), color = col_mod, width = "2.0",
           cls = "opt awsom", opacity = 0.55),
      list(x = x_day, y = ap_of(outlook$y_boost), color = "#0e7490", width = "2.0",
           cls = "opt boost", opacity = 0.8)
    ))
    ap_svg = forecast_svg_chart(
      ap_series, oxmin, oxmax, 0, ap_ymax, ap_ticks,
      as.numeric(oticks), trimws(format(oticks, "%b %e", tz = "UTC")),
      "ap", "Date (UTC)", "27 day ap outlook")
    bin_start = day_end - 3 * 3600
    day_key = format(bin_start, "%Y-%m-%d", tz = "UTC")
    daily_rows = paste(vapply(unique(day_key), function(d) {
      ix = which(day_key == d)
      sprintf("<tr><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
              html_escape(d), length(ix),
              fmt_num(mean(ap_r[ix]), 2), fmt_num(max(ap_r[ix]), 2),
              fmt_num(mean(kp_r[ix]), 3), fmt_num(max(kp_r[ix]), 3))
    }, ""), collapse = "\n")
    axis_xmin = oxmin
    axis_xmax = oxmax
    axis_xtick = as.numeric(oticks)
    axis_xlab = trimws(format(oticks, "%b %e", tz = "UTC"))
    meas_ap_by_end = if (length(meas_t)) stats::setNames(meas_ap, sprintf("%.0f", meas_t)) else numeric()
    meas_kp_by_end = if (length(meas_t)) stats::setNames(meas_kp, sprintf("%.0f", meas_t)) else numeric()
    outlook_rows = paste(vapply(seq_len(nrow(outlook)), function(i) {
      key = sprintf("%.0f", as.numeric(parse_utc(outlook$day_end[i])))
      got = key %in% names(meas_kp_by_end)
      sprintf("<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
              outlook$day[i], html_escape(outlook$day_end[i]),
              fmt_num(outlook$y_insitu[i], 3), fmt_num(outlook$y_awsom[i], 3),
              fmt_num(y_line27[i], 3), fmt_num(outlook$y_boost[i], 3),
              if (got) fmt_num(meas_kp_by_end[[key]], 3) else "",
              fmt_num(outlook$y_matched[i], 3))
    }, ""), collapse = "\n")
    ap_rows = paste(vapply(seq_len(nrow(outlook)), function(i) {
      key = sprintf("%.0f", as.numeric(parse_utc(outlook$day_end[i])))
      got = key %in% names(meas_kp_by_end)
      sprintf("<tr><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
              outlook$day[i], html_escape(outlook$day_end[i]),
              fmt_num(ap_r[i], 2), fmt_num(ap_of(outlook$y_boost[i]), 2),
              if (got) fmt_num(meas_ap_by_end[[key]], 2) else "",
              fmt_num(ap_of(outlook$y_insitu[i]), 2),
              fmt_num(ap_of(outlook$y_awsom[i]), 2))
    }, ""), collapse = "\n")
    obs_note = ""
    if (has_obs) {
      obs_span = if (length(kp_obs_time) == 1L) {
        paste0("The latest reading is ", format(kp_obs_time[1], "%d %b %H:%M", tz = "UTC"), " UTC. ")
      } else {
        paste0(format(min(kp_obs_time), "%d %b %H:%M", tz = "UTC"), " through ",
               format(max(kp_obs_time), "%d %b %H:%M", tz = "UTC"), " UTC. ")
      }
      obs_note = paste0(
        '<p>Gray is observed Kp from the day the forecast was issued. ',
        obs_span)
      if (identical(kp_obs_source, "NOAA")) {
        obs_note = paste0(obs_note,
          'These readings are NOAA planetary Kp, rounded to the nearest third. ')
      }
      obs_note = paste0(obs_note, 'Later periods have not ended, so they have no observed Kp yet.</p>')
    }
    legend_items = list(
      list(color = col_boost, label = if (has_issued) "Issued" else "Boosted"),
      list(color = col_match, label = "Best matched cycle", opacity = 0.7)
    )
    if (has_obs) {
      legend_items = c(list(list(color = col_seen, label = "Observed Kp")), legend_items)
    }
    ap_legend = list(list(color = col_boost, label = if (has_issued) "Issued ap" else "Boosted ap"))
    if (any(on_chart)) {
      ap_legend = c(list(list(color = col_seen, label = "Measured ap")), ap_legend)
    }
    meas_ap_at = stats::setNames(meas_ap, sprintf("%.0f", meas_t))
    issued_ap_note = ""
    if (length(meas_t)) {
      issued_ap_note = paste0(
        "Gray is measured ap, the Bartels value of Kp already on record. At issue time the last Kp was ",
        fmt_num(kp_bin_value[length(kp_bin_value)], 3),
        ", so the measured ap was ",
        fmt_num(ap_of(kp_bin_value[length(kp_bin_value)]), 2), ". ")
      later = meas_t > as.numeric(issue_time)
      if (any(later)) {
        j = which(later)[which.max(meas_t[later])]
        issued_ap_note = paste0(
          issued_ap_note, "The newest ended bin is ",
          fmt_utc(as.POSIXct(meas_t[j], origin = "1970-01-01", tz = "UTC")),
          ", Kp ", fmt_num(meas_kp[j], 3), ", ap ", fmt_num(meas_ap[j], 2), ". ")
      }
    }
    kp27_html = paste0(
      '<h2>Next 27 days</h2>',
      '<p>Each point is one 3-hour Kp, at the end of that period. The first point is 3 hours ahead. The last is 27 days ahead.</p>',
      if (has_issued) {
        paste0(
          '<p>Green is the issued forecast. The first 24 hours are the same XGBoost values as the 24-hour Kp view. ',
          'After that, each step is the in-situ regression moved forward one day at a time. ',
          'Later days use the solar wind that followed the closest aligned OMNI cycle, from 1995 through January 2026, and the predicted Kp becomes the next day&rsquo;s recent Kp. ',
          'The buttons add In-situ, AWSoM, and Boosted for every step. Violet is the Kp that followed the closest cycle, for comparison.</p>')
      } else {
        paste0(
          '<p>Green is Boosted, the issued forecast. It is the 24-hour model stepped forward one day at a time. ',
          'The first day is that 24-hour forecast. Later days use the solar wind that followed the closest aligned OMNI cycle, from 1995 through January 2026, and the predicted Kp becomes the next day&rsquo;s recent Kp. ',
          'Boosted uses AWSoM only when the latest observed Kp is high enough that earlier storms at that lead favored the blend. Otherwise it stays on the in-situ forecast. ',
          'The buttons add In-situ, AWSoM, and Boosted. Violet is the Kp that followed the closest cycle, for comparison.</p>')
      },
      obs_note,
      '<figure>',
      legend_html(legend_items),
      toggle_btns(TRUE),
      outlook_svg,
      '</figure>',
      '<div class="wrap scroll"><table><thead><tr>',
      '<th>Day</th><th>Forecast time (UTC)</th><th>In-situ</th><th>AWSoM</th><th>',
      if (has_issued) "Issued" else "Boosted",
      '</th><th>Boosted</th><th>Measured Kp</th><th>Best matched cycle</th>',
      '</tr></thead><tbody>', outlook_rows, '</tbody></table></div>'
    )
    ap27_html = paste0(
      '<h2>Next 27 days</h2>',
      '<p>ap is the 3-hour planetary amplitude of the issued Kp, on the Bartels scale used by GFZ and NOAA. ',
      'Between scale points it is interpolated in log(ap). It is not a separate forecast model. ',
      issued_ap_note,
      'Green is the ap of the issued Kp. The buttons add the ap of the in-situ, AWSoM, and boosted Kp. ',
      'The daily table gives the mean and the maximum of those 3-hour ap values whose period starts on that UT day, and the same two numbers for Kp.</p>',
      '<figure>',
      legend_html(ap_legend),
      toggle_btns(TRUE),
      ap_svg,
      '</figure>',
      '<div class="wrap"><table><thead><tr>',
      '<th>Date (UTC)</th><th>Bins</th><th>Ap mean</th><th>Ap max</th><th>Kp mean</th><th>Kp max</th>',
      '</tr></thead><tbody>', daily_rows, '</tbody></table></div>',
      '<div class="wrap scroll"><table><thead><tr>',
      '<th>Day</th><th>Forecast time (UTC)</th><th>Issued ap</th><th>Boosted ap</th><th>Measured ap</th><th>In-situ ap</th><th>AWSoM ap</th>',
      '</tr></thead><tbody>', ap_rows, '</tbody></table></div>'
    )
  }

  kp_end = kp_bin_start + kp_step_hr * 3600
  kp_rows_html = paste(vapply(seq_along(kp_bin_value), function(i) {
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>",
            html_escape(fmt_utc(kp_bin_start[i])),
            html_escape(fmt_utc(kp_end[i])),
            fmt_num(kp_bin_value[i], 3))
  }, ""), collapse = "\n")
  ap_rows_html = paste(vapply(seq_along(kp_bin_value), function(i) {
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>",
            html_escape(fmt_utc(kp_bin_start[i])),
            html_escape(fmt_utc(kp_end[i])),
            fmt_num(ap_of(kp_bin_value[i]), 2))
  }, ""), collapse = "\n")

  feature_meta = list(
    B = list(label = "Magnetic field strength", unit = "nT", digits = 2,
             note = "How strong the magnetic field in the solar wind is."),
    Bz = list(label = "Magnetic field, north-south", unit = "nT", digits = 2,
              note = "Below zero means the field points south. Southward field disturbs Earth more."),
    Vs = list(label = "Solar-wind speed", unit = "km/s", digits = 1,
              note = "How fast the solar wind is moving."),
    Pd = list(label = "Proton density", unit = "per cm3", digits = 2,
              note = "How many solar-wind protons are packed into a cubic centimeter."),
    Pt = list(label = "Temperature", unit = "K", digits = 0,
              note = "How hot the solar wind is.")
  )
  obs_d = drivers[drivers$input == "Earth_obs", ]
  mod_d = drivers[drivers$input == "model_affine", ]
  obs_d = obs_d[order(obs_d$hour_start), ]
  mod_d = mod_d[order(mod_d$hour_start), ]
  hours_x = seq_along(obs_d$hour_start)
  hour_lab = format(parse_utc(obs_d$hour_start), "%m-%d %H:%M", tz = "UTC")
  hour_tick = format(parse_utc(obs_d$hour_start), "%H:%M", tz = "UTC")
  charts = vapply(names(feature_meta), function(nm) {
    meta = feature_meta[[nm]]
    yo = obs_d[[nm]]
    ym = mod_d[[nm]]
    rng = range(c(yo, ym))
    pad = diff(rng)
    if (!is.finite(pad) || pad < 1e-8) pad = max(1, abs(rng[1]) * 0.05)
    else pad = pad * 0.1
    ylo = rng[1] - pad
    yhi = rng[2] + pad
    cw = 520; ch = 220; cx0 = 78; cx1 = 500; cy0 = 16; cy1 = 168
    xt = c(1L, 8L, 16L, length(hours_x))
    xlabels = paste(vapply(xt, function(i) {
      sprintf('<text x="%.1f" y="190" text-anchor="middle" font-size="11" fill="#57534e">%s</text>',
              svg_x(i, 1, length(hours_x), cx0, cx1), html_escape(hour_tick[i]))
    }, ""), collapse = "")
    yticks2 = pretty(c(ylo, yhi), n = 4)
    yticks2 = yticks2[yticks2 >= ylo & yticks2 <= yhi]
    ylabels = paste(vapply(yticks2, function(g) {
      sprintf('<text x="70" y="%.1f" text-anchor="end" font-size="11" fill="#57534e">%s</text>',
              svg_y(g, ylo, yhi, cy0, cy1) + 4, fmt_num(g, meta$digits))
    }, ""), collapse = "")
    grids = paste(vapply(yticks2, function(g) {
      yy = svg_y(g, ylo, yhi, cy0, cy1)
      sprintf('<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="#e7e5e4"/>', cx0, yy, cx1, yy)
    }, ""), collapse = "")
    sw_legend = legend_html(list(
      list(color = col_obs, label = "Measured near Earth"),
      list(color = col_mod, label = "AWSoM model")
    ))
    sprintf(paste0(
      '<figure><figcaption>%s <span>(%s)</span></figcaption>',
      '<p>%s</p>',
      '%s',
      '<svg viewBox="0 0 %d %d" role="img" aria-label="%s">',
      grids, ylabels, xlabels,
      '<polyline fill="none" stroke="%s" stroke-width="2.2" points="%s"/>',
      '<polyline fill="none" stroke="%s" stroke-width="2.2" points="%s"/>',
      '</svg></figure>'),
      html_escape(meta$label), html_escape(meta$unit),
      html_escape(meta$note),
      sw_legend,
      cw, ch, html_escape(meta$label),
      col_obs, svg_pts(hours_x, yo, 1, length(hours_x), ylo, yhi, cx0, cx1, cy0, cy1),
      col_mod, svg_pts(hours_x, ym, 1, length(hours_x), ylo, yhi, cx0, cx1, cy0, cy1))
  }, "")

  sw_rows = paste(vapply(seq_along(hours_x), function(i) {
    sprintf("<tr><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
            html_escape(hour_lab[i]), obs_d$n_native[i],
            fmt_num(obs_d$B[i], 2), fmt_num(obs_d$Bz[i], 2), fmt_num(obs_d$Vs[i], 1),
            fmt_num(obs_d$Pd[i], 2), fmt_num(obs_d$Pt[i], 0),
            mod_d$n_native[i],
            fmt_num(mod_d$B[i], 2), fmt_num(mod_d$Bz[i], 2), fmt_num(mod_d$Vs[i], 1),
            fmt_num(mod_d$Pd[i], 2), fmt_num(mod_d$Pt[i], 0))
  }, ""), collapse = "\n")

  f107_24_html = "<p>No F10.7 forecast falls in the next 24 hours.</p>"
  f107_27_html = "<p>No F10.7 forecast is available for this issue.</p>"
  has_f107 = !is.null(f107) && nrow(f107)
  if (has_f107) {
    show = f107
    if (!is.null(outlook) && nrow(outlook)) {
      ends = as.Date(substr(outlook$day_end, 1, 10))
      show = f107[f107$date %in% ends | f107$date %in% (ends - 1), ]
    }
    starts24 = parse_utc(f24o$target_bin_end) - kp_step_hr * 3600
    dates24 = unique(as.Date(format(starts24, "%Y-%m-%d", tz = "UTC")))
    flux_view = function(keep, hist_days, tick_by, lead_note) {
      part = show[show$date %in% keep, ]
      if (!nrow(part)) return(sprintf("<p>%s</p>", html_escape(lead_note)))
      flux_rows = paste(vapply(seq_len(nrow(part)), function(i) {
        got = !is.null(f107_measured) && nrow(f107_measured) && part$date[i] %in% f107_measured$date
        meas_row = if (got) f107_measured[match(part$date[i], f107_measured$date), ] else NULL
        sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
                html_escape(format(part$date[i])),
                if (got) fmt_num(meas_row$f107_obs, 1) else "",
                if (got) fmt_num(meas_row$f107_adj, 1) else "",
                fmt_num(part$f107_obs[i], 1), fmt_num(part$f107_adj[i], 1))
      }, ""), collapse = "\n")
      day_x = function(d) as.numeric(as.POSIXct(paste(format(d), "00:00:00"), tz = "UTC"))
      flux_x = day_x(part$date)
      meas_flux = f107_measured
      if (!is.null(meas_flux) && nrow(meas_flux)) {
        meas_flux = meas_flux[is.finite(meas_flux$f107_obs) & meas_flux$date < min(part$date), ]
        if (nrow(meas_flux)) {
          meas_flux = if (hist_days > 0) {
            meas_flux[meas_flux$date >= min(part$date) - hist_days, ]
          } else {
            meas_flux[which.max(meas_flux$date), ]
          }
        }
      }
      meas_x = if (!is.null(meas_flux) && nrow(meas_flux)) day_x(meas_flux$date) else numeric()
      y_obs = part$f107_obs
      y_adj = part$f107_adj
      x_line = flux_x
      if (length(meas_x)) {
        last = which.max(meas_flux$date)
        x_line = c(meas_x[last], flux_x)
        y_obs = c(meas_flux$f107_obs[last], part$f107_obs)
        y_adj = c(meas_flux$f107_adj[last], part$f107_adj)
      }
      fx_all = c(meas_x, flux_x)
      fxmin = min(fx_all) - 18 * 3600
      fxmax = max(fx_all) + 18 * 3600
      ft = seq(as.POSIXct(format(as.POSIXct(fxmin, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d"), tz = "UTC"),
               as.POSIXct(format(as.POSIXct(fxmax, origin = "1970-01-01", tz = "UTC"), "%Y-%m-%d"), tz = "UTC"),
               by = tick_by)
      ft = ft[as.numeric(ft) >= fxmin & as.numeric(ft) <= fxmax]
      if (!length(ft)) ft = as.POSIXct(format(part$date[1], "%Y-%m-%d"), tz = "UTC")
      fxtick = as.numeric(ft)
      fxlab = trimws(format(ft, "%b %e", tz = "UTC"))
      flux_vals = c(part$f107_obs, part$f107_adj)
      if (length(meas_x)) flux_vals = c(flux_vals, meas_flux$f107_obs, meas_flux$f107_adj)
      flux_rng = range(flux_vals, na.rm = TRUE)
      flux_pad = diff(flux_rng)
      if (!is.finite(flux_pad) || flux_pad < 1) flux_pad = 1
      flux_ticks = pretty(c(flux_rng[1] - 0.08 * flux_pad, flux_rng[2] + 0.08 * flux_pad), n = 4)
      flux_series = list()
      if (length(meas_x)) {
        flux_series = list(
          list(x = meas_x, y = meas_flux$f107_obs, color = "#44403c", width = "2.2"),
          list(x = meas_x, y = meas_flux$f107_adj, color = "#a8a29e", width = "1.8", dash = TRUE)
        )
      }
      flux_series = c(flux_series, list(
        list(x = x_line, y = y_obs, color = "#047857", width = "2.4"),
        list(x = x_line, y = y_adj, color = "#0369a1", width = "2.2")
      ))
      flux_svg = forecast_svg_chart(
        flux_series, fxmin, fxmax, min(flux_ticks), max(flux_ticks), flux_ticks,
        fxtick, fxlab, "sfu", "Date (UTC)", "F10.7 forecast")
      meas_note = "Penticton measures the flux at 20:00 UT, so the value for the issue date is still a forecast."
      meas_table = ""
      if (length(meas_x)) {
        last_i = which.max(meas_flux$date)
        meas_note = paste0(
          if (hist_days <= 0) "The gray point is only the last measurement, " else "Gray is the measured flux. The last measurement is ",
          format(meas_flux$date[last_i]), ": ",
          fmt_num(meas_flux$f107_obs[last_i], 1), " sfu, adjusted to 1 AU ",
          fmt_num(meas_flux$f107_adj[last_i], 1),
          ". Penticton measures at 20:00 UT, so the issue date is not missing: it has not been measured, and the green point on that date is the forecast.")
        recent_m = if (hist_days <= 0) meas_flux[0, ] else if (hist_days <= 7) meas_flux else meas_flux[meas_flux$date >= max(meas_flux$date) - 13, ]
        meas_table = if (!nrow(recent_m)) "" else paste0(
          '<div class="wrap"><table><thead><tr>',
          '<th>Date (UTC)</th><th>Measured F10.7</th><th>Measured, adjusted to 1 AU</th>',
          '</tr></thead><tbody>',
          paste(vapply(seq_len(nrow(recent_m)), function(i) {
            sprintf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>",
                    html_escape(format(recent_m$date[i])),
                    fmt_num(recent_m$f107_obs[i], 1), fmt_num(recent_m$f107_adj[i], 1))
          }, ""), collapse = "\n"),
          '</tbody></table></div>')
      }
      paste0(
        '<h2>F10.7</h2>',
        '<p>This is the observed 10.7 cm solar flux, in solar flux units. ', lead_note, ' ', meas_note, ' ',
        'Each day ahead has its own ridge regression, not the Kp XGBoost, fit on the measured flux from 1957 through the last day already on record. ',
        'Green is that forecast, drawn from the last measurement so the step into the forecast is not a gap. Blue is the same value adjusted to 1 AU. The pale line is the measured flux adjusted to 1 AU. Each point sits on its UT date.</p>',
        '<figure>',
        legend_html(c(
          if (length(meas_x)) list(
            list(color = "#44403c", label = "Measured F10.7"),
            list(color = "#a8a29e", label = "Measured, adjusted to 1 AU", dash = TRUE)
          ) else list(),
          list(
            list(color = "#047857", label = "Forecast F10.7"),
            list(color = "#0369a1", label = "Forecast, adjusted to 1 AU")
          )
        )),
        flux_svg,
        '</figure>',
        meas_table,
        '<div class="wrap"><table><thead><tr>',
        '<th>Date (UTC)</th><th>Measured F10.7</th><th>Measured, adjusted to 1 AU</th><th>Forecast F10.7</th><th>Forecast, adjusted to 1 AU</th>',
        '</tr></thead><tbody>', flux_rows, '</tbody></table></div>'
      )
    }
    f107_24_html = flux_view(
      dates24, 0, "1 day",
      "This chart is the forecast for the next 24 hours, not the past week. F10.7 is one value per day, so those hours cover the UT dates below.")
    f107_27_html = flux_view(
      show$date, 27, "7 days",
      "This view is the daily flux through the end of the 27-day Kp outlook.")
  }
  f107_past_html = "<p>No measured F10.7 is on record before this issue.</p>"
  if (has_f107 && !is.null(f107_measured) && nrow(f107_measured)) {
    issue_day = as.Date(format(issue_time, "%Y-%m-%d", tz = "UTC"))
    histf = f107_measured[is.finite(f107_measured$f107_obs) & f107_measured$date < issue_day, ]
    if (nrow(histf)) {
      histf = histf[order(histf$date), ]
      histf = histf[histf$date >= max(histf$date) - 13, ]
      f107_past_rows = paste(vapply(seq_len(nrow(histf)), function(i) {
        sprintf("<tr><td>%s</td><td>%s</td><td>%s</td></tr>",
                html_escape(format(histf$date[i])),
                fmt_num(histf$f107_obs[i], 1), fmt_num(histf$f107_adj[i], 1))
      }, ""), collapse = "\n")
      f107_past_html = paste0(
        '<h2>F10.7 already measured</h2>',
        '<p>These are the daily flux values on record when the forecast was issued. ',
        'The forecast chart starts from the last of them. ',
        'Penticton measures at 20:00 UT, so the issue date is not in this table.</p>',
        '<div class="wrap"><table><thead><tr>',
        '<th>Date (UTC)</th><th>Measured F10.7</th><th>Measured, adjusted to 1 AU</th>',
        '</tr></thead><tbody>', f107_past_rows, '</tbody></table></div>'
      )
    }
  }

  page = paste0(
    '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>Kp, ap, and F10.7 forecast</title><style>',
    page_css(paste0(
      '.charts{display:grid;grid-template-columns:1fr;gap:8px;}',
      '@media(min-width:800px){.charts{grid-template-columns:1fr 1fr;}}',
      'figure{margin:0 0 12px;background:#fff;border:1px solid #e7e5e4;border-radius:10px;padding:8px 10px 4px;min-width:0;}',
      'figcaption{font-size:14px;margin-bottom:2px;} figcaption span{color:#a8a29e;}',
      'figure p{margin:4px 0 6px;font-size:14px;}',
      'details{margin-top:12px;} summary{cursor:pointer;color:#57534e;}',
      '.metric-panel[hidden]{display:none;}',
      '.wind-block[hidden]{display:none;}',
      '.meta-issued{display:grid;grid-template-columns:1fr 1fr;gap:12px;}',
      '.meta-issued>div{min-width:0;}',
      '.meta-measured{margin-top:0;}'
    )),
    '</style></head><body><main>',
    '<nav><span class="here">Realtime</span><a href="VARX_Kp_baseline.html">Baseline</a></nav>',
    '<h1>', if (is.null(f107) || !nrow(f107)) "Kp and ap forecast" else "Kp, ap, and F10.7 forecast", '</h1>',
    '<p>Kp is a number from 0 to 9 that says how disturbed Earth&rsquo;s magnetic field is. ',
    'Higher means a stronger disturbance. As a rough guide, 0 to 2 is quiet, 3 is unsettled, 4 is active, 5 is a minor storm, and 6 or more is a stronger storm. ',
    'ap is the planetary amplitude of that Kp. ',
    'Times on this page are UTC.</p>',
    '<p>',
    if (is.null(f107) || !nrow(f107)) {
      'The charts are the issued Kp and the ap converted from that Kp. '
    } else {
      'The charts are the issued forecast for Kp, for ap, and for the daily 10.7 cm solar flux. '
    },
    'Choose Kp, ap, or F10.7, then next 24 hours or next 27 days. ',
    'The forecast was issued at the time shown below, using only information that had already arrived by then.</p>',
    issued_models_note(),
    if (is.null(kp_history_note) || !nzchar(kp_history_note)) "" else paste0(
      '<p>', html_escape(kp_history_note), '</p>'),
    '<div class="meta meta-issued">',
    '<div>Forecast issued (UTC)<b>', html_escape(fmt_utc(issue_time)), '</b></div>',
    '<div>Page made (UTC)<b>', html_escape(generated_at), '</b></div>',
    '</div>',
    '<div class="meta meta-measured">',
    '<div>Latest measured Kp<b>', fmt_num(kp_bin_value[length(kp_bin_value)], 3), '</b></div>',
    '<div>ap of that Kp<b>', fmt_num(ap_of(kp_bin_value[length(kp_bin_value)]), 2), '</b></div>',
    if (is.null(f107_measured) || !nrow(f107_measured)) "" else paste0(
      '<div>Latest measured F10.7<b>',
      fmt_num(f107_measured$f107_obs[which.max(f107_measured$date)], 1),
      '</b></div>'),
    '</div>',
    '<p class="downloads">Download this forecast: ',
    '<a href="issued_forecast_3h.csv" download>3-hour Kp, ap, and F10.7</a> ',
    '<a href="issued_forecast_27d.csv" download>next 27 days</a></p>',
    '<div class="togs metric-switch">',
    '<button type="button" class="tog metric-btn" data-metric="kp" aria-pressed="true">Kp</button>',
    '<button type="button" class="tog metric-btn" data-metric="ap" aria-pressed="false">ap</button>',
    if (has_f107) '<button type="button" class="tog metric-btn" data-metric="f107" aria-pressed="false">F10.7</button>' else "",
    '<button type="button" class="tog wind-btn" aria-pressed="false">Solar wind over the past day</button>',
    '</div>',
    '<div class="wind-block" hidden>',
    '<p>These charts are not a forecast. They are the solar wind during the 24 hours before the forecast was issued, and they are the solar-wind input to the 24-hour Kp chart. ',
    'Blue is what was measured near Earth. Orange is the AWSoM model, adjusted so its usual level matches the long-term solar-wind record (',
    html_escape(fmt_utc(ab_live$start)), ' through ', html_escape(fmt_utc(ab_live$end)),
    '). The hours being forecast are not included. ',
    'The measured record for this forecast looks back ', fmt_num(lag_obs_hr / 24, 0),
    ' days for a similar solar-wind pattern, and the model record looks back ',
    fmt_num(lag_mod_hr / 24, 0), ' days.</p>',
    '<div class="charts">', paste(charts, collapse = ""), '</div>',
    '<details><summary>Show the 24 hourly values</summary>',
    '<div class="wrap scroll"><table><thead><tr>',
    '<th>Hour start UTC</th><th>Measured samples</th><th>Measured |B|</th><th>Measured Bz</th><th>Measured Vs</th><th>Measured Pd</th><th>Measured Pt</th>',
    '<th>AWSoM samples</th><th>AWSoM |B|</th><th>AWSoM Bz</th><th>AWSoM Vs</th><th>AWSoM Pd</th><th>AWSoM Pt</th>',
    '</tr></thead><tbody>', sw_rows, '</tbody></table></div></details>',
    '</div>',
    range_switch(),
    '<p>One quantity and one range at a time.</p>',
    view_panel("kp", "24", paste0(
    '<h2>Next 24 hours</h2>',
    '<p>Left is 3 hours ahead. Right is 24 hours ahead. The height is Kp at the end of that step.</p>',
    if (has_issued) {
      paste0(
        '<p>Green is the issued forecast. Each of the eight steps has its own XGBoost model, trained on Kp and near-Earth solar wind through January 2026. ',
        'It uses the last 24 hours of Kp, the last 24 hours of solar wind measured near Earth, and the closest solar-wind pattern 25 to 29 days earlier. ',
        'The eight steps are forecast separately, and no predicted Kp is fed back in. ',
        'The buttons add the in-situ and AWSoM regressions. The dashed line repeats the latest Kp. It is not a forecast.</p>')
    } else {
      paste0(
        '<p>Green is Boosted, the forecast this page issues. It uses solar wind measured near Earth. ',
        'When the latest Kp is already high, and earlier storms at that lead were closer to AWSoM, it blends the two. Otherwise it stays on the in-situ forecast. ',
        'The buttons add In-situ and AWSoM. The dashed line repeats the latest Kp. It is not a forecast.</p>')
    },
    '<p>Only solar wind that has already reached Earth is used. The forecast does not see solar wind that has not arrived, so it will not rise ahead of a storm that is not yet in the measurements.</p>',
    '<figure>',
    legend_html(list(
      list(color = col_boost, label = issued_label),
      list(color = col_pers, label = "Latest Kp, held flat", dash = TRUE)
    )),
    toggle_btns(),
    forecast_svg,
    '</figure>',
    '<div class="wrap"><table><thead><tr>',
    '<th>Forecast time (UTC)</th><th>Hours ahead</th><th>', issued_label, ' Kp</th>',
    '<th>Measured Kp</th><th>In-situ</th><th>AWSoM</th><th>Latest Kp</th>',
    '</tr></thead><tbody>', forecast_rows_html, '</tbody></table></div>'
    )),
    view_panel("kp", "27", kp27_html, hidden = TRUE),
    view_panel("ap", "24", paste0(
      '<h2>Next 24 hours</h2>',
      '<p>ap is the Bartels amplitude of the Kp in the Kp view. It is not a separate model. ',
      'Left is 3 hours ahead. Right is 24 hours ahead. ',
      'Green is the issued ap. The dashed line holds the ap of the latest measured Kp. ',
      'The buttons add the ap of the in-situ and AWSoM Kp.</p>',
      '<figure>',
      legend_html(c(
        if (any(got24)) list(list(color = "#44403c", label = "Measured ap")) else list(),
        list(
          list(color = col_boost, label = paste(issued_label, "ap")),
          list(color = col_pers, label = "Latest ap, held flat", dash = TRUE)
        )
      )),
      toggle_btns(),
      ap24_svg,
      '</figure>',
      '<div class="wrap"><table><thead><tr>',
      '<th>Forecast time (UTC)</th><th>Hours ahead</th><th>', issued_label, ' ap</th>',
      '<th>Measured ap</th><th>In-situ ap</th><th>AWSoM ap</th><th>Latest ap</th>',
      '</tr></thead><tbody>', ap24_rows_html, '</tbody></table></div>'
    ), hidden = TRUE),
    view_panel("ap", "27", ap27_html, hidden = TRUE),
    if (has_f107) view_panel("f107", "24", f107_24_html, hidden = TRUE) else "",
    if (has_f107) view_panel("f107", "27", f107_27_html, hidden = TRUE) else "",
    view_panel("kp", NULL, paste0(
    '<h2>Kp over the past day</h2>',
    '<p>These eight values are the Kp that had already been recorded when the forecast was issued. Both the 24-hour forecast and the 27-day outlook start from this history.</p>',
    '<div class="wrap"><table><thead><tr><th>Period starts (UTC)</th><th>Period ends (UTC)</th><th>Kp</th></tr></thead><tbody>',
    kp_rows_html, '</tbody></table></div>'
    )),
    view_panel("ap", NULL, paste0(
      '<h2>ap over the past day</h2>',
      '<p>These are the Bartels amplitudes of the eight Kp values recorded before the forecast was issued. ',
      'ap uses that same history. It is not a separate measurement.</p>',
      '<div class="wrap"><table><thead><tr><th>Period starts (UTC)</th><th>Period ends (UTC)</th><th>ap</th></tr></thead><tbody>',
      ap_rows_html, '</tbody></table></div>'
    ), hidden = TRUE),
    if (has_f107) view_panel("f107", NULL, f107_past_html, hidden = TRUE) else "",
    chart_script(),
    '</main></body></html>'
  )
  write_html(path, page)
}

issued_quantity_baseline = function(xgb = NULL, roll = NULL, flux = NULL) {
  score_one = function(y, yhat, ref) {
    ok = is.finite(y) & is.finite(yhat) & is.finite(ref)
    n = sum(ok)
    if (!n) return(data.frame(n = 0L, RMSE = NA_real_, MAE = NA_real_, Correlation = NA_real_, skill = NA_real_))
    rm = sqrt(mean((y[ok] - yhat[ok])^2))
    rr = sqrt(mean((y[ok] - ref[ok])^2))
    data.frame(
      n = n, RMSE = rm, MAE = mean(abs(y[ok] - yhat[ok])),
      Correlation = if (n > 1L) suppressWarnings(cor(y[ok], yhat[ok])) else NA_real_,
      skill = if (!is.finite(rr) || rr == 0) NA_real_ else 1 - rm / rr
    )
  }
  skill_td = function(x) {
    if (!is.finite(x)) return("<td></td>")
    color = if (x > 0) "#087F23" else if (x < 0) "#b91c1c" else "#57534e"
    sprintf('<td style="color:%s">%+.3f</td>', color, x)
  }
  row_of = function(lead, input, sc, lead_unit = "h") {
    sprintf("<tr><td>%d%s</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td>%s</tr>",
            lead, lead_unit, html_escape(input), sc$n,
            fmt_num(sc$RMSE, 3), fmt_num(sc$MAE, 3),
            if (!is.finite(sc$Correlation)) "" else fmt_num(sc$Correlation, 3),
            skill_td(sc$skill))
  }
  head_row = paste0(
    '<div class="wrap"><table><thead><tr>',
    '<th>Lead</th><th>Input</th><th>n</th><th>RMSE</th><th>MAE</th><th>Correlation</th><th>Skill</th>',
    '</tr></thead><tbody>')
  pair_chart = function(x, y_true, y_hat, y_ref, ylab, aria) {
    ok = is.finite(x) & is.finite(y_true) & is.finite(y_hat)
    x = x[ok]; y_true = y_true[ok]; y_hat = y_hat[ok]; y_ref = y_ref[ok]
    if (length(x) < 2L) return("")
    recent = x >= max(x) - 14 * 86400
    if (sum(recent) < 2L) recent = rep(TRUE, length(x))
    x = x[recent]; y_true = y_true[recent]; y_hat = y_hat[recent]; y_ref = y_ref[recent]
    ymax = max(c(y_true, y_hat, y_ref), na.rm = TRUE)
    ymin = 0
    ticks = pretty(c(ymin, ymax), n = 4)
    ticks = ticks[ticks >= 0]
    ymax = max(ticks)
    t0 = as.POSIXct(min(x), origin = "1970-01-01", tz = "UTC")
    t1 = as.POSIXct(max(x), origin = "1970-01-01", tz = "UTC")
    xt = seq(as.POSIXct(format(t0, "%Y-%m-%d"), tz = "UTC"), t1, by = "2 days")
    xt = xt[as.numeric(xt) >= min(x) & as.numeric(xt) <= max(x)]
    paste0(
      '<figure>',
      legend_html(list(
        list(color = "#a8a29e", label = "Measured"),
        list(color = "#047857", label = "Issued"),
        list(color = "#6b7280", label = "Reference", dash = TRUE)
      )),
      forecast_svg_chart(
        list(
          list(x = x, y = y_true, color = "#a8a29e", width = "1.6"),
          list(x = x, y = y_hat, color = "#047857", width = "2.2"),
          list(x = x, y = y_ref, color = "#6b7280", width = "1.4", dash = TRUE)
        ),
        min(x), max(x), ymin, ymax, ticks,
        as.numeric(xt), trimws(format(xt, "%b %e", tz = "UTC")),
        ylab, "Target time (UTC)", aria),
      '</figure>')
  }
  kp_parts = character()
  ap_parts = character()
  f107_html = ""
  if (!is.null(xgb) && nrow(xgb)) {
    leads = seq(3L, 24L, by = 3L)
    kp_rows = paste(vapply(leads, function(h) {
      y = xgb[[sprintf("y_true_%02dh", h)]]
      yhat = xgb[[sprintf("pred_%02dh", h)]]
      ref = xgb$last_observed_Kp
      paste0(
        row_of(h, "Issued XGBoost", score_one(y, yhat, ref)),
        row_of(h, "Persistence", score_one(y, ref, ref)))
    }, ""), collapse = "")
    ap_rows = paste(vapply(leads, function(h) {
      y = kp_to_ap(xgb[[sprintf("y_true_%02dh", h)]])
      yhat = kp_to_ap(xgb[[sprintf("pred_%02dh", h)]])
      ref = kp_to_ap(xgb$last_observed_Kp)
      paste0(
        row_of(h, "Issued XGBoost", score_one(y, yhat, ref)),
        row_of(h, "Persistence", score_one(y, ref, ref)))
    }, ""), collapse = "")
    y24 = kp_to_ap(xgb$y_true_24h)
    p24 = kp_to_ap(xgb$pred_24h)
    r24 = kp_to_ap(xgb$last_observed_Kp)
    x24 = as.numeric(parse_utc(xgb$issue_time)) + 24 * 3600
    span = paste(html_escape(xgb$issue_time[1]), "through", html_escape(tail(xgb$issue_time, 1)))
    kp_parts = c(kp_parts, paste0(
      '<h2>Issued Kp, next 24 hours</h2>',
      '<p>These are the vanilla XGBoost models used for the first 24 hours of the issued forecast, scored where the target Kp had already been observed. ',
      'Skill is against repeating the Kp known at issue time. ', span, '.</p>',
      head_row, kp_rows, '</tbody></table></div>'
    ))
    ap_parts = c(ap_parts, paste0(
      '<h2>Issued ap, next 24 hours</h2>',
      '<p>ap is the Bartels amplitude of that Kp. The same issues are scored in ap. ',
      'Gray is the measured ap. Green is the issued forecast. The dashed line repeats the ap known at issue time. ',
      span, '.</p>',
      head_row, ap_rows, '</tbody></table></div>',
      pair_chart(x24, y24, p24, r24, "ap", "24 hour ap versus measured ap")
    ))
  }
  if (!is.null(roll) && nrow(roll)) {
    leads = c(27L, 72L, 168L, 336L, 648L)
    kp_rows = paste(vapply(leads, function(h) {
      sub = roll[roll$lead_end_hr == h, ]
      y = sub$y_true
      yhat = sub$y_insitu
      yaw = sub$y_awsom
      ref = sub$y_matched
      paste0(
        row_of(h, "Issued in-situ", score_one(y, yhat, ref)),
        row_of(h, "Issued AWSoM", score_one(y, yaw, ref)),
        row_of(h, "Best matched cycle", score_one(y, ref, ref)))
    }, ""), collapse = "")
    ap_rows = paste(vapply(leads, function(h) {
      sub = roll[roll$lead_end_hr == h, ]
      ap_keep = function(v) {
        out = rep(NA_real_, length(v))
        ok = is.finite(v)
        if (any(ok)) out[ok] = kp_to_ap(v[ok])
        out
      }
      y = ap_keep(sub$y_true)
      yhat = ap_keep(sub$y_insitu)
      yaw = ap_keep(sub$y_awsom)
      ref = ap_keep(sub$y_matched)
      paste0(
        row_of(h, "Issued in-situ", score_one(y, yhat, ref)),
        row_of(h, "Issued AWSoM", score_one(y, yaw, ref)),
        row_of(h, "Best matched cycle", score_one(y, ref, ref)))
    }, ""), collapse = "")
    sub = roll[roll$lead_end_hr == 648L, ]
    finite_issue = roll$issue_time[is.finite(roll$y_insitu) & is.finite(roll$y_true)]
    roll_span = if (length(finite_issue)) {
      paste(html_escape(min(finite_issue)), "through", html_escape(max(finite_issue)))
    } else ""
    kp_parts = c(kp_parts, paste0(
      '<h2>Issued Kp, after 24 hours</h2>',
      '<p>From 27 hours through 27 days the issued forecast is the in-situ regression stepped forward on the closest OMNI cycle. ',
      'AWSoM is the same outlook with the AWSoM solar-wind block. It is not the issued forecast. ',
      'The reference for both is the Kp that followed that cycle. Skill = 1 - RMSE / RMSE of that cycle. ',
      'Issues with both a forecast and a later Kp run ', roll_span, '.</p>',
      head_row, kp_rows, '</tbody></table></div>'
    ))
    ap_parts = c(ap_parts, paste0(
      '<h2>Issued ap, after 24 hours</h2>',
      '<p>The same in-situ Kp, converted to ap, for issues ', roll_span, '. ',
      'AWSoM is that outlook converted to ap. It is not the issued forecast. ',
      'The reference is the ap of the Kp that followed the closest cycle. ',
      'The chart is the bin that ends 27 days after issue. Gray is the ap of the Kp observed later.</p>',
      head_row, ap_rows, '</tbody></table></div>',
      if (!nrow(sub)) "" else pair_chart(
        as.numeric(parse_utc(sub$target_end)),
        kp_to_ap(sub$y_true), kp_to_ap(sub$y_insitu), kp_to_ap(sub$y_matched),
        "ap", "27 day ap versus measured ap")
    ))
  }
  if (!is.null(flux) && nrow(flux)) {
    show_leads = c(1L, 3L, 7L, 14L, 27L)
    use = flux$model %in% c("persistence", "ridge") & flux$lead_day %in% show_leads
    shown = flux[use, ]
    shown = shown[order(shown$split, shown$lead_day, match(shown$model, c("ridge", "persistence"))), ]
    flux_rows = paste(vapply(seq_len(nrow(shown)), function(i) {
      sprintf("<tr><td>%s</td><td>%d day</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td></tr>",
              html_escape(shown$split[i]), shown$lead_day[i], html_escape(shown$model[i]),
              shown$n[i], fmt_num(shown$RMSE[i], 2),
              if (shown$model[i] == "persistence") "+0.000" else sprintf("%+.3f", shown$skill_vs_persistence[i]))
    }, ""), collapse = "\n")
    mean_rows = paste(vapply(unique(flux$split), function(sp) {
      sub = flux[flux$split == sp & flux$lead_day %in% 1:27, ]
      paste(vapply(c("ridge", "persistence"), function(model) {
        s = sub[sub$model == model, ]
        sprintf("<tr><td>%s</td><td>1–27 day mean</td><td>%s</td><td>%d</td><td>%s</td><td></td></tr>",
                html_escape(sp), model, s$n[1], fmt_num(mean(s$RMSE), 2))
      }, ""), collapse = "")
    }, ""), collapse = "")
    window = if (nrow(shown)) shown$split[1] else ""
    f107_html = paste0(
      '<p>The table scores the ridge from the note above, not the Kp XGBoost. ',
      'That ridge was chosen on 2019&ndash;2022, where it beat repeating the latest flux, repeating the flux from one rotation earlier, and a blend of those two. ',
      'The rows are that same ridge on every target day from ', html_escape(window),
      '. Coefficients are fit only on issues that ended before this window. ',
      'Skill is against repeating the latest measured flux. The mean is the average of the daily RMSE values from lead 1 through lead 27.</p>',
      '<div class="wrap"><table><thead><tr>',
      '<th>Period</th><th>Lead</th><th>Model</th><th>n</th><th>RMSE (sfu)</th><th>Skill</th>',
      '</tr></thead><tbody>', flux_rows, mean_rows, '</tbody></table></div>'
    )
  }
  list(
    kp = paste(kp_parts, collapse = ""),
    ap = paste(ap_parts, collapse = ""),
    f107 = f107_html
  )
}

write_baseline_page = function(path, base, generated_at, outlook_eval = NULL, issued_section = "") {
  col_obs = "#1d4ed8"
  col_mod = "#c2410c"
  col_boost = "#047857"
  col_pers = "#6b7280"
  col_true = "#a8a29e"
  metrics = base$metrics
  tau = base$tau
  y_true = base$y_true
  pers = base$pers
  order_in = c("Persistence", lab_meas, lab_awso, lab_boost)
  input_col = c(col_pers, col_obs, col_mod, col_boost)
  names(input_col) = order_in

  fmt_cor = function(x) if (!is.finite(x)) "" else fmt_num(x, 3)
  skill_td = function(x) {
    if (!is.finite(x)) return("<td></td>")
    color = if (x > 0) "#087F23" else if (x < 0) "#b91c1c" else "#57534e"
    sprintf('<td style="color:%s">%+.3f</td>', color, x)
  }
  metric_rows = function(frame) {
    paste(vapply(seq_len(nrow(frame)), function(i) {
      sprintf("<tr><td>%dh</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td>%s</tr>",
              frame$lead_end_hr[i], html_escape(frame$input[i]), frame$n[i],
              fmt_num(frame$RMSE[i], 3), fmt_num(frame$MAE[i], 3), fmt_cor(frame$Correlation[i]),
              skill_td(frame$skill_RMSE[i]))
    }, ""), collapse = "\n")
  }
  endpoint = metrics[metrics$lead_end_hr == metrics$horizon_hr, ]
  endpoint = endpoint[order(match(endpoint$input, order_in)), ]

  rmse_at = function(lead, input) {
    hit = metrics[metrics$horizon_hr == 24 & metrics$lead_end_hr == lead & metrics$input == input, ]
    if (!nrow(hit)) NA_real_ else hit$RMSE[1]
  }
  bar_leads = c(6L, 24L)
  bar_groups = c("6 hours", "24 hours")
  bar_spec = do.call(rbind, lapply(seq_along(bar_leads), function(g) {
    data.frame(
      group = bar_groups[g], input = order_in,
      value = vapply(order_in, function(inp) rmse_at(bar_leads[g], inp), 0),
      stringsAsFactors = FALSE)
  }))
  W = 760; H = 280; x0 = 52; x1 = 730; y0 = 28; y1 = 196
  ymax = max(bar_spec$value) * 1.28
  yticks = pretty(c(0, max(bar_spec$value)), n = 4)
  yticks = yticks[yticks >= 0 & yticks <= ymax]
  grid = paste(vapply(yticks, function(g) {
    yy = svg_y(g, 0, ymax, y0, y1)
    sprintf('<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="#e7e5e4"/>', x0, yy, x1, yy)
  }, ""), collapse = "")
  ylab = paste(vapply(yticks, function(g) {
    sprintf('<text x="44" y="%.1f" text-anchor="end" font-size="12" fill="#57534e">%s</text>',
            svg_y(g, 0, ymax, y0, y1) + 4, fmt_num(g, 2))
  }, ""), collapse = "")
  bar_w = 34
  gap = 8
  group_gap = 36
  cluster = length(order_in) * bar_w + (length(order_in) - 1) * gap
  total_w = length(bar_leads) * cluster + group_gap
  x_cursor = x0 + ((x1 - x0) - total_w) / 2
  bars = character(0)
  for (g in seq_along(bar_leads)) {
    group_x0 = x_cursor
    rows = bar_spec[bar_spec$group == bar_groups[g], ]
    for (i in seq_len(nrow(rows))) {
      bh = (rows$value[i] / ymax) * (y1 - y0)
      yy = y1 - bh
      bars = c(bars, sprintf(paste0(
        '<rect x="%.1f" y="%.1f" width="%g" height="%.1f" fill="%s" rx="3"/>',
        '<text x="%.1f" y="%.1f" text-anchor="middle" font-size="11" fill="#1c1917">%s</text>'),
        x_cursor, yy, bar_w, bh, input_col[[rows$input[i]]],
        x_cursor + bar_w / 2, max(y0 + 12, yy - 6), fmt_num(rows$value[i], 3)))
      x_cursor = x_cursor + bar_w + gap
    }
    x_cursor = x_cursor - gap
    bars = c(bars, sprintf(
      '<text x="%.1f" y="220" text-anchor="middle" font-size="13" fill="#1c1917">%s</text>',
      (group_x0 + x_cursor) / 2, bar_groups[g]))
    x_cursor = x_cursor + group_gap
  }
  bar_svg = sprintf(paste0(
    '<figure><figcaption>RMSE at 6 hours and 24 hours</figcaption>',
    '<p>Each group is one lead. All four inputs are shown.</p>',
    legend_html(list(
      list(color = col_pers, label = "Persistence"),
      list(color = col_obs, label = lab_meas),
      list(color = col_mod, label = lab_awso),
      list(color = col_boost, label = lab_boost)
    )),
    '<svg viewBox="0 0 %d %d" role="img" aria-label="RMSE at 6 and 24 hours">',
    grid, ylab, paste(bars, collapse = ""),
    '<text x="%g" y="252" font-size="12" fill="#57534e">RMSE of that 3-hour Kp bin</text>',
    '</svg></figure>'), W, H, x0)

  series_svg = function(x, y_true_v, y_obs_v, y_mod_v, y_boost_v, tick_by) {
    xmin = min(x); xmax = max(x)
    if (xmin == xmax) xmax = xmin + 1
    ymax_s = max(9, ceiling(max(c(y_true_v, y_obs_v, y_mod_v, y_boost_v), na.rm = TRUE)))
    ymin_s = 0
    yticks_s = if (ymax_s <= 9) c(0, 3, 6, 9) else pretty(c(0, ymax_s), n = 4)
    sw = 760; sh = 270; sx0 = 44; sx1 = 730; sy0 = 16; sy1 = 214
    grid_s = paste(vapply(yticks_s, function(g) {
      yy = svg_y(g, ymin_s, ymax_s, sy0, sy1)
      sprintf('<line x1="%g" y1="%.1f" x2="%g" y2="%.1f" stroke="#e7e5e4"/>', sx0, yy, sx1, yy)
    }, ""), collapse = "")
    ylab_s = paste(vapply(yticks_s, function(g) {
      sprintf('<text x="36" y="%.1f" text-anchor="end" font-size="12" fill="#57534e">%g</text>',
              svg_y(g, ymin_s, ymax_s, sy0, sy1) + 4, g)
    }, ""), collapse = "")
    t_min = as.POSIXct(xmin, origin = "1970-01-01", tz = "UTC")
    t_max = as.POSIXct(xmax, origin = "1970-01-01", tz = "UTC")
    if (tick_by == "month") {
      ticks = seq(as.POSIXct(format(t_min, "%Y-%m-01"), tz = "UTC"), t_max, by = "month")
      lab = format(ticks, "%b", tz = "UTC")
    } else {
      ticks = seq(as.POSIXct(format(t_min, "%Y-%m-%d"), tz = "UTC"), t_max, by = "2 days")
      lab = format(ticks, "%m-%d", tz = "UTC")
    }
    keep_tick = ticks >= t_min & ticks <= t_max
    ticks = ticks[keep_tick]
    lab = lab[keep_tick]
    xlab = if (!length(ticks)) "" else paste(vapply(seq_along(ticks), function(i) {
      sprintf('<text x="%.1f" y="236" text-anchor="middle" font-size="11" fill="#57534e">%s</text>',
              svg_x(as.numeric(ticks[i]), xmin, xmax, sx0, sx1), lab[i])
    }, ""), collapse = "")
    sprintf(paste0(
      '<svg viewBox="0 0 %d %d" role="img">',
      grid_s, ylab_s, xlab,
      '%s%s%s%s',
      '</svg>'),
      sw, sh,
      svg_runs(x, y_true_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_true, "1.2"),
      svg_runs(x, y_boost_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_boost, "2.2"),
      svg_runs(x, y_obs_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_obs, "1.8", cls = "opt insitu", opacity = 0.55),
      svg_runs(x, y_mod_v, xmin, xmax, ymin_s, ymax_s, sx0, sx1, sy0, sy1, col_mod, "1.8", cls = "opt awsom", opacity = 0.55))
  }

  x_all = as.numeric(tau)
  recent = tau >= (max(tau) - 14 * 86400)
  series_legend = paste0(legend_html(list(
    list(color = col_true, label = "Observed Kp"),
    list(color = col_boost, label = "Boosted")
  )), toggle_btns())
  chart_block = function(caption, note, svg) {
    paste0("<figure><figcaption>", html_escape(caption),
           "</figcaption><p>", html_escape(note), "</p>", series_legend, svg, "</figure>")
  }
  charts = paste0(
    chart_block("6 hours ahead, whole test period",
                "Gray is Kp observed later. Green is Boosted. The buttons add In-situ and AWSoM.",
                series_svg(x_all, y_true[, 2], base$y24o[, 2], base$y24m[, 2], base$yboost[, 2], "month")),
    chart_block("6 hours ahead, last 14 days",
                "The same 6-hour chart, kept to the latest 14 days so the rises and falls are easier to see.",
                series_svg(as.numeric(tau[recent]), y_true[recent, 2], base$y24o[recent, 2], base$y24m[recent, 2], base$yboost[recent, 2], "day")),
    chart_block("24 hours ahead, whole test period",
                "Each point is a forecast issued in the past, for the 3-hour period that ends a day later. The forecast usually peaks about 18 hours after the observed Kp, because it only knows conditions at the moment it is issued, not the solar wind that arrives the next day.",
                series_svg(x_all, y_true[, 8], base$y24o[, 8], base$y24m[, 8], base$yboost[, 8], "month")),
    chart_block("24 hours ahead, last 14 days",
                "The same 24-hour chart, kept to the latest 14 days.",
                series_svg(as.numeric(tau[recent]), y_true[recent, 8], base$y24o[recent, 8], base$y24m[recent, 8], base$yboost[recent, 8], "day"))
  )

  ord = order(tau, decreasing = TRUE)
  take = ord[recent[ord]]
  fmt_cell = function(x) if (!is.finite(x)) "" else fmt_num(x, 3)
  recent_rows = paste(vapply(take, function(i) {
    sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
            html_escape(fmt_utc(tau[i])),
            html_escape(fmt_utc(tau[i] + 24 * 3600)),
            fmt_num(y_true[i, 8], 3), fmt_cell(base$y24o[i, 8]),
            fmt_cell(base$y24m[i, 8]), fmt_cell(base$yboost[i, 8]),
            fmt_num(pers[i], 3))
  }, ""), collapse = "\n")

  by_lead = metrics[order(metrics$lead_end_hr, match(metrics$input, order_in)), ]
  n_iss = length(tau)
  page = paste0(
    '<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1">',
    '<title>Kp, ap, and F10.7 baseline</title><style>',
    page_css(paste0(
      'th:nth-child(2),td:nth-child(2){text-align:left;}',
      'figure{margin:0 0 14px;background:#fff;border:1px solid #e7e5e4;border-radius:10px;padding:8px 10px 4px;min-width:0;}',
      'figcaption{font-size:15px;color:#1c1917;font-weight:600;margin-bottom:2px;}',
      'figure p{margin:0 0 4px;font-size:13px;}',
      'h3{font-size:16px;margin:18px 0 6px;}',
      '.metric-panel[hidden]{display:none;}',
      '.metric-switch{margin:0 0 8px;}'
    )),
    '</style></head><body><main>',
    '<nav><a href="VARX_Kp_realtime.html">Realtime</a><span class="here">Baseline</span></nav>',
    metric_switch(),
    '<div class="metric-panel" data-metric="kp">',
    '<h1>Kp baseline</h1>',
    '<p>This is the Kp model the forecast page issues, scored against Kp that was observed later. ',
    'The first 24 hours are the vanilla XGBoost. From 27 hours through 27 days it is the in-situ regression stepped forward on the closest OMNI cycle. ',
    'The XGBoost record starts on 2026-02-25, after its training cutoff of 2026-01-25. ',
    'The longer leads use every 2026 issue whose target bin has already ended.</p>',
    if (is.list(issued_section)) issued_section$kp else issued_section,
    '</div>',
    '<div class="metric-panel" data-metric="ap" hidden>',
    '<h1>ap baseline</h1>',
    '<p>ap is the Bartels amplitude of the issued Kp. The issues and the dates are the same as on the Kp button.</p>',
    if (is.list(issued_section)) issued_section$ap else "",
    '</div>',
    '<div class="metric-panel" data-metric="f107" hidden>',
    '<h1>F10.7 baseline</h1>',
    if (is.list(issued_section)) issued_section$f107 else "",
    '</div>',
    chart_script(),
    '</main></body></html>'
  )
  write_html(path, page)
}

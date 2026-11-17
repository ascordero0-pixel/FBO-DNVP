# =============================================================================
# deviation_weekend.R
# -----------------------------------------------------------------------------
# Temperature-deviation (residual) figures for the weekend free-running windows
# (Friday 22:00 -> Monday 07:00, ~57 h). For each window and for all windows
# overlaid, it plots how v3 (single zone) and v6 (two zone) deviate from the
# MONITORED indoor temperature, plus the v6-v3 difference.
#
# Residuals (sign: negative = model colder than measured):
#     r_v3(t) = T_v3(t)  - T_meas(t)
#     r_v6(t) = T_v6(t)  - T_meas(t)
#     Delta(t)= T_v6(t)  - T_v3(t)   = r_v6 - r_v3   (pure model difference)
#
# Smoothing = centred moving average of window w points (w=12 -> 1 h at 5 min),
# a rectangular FIR low-pass that suppresses sub-hourly sensor noise while
# preserving the slow bias trend. Raw series drawn faint underneath.
#
# DATA SOURCE: D_win is built by load_calibrated_series.R from the CALIBRATED
# Track A runs (runs/v3_trackA, runs/v6_trackA), anchored to their validated
# metrics (v3 -5.41% / v6 +6.05% nocturnal). It is NOT the uncalibrated
# "Option A" series. Columns: k (POSIXct 5-min), T_out, v3, v6, T_meas, win.
# Output: PNG + PDF, one figure per window + one overlay, in OUT_DIR/deviation_figs.
# =============================================================================

suppressWarnings({library(lubridate)})
# Ensure D_win comes from the CALIBRATED Track A series (not the Option A run).
if(!exists(".CALIB_SERIES_LOADED") || !isTRUE(.CALIB_SERIES_LOADED))
  source("load_calibrated_series.R")
stopifnot(exists("D_win"))
if(nrow(D_win) == 0) stop("D_win is empty - check load_calibrated_series.R.")

OUTFIG <- if(exists("OUT_DIR")) file.path(OUT_DIR, "deviation_figs") else
                                 file.path(getwd(), "deviation_figs")
dir.create(OUTFIG, recursive = TRUE, showWarnings = FALSE)

SMOOTH_W <- 12                                   # 12 x 5 min = 1 h moving-average
cma <- function(x, w = SMOOTH_W)                 # centred moving average (NA at edges)
  as.numeric(stats::filter(x, rep(1/w, w), sides = 2))

# palette (colour-blind-safe blues/reds)
col_meas <- "black"; col_v3 <- "#1f6feb"; col_v6 <- "#d1242f"; col_out <- "grey55"
col_d    <- "#6f42c1"
faint <- function(col, a = 0.22) adjustcolor(col, alpha.f = a)

# plotmath axis labels (degree sign renders in BOTH png and pdf devices)
lab_Tin  <- expression("Indoor T ("*degree*"C)")
lab_Tout <- expression("Outdoor T ("*degree*"C)")
lab_res  <- expression("Residual  (sim "-" meas)  ("*degree*"C)")
lab_del  <- expression(Delta*" = v6 "-" v3  ("*degree*"C)")

wins <- sort(unique(D_win$win))

# helper: save the same drawing to PNG and PDF ------------------------------
save_both <- function(name, draw_fn, w = 9, h = 8.5){
  png(file.path(OUTFIG, paste0(name, ".png")), width = w, height = h,
      units = "in", res = 150); draw_fn(); dev.off()
  pdf(file.path(OUTFIG, paste0(name, ".pdf")), width = w, height = h,
      encoding = "ISOLatin1", useDingbats = FALSE); draw_fn(); dev.off()
}

# time axis: a tick every 12 h, labelled 'Fri 22h' ...
time_axis <- function(k)
  axis.POSIXct(1, at = seq(min(k), max(k), by = "12 hours"), format = "%a %Hh")

# ---------------------------------------------------------------------------
# 1. PER-WINDOW 3-panel figure
#    A: absolute temperatures (+ outdoor on right axis)
#    B: residuals r_v3, r_v6 (raw faint + 1-h moving average bold)
#    C: model difference Delta = v6 - v3 (raw faint + MA bold)
# ---------------------------------------------------------------------------
draw_window <- function(g){
  g   <- g[order(g$k), ]
  rv3 <- g$v3 - g$T_meas; rv6 <- g$v6 - g$T_meas; dd <- g$v6 - g$v3
  mbar <- mean(g$T_meas)
  nmbe3 <- 100 * mean(rv3) / mbar; nmbe6 <- 100 * mean(rv6) / mbar
  d0 <- format(min(g$k), "%a %d %b %Y")

  layout(matrix(1:3, ncol = 1), heights = c(2.2, 1.7, 1.3))

  # -- Panel A --
  par(mar = c(2.4, 4.6, 2.6, 4.6))
  yl <- range(c(g$T_meas, g$v3, g$v6), na.rm = TRUE)
  plot(g$k, g$T_meas, type = "l", lwd = 2.2, col = col_meas, ylim = yl,
       xlab = "", ylab = lab_Tin, xaxt = "n",
       main = sprintf("Free-running weekend  -  %s     (NMBE: v3 %+.1f%%,  v6 %+.1f%%)",
                      d0, nmbe3, nmbe6))
  lines(g$k, g$v3, lwd = 2, col = col_v3)
  lines(g$k, g$v6, lwd = 2, col = col_v6)
  time_axis(g$k)
  par(new = TRUE)
  plot(g$k, g$T_out, type = "l", lty = 2, col = col_out, axes = FALSE,
       xlab = "", ylab = "")
  axis(4, col = col_out, col.axis = col_out)
  mtext(lab_Tout, side = 4, line = 2.8, col = col_out)
  legend("topright", bty = "n", lwd = c(2.2, 2, 2, 1.4), lty = c(1, 1, 1, 2),
         col = c(col_meas, col_v3, col_v6, col_out),
         legend = c("Measured", "v3 (1-zone)", "v6 (2-zone)", "Outdoor"))

  # -- Panel B: residuals --
  par(mar = c(2.4, 4.6, 1.0, 4.6))
  yb <- range(c(rv3, rv6), na.rm = TRUE)
  plot(g$k, rv3, type = "l", col = faint(col_v3), lwd = 1, ylim = yb,
       xlab = "", ylab = lab_res, xaxt = "n")
  lines(g$k, rv6, col = faint(col_v6), lwd = 1)
  lines(g$k, cma(rv3), col = col_v3, lwd = 2.6)
  lines(g$k, cma(rv6), col = col_v6, lwd = 2.6)
  abline(h = 0, col = "grey30", lty = 3)
  time_axis(g$k)
  legend("bottomleft", bty = "n", ncol = 2, lwd = c(1, 2.6, 1, 2.6),
         col = c(faint(col_v3), col_v3, faint(col_v6), col_v6),
         legend = c("v3 raw", "v3 1-h MA", "v6 raw", "v6 1-h MA"))

  # -- Panel C: model difference --
  par(mar = c(3.6, 4.6, 1.0, 4.6))
  plot(g$k, dd, type = "l", col = faint(col_d), lwd = 1, xlab = "", ylab = lab_del)
  lines(g$k, cma(dd), col = col_d, lwd = 2.6)
  abline(h = 0, col = "grey30", lty = 3)
  abline(h = mean(dd), col = col_d, lty = 2)
  time_axis(g$k)
  legend("topright", bty = "n", lty = c(1, 2), lwd = c(2.6, 1),
         col = c(col_d, col_d),
         legend = c(expression(Delta~"(1-h MA)"),
                    as.expression(bquote("mean "*Delta*" = "*.(sprintf("%+.2f", mean(dd)))*degree*"C"))))
}

for(wv in wins){
  g  <- D_win[D_win$win == wv, ]
  d0 <- format(min(g$k), "%Y%m%d")
  save_both(paste0("deviation_weekend_", d0), function() draw_window(g))
}

# ---------------------------------------------------------------------------
# 2. OVERLAY: residuals of ALL windows on a common normalised time axis
#    (hours since Friday 22:00), each window a 1-h-MA line coloured by its mean
#    outdoor temperature. This visualises the key finding: colder windows carry
#    a larger (more negative) residual -> the bias ~ T_out dependence.
# ---------------------------------------------------------------------------
win_tout <- sapply(wins, function(wv) mean(D_win$T_out[D_win$win == wv]))
trng <- range(win_tout)
pal  <- colorRampPalette(c("#2166ac", "grey65", "#b2182b"))(101)   # cold -> warm
col_for <- function(tt)
  pal[pmax(1, pmin(101, round((tt - trng[1]) / (diff(trng) + 1e-9) * 100) + 1))]

draw_overlay <- function(){
  layout(matrix(1:2, ncol = 1), heights = c(1, 1))
  yr <- range(c(D_win$v3 - D_win$T_meas, D_win$v6 - D_win$T_meas), na.rm = TRUE)
  for(mdl in c("v3", "v6")){
    par(mar = c(4, 4.6, 2.2, 2))
    plot(NA, xlim = c(0, 57), ylim = yr,
         xlab = "Hours since Friday 22:00", ylab = lab_res,
         main = sprintf("Weekend residuals - %s  (each line = one window, coloured by mean outdoor T)", mdl))
    abline(h = 0, col = "grey30", lty = 3)
    for(wv in wins){
      g <- D_win[D_win$win == wv, ]; g <- g[order(g$k), ]
      hrs <- as.numeric(g$k - min(g$k), units = "hours")
      r <- if(mdl == "v3") g$v3 - g$T_meas else g$v6 - g$T_meas
      lines(hrs, cma(r), col = col_for(mean(g$T_out)), lwd = 2.3)
    }
    if(mdl == "v3"){
      labs <- sapply(wins, function(wv){ g <- D_win[D_win$win == wv, ]
        sprintf("%s  (%.1f)", format(min(g$k), "%d-%b"), mean(g$T_out)) })
      legend("bottomleft", bty = "n", lwd = 2.3, col = col_for(win_tout),
             legend = labs, cex = 0.8,
             title = expression(bold("Window (mean "*T[out]*", "*degree*"C)")))
    }
  }
}
save_both("deviation_weekend_OVERLAY", draw_overlay, w = 9, h = 8)

cat(sprintf("Saved %d per-window figures + 1 overlay (PNG + PDF) to:\n  %s\n",
            length(wins), OUTFIG))
print(list.files(OUTFIG, full.names = TRUE))

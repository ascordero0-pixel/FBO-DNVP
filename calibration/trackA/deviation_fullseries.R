# =============================================================================
# deviation_fullseries.R
# -----------------------------------------------------------------------------
# Deviation (residual) of the v3 (single-zone) and v6 (two-zone) models with
# respect to the MEASURED indoor temperature, for EVERY 5-min step of the whole
# RunPeriod (1 Nov -> 20 Dec 2021, 50 days).
#
# Residuals (sign: negative = model colder than measured):
#     r_v3(t) = T_v3(t) - T_meas(t)
#     r_v6(t) = T_v6(t) - T_meas(t)
#
# Global fit metrics annotated on the figure (mbar = mean measured indoor T):
#     NMBE    = 100/mbar * mean(r)                       [%]
#     CV-RMSE = 100/mbar * sqrt(mean(r^2))               [%]
#
# The bold overlay is a CENTRED 24-h moving average (w = 288 points): by
# averaging over full days it removes the diurnal cycle (weekday occupancy /
# heating) and exposes the day-to-day trend of the bias. Weekends (Fri 22:00 ->
# Mon 07:00, free-running) are lightly shaded so the systematic weekend deficit
# is visible against the calibrated weekday nights.
#
# DATA SOURCE: df (k, T_out, v3, v6) and meas (k, T_meas, ...) are built by
# load_calibrated_series.R from the CALIBRATED Track A runs, anchored to their
# validated metrics. This script sources that loader automatically if needed.
# Output: PNG + PDF in OUT_DIR/deviation_figs.
# =============================================================================

suppressWarnings({library(lubridate)})
# Ensure df/meas come from the CALIBRATED Track A series (not the Option A run).
if(!exists(".CALIB_SERIES_LOADED") || !isTRUE(.CALIB_SERIES_LOADED))
  source("load_calibrated_series.R")
stopifnot(exists("df"), exists("meas"))

# --- build the full merged 5-min series over the RunPeriod -------------------
DF <- merge(df, meas[, c("k","T_meas")], by = "k")
DF <- DF[order(DF$k), ]
t0 <- as.POSIXct("2021-11-01 00:00:00", tz = "UTC")
t1 <- as.POSIXct("2021-12-21 00:00:00", tz = "UTC")   # through 20 Dec inclusive
DF <- DF[DF$k >= t0 & DF$k < t1, ]
stopifnot(nrow(DF) > 0)

# --- residuals and global metrics -------------------------------------------
DF$r_v3 <- DF$v3 - DF$T_meas
DF$r_v6 <- DF$v6 - DF$T_meas
mbar    <- mean(DF$T_meas)
nmbe    <- function(r) 100 * mean(r) / mbar
cvrmse  <- function(r) 100 * sqrt(mean(r^2)) / mbar
NMBE_v3 <- nmbe(DF$r_v3);  CV_v3 <- cvrmse(DF$r_v3)
NMBE_v6 <- nmbe(DF$r_v6);  CV_v6 <- cvrmse(DF$r_v6)
cat(sprintf("Full series 1-Nov..20-Dec  N=%d  mean measured=%.2f C\n", nrow(DF), mbar))
cat(sprintf("  v3:  NMBE=%+.2f %%   CV-RMSE=%.2f %%\n", NMBE_v3, CV_v3))
cat(sprintf("  v6:  NMBE=%+.2f %%   CV-RMSE=%.2f %%\n", NMBE_v6, CV_v6))

# --- centred 24-h moving average --------------------------------------------
SMOOTH_W <- 288                                    # 288 x 5 min = 24 h
cma <- function(x, w = SMOOTH_W)
  as.numeric(stats::filter(x, rep(1/w, w), sides = 2))

# --- weekend spans (Fri 22:00 -> Mon 07:00) for background shading ----------
wd <- wday(DF$k); hh <- hour(DF$k) + minute(DF$k)/60
inwk <- (wd==6 & hh>=22) | (wd==7) | (wd==1) | (wd==2 & hh<7)
rr <- rle(inwk); ce <- cumsum(rr$lengths); cs <- c(1, head(ce,-1)+1)
wk_spans <- lapply(which(rr$values), function(i) c(DF$k[cs[i]], DF$k[ce[i]]))

# --- palette / labels -------------------------------------------------------
col_meas <- "black"; col_v3 <- "#1f6feb"; col_v6 <- "#d1242f"; col_out <- "grey55"
faint <- function(col, a = 0.20) adjustcolor(col, alpha.f = a)
lab_Tin  <- expression("Indoor T ("*degree*"C)")
lab_res  <- expression("Deviation  (model "-" measured)  ("*degree*"C)")

OUTFIG <- if(exists("OUT_DIR")) file.path(OUT_DIR, "deviation_figs") else
                                 file.path(getwd(), "deviation_figs")
dir.create(OUTFIG, recursive = TRUE, showWarnings = FALSE)

tick_axis <- function()
  axis.POSIXct(1, at = seq(t0, t1, by = "5 days"), format = "%d-%b")

shade_weekends <- function(yb, yt)
  for(s in wk_spans) rect(s[1], yb, s[2], yt,
                          col = faint("grey40", 0.12), border = NA)

# --- the figure: (A) absolute temperatures, (B) residuals -------------------
draw <- function(){
  layout(matrix(1:2, ncol = 1), heights = c(1.5, 2.0))

  # (A) absolute temperatures
  par(mar = c(2.4, 4.8, 2.6, 1.2))
  yl <- range(c(DF$T_meas, DF$v3, DF$v6), na.rm = TRUE)
  plot(DF$k, DF$T_meas, type = "n", ylim = yl, xaxt = "n", xlab = "", ylab = lab_Tin,
       main = "Indoor temperature: measured vs models (v3, v6) - full RunPeriod 1 Nov - 20 Dec 2021")
  shade_weekends(yl[1], yl[2])
  lines(DF$k, DF$T_meas, col = col_meas, lwd = 0.8)
  lines(DF$k, DF$v3,     col = faint(col_v3, 0.85), lwd = 0.7)
  lines(DF$k, DF$v6,     col = faint(col_v6, 0.85), lwd = 0.7)
  tick_axis()
  legend("topright", bty = "n", lwd = 2, col = c(col_meas, col_v3, col_v6),
         legend = c("Measured", "v3 (1-zone)", "v6 (2-zone)"), horiz = TRUE, cex = 0.9)

  # (B) residuals (deviation) + 24-h moving average
  par(mar = c(3.6, 4.8, 1.2, 1.2))
  yb <- range(c(DF$r_v3, DF$r_v6), na.rm = TRUE)
  plot(DF$k, DF$r_v3, type = "n", ylim = yb, xaxt = "n", xlab = "", ylab = lab_res)
  shade_weekends(yb[1], yb[2])
  abline(h = 0, col = "grey20", lty = 3)
  lines(DF$k, DF$r_v3, col = faint(col_v3, 0.18), lwd = 0.5)   # raw 5-min
  lines(DF$k, DF$r_v6, col = faint(col_v6, 0.18), lwd = 0.5)
  lines(DF$k, cma(DF$r_v3), col = col_v3, lwd = 2.4)           # 24-h MA
  lines(DF$k, cma(DF$r_v6), col = col_v6, lwd = 2.4)
  tick_axis()
  mtext("Date (2021)", side = 1, line = 2.3, cex = 0.95)
  legend("bottomleft", bty = "n", ncol = 2, lwd = c(0.8,2.4,0.8,2.4),
         col = c(faint(col_v3,0.5), col_v3, faint(col_v6,0.5), col_v6),
         legend = c("v3 raw 5-min", "v3 24-h MA", "v6 raw 5-min", "v6 24-h MA"))
  legend("topright", bty = "n", cex = 0.9,
         legend = c(sprintf("v3:  NMBE %+.1f%%   CV-RMSE %.1f%%", NMBE_v3, CV_v3),
                    sprintf("v6:  NMBE %+.1f%%   CV-RMSE %.1f%%", NMBE_v6, CV_v6),
                    "shaded = weekend (free-running)"))
}

# --- save PNG + PDF (wide, to resolve the 50-day 5-min detail) ---------------
png(file.path(OUTFIG, "deviation_fullseries.png"),
    width = 15, height = 9, units = "in", res = 150); draw(); dev.off()
pdf(file.path(OUTFIG, "deviation_fullseries.pdf"),
    width = 15, height = 9, encoding = "ISOLatin1", useDingbats = FALSE); draw(); dev.off()

cat("Saved deviation_fullseries.{png,pdf} to:\n  ", OUTFIG, "\n")

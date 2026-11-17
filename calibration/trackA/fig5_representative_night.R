# =============================================================================
# fig5_representative_night.R
# -----------------------------------------------------------------------------
# Figure 5 (Section 3.1.2): measured vs simulated (calibrated v3) indoor
# temperature over a REPRESENTATIVE free-running weekday night (22:00 -> 07:00),
# with the residual. The representative night is chosen objectively as the
# complete weekday night whose per-night NMBE is closest to the headline
# -5.41 %, so it illustrates typical calibration performance (not a hand-picked
# best case). Override with FIG5_NIGHT <- "YYYY-MM-DD" before sourcing.
#
# Uses the CALIBRATED series via load_calibrated_series.R (anchored to -5.41 %).
# Output: PNG + PDF in OUT_DIR/deviation_figs.
# =============================================================================

suppressWarnings({library(lubridate)})
if(!exists(".CALIB_SERIES_LOADED") || !isTRUE(.CALIB_SERIES_LOADED))
  source("load_calibrated_series.R")
stopifnot(exists("D"))

# --- nocturnal subset + night labels (evening date; early morning -> D-1) ----
hh   <- hour(D$k) + minute(D$k)/60
Dn   <- D[hh >= 22 | hh < 7, ]
hhn  <- hour(Dn$k) + minute(Dn$k)/60
night <- as.Date(Dn$k, tz = "UTC"); night[hhn < 7] <- night[hhn < 7] - 1
Dn$night <- night
Dn$wend  <- wday(night) %in% c(6, 7, 1)

# --- per-night NMBE (v3) and completeness ------------------------------------
agg <- do.call(rbind, lapply(split(Dn, Dn$night), function(g){
  data.frame(night = g$night[1], n = nrow(g), wend = g$wend[1],
             nmbe = 100*mean(g$v3 - g$T_meas)/mean(g$T_meas),
             Tout = mean(g$T_out))
}))
# representative = complete (~108 pts) WEEKDAY night closest to headline -5.41 %
cand <- agg[agg$n >= 104 & !agg$wend, ]
cand$dist <- abs(cand$nmbe - (-5.41))
cand <- cand[order(cand$dist), ]
cat("Representative weekday-night candidates (per-night NMBE closest to -5.41%):\n")
print(head(cand[, c("night","n","nmbe","Tout")], 5), row.names = FALSE)

sel <- if(exists("FIG5_NIGHT")) as.Date(FIG5_NIGHT) else cand$night[1]
cat(sprintf("\nSelected representative night: %s\n", format(sel)))

# --- extract that night and its metrics --------------------------------------
g <- Dn[Dn$night == sel, ]; g <- g[order(g$k), ]
r    <- g$v3 - g$T_meas
mbar <- mean(g$T_meas)
NMBE <- 100 * mean(r) / mbar
CVR  <- 100 * sqrt(mean(r^2)) / mbar
cat(sprintf("  night %s: n=%d  NMBE=%+.2f%%  CV-RMSE=%.2f%%  mean meas=%.2f C\n",
            format(sel), nrow(g), NMBE, CVR, mbar))

# --- styling -----------------------------------------------------------------
col_meas <- "black"; col_v3 <- "#1f6feb"; col_out <- "grey55"
lab_T   <- expression("Indoor / outdoor T ("*degree*"C)")
lab_res <- expression("Residual  (sim "-" meas)  ("*degree*"C)")
OUTFIG  <- if(exists("OUT_DIR")) file.path(OUT_DIR,"deviation_figs") else
                                  file.path(getwd(),"deviation_figs")
dir.create(OUTFIG, recursive = TRUE, showWarnings = FALSE)
t_ax <- function(k) axis.POSIXct(1, at = seq(min(k), max(k), by = "1 hour"),
                                 format = "%H:%M")

draw <- function(){
  layout(matrix(1:2, ncol = 1), heights = c(2.1, 1.0))
  # (A) measured vs simulated indoor T (+ outdoor reference)
  par(mar = c(2.2, 4.6, 2.4, 1.0))
  yl <- range(c(g$T_meas, g$v3, g$T_out), na.rm = TRUE)
  plot(g$k, g$T_meas, type = "l", lwd = 2.2, col = col_meas, ylim = yl,
       xaxt = "n", xlab = "", ylab = lab_T,
       main = "")
  lines(g$k, g$v3,    lwd = 2.2, col = col_v3)
  lines(g$k, g$T_out, lwd = 1.4, col = col_out, lty = 2)
  .en_mon <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
  .en_dow <- c("Sun","Mon","Tue","Wed","Thu","Fri","Sat")
  sel_lbl <- sprintf("%s %d %s %d", .en_dow[lubridate::wday(sel)],
                     lubridate::mday(sel), .en_mon[lubridate::month(sel)], lubridate::year(sel))
  title(main = sprintf("Representative free-running night — %s", sel_lbl), cex.main = 1.15)
  title(main = sprintf("\n\nv3 nocturnal NMBE %+.2f %%,  CV-RMSE %.2f %%", NMBE, CVR),
        cex.main = 0.92, font.main = 1)
  t_ax(g$k)
  legend("topright", bty = "n", lwd = c(2.2,2.2,1.4), lty = c(1,1,2),
         col = c(col_meas, col_v3, col_out),
         legend = c("Measured", "Simulated (v3)", "Outdoor"))
  # (B) residual
  par(mar = c(3.4, 4.6, 0.8, 1.0))
  plot(g$k, r, type = "l", lwd = 2, col = col_v3, xaxt = "n",
       xlab = "", ylab = lab_res)
  abline(h = 0, col = "grey30", lty = 3)
  t_ax(g$k); mtext("Time (22:00 -> 07:00)", side = 1, line = 2.2, cex = 0.95)
}

png(file.path(OUTFIG,"Fig5_representative_night.png"),
    width = 7, height = 6, units = "in", res = 200); draw(); dev.off()
pdf(file.path(OUTFIG,"Fig5_representative_night.pdf"),
    width = 7, height = 6, encoding = "ISOLatin1", useDingbats = FALSE); draw(); dev.off()
cat("\nSaved Fig5_representative_night.{png,pdf} to:\n  ", OUTFIG, "\n")

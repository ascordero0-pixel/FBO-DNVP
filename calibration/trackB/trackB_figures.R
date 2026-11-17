# =============================================================================
# trackB_figures.R  --  Calibration fit figures (Track B)
# Thesis T3 -- FBO framework
# -----------------------------------------------------------------------------
# Figure B.1 : CO2 simulated vs measured over the anchor window (09:00-14:00)
#              for the three campaign days, with WWeqR (switching) overlaid.
#              Shows the DNVP "loops" (T1 Fig. 8 analogue) reproduced by the
#              calibrated model, and the SNVP quasi-steady behaviour (Dec 17).
# Figure B.2 : Measured-vs-simulated scatter with 1:1 line and +/-40 ppm
#              sensor-accuracy band (TEMTOP M2000C, T1 Table 3).
#
# Uses base-R graphics (no external plotting deps) -> PNG output.
# =============================================================================

source("load_co2_campaign.R")
source("trackB_massbalance.R")
source("trackB_calibrate.R")

G_HAT <- 4.0943e-06     # m3/s (calibrated, Option A headline)

make_figures <- function(outdir = "runs") {
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  tb <- build_trackB(verbose = FALSE)
  dl <- tb$days
  day_ord <- c("18","26","17")
  day_lab <- c("18"="18 Nov 2021 (DNVP, storm)",
               "26"="26 Nov 2021 (DNVP, moderate)",
               "17"="17 Dec 2021 (SNVP / S2)")

  sims <- lapply(day_ord, function(dd) simulate_day(dl[[dd]], G_HAT))
  names(sims) <- day_ord

  # ---- Figure B.1: time-series panels ----
  png(file.path(outdir,"FigB1_trackB_fit.png"), width=1500, height=1500, res=150)
  par(mfrow=c(3,1), mar=c(4,4.5,3,4.5))
  for (dd in day_ord) {
    s <- sims[[dd]]
    tt <- as.POSIXct(s$ts)
    ymax <- max(1300, max(s$Csim, s$Cmeas, na.rm=TRUE)*1.05)
    plot(tt, s$Cmeas, type="n",
         ylim=c(400, ymax), xlab="", ylab=expression(CO[2]~"(ppm)"),
         main=day_lab[[dd]], las=1)
    # +/-40 ppm sensor-accuracy band around the measured series (drawn first)
    okm <- !is.na(s$Cmeas)
    if (sum(okm) > 1) {
      xb <- tt[okm]; yb <- s$Cmeas[okm]
      polygon(c(xb, rev(xb)), c(yb + 40, rev(yb - 40)),
              col = "#1b6ca833", border = NA)
    }
    points(tt, s$Cmeas, pch=19, col="#1b6ca8", cex=0.7)
    lines(tt, s$Csim, col="#c8102e", lwd=2)
    abline(h=800, lty=2, col="grey40")
    # WWeqR on secondary axis (switching)
    par(new=TRUE)
    plot(tt, s$WWeqR*100, type="s", col="#2e8b57", lwd=1.3, axes=FALSE,
         xlab="", ylab="", ylim=c(0, 6))
    axis(4, las=1); mtext("WWeqR (%)", side=4, line=3, cex=0.8)
    legend("topright", bty="n", cex=0.8,
           legend=c("CO2 measured","+/-40 ppm sensor band","CO2 simulated",
                    "WWeqR (switched)","800 ppm limit"),
           col=c("#1b6ca8","#1b6ca855","#c8102e","#2e8b57","grey40"),
           pch=c(19,15,NA,NA,NA), lty=c(NA,NA,1,1,2), lwd=c(NA,NA,2,1.3,1))
  }
  dev.off()

  # ---- Figure B.2: measured vs simulated scatter ----
  png(file.path(outdir,"FigB2_trackB_scatter.png"), width=1000, height=1000, res=150)
  par(mar=c(4.5,4.5,3,1))
  allC <- do.call(rbind, sims)
  ok <- !is.na(allC$Cmeas)
  cols <- c("18"="#c8102e","26"="#1b6ca8","17"="#2e8b57")
  daycol <- unlist(lapply(day_ord, function(dd) rep(cols[[dd]], nrow(sims[[dd]]))))
  lim <- c(400, max(allC$Csim, allC$Cmeas, na.rm=TRUE)*1.03)
  plot(allC$Cmeas[ok], allC$Csim[ok], col=daycol[ok], pch=19, cex=0.7,
       xlim=lim, ylim=lim, xlab="Measured CO2 (ppm)", ylab="Simulated CO2 (ppm)",
       main="Track B: simulated vs measured (anchor window)", las=1)
  abline(0,1,lwd=2)
  abline(40, 1, lty=3, col="grey50"); abline(-40, 1, lty=3, col="grey50")
  legend("topleft", bty="n",
         legend=c("18 Nov (storm)","26 Nov (moderate)","17 Dec (SNVP)",
                  "1:1 line","+/-40 ppm sensor"),
         col=c("#c8102e","#1b6ca8","#2e8b57","black","grey50"),
         pch=c(19,19,19,NA,NA), lty=c(NA,NA,NA,1,3), lwd=c(NA,NA,NA,2,1))
  dev.off()

  cat("Figures written: FigB1_trackB_fit.png, FigB2_trackB_scatter.png\n")
}

if (sys.nframe() == 0) make_figures()

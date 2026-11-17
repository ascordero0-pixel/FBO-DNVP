# =============================================================================
# figA3_tau.R   (Annex A - per-night thermal time constant)
# -----------------------------------------------------------------------------
# Figure A.3 (Appendix A.1.4): distribution of the per-night thermal time
# constant tau, estimated by an independent first-order exponential decay fit
# to each free-running night. Supports the choice of the per-night median
# (38.1 h) as the citable value and contrasts it with the pooled lumped-
# capacitance regression (57.4 h, R2 = 0.04), which the histogram range
# deliberately excludes to show it is an outlier estimator.
#
# AUTONOMOUS: reads only runs/tau_per_night.csv (produced by tau_per_night.R).
# Output: PNG + PDF in runs/deviation_figs.
# =============================================================================

CSV     <- "runs/tau_per_night.csv"
POOLED  <- 57.38          # pooled lumped-capacitance tau (R2 = 0.04) - contrast
OUTFIG  <- "runs/deviation_figs"
dir.create(OUTFIG, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(CSV))
d <- read.csv(CSV)
tau <- d$tau[is.finite(d$tau)]
med <- median(tau); mn <- mean(tau); sdv <- sd(tau)
q1  <- quantile(tau, 0.25); q3 <- quantile(tau, 0.75)
r2med <- median(d$r2)
cat(sprintf("[figA3] n=%d nights | median tau=%.2f h | mean=%.2f | SD=%.2f | IQR=[%.2f, %.2f] | median R2=%.3f\n",
            length(tau), med, mn, sdv, q1, q3, r2med))

col_bar <- "#9ecae1"; col_med <- "#08519c"; col_iqr <- "#c6dbef"; col_pool <- "#cb181d"

draw <- function(){
  par(mar = c(4.4, 4.6, 2.6, 1.2))
  # histogram of per-night tau, x-range [0, 48] (the fit-quality filter cap)
  br <- seq(0, 48, by = 3)
  h  <- hist(tau, breaks = br, plot = FALSE)
  plot(h, col = col_bar, border = "white",
       xlim = c(0, 60), ylim = c(0, max(h$counts) + 1),
       main = "Per-night thermal time constant \u03C4  (43 free-running nights)",
       xlab = expression("Time constant  "*tau*"  (h)"),
       ylab = "Number of nights", cex.main = 1.1)
  # IQR band
  rect(q1, 0, q3, max(h$counts) + 1, col = adjustcolor(col_iqr, 0.35), border = NA)
  plot(h, col = adjustcolor(col_bar, 0.9), border = "white", add = TRUE)
  # median line
  abline(v = med, col = col_med, lwd = 2.6)
  text(med, max(h$counts) + 0.6, sprintf("median = %.1f h", med),
       col = col_med, pos = 4, font = 2, cex = 0.95)
  # pooled lumped estimator (out of the per-night range) - contrast
  abline(v = POOLED, col = col_pool, lwd = 2.2, lty = 2)
  arrows(POOLED - 4, max(h$counts) * 0.55, POOLED, max(h$counts) * 0.55,
         length = 0.08, col = col_pool)
  text(POOLED, max(h$counts) * 0.62,
       sprintf("pooled lumped fit\n%.1f h  (R\u00B2 = 0.04)", POOLED),
       col = col_pool, pos = 2, cex = 0.85)
  # summary box
  legend("topleft", bty = "n", cex = 0.9,
    legend = c(sprintf("median  = %.1f h", med),
               sprintf("mean    = %.1f h", mn),
               sprintf("SD      = %.1f h", sdv),
               sprintf("IQR     = [%.1f, %.1f] h", q1, q3),
               sprintf("median R\u00B2 = %.3f", r2med)))
}

png(file.path(OUTFIG, "FigA3_tau_per_night.png"),
    width = 7.2, height = 4.6, units = "in", res = 200); draw(); dev.off()
pdf(file.path(OUTFIG, "FigA3_tau_per_night.pdf"),
    width = 7.2, height = 4.6, useDingbats = FALSE); draw(); dev.off()
cat("[OK] Saved FigA3_tau_per_night.{png,pdf} in", OUTFIG, "\n")

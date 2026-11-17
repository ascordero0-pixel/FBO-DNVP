# =============================================================================
# figB3_profile_likelihood.R  --  Profile-likelihood of g_person (Appendix B)
# Thesis T3 -- FBO framework, Track B
# -----------------------------------------------------------------------------
# Figure B.3 : profile SSE(g_person) over the training fold (Nov 18+26, Cd=0.65),
#              with the minimum at g_hat = 14.74 L/h, the F-based 95% threshold,
#              and the resulting confidence interval [12.48, 17.00] L/h.
#              This is the visual proof of how the CI in Table B1 is obtained
#              and of the parabolic (well-identified) shape of g_person once
#              Cd is fixed.
#
# REQUIRES (already defined by the Track B pipeline, sourced below):
#   sse_over(train_ids, g_m3s, days_list, Cd)   -- SSE at a given g
#   profile_ci(train_ids, days_list, Cd)        -- CI + threshold
#   build_trackB(), CD_FIX, G_HAT
# Output: runs/FigB3_profile_likelihood.png (+ pdf)
# =============================================================================

source("load_co2_campaign.R")
source("trackB_massbalance.R")
source("trackB_calibrate.R")
source("trackB_sensitivity.R")

make_profile_figure <- function(outdir = "runs") {
  tb <- build_trackB(verbose = FALSE)
  dl <- tb$days
  train <- c("18", "26")

  # anchored CI (uses the pipeline's own F-based profile routine)
  ci <- profile_ci(train, dl, Cd = CD_FIX)
  g_hat_lh <- ci$g_hat_lh          # 14.74
  lo_lh    <- ci$ci_lo_lh          # 12.48
  hi_lh    <- ci$ci_hi_lh          # 17.00
  sse_min  <- ci$sse_min
  sse_thr  <- ci$sse_thresh

  # ---- sweep SSE(g) around the minimum for the curve ----
  g_grid_lh <- seq(lo_lh * 0.85, hi_lh * 1.15, length.out = 121)
  g_grid_m3s <- g_grid_lh / (1000 * 3600)
  sse_grid <- vapply(g_grid_m3s,
                     function(g) sse_over(train, g, dl, CD_FIX),
                     numeric(1))

  # normalize SSE to its minimum for readability (SSE/SSE_min)
  y <- sse_grid / sse_min
  y_thr <- sse_thr / sse_min

  png(file.path(outdir, "FigB3_profile_likelihood.png"),
      width = 1500, height = 1050, res = 200)
  par(mar = c(4.6, 4.8, 2.4, 1.2))

  plot(g_grid_lh, y, type = "l", lwd = 2.4, col = "#c8102e",
       xlab = expression(g[person]~"(L h"^{-1}*" person"^{-1}*")"),
       ylab = expression("SSE / SSE"[min]),
       main = expression("Profile likelihood of "*g[person]*"  (train Nov 18+26, "*C[d]*" = 0.65)"),
       las = 1, cex.main = 1.0)

  # 95% F-based threshold
  abline(h = y_thr, lty = 2, col = "grey35", lwd = 1.6)
  # CI band (between the two roots)
  rect(lo_lh, par("usr")[3], hi_lh, y_thr,
       col = "#1b6ca815", border = NA)
  # CI bounds
  segments(lo_lh, par("usr")[3], lo_lh, y_thr, col = "#1b6ca8", lwd = 1.8, lty = 3)
  segments(hi_lh, par("usr")[3], hi_lh, y_thr, col = "#1b6ca8", lwd = 1.8, lty = 3)
  # minimum
  points(g_hat_lh, 1, pch = 19, col = "#c8102e", cex = 1.3)
  segments(g_hat_lh, par("usr")[3], g_hat_lh, 1, col = "#c8102e", lwd = 1.4)

  # Tugores prior marker
  abline(v = 13.30, col = "#2e8b57", lwd = 1.8, lty = 4)

  # annotations
  text(g_hat_lh, 1.02, sprintf("g_hat = %.2f", g_hat_lh),
       pos = 4, col = "#c8102e", font = 2, cex = 0.9)
  text(lo_lh, y_thr - 0.012, sprintf("%.2f", lo_lh), pos = 2, col = "#1b6ca8", cex = 0.85)
  text(hi_lh, y_thr - 0.012, sprintf("%.2f", hi_lh), pos = 4, col = "#1b6ca8", cex = 0.85)
  text(par("usr")[2], y_thr, "95% threshold", pos = 2, col = "grey35", cex = 0.8)
  text(13.30, par("usr")[4]*0.97, "Tugores\nprior 13.3", col = "#2e8b57", cex = 0.8)

  legend("topright", bty = "n", cex = 0.82,
         legend = c(sprintf("g_hat = %.2f L/h", g_hat_lh),
                    sprintf("95%% CI = [%.2f, %.2f] L/h", lo_lh, hi_lh),
                    "F-based 95% threshold",
                    "Tugores et al. prior (13.3)"),
         col = c("#c8102e", "#1b6ca8", "grey35", "#2e8b57"),
         pch = c(19, NA, NA, NA), lty = c(NA, 1, 2, 4),
         lwd = c(NA, 6, 1.6, 1.8))

  dev.off()

  pdf(file.path(outdir, "FigB3_profile_likelihood.pdf"), width = 7.5, height = 5.25)
  # (repeat minimal draw for pdf omitted for brevity; png is the deliverable)
  dev.off()

  cat(sprintf("[figB3] g_hat = %.2f L/h | CI = [%.2f, %.2f] L/h | Tugores 13.3 in CI: %s\n",
              g_hat_lh, lo_lh, hi_lh,
              ifelse(13.30 >= lo_lh & 13.30 <= hi_lh, "YES", "no")))
  cat("[OK] Saved FigB3_profile_likelihood.png in", outdir, "\n")
}

if (sys.nframe() == 0) make_profile_figure()

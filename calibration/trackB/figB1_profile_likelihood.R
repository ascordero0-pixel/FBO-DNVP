# =============================================================================
# figB1_profile_likelihood.R  --  Profile likelihood of G (Appendix B, Fig. B1)
# Thesis T3 -- FBO framework, Track B
# -----------------------------------------------------------------------------
# Figure B1 : profile SSE(G) over the training fold (18 + 26 November, Cd = 0.65),
#             with the minimum at G_hat = 14.74 L/h, the F-based 95 % threshold
#             and the resulting confidence interval [12.48, 17.00] L/h
#             (values reported in Table B4 and Section 3.1).
#
# Changes with respect to figB3_profile_likelihood.R (manuscript REV17):
#   * Notation C2: g_person -> G, g_hat -> G-hat, prior -> G0 (as in the text
#     and in Table B4). In the code, the variable keeps the name g_person.
#   * Units written as L h^-1 (as in the text).
#   * Title inside the plot removed: the caption already carries it (Elsevier style).
#   * The PDF was being saved blank: drawing moved into a function that is
#     called for both devices.
#   * Output renamed to FigB1_* (new appendix numbering).
#
# REQUIRES (already defined by the Track B pipeline, sourced below):
#   sse_over(train_ids, g_m3s, days_list, Cd)   -- SSE at a given g
#   profile_ci(train_ids, days_list, Cd)        -- CI + threshold
#   build_trackB(), CD_FIX
# Output: runs/FigB1_profile_likelihood.png (+ pdf)
# =============================================================================

# ---- working directory = folder of this script (Track B pipeline) ----------
# The pipeline files are sourced with relative paths, so R must run from the
# TRACK B folder. The block below tries to detect it automatically
# (source() or RStudio); otherwise it uses TRACKB_DIR.
TRACKB_DIR <- getwd()   # run from calibration/trackB [repo patch P3]
.script_dir <- tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) NULL)
if (is.null(.script_dir) && requireNamespace("rstudioapi", quietly = TRUE) &&
    rstudioapi::isAvailable()) {
  .p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
  if (nzchar(.p)) .script_dir <- dirname(.p)
}
if (is.null(.script_dir) || !file.exists(file.path(.script_dir, "load_co2_campaign.R")))
  .script_dir <- TRACKB_DIR
setwd(.script_dir)
message("[figB1] working directory: ", getwd())

source("load_co2_campaign.R")
source("trackB_massbalance.R")
source("trackB_calibrate.R")
source("trackB_sensitivity.R")

G0_PRIOR <- 13.30   # Tugores et al. prior (L/h)

draw_profile <- function(g_grid_lh, y, y_thr, g_hat_lh, lo_lh, hi_lh) {
  par(mar = c(4.6, 4.8, 1.2, 1.2))

  plot(g_grid_lh, y, type = "l", lwd = 2.4, col = "#c8102e",
       xlab = expression(italic(G) ~ "(L h"^{-1} * ")"),
       ylab = expression("SSE / SSE"[min]),
       main = NULL, las = 1)

  # 95 % F-based threshold
  abline(h = y_thr, lty = 2, col = "grey35", lwd = 1.6)
  # CI band (between the two roots)
  rect(lo_lh, par("usr")[3], hi_lh, y_thr, col = "#1b6ca815", border = NA)
  # CI bounds
  segments(lo_lh, par("usr")[3], lo_lh, y_thr, col = "#1b6ca8", lwd = 1.8, lty = 3)
  segments(hi_lh, par("usr")[3], hi_lh, y_thr, col = "#1b6ca8", lwd = 1.8, lty = 3)
  # minimum
  points(g_hat_lh, 1, pch = 19, col = "#c8102e", cex = 1.3)
  segments(g_hat_lh, par("usr")[3], g_hat_lh, 1, col = "#c8102e", lwd = 1.4)

  # prior G0 (Tugores et al.)
  abline(v = G0_PRIOR, col = "#2e8b57", lwd = 1.8, lty = 4)

  # annotations
  text(g_hat_lh, 1.02,
       bquote(hat(italic(G)) == .(sprintf("%.2f", g_hat_lh))),
       pos = 4, col = "#c8102e", font = 2, cex = 0.9)
  text(lo_lh, y_thr - 0.012, sprintf("%.2f", lo_lh), pos = 2, col = "#1b6ca8", cex = 0.85)
  text(hi_lh, y_thr - 0.012, sprintf("%.2f", hi_lh), pos = 4, col = "#1b6ca8", cex = 0.85)
  text(par("usr")[2], y_thr, "95 % threshold", pos = 2, col = "grey35", cex = 0.8)
  text(G0_PRIOR, par("usr")[4] * 0.97,
       bquote(italic(G)[0] == .(sprintf("%.1f", G0_PRIOR))),
       col = "#2e8b57", cex = 0.8)

  legend("topright", bty = "n", cex = 0.82,
         legend = c(
           as.expression(bquote(hat(italic(G)) == .(sprintf("%.2f", g_hat_lh)) ~ "L h"^{-1})),
           as.expression(bquote("95 % CI = [" * .(sprintf("%.2f", lo_lh)) * ", " *
                                  .(sprintf("%.2f", hi_lh)) * "] L h"^{-1})),
           "F-based 95 % threshold",
           as.expression(bquote("Prior " * italic(G)[0] == .(sprintf("%.1f", G0_PRIOR)) ~
                                  "L h"^{-1} ~ "(Tugores et al.)"))),
         col = c("#c8102e", "#1b6ca8", "grey35", "#2e8b57"),
         pch = c(19, NA, NA, NA), lty = c(NA, 1, 2, 4),
         lwd = c(NA, 6, 1.6, 1.8))
}

make_profile_figure <- function(outdir = "runs") {
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
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

  # ---- sweep SSE(G) around the minimum for the curve ----
  g_grid_lh  <- seq(lo_lh * 0.85, hi_lh * 1.15, length.out = 121)
  g_grid_m3s <- g_grid_lh / (1000 * 3600)
  sse_grid <- vapply(g_grid_m3s,
                     function(g) sse_over(train, g, dl, CD_FIX),
                     numeric(1))

  # normalize SSE to its minimum for readability (SSE / SSE_min)
  y     <- sse_grid / sse_min
  y_thr <- sse_thr / sse_min

  png(file.path(outdir, "FigB1_profile_likelihood.png"),
      width = 1500, height = 1000, res = 200)
  draw_profile(g_grid_lh, y, y_thr, g_hat_lh, lo_lh, hi_lh)
  dev.off()

  pdf(file.path(outdir, "FigB1_profile_likelihood.pdf"), width = 7.5, height = 5.0)
  draw_profile(g_grid_lh, y, y_thr, g_hat_lh, lo_lh, hi_lh)
  dev.off()

  cat(sprintf("[figB1] G_hat = %.2f L/h | CI = [%.2f, %.2f] L/h | prior %.1f in CI: %s\n",
              g_hat_lh, lo_lh, hi_lh, G0_PRIOR,
              ifelse(G0_PRIOR >= lo_lh & G0_PRIOR <= hi_lh, "YES", "no")))
  cat("[OK] Saved FigB1_profile_likelihood.png/.pdf in", outdir, "\n")
}

# Runs both with Rscript and with source() (the original only ran with Rscript).
make_profile_figure()

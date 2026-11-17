# =============================================================================
# F2/figs.R  --  T3 FBO Phase 2: benchmark figures (vector PDF + 300 dpi PNG).
#   Fig 1  SAMPLE EFFICIENCY: best-U vs number of E+ evaluations (FBO vs NSGA-II)
#          -> makes the ratio N_NSGA / N_FBO visible (metric a).
#   Fig 2  Real PARETO FRONT (f1=CO2 vs f3=Cost, colour = f2=comfort) with the
#          FBO optimum marked -> illustrates (b) on the front and the trade-off.
#
#   Base R (no ggplot2). Uses the CURRENT FIS (re-scores the logs) -> consistent with SC7.
#   Data: F2/out/de_{log,front}_<scn>.csv (from nsga2_direct_ep.R). Output: F2/fig/.
#
# USAGE:
#   source(here::here("F2","figs.R"))
#   make_figs()                                  # both figures, the 2 scenarios (language = LANG, default "en")
#   make_figs(lang = "es")                       # Spanish version (e.g. for the team presentation)
#   # or individually:  fig_efficiency(lang="en");  fig_front("N_w410", lang="es")
#   # Figure-text language: global variable LANG ("en"|"es") or the lang= argument.
# =============================================================================

source(here::here("F2", "metrics_bc.R"))   # f2_load_training, f2_fbo_optimum, cost_transparent, .nd_mask, SF, VARS, EPS_U
DE_OUT  <- here::here("F2", "out")
FIG_DIR <- here::here("F2", "fig"); if (!dir.exists(FIG_DIR)) dir.create(FIG_DIR, recursive = TRUE)

COL_FBO  <- "#D8880F"   # amber  (FBO)
COL_NSGA <- "#1C7293"   # teal   (NSGA-II)
GREYL    <- "#9AA7AB"

# --- i18n: language of the FIGURE TEXT --------------------------------------
#   "en" -> paper ;  "es" -> team presentation.
#   Change LANG here, or pass it as an argument: make_figs(lang = "es").
LANG <- "en"
.LABS <- list(
  en = list(
    eff_x   = "EnergyPlus evaluations",
    eff_y   = "Best utility U attained",
    eff_tgt = "FBO optimum quality (U*−ε)",   # U*-epsilon
    fro_x   = "f1 = CO2 excess  (P95, ppm)  →  worse IAQ",
    fro_y   = "f3 = Operational cost  →  worse",
    fro_main= "Real Pareto front — %s",
    fro_opt = "FBO optimum",
    fro_bar = "f2 = p_III (% Cat III)"
  ),
  es = list(
    eff_x   = "Evaluaciones EnergyPlus",
    eff_y   = "Mejor utilidad U alcanzada",
    eff_tgt = "calidad del óptimo FBO (U*−ε)",
    fro_x   = "f1 = exceso de CO2  (P95, ppm)  →  peor IAQ",
    fro_y   = "f3 = Coste operativo  →  peor",
    fro_main= "Frente de Pareto real — %s",
    fro_opt = "óptimo FBO",
    fro_bar = "f2 = p_III (% Cat III)"
  )
)
.lab <- function(lang) { if (is.null(.LABS[[lang]])) stop("LANG must be 'en' or 'es'"); .LABS[[lang]] }

# --- load a de_log and re-score U with the current FIS -----------------------
.load_de <- function(scn) {
  d <- read.csv(file.path(DE_OUT, paste0("de_log_", scn, ".csv")), stringsAsFactors = FALSE)
  d <- d[order(d$eval), ]
  d$U <- vapply(seq_len(nrow(d)), function(i)
    if (!is.finite(d$g1[i])) NA_real_ else
      scoring_utility(SF, g1 = d$g1[i], g2 = c(d$g2_pII[i], d$g2_pIII[i]),
                      g3 = d$g3[i], g4 = d$g4[i], g5 = d$g5[i]), numeric(1))
  d
}

# helper: open PDF and PNG, draw with draw(), close both
.render <- function(base, w, h, draw) {
  pdf(file.path(FIG_DIR, paste0(base, ".pdf")), width = w, height = h, pointsize = 11); draw(); dev.off()
  png(file.path(FIG_DIR, paste0(base, ".png")), width = w, height = h, units = "in", res = 300, pointsize = 11); draw(); dev.off()
  cat("  ->", file.path(FIG_DIR, paste0(base, ".{pdf,png}")), "\n")
}

# ---- FIG 1: sample efficiency (one panel per scenario) ---------------------
.eff_panel <- function(scn, L = .lab(LANG)) {
  d <- .load_de(scn)
  nsga <- cummax(ifelse(is.finite(d$U), d$U, -Inf))          # NSGA-II incumbent
  bo   <- f2_load_training(scn)                              # seed + BO in evaluation order
  fbo  <- cummax(bo$U)
  opt  <- max(bo$U); tgt <- opt - EPS_U
  Nf   <- which(fbo  >= tgt)[1]
  Nn   <- which(nsga >= tgt)[1]
  xmax <- nrow(d); ymin <- min(c(nsga[is.finite(nsga)], fbo), na.rm = TRUE)

  plot(NA, xlim = c(1, xmax), ylim = c(ymin, opt + 0.006), xlab = L$eff_x,
       ylab = L$eff_y, main = scn, cex.main = 1.15, font.main = 1, las = 1, bty = "n")
  abline(h = tgt, lty = 3, col = GREYL)
  text(xmax, tgt, L$eff_tgt, col = GREYL, cex = 0.72, adj = c(1, 1.5))  # BELOW the line (avoids the ratio)
  lines(seq_along(fbo),  fbo,  col = COL_FBO,  lwd = 3)
  lines(seq_len(xmax),   nsga[seq_len(xmax)], col = COL_NSGA, lwd = 3)
  # N_FBO and N_NSGA markers
  if (!is.na(Nf)) { segments(Nf, ymin, Nf, tgt, col = COL_FBO, lty = 2);  points(Nf, tgt, pch = 19, col = COL_FBO, cex = 1.1) }
  if (!is.na(Nn)) { segments(Nn, ymin, Nn, tgt, col = COL_NSGA, lty = 2); points(Nn, tgt, pch = 19, col = COL_NSGA, cex = 1.1) }
  ratio <- if (!is.na(Nn) && !is.na(Nf)) sprintf("%.1f×", Nn / Nf) else ">budget"
  legend("bottomright", bty = "n", cex = 0.9, lwd = 3, col = c(COL_FBO, COL_NSGA),
         legend = c(sprintf("FBO  (N = %s)", ifelse(is.na(Nf), "—", Nf)),
                    sprintf("NSGA-II  (N = %s)", ifelse(is.na(Nn), ">bud", Nn))))
  mtext(sprintf("ratio N_NSGA / N_FBO = %s", ratio), side = 3, line = -1.4, adj = 0.98, cex = 0.9, col = "#333333", font = 2)
  invisible(list(scn = scn, Nf = Nf, Nn = Nn))
}
fig_efficiency <- function(scns = c("N_w410", "NE_w216"), base = "Fig_efficiency", lang = LANG) {
  L <- .lab(lang)
  cat("[fig] sample efficiency ...\n")
  .render(base, w = 4.7 * length(scns), h = 4.4, function() {
    op <- par(mfrow = c(1, length(scns)), mar = c(4.2, 4.2, 2.4, 1.2), mgp = c(2.4, 0.7, 0)); on.exit(par(op))
    for (s in scns) .eff_panel(s, L)
  })
}

# ---- FIG 2: Pareto front (f1 vs f3, colour = f2) + FBO optimum --------------
fig_front <- function(scn, base = NULL, lang = LANG) {
  L <- .lab(lang)
  if (is.null(base)) base <- paste0("Fig_front_", scn)
  d <- .load_de(scn); real <- d[is.finite(d$U), , drop = FALSE]
  F <- as.matrix(real[, c("f1", "f2", "f3")]); nd <- .nd_mask(F); fr <- real[nd, , drop = FALSE]
  opt <- f2_fbo_optimum(scn)
  # colour by f2 (p_III, comfort): viridis (colourblind-safe, R>=3.6)
  pal <- hcl.colors(100, "viridis")
  f2rng <- range(fr$g2_pIII); to_col <- function(v) pal[pmax(1, pmin(100, round((v - f2rng[1]) / diff(f2rng) * 99) + 1))]
  cat("[fig] Pareto front", scn, "...\n")
  .render(base, w = 6.2, h = 5.0, function() {
    op <- par(mar = c(4.2, 4.4, 2.6, 5.6), mgp = c(2.5, 0.7, 0)); on.exit(par(op))
    plot(fr$f1, fr$f3, pch = 19, col = to_col(fr$g2_pIII), cex = 1.05, las = 1, bty = "n",
         xlab = L$fro_x,
         ylab = L$fro_y,
         main = sprintf(L$fro_main, scn), font.main = 1, cex.main = 1.1)
    # FBO optimum as an amber star
    points(opt$f["f1"], opt$f["f3"], pch = 8, col = COL_FBO, cex = 2.2, lwd = 3)
    points(opt$f["f1"], opt$f["f3"], pch = 21, bg = NA, col = COL_FBO, cex = 3.0, lwd = 1.5)
    text(opt$f["f1"], opt$f["f3"], L$fro_opt, pos = 3, col = COL_FBO, font = 2, cex = 0.9, offset = 0.9)
    # colour bar for f2 (p_III)
    xr <- par("usr"); bx <- xr[2] + diff(xr[1:2]) * 0.03; bw <- diff(xr[1:2]) * 0.03
    ys <- seq(xr[3], xr[4], length.out = 101)
    rect(bx, ys[-101], bx + bw, ys[-1], col = pal, border = NA, xpd = TRUE)
    text(bx + bw, xr[4], sprintf(" %.1f", f2rng[2]), xpd = TRUE, adj = c(0, 1), cex = 0.8)
    text(bx + bw, xr[3], sprintf(" %.1f", f2rng[1]), xpd = TRUE, adj = c(0, 0), cex = 0.8)
    text(bx + bw*1.6, mean(xr[3:4]), L$fro_bar, xpd = TRUE, srt = 90, cex = 0.85)
  })
  invisible(fr)
}

# ---- everything ------------------------------------------------------------
make_figs <- function(scns = c("N_w410", "NE_w216"), lang = LANG) {
  fig_efficiency(scns, lang = lang)
  for (s in scns) fig_front(s, lang = lang)
  cat("Figures written to:", FIG_DIR, "\n")
}

if (sys.nframe() == 0L) make_figs()

# ============================================================================
# indicators.R  --  T3 FBO framework, Phase 0: performance-indicator computation.
#
# Maps the 5-min simulation series of ONE candidate/scenario to the direct
# performance indicators g1..g4 (g5 = outdoor-import is Tanda 2, separate).
# Native R, no EnergyPlus dependency: pure time-series post-processing.
# Feeds fis_scoring.R:  scoring_utility(sf, g1, g2=c(pII,pIII), g3, g4, g5).
#
# Closed definitions (see manuscript S2, revised):
#   Occupancy mask T_occ : passed as a logical column `occ` (L-Fri 09:00-14:00
#                          within the school calendar; summer closure excluded).
#   Adaptive comfort     : theta_comf(d) = 0.33*clamp(theta_rm(d),10,30) + 18.8
#                          theta_rm(d) = (1-a)*theta_ed(d-1) + a*theta_rm(d-1), a=0.8
#                          (recursive EN 16798-1 running mean; computed over the
#                           WHOLE year so summer memory reaches September).
#
#   PALANCA 2 (air-speed comfort extension, EN 16798-1 sec.6 / ASHRAE 55 sec.5.3.3):
#     when the window is OPEN (s=1) the elevated indoor air speed raises the
#     acceptable UPPER operative-temperature limit by `dtheta_airspeed` [C].
#     Applied per-timestep to the UPPER bounds of g2 only (air movement does not
#     offset cold, so the lower bounds are unchanged). Rewards ventilation
#     exactly when it helps (warm hours with the window open).
#
#   g1  = P95{ C_in - C_out : T_occ }                                    [ppm]
#   g2  = ( p_II , p_III ),  p_* = 100*|theta_op outside Cat II(+-3)/III(+-4)| / |T_occ|
#         upper band widened by dtheta_airspeed*s when the window is open
#   g3  = mean_d  sum_{T_occ(d)} max(0,(theta_comf-3)-theta_op)*dt        [C.min/day]
#   g4  = mean_d  #{ s(t) != s(t-dt) within T_occ(d) }                    [actions/session]
#
# Input data.frame columns (5-min, time-ordered within day):
#   day, theta_op, C_in, C_out, theta_out, s (0/1 window state), occ (logical)
# ============================================================================


# ---------------------------------------------------------------------------
# adaptive comfort temperature by day (recursive running mean + clamp)
# ---------------------------------------------------------------------------
comfort_by_day <- function(df, alpha = 0.8, clamp_rm = c(10, 30)) {
  days     <- sort(unique(df$day))
  theta_ed <- tapply(df$theta_out, df$day, mean)[as.character(days)]   # daily mean outdoor
  theta_rm <- numeric(length(days)); names(theta_rm) <- as.character(days)
  theta_rm[1] <- theta_ed[1]                                          # initialisation
  if (length(days) > 1)
    for (i in 2:length(days))
      theta_rm[i] <- (1 - alpha) * theta_ed[i - 1] + alpha * theta_rm[i - 1]
  rm_cl <- pmin(pmax(theta_rm, clamp_rm[1]), clamp_rm[2])
  0.33 * rm_cl + 18.8                                                  # named by day
}


# ---------------------------------------------------------------------------
# main: compute g1..g4
# ---------------------------------------------------------------------------
compute_indicators <- function(df, alpha = 0.8, dt_min = 5,
                               clamp_rm = c(10, 30), p95_type = 7,
                               dtheta_airspeed = 1.5) {
  stopifnot(all(c("day","theta_op","C_in","C_out","theta_out","s","occ") %in% names(df)))
  occ  <- as.logical(df$occ)
  tcd  <- comfort_by_day(df, alpha, clamp_rm)
  tc   <- tcd[as.character(df$day)]                 # comfort temp broadcast to timesteps

  # g1 -- P95 of CO2 excess over occupied hours
  exc <- (df$C_in - df$C_out)[occ]
  g1  <- unname(stats::quantile(exc, 0.95, type = p95_type))

  # g2 -- nested-category exceedance percentages (pair)
  # Palanca 2: elevated-air-speed extension raises the UPPER comfort limit by
  # dtheta_airspeed when the window is open (s=1). Lower limits unchanged.
  nocc   <- sum(occ)
  up_off <- dtheta_airspeed * df$s
  outII  <- occ & (df$theta_op < tc - 3 | df$theta_op > tc + 3 + up_off)
  outIII <- occ & (df$theta_op < tc - 4 | df$theta_op > tc + 4 + up_off)
  pII    <- 100 * sum(outII)  / nocc
  pIII   <- 100 * sum(outIII) / nocc

  # g3 -- overcooling indicator (mean daily overcooling below lower Cat II bound)
  deficit  <- pmax(0, (tc - 3) - df$theta_op)
  ovc_step <- ifelse(occ, deficit * dt_min, 0)
  ovc_day  <- tapply(ovc_step, df$day, sum)
  occ_days <- sort(unique(df$day[occ]))
  g3       <- mean(ovc_day[as.character(occ_days)])

  # g4 -- control effort (mean window-state transitions per occupied session)
  s_occ <- df$s[occ]; d_occ <- df$day[occ]
  Nd    <- tapply(s_occ, d_occ, function(sv) if (length(sv) < 2) 0 else sum(diff(sv) != 0))
  g4    <- mean(Nd)

  list(g1 = g1, g2 = c(p_II = pII, p_III = pIII), g3 = g3, g4 = g4,
       theta_comf_by_day = tcd)
}


# ---------------------------------------------------------------------------
# validation harness (synthetic series with analytic answers)
# ---------------------------------------------------------------------------
# row builder: n occupied (or not) steps for one day
.mk <- function(day, theta_op, excess, s, occ = TRUE, theta_out = 14, C_out = 420) {
  n <- max(length(theta_op), length(excess), length(s))
  data.frame(day = day,
             theta_op  = rep_len(theta_op, n),
             C_out     = C_out,
             C_in      = C_out + rep_len(excess, n),
             theta_out = theta_out,
             s         = rep_len(s, n),
             occ       = occ)
}

run_tests_indicators <- function() {
  e <- new.env(); e$ok <- 0L; e$fail <- 0L
  chk <- function(cond, msg) {
    if (isTRUE(cond)) e$ok <- e$ok + 1L
    else { e$fail <- e$fail + 1L; cat("  FAIL:", msg, "\n") }
  }
  approx <- function(a, b, tol = 1e-6) abs(a - b) < tol

  TC14 <- 0.33 * 14 + 18.8            # = 23.42  (theta_out=14 -> theta_rm=14, unclamped)

  cat("[C] comfort temperature + clamp\n")
  chk(approx(unname(comfort_by_day(.mk(1, 22, 0, 0))["1"]), TC14), "theta_comf(14) wrong")
  chk(approx(unname(comfort_by_day(.mk(1, 22, 0, 0, theta_out = 5))["1"]), 0.33*10 + 18.8),
      "cold clamp (rm<10) not applied")
  chk(approx(unname(comfort_by_day(.mk(1, 22, 0, 0, theta_out = 35))["1"]), 0.33*30 + 18.8),
      "hot clamp (rm>30) not applied")

  cat("[g1] P95 of CO2 excess, occupancy masking\n")
  df1 <- rbind(.mk(1, 22, rep(500, 60), 0),                      # 60 occ steps, excess 500
               .mk(1, 22, rep(9999, 20), 0, occ = FALSE))        # junk non-occupied
  r1 <- compute_indicators(df1)
  chk(approx(r1$g1, 500), paste("g1 =", round(r1$g1, 3), "expected 500 (masking?)"))

  cat("[g2] nested-ladder exceedance percentages (windows closed, s=0)\n")
  # theta_comf=23.42 -> Cat II [20.42,26.42], Cat III [19.42,27.42]
  top <- c(rep(23.42, 40), rep(27.0, 10), rep(28.0, 10))  # 40 inside, 10 out-II-in-III, 10 out-III
  r2  <- compute_indicators(.mk(1, top, 0, 0))            # s=0 -> Palanca 2 inactive
  chk(approx(r2$g2["p_II"],  100 * 20 / 60, 1e-4), paste("p_II =", round(r2$g2["p_II"], 3)))
  chk(approx(r2$g2["p_III"], 100 * 10 / 60, 1e-4), paste("p_III =", round(r2$g2["p_III"], 3)))
  chk(r2$g2["p_III"] <= r2$g2["p_II"], "nesting p_III <= p_II violated")

  cat("[g2b] air-speed extension: open window widens the UPPER band only\n")
  # same hot series but window OPEN (s=1), dtheta=1.5 -> Cat II upper 26.42->27.92,
  # Cat III upper 27.42->28.92. Now the 27.0 block is inside Cat II, the 28.0 inside Cat III.
  r2o <- compute_indicators(.mk(1, top, 0, 1))            # s=1
  chk(approx(r2o$g2["p_II"],  100 * 10 / 60, 1e-4), paste("open p_II =", round(r2o$g2["p_II"],3), "expected ~16.7"))
  chk(approx(r2o$g2["p_III"], 0, 1e-4),             paste("open p_III =", round(r2o$g2["p_III"],3), "expected 0"))
  # cold is NOT offset by air speed: an open-window cold series still counts
  cold <- compute_indicators(.mk(1, rep(15, 60), 0, 1))  # 15 << 19.42 -> out-III cold
  chk(approx(cold$g2["p_III"], 100), "air speed must NOT offset cold")

  cat("[g3] overcooling indicator below lower bound\n")
  # lower bound = 23.42 - 3 = 20.42 ; theta_op = 18.42 -> deficit 2 C over 60 steps
  r3a <- compute_indicators(.mk(1, rep(18.42, 60), 0, 0))
  chk(approx(r3a$g3, 2 * 60 * 5), paste("g3 all-deficit =", round(r3a$g3, 2), "expected 600"))
  # half deficit-2, half inside band (no penalty) -> 2*30*5 = 300
  r3b <- compute_indicators(.mk(1, c(rep(18.42, 30), rep(24, 30)), 0, 0))
  chk(approx(r3b$g3, 2 * 30 * 5), paste("g3 half-deficit =", round(r3b$g3, 2), "expected 300"))
  # two identical all-deficit days -> mean stays 600
  r3c <- compute_indicators(rbind(.mk(1, rep(18.42, 60), 0, 0), .mk(2, rep(18.42, 60), 0, 0)))
  chk(approx(r3c$g3, 600), paste("g3 two-day mean =", round(r3c$g3, 2)))

  cat("[g4] control-effort transitions per session\n")
  s_seq <- c(0,0,1,1,0,1,1)      # transitions at pos 3,5,6 -> 3 toggles
  r4 <- compute_indicators(.mk(1, 22, 0, s_seq))
  chk(r4$g4 == 3, paste("g4 single-day =", r4$g4, "expected 3"))
  r4b <- compute_indicators(rbind(.mk(1, 22, 0, c(0,1,0,1)),      # 3 toggles
                                  .mk(2, 22, 0, c(0,1,0,1,0))))   # 4 toggles
  chk(approx(r4b$g4, (3 + 4) / 2), paste("g4 two-day mean =", r4b$g4, "expected 3.5"))

  cat("[I] end-to-end integration with the scoring FIS (g5 placeholder)\n")
  # location-agnostic: use the scoring FIS if it is already loaded (MASTER_FIS sources
  # core/fis_scoring.R before this), instead of a bare-path source() that breaks under core/.
  if (exists("build_scoring_fis", mode = "function") &&
      exists("scoring_utility",   mode = "function")) {
    sf <- build_scoring_fis()
    # a clean-comfort, low-CO2, cheap config -> expect high U
    good <- .mk(1, 23.42, rep(300, 60), c(rep(0, 30), rep(1, 30)))   # inside band, low excess, 1 toggle
    ri   <- compute_indicators(good)
    U    <- scoring_utility(sf, g1 = ri$g1, g2 = ri$g2, g3 = ri$g3, g4 = ri$g4, g5 = 1.5)
    chk(U > 0.70, paste("integration U =", round(U, 3), "expected Good (>0.70)"))
    cat(sprintf("      indicators: g1=%.1f  g2=(%.1f, %.1f)  g3=%.1f  g4=%.1f  -> U=%.3f\n",
                ri$g1, ri$g2["p_II"], ri$g2["p_III"], ri$g3, ri$g4, U))
  } else cat("      (scoring FIS not loaded; integration test skipped -- run via MASTER_FIS)\n")

  cat(sprintf("\nRESULT: %d passed, %d failed\n", e$ok, e$fail))
  if (e$fail > 0) stop("indicator validation failed")
  invisible(e$ok)
}

if (sys.nframe() == 0L) run_tests_indicators()

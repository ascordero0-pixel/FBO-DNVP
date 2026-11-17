# ============================================================================
# control_law.R  --  T3 FBO framework, Phase 0: crisp DNVP opening law (Eq. 9).
#
# The per-time-step ventilation controller that runs INSIDE EnergyPlus. Given the
# instantaneous state (indoor CO2, outdoor air-quality sub-index) and the five
# operational set-points x, it returns the binary window state s(t). This is the
# ONLY object the optimizer moves; the scoring FIS (fis_scoring.R) is separate.
#
# Decision (evaluated at every actuation-slot boundary, then held):
#   s = [ DNVP_req AND (AQI_ext <= tau_AQ) ]  OR  [ C_CO2 > C_ceiling ]
#
# Design (confirmed defaults):
#   * DNVP_req is THERMALLY BLIND: a Schmitt trigger on CO2 only --
#       open  when C_CO2 > tau_CO2 ; close when C_CO2 < tau_CO2 - delta_hys ;
#       hold in between (delta_hys = closing hysteresis -> no chattering).
#     Thermal protection lives entirely in the scoring (g2, g3), not the control.
#   * s(t) in {0,1} is the sole time signal (what g4 counts, what E+ reads).
#     f_open is the STATIC opening-area fraction when s=1 (a set-point in x); it
#     scales Q but does not appear in s(t).  Returned as `area = s * f_open`.
#   * d_slot is the DECISION GRID: Eq. 9 is (re)evaluated only at slot boundaries
#     (every d_slot minutes); s is held between them. Gives N_max = T_occ/d_slot.
#   * C_ceiling is a FIXED health override (1500 ppm; not in x), above tau_CO2's
#     range; it forces opening even when the outdoor-quality gate is shut.
#   * AQI_ext is the worst sub-index of the OUTDOOR air (boundary condition).
#
# The R module implements the decision LOGIC and is validated on synthetic input
# trajectories; the closed loop (C_CO2 depends on s) is closed by EnergyPlus.
# Base R only. Location-agnostic: the g4-integration self-check uses
# compute_indicators() only if a driver (MASTER_FIS.R / smoke_test.R) has already
# sourced indicators.R -- this module never sources a sibling by path.
# ============================================================================


# ---------------------------------------------------------------------------
# 1. Control law
# ---------------------------------------------------------------------------
# C_CO2   : indoor CO2 series [ppm]           (5-min grid)
# AQI_ext : outdoor worst sub-index series    (1..6, same grid)
# x set-points: tau_CO2 [ppm], delta_hys [ppm], d_slot [min], f_open [-], tau_AQ [1..6]
control_law <- function(C_CO2, AQI_ext,
                        tau_CO2, delta_hys, d_slot, f_open, tau_AQ,
                        C_ceiling = 1500, dt_min = 5, s0 = 0L) {
  n <- length(C_CO2)
  stopifnot(length(AQI_ext) == n, f_open >= 0, f_open <= 1,
            delta_hys >= 0, d_slot >= dt_min)
  slot_steps <- max(1L, as.integer(round(d_slot / dt_min)))

  s      <- integer(n)
  dnvp   <- as.integer(s0)   # Schmitt CO2-demand state (carried across slots)
  cur_s  <- as.integer(s0)   # held window state between slot boundaries

  for (t in seq_len(n)) {
    if (((t - 1L) %% slot_steps) == 0L) {          # actuation-slot boundary
      # --- Schmitt trigger on CO2 (closing hysteresis) ---
      if (C_CO2[t] > tau_CO2)                 dnvp <- 1L
      else if (C_CO2[t] < tau_CO2 - delta_hys) dnvp <- 0L
      # else: hold dnvp

      # --- Eq. 9: quality-gated demand OR health override ---
      open_gated <- (dnvp == 1L) && (AQI_ext[t] <= tau_AQ)
      override   <- C_CO2[t] > C_ceiling
      cur_s      <- as.integer(open_gated || override)
    }
    s[t] <- cur_s
  }
  list(s = s, area = s * f_open, slot_steps = slot_steps,
       n_switch = sum(abs(diff(s)) != 0))
}


# ---------------------------------------------------------------------------
# 2. Validation harness (synthetic trajectories)
# ---------------------------------------------------------------------------
run_tests_control <- function() {
  e <- new.env(); e$ok <- 0L; e$fail <- 0L
  chk <- function(c, m) if (isTRUE(c)) e$ok <- e$ok + 1L
                        else { e$fail <- e$fail + 1L; cat("  FAIL:", m, "\n") }
  clean <- 1L; dirty <- 6L                       # outdoor AQI: clean vs. dirty

  cat("[B] baseline demand logic\n")
  # low CO2 (below close threshold) -> never opens
  r <- control_law(rep(700, 60), rep(clean, 60), 1000, 100, 5, 0.5, 3)
  chk(all(r$s == 0L), "low CO2 should keep window shut")
  # high CO2, clean air -> opens and stays
  r <- control_law(rep(1100, 60), rep(clean, 60), 1000, 100, 5, 0.5, 3)
  chk(all(r$s[-1] == 1L), "high CO2 + clean air should open")
  chk(all(r$area[-1] == 0.5), "area should equal f_open when open")

  cat("[H] closing hysteresis prevents chattering\n")
  # CO2 oscillates in [940,1060]: crosses tau_CO2=1000 many times but never drops
  # below the close threshold 900 -> exactly one 0->1 transition, then held.
  co2 <- 1000 + 60 * sin(seq(0, 20 * pi, length.out = 240))
  r <- control_law(co2, rep(clean, 240), 1000, 100, 5, 0.5, 3)
  chk(r$n_switch <= 1, paste("hysteresis chattering: n_switch =", r$n_switch))
  # sanity: WITHOUT hysteresis (delta_hys=0) the same signal chatters a lot
  r0 <- control_law(co2, rep(clean, 240), 1000, 0, 5, 0.5, 3)
  chk(r0$n_switch > 10, paste("no-hysteresis control should chatter; got", r0$n_switch))

  cat("[G] outdoor-quality gate closes the window\n")
  # demand is on (high CO2) but outdoor air is dirty and > tau_AQ -> stays shut
  r <- control_law(rep(1100, 60), rep(dirty, 60), 1000, 100, 5, 0.5, 3)
  chk(all(r$s == 0L), "gate should keep window shut under dirty outdoor air")
  # gate exactly at threshold (AQI == tau_AQ) -> allowed open (<=)
  r <- control_law(rep(1100, 60), rep(3, 60), 1000, 100, 5, 0.5, 3)
  chk(all(r$s[-1] == 1L), "AQI == tau_AQ should still permit opening")

  cat("[O] health override prevails over the gate\n")
  # CO2 above ceiling (1600) with dirtiest outdoor air -> forced open
  r <- control_law(rep(1600, 60), rep(dirty, 60), 1000, 100, 5, 0.5, 3)
  chk(all(r$s[-1] == 1L), "health override should force opening despite dirty air")

  cat("[D] actuation slot limits switching frequency\n")
  # CO2 clears the hysteresis band both ways every step (would toggle each step)
  co2 <- rep(c(1200, 800), length.out = 120)      # 800 < 900 (close), 1200 > 1000 (open)
  r5  <- control_law(co2, rep(clean, 120), 1000, 100,  5, 0.5, 3)   # decide every step
  r30 <- control_law(co2, rep(clean, 120), 1000, 100, 30, 0.5, 3)   # decide every 6 steps
  chk(r30$slot_steps == 6L, "d_slot=30min should give 6 steps/slot")
  chk(r30$n_switch < r5$n_switch, paste("slot should reduce switching:",
                                        r30$n_switch, "vs", r5$n_switch))
  chk(r30$n_switch <= ceiling(120 / 6), "switches must not exceed number of slots")

  cat("[I] s(t) feeds a coherent g4 through indicators.R\n")
  if (exists("compute_indicators", mode = "function")) {
    co2 <- rep(c(1200, 800), length.out = 60)
    r <- control_law(co2, rep(clean, 60), 1000, 100, 5, 0.5, 3)
    df <- data.frame(day = 1, theta_op = 22, C_in = 420, C_out = 420,
                     theta_out = 14, s = r$s, occ = TRUE)
    g4 <- compute_indicators(df)$g4
    chk(g4 == r$n_switch, paste("g4 =", g4, "should equal control n_switch =", r$n_switch))
    cat(sprintf("      slot=5min -> %d switches over the session (N_max=60)\n", r$n_switch))
  } else cat("      (indicators.R not loaded; g4 integration skipped)\n")

  cat(sprintf("\nRESULT: %d passed, %d failed\n", e$ok, e$fail))
  if (e$fail > 0) stop("control_law validation failed")
  invisible(e$ok)
}

if (sys.nframe() == 0L) run_tests_control()

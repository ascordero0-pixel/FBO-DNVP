# =============================================================================
# trackB_massbalance.R  --  CO2 transient mass-balance engine (Track B)
# Thesis T3 -- FBO framework
# -----------------------------------------------------------------------------
# Implements the EnergyPlus Carbon Dioxide Predictor-Corrector used in T1
# (Appl. Sci. 2024, 14, 4651, Eq. 2 / Appendix A), single well-mixed zone:
#
#     rho * V * dC/dt = G(t) + m_inf(t) * (C_ext - C(t))
#
# with C in ppm, integrated by BACKWARD EULER (implicit, unconditionally
# stable; Courant number reaches ~0.77-9 under gale conditions, so explicit
# Euler is unsafe):
#
#     C_{k+1} = [ C_k + (dt/V)*(Gvol_k + Vext_k*C_ext) ] / [ 1 + (dt/V)*Vext_k ]
#
# KEY DESIGN POINT (corrects earlier error): the DNVP is a SWITCHING controller.
# The WWeqR column already contains the human-executed switching between the
# base opening and the 2.798% flush (confirmed by data audit: flushes coincide
# with CO2 >= 800 ppm, with realistic human-execution lag). Therefore the engine
# CONSUMES WWeqR(t) directly as the effective opening forcing -- it does NOT
# re-implement any control logic. The human-execution margin of error is thereby
# automatically embedded in the driving data.
#
# Ventilation flow (cross-ventilation, T1 model, Eq. 1 for A_eq):
#     A_eq(t) = WWeqR(t) * A_WALL           (T1 single-wall definition)
#     Q_wind  = C_d * A_eq * u_site * sqrt(dCp)
#     Q_stack = C_d * A_eq * sqrt(2 g dTheta H / Tbar)
#     Q_vent  = sqrt(Q_wind^2 + Q_stack^2)
#     u_site  = WIND_FACTOR * min(u_10, WIND_CAP)   (EN 16798-7 log profile,
#                                                    urban, z~11 m; winsor 20 m/s)
# =============================================================================

# ---- FIXED PARAMETERS (with provenance) ------------------------------------
V         <- 126.7031      # m3   zone volume            [Track A / IDF v3]
Q_INF     <- 0.64          # ren/h infiltration           [Track A, fixed interface]
V_INF     <- Q_INF * V / 3600            # m3/s = 0.0225169
DT        <- 300.0         # s    timestep (5 min)
CEXT      <- 420.0         # ppm  outdoor reference        [default; not measured]
A_WALL    <- 15.67         # m2   single-wall ref          [T1 WWeqR denominator]
H_STACK   <- 0.40          # m    N/S opening mid-height diff [window elevations]
DELTA_CP  <- 0.75          # -    windward-leeward Cp diff  [EN 16798-7 B.3.3.4]
G_GRAV    <- 9.81          # m/s2
CD_FIX    <- 0.65          # -    discharge coefficient     [T1, FIXED - see identifiability]
WIND_FACTOR <- 0.52141     # -    u_site/u_10 log profile   [EN 16798-7, urban, z=11 m]
WIND_CAP    <- 20.0        # m/s  winsor threshold on u_10  [Decision 2b]

# ---- ventilation flow -------------------------------------------------------
u_site <- function(u10) {
  if (is.na(u10)) u10 <- 0
  WIND_FACTOR * min(u10, WIND_CAP)
}

vent_flow <- function(wweqr, u10, Ti, Ot, Cd = CD_FIX) {
  if (is.na(wweqr) || wweqr <= 0) return(0.0)
  A_eq   <- wweqr * A_WALL
  U      <- u_site(u10)
  dTheta <- abs(Ti - Ot)
  Tbar   <- 0.5 * (Ti + Ot) + 273.15
  Q_wind  <- Cd * A_eq * U * sqrt(DELTA_CP)
  Q_stack <- if (dTheta > 0) Cd * A_eq * sqrt(2 * G_GRAV * dTheta * H_STACK / Tbar) else 0.0
  sqrt(Q_wind^2 + Q_stack^2)
}

# ---- forward integration over one day's anchor window ----------------------
# steps: data frame with columns Ti, Ot, WindS, CO2i, Occ, WWeqR (chronological)
simulate_day <- function(steps, g_person, Cd = CD_FIX) {
  n <- nrow(steps)
  Csim <- numeric(n)
  # initial condition = first measured CO2i in the window (fallback C_ext)
  idx0 <- which(!is.na(steps$CO2i))[1]
  C <- if (length(idx0)) steps$CO2i[idx0] else CEXT
  for (k in seq_len(n)) {
    s <- steps[k, ]
    Vext <- vent_flow(s$WWeqR, s$WindS, s$Ti, s$Ot, Cd) + V_INF
    occ  <- if (is.na(s$Occ)) 0 else s$Occ
    Gvol <- occ * g_person * 1e6          # m3 * ppm / s
    C <- (C + (DT / V) * (Gvol + Vext * CEXT)) / (1 + (DT / V) * Vext)
    Csim[k] <- C
  }
  data.frame(ts = steps$ts, Csim = Csim, Cmeas = steps$CO2i, WWeqR = steps$WWeqR)
}

# ---- goodness-of-fit metrics (ASHRAE Guideline 14) -------------------------
gof <- function(sim_df) {
  ok <- !is.na(sim_df$Cmeas)
  c_ <- sim_df$Csim[ok]; m_ <- sim_df$Cmeas[ok]
  n <- length(m_); mbar <- mean(m_)
  nmbe   <- sum(c_ - m_) / (n * mbar) * 100
  rmse   <- sqrt(mean((c_ - m_)^2))
  cvrmse <- rmse / mbar * 100
  list(n = n, mean_meas = mbar, NMBE = nmbe, CVRMSE = cvrmse, RMSE = rmse)
}

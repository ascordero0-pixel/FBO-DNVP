# ============================================================================
# g5.R  --  T3 FBO framework, Phase 0 (Tanda 2): outdoor-pollution-import indicator.
#
# g5 = P95 over occupied hours of the worst indoor air-quality sub-index among
# PM2.5, NO2 and O3, from a per-species single-compartment dynamic balance:
#
#     dC_in,s/dt = (Q/V)*(P_s*C_out,s - C_in,s) - k_s*C_in,s
#
# solved by IMPLICIT EULER (unconditionally stable, positivity-safe):
#     C_in[t+1] = ( C_in[t] + dt*a[t+1]*P*C_out[t+1] ) / ( 1 + dt*(a[t+1]+k) ),
#     a(t) = Q(t)/V  [h^-1]  (total air exchange: ventilation + infiltration)
#
# Open-window natural ventilation => P ~ 1; species differ in loss rate k_s
# (O3 >> PM2.5 > NO2; NO2 is the least-removed / highest I/O species).
#
# Sub-index mapping: SPANISH NATIONAL ICA, Orden TEC/351/2019 (Table 1, in force).
# The sub-index is applied to the pollutant's REFERENCE-TIME RUNNING MEAN, not the
# instantaneous 5-min value:  PM2.5 -> 24 h,  O3 -> 8 h,  NO2 -> 1 h.
# Worst-sub-index rule = official ICA methodology. Equivalent to current EAQI (EEA 2017).
#
# Anchored (assumed, not calibrated) loss rates:
#   PM2.5: k=0.20 h^-1 (Chen & Zhao 2011); NO2: v_d=0.002 cm/s (Grontoft&Raychaudhuri 2004);
#   O3: v_d=0.04 cm/s (Weschler 2000). A/V ~ 1.31 m^-1 (V=126.70 m3, A~166 m2 -- confirm IDF).
# INJECTABLE: C_out,s(t) (Torneo 2024 hourly, held to 5-min) and the ICA breakpoints.
# ============================================================================


# ---------------------------------------------------------------------------
# 1. Species parameters and ICA sub-index breakpoints
# ---------------------------------------------------------------------------
k_from_vd <- function(v_d_cm_s, A_over_V) v_d_cm_s * 0.01 * 3600 * A_over_V   # cm/s -> h^-1

default_species_params <- function(A_over_V = 1.31) {
  list(
    PM25 = list(k = 0.20,                       P = 1.0, ref_h = 24),  # 24 h running mean
    NO2  = list(k = k_from_vd(0.002, A_over_V), P = 1.0, ref_h = 1),   # 1 h
    O3   = list(k = k_from_vd(0.040, A_over_V), P = 1.0, ref_h = 8)    # 8 h
  )
}

# Orden TEC/351/2019, Tabla 1 (ug/m3): upper bound of bands 1..5; band = findInterval+1 in 1..6
ica_2019_breaks <- function() {
  list(
    PM25 = c(10, 20, 25, 50, 75),      # 24 h
    O3   = c(50, 100, 130, 240, 380),  # 8 h
    NO2  = c(40, 90, 120, 230, 340)    # 1 h
  )
}
default_subindex_breaks <- ica_2019_breaks


# ---------------------------------------------------------------------------
# 2. Numerics: implicit-Euler balance, reference-time running mean, sub-index
# ---------------------------------------------------------------------------
solve_species <- function(C_out, a, k, P = 1.0, dt_h = 5/60, C0 = NULL) {
  n <- length(C_out); a <- rep_len(a, n); C <- numeric(n)
  C[1] <- if (is.null(C0)) C_out[1] * a[1] * P / (a[1] + k + 1e-12) else C0
  for (t in 2:n)
    C[t] <- (C[t-1] + dt_h * a[t] * P * C_out[t]) / (1 + dt_h * (a[t] + k))
  C
}

# trailing running mean over w steps (expanding window at the start; no NA)
running_mean <- function(x, w) {
  n <- length(x); if (w <= 1L) return(as.numeric(x))
  cs <- c(0, cumsum(x)); idx <- seq_len(n); lo <- pmax(1L, idx - w + 1L)
  (cs[idx + 1L] - cs[lo]) / (idx - lo + 1L)
}

subindex <- function(conc, breaks) findInterval(conc, breaks) + 1L        # -> 1..6


# ---------------------------------------------------------------------------
# 3. Main: compute g5
# ---------------------------------------------------------------------------
# df columns: Q (m3/h), Cout_PM25, Cout_NO2, Cout_O3, occ (logical)
compute_g5 <- function(df, V = 126.70, params = default_species_params(),
                       breaks = default_subindex_breaks(), dt_h = 5/60, p95_type = 7) {
  stopifnot(all(c("Q","Cout_PM25","Cout_NO2","Cout_O3","occ") %in% names(df)))
  a <- df$Q / V; occ <- as.logical(df$occ)
  wsteps <- function(ref_h) max(1L, as.integer(round(ref_h / dt_h)))

  Cin <- list(
    PM25 = solve_species(df$Cout_PM25, a, params$PM25$k, params$PM25$P, dt_h),
    NO2  = solve_species(df$Cout_NO2,  a, params$NO2$k,  params$NO2$P,  dt_h),
    O3   = solve_species(df$Cout_O3,   a, params$O3$k,   params$O3$P,   dt_h)
  )
  # reference-time running means, then sub-index per species
  Cavg <- list(
    PM25 = running_mean(Cin$PM25, wsteps(params$PM25$ref_h)),
    NO2  = running_mean(Cin$NO2,  wsteps(params$NO2$ref_h)),
    O3   = running_mean(Cin$O3,   wsteps(params$O3$ref_h))
  )
  idx <- cbind(subindex(Cavg$PM25, breaks$PM25),
               subindex(Cavg$NO2,  breaks$NO2),
               subindex(Cavg$O3,   breaks$O3))
  worst <- apply(idx, 1, max)
  g5 <- unname(stats::quantile(worst[occ], 0.95, type = p95_type))
  list(g5 = g5, Cin = Cin, Cavg = Cavg, worst = worst)
}


# ---------------------------------------------------------------------------
# 4. Validation harness
# ---------------------------------------------------------------------------
run_tests_g5 <- function() {
  e <- new.env(); e$ok <- 0L; e$fail <- 0L
  chk <- function(cond, msg) { if (isTRUE(cond)) e$ok <- e$ok + 1L
                               else { e$fail <- e$fail + 1L; cat("  FAIL:", msg, "\n") } }
  apx <- function(x, y, tol = 1e-2) abs(x - y) < tol

  cat("[S] steady-state  C_ss = C_out*a*P/(a+k)\n")
  chk(apx(tail(solve_species(rep(100, 5000), 2, 1, 1), 1), 100*2/3), "steady state wrong")
  chk(apx(tail(solve_species(rep(100, 5000), 2, 1, 0.5), 1), 0.5*100*2/3), "penetration P not applied")

  cat("[H] loss hierarchy  O3 >> PM2.5 > NO2\n")
  pr <- default_species_params(); css <- function(k) 100*3/(3+k)
  chk(css(pr$O3$k) < css(pr$PM25$k) && css(pr$PM25$k) < css(pr$NO2$k), "steady-state ordering wrong")
  chk(pr$O3$k > pr$PM25$k && pr$PM25$k > pr$NO2$k, "k ordering wrong")

  cat("[V] ventilation increases import\n")
  chk(tail(solve_species(rep(100,5000), 6, 0.2),1) > tail(solve_species(rep(100,5000), 0.5, 0.2),1),
      "import not increasing with ventilation")
  chk(tail(solve_species(rep(100,3000), 0, 0.5),1) < 1e-6, "closed (a=0) should decay to ~0")

  cat("[R] reference-time running mean\n")
  chk(all(abs(running_mean(c(0,0,0,12), 4) - c(0,0,0,3)) < 1e-9), "running_mean window wrong")
  chk(apx(running_mean(rep(10,300), 288)[300], 10), "running_mean of constant != constant")

  cat("[X] sub-index mapping (ICA 2019) + worst rule\n")
  br <- ica_2019_breaks()
  chk(subindex(30, br$PM25) == 4L, "PM2.5=30 -> band 4")
  chk(subindex(5,  br$PM25) == 1L, "PM2.5=5 -> band 1")
  chk(subindex(500, br$O3) == 6L,  "O3=500 -> band 6")
  chk(max(subindex(5,br$PM25), subindex(300,br$NO2), subindex(5,br$O3)) == 5L, "worst rule wrong")

  cat("[g5] worst sub-index over occupied hours (with reference-time averaging)\n")
  n <- 700                                   # ~58 h: fills 24 h PM window, reaches steady state
  base <- data.frame(Q = 200, Cout_PM25 = 40, Cout_NO2 = 30, Cout_O3 = 40, occ = TRUE)[rep(1,n),]
  r <- compute_g5(base)
  chk(r$g5 == 4, paste("steady-state g5 =", r$g5, "expected 4 (PM2.5 ~35.5 -> band 4)"))
  # a 30-min PM2.5 dust spike must be SMOOTHED by the 24 h window (not blow up to band 6)
  clean <- data.frame(Q = 200, Cout_PM25 = 5, Cout_NO2 = 10, Cout_O3 = 20, occ = TRUE)[rep(1,n),]
  clean$Cout_PM25[400:405] <- 800
  rs <- compute_g5(clean)
  chk(rs$g5 <= 2, paste("24h averaging should smooth a 30-min PM spike; g5 =", rs$g5))
  # occupancy masking: a non-occupied tail spike must not enter g5
  md <- base; md$occ <- c(rep(TRUE, n-6), rep(FALSE, 6)); md$Cout_PM25[(n-5):n] <- 999
  chk(compute_g5(md)$g5 == 4, "masking failed")

  cat(sprintf("\n  example: steady g5 = %.0f | smoothed-spike g5 = %.0f\n", r$g5, rs$g5))
  cat(sprintf("RESULT: %d passed, %d failed\n", e$ok, e$fail))
  if (e$fail > 0) stop("g5 validation failed")
  invisible(e$ok)
}

if (sys.nframe() == 0L) run_tests_g5()

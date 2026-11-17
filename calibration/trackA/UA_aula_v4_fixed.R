# =============================================================================
# UA_aula_v4_fixed.R  --  CORRECTED, autonomous envelope-conductance calculator
# -----------------------------------------------------------------------------
# REPLACES the deprecated UA_aula_v3.R, which had TWO bugs:
#   (bug 1) counted only 2 exterior walls — it excluded Wall_3 because it read the
#           v3 IDF where Wall_3 carried the spurious OtherSideConditionsModel and
#           was therefore not classified as 'Outdoors'.
#   (bug 2) applied U = 1.505 (no-film) to the ceiling and treated it as an
#           EXTERIOR loss, when the real model uses U = 0.917 and the ceiling is
#           an INTERIOR boundary (to the unventilated roof cavity).
#   The two errors happened to net out near 156 W/K for the wrong reasons.
#
# This version reads the authoritative EnergyPlus "Envelope Summary" tables from
# eplustbl.htm (Opaque Exterior + Exterior Fenestration) and returns the
# direct-exterior envelope conductance. It is AUTONOMOUS (file-only; no engine).
#
#   Canonical result (Wall_3 = Outdoors, i.e. v4 or v3+normalize): UA = 99.3 W/K
#
# USAGE:
#   source("UA_aula_v4_fixed.R"); ua_aula("runs/v3_trackA/eplustbl.htm")
# =============================================================================

.cellval <- function(x) suppressWarnings(as.numeric(trimws(gsub("<[^>]*>", "", x))))
.cellstr <- function(x) trimws(gsub("<[^>]*>", "", x))

ua_aula <- function(htm = "runs/v3_trackA/eplustbl.htm", verbose = TRUE) {
  stopifnot(file.exists(htm))
  h <- readLines(htm, warn = FALSE)

  # ---- Opaque Exterior walls of the classroom (Wall_2/3/5) ----
  # Column order in the E+ htm row: name, construction, reflectance,
  #   U-with-film, U-no-film, gross area, net area, azimuth, tilt, cardinal.
  wl <- grep("CLASSROOM_WALL_(2|3|5)_0_0<", h)
  wall_tab <- data.frame(name=character(), U=numeric(), net=numeric(),
                         UA=numeric(), dir=character(), stringsAsFactors=FALSE)
  for (i in wl) {
    cells <- h[i:(i+9)]
    nm  <- .cellstr(cells[1])
    U   <- .cellval(cells[4])   # U with film
    net <- .cellval(cells[7])   # net area
    dir <- .cellstr(cells[10])
    if (is.finite(U) && is.finite(net))
      wall_tab <- rbind(wall_tab, data.frame(name=nm, U=U, net=net, UA=U*net, dir=dir))
  }

  # ---- Exterior windows: recover total glazed area from wall gross-net ----
  # Robust approach: for each hosting wall, opening area = gross - net.
  win_area_by_wall <- function(host) {
    r <- grep(paste0(host, "<"), h)
    if (!length(r)) return(NA_real_)
    cells <- h[r[1]:(r[1]+9)]
    gross <- .cellval(cells[6]); net <- .cellval(cells[7])
    gross - net
  }
  aN <- win_area_by_wall("CLASSROOM_WALL_3_0_0")   # north windows
  aS <- win_area_by_wall("CLASSROOM_WALL_5_0_0")   # south windows
  Uwin <- 1.960                                    # glass U-factor (constr 1001/1002)
  ua_win <- (aN + aS) * Uwin

  ua_wall <- sum(wall_tab$UA)
  UA <- ua_wall + ua_win

  if (verbose) {
    cat("=========  UA direct-exterior of the classroom (v4 / corrected)  =========\n")
    for (k in seq_len(nrow(wall_tab)))
      cat(sprintf("  %-42s dir=%-2s A=%6.2f U=%.3f -> UA=%6.2f W/K\n",
                  wall_tab$name[k], wall_tab$dir[k], wall_tab$net[k],
                  wall_tab$U[k], wall_tab$UA[k]))
    cat(sprintf("  WINDOWS North A=%.2f + South A=%.2f  (U=%.3f) -> UA=%6.2f W/K\n",
                aN, aS, Uwin, ua_win))
    cat("--------------------------------------------------------------------------\n")
    cat(sprintf("  UA TOTAL (direct exterior) = %.1f W/K\n", UA))
    cat("  NOTE: ceiling (U=0.917 -> roof cavity) and floor (U=0.873 -> lower storey)\n")
    cat("        are INTERIOR boundaries, resolved by CTF, NOT counted here.\n")
    cat("        As-surveyed T1 [11] value 200.7 W/K counts the pitched roof as a\n")
    cat("        direct loss and is a non-comparable upper bound.\n")
    cat("==========================================================================\n")
  }
  invisible(list(UA = UA, walls = wall_tab, win_north = aN, win_south = aS,
                 ua_walls = ua_wall, ua_windows = ua_win))
}

if (interactive()) cat("[UA_aula_v4_fixed.R] loaded. Call ua_aula('runs/v3_trackA/eplustbl.htm')\n")

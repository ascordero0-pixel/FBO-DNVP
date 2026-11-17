# =============================================================================
# load_co2_campaign.R  --  Track B (CO2 mass-balance) data loader
# Thesis T3 -- FBO framework -- 1946 uninsulated corner classroom, Seville
# -----------------------------------------------------------------------------
# Mirrors the Track A architecture: builds the campaign data frame from the
# audited co2_clean.csv (exported via openpyxl to avoid the Excel serial-number
# bug) and performs AUTOMATIC ANCHOR VERIFICATION with fail-fast on mismatch.
#
# Anchor-first rule (non-negotiable, per project methodology):
#   No derived analysis proceeds until integrity anchors are reproduced.
#   Track B anchors verified here:
#     (a) 863 total 5-min rows over 3 days (Nov 18, Nov 26, Dec 17 2021)
#     (b) 183 rows inside the anchor window 09:00-14:00
#     (c) per-day occupancy peaks 26 / 26 / 21
#     (d) ventilation labels: Nov18 DNVP, Nov26 DNVP, Dec17 SNVP
#     (e) WWeqR flush level 2.798% present on DNVP days (switching already baked in)
#     (f) C_ext constant 420 ppm (default, not measured)
#
# No-fabrication policy: this script only reads verified data; it computes no
# calibration result. Parameters carry explicit provenance labels.
# =============================================================================

# base R only; co2_clean.csv already audited upstream (openpyxl export)

# ---- configuration ----------------------------------------------------------
CO2_CSV <- "co2_clean.csv"

# Anchor window (Track B anchor definition, per user): 09:00-14:00 inclusive
ANCHOR_START_SEC <- 9 * 3600
ANCHOR_END_SEC   <- 14 * 3600

# Expected integrity anchors
EXP_TOTAL_ROWS   <- 863L
EXP_ANCHOR_ROWS  <- 183L
EXP_OCC_PEAK     <- c("18" = 26, "26" = 26, "17" = 21)
EXP_VENT         <- c("18" = "DNVP", "26" = "DNVP", "17" = "SNVP")
EXP_FLUSH_WWEQR  <- 0.027984            # 2.798% flush level (DNVP days)
EXP_CEXT         <- 420.0               # ppm, default outdoor reference
TOL              <- 1e-6

# ---- load -------------------------------------------------------------------
load_co2 <- function(path = CO2_CSV) {
  if (!file.exists(path))
    stop(sprintf("[FAIL] co2_clean.csv not found at '%s'. Run the openpyxl export first.", path))

  d <- read.csv(path, stringsAsFactors = FALSE)

  # parse timestamp -> seconds-of-day for anchor masking
  ts <- as.POSIXct(d$datetime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
  d$sod <- as.integer(format(ts, "%H")) * 3600L +
           as.integer(format(ts, "%M")) * 60L +
           as.integer(format(ts, "%S"))
  d$ts  <- ts
  d$day <- as.character(d$day)

  # numeric coercions (empty strings -> NA)
  num <- function(x) suppressWarnings(as.numeric(x))
  for (col in c("Ti","Ot","WindS","CO2i","CO2o","Occ","WWeqR"))
    if (col %in% names(d)) d[[col]] <- num(d[[col]])

  d
}

# ---- anchor verification (fail-fast) ---------------------------------------
verify_anchors <- function(d) {
  msg <- function(ok, txt) cat(sprintf("  [%s] %s\n", ifelse(ok,"OK","FAIL"), txt))
  fail <- FALSE

  # (a) total rows
  ok <- nrow(d) == EXP_TOTAL_ROWS
  msg(ok, sprintf("total rows = %d (expected %d)", nrow(d), EXP_TOTAL_ROWS)); fail <- fail || !ok

  # (b) anchor-window rows
  anc <- d[d$sod >= ANCHOR_START_SEC & d$sod <= ANCHOR_END_SEC, ]
  ok <- nrow(anc) == EXP_ANCHOR_ROWS
  msg(ok, sprintf("anchor-window rows = %d (expected %d)", nrow(anc), EXP_ANCHOR_ROWS)); fail <- fail || !ok

  # (c) per-day occupancy peaks
  for (dd in names(EXP_OCC_PEAK)) {
    pk <- suppressWarnings(max(d$Occ[d$day == dd], na.rm = TRUE))
    ok <- abs(pk - EXP_OCC_PEAK[[dd]]) < TOL
    msg(ok, sprintf("day %s occupancy peak = %g (expected %g)", dd, pk, EXP_OCC_PEAK[[dd]]))
    fail <- fail || !ok
  }

  # (d) ventilation labels
  for (dd in names(EXP_VENT)) {
    labs <- unique(d$Vent[d$day == dd & d$Vent %in% c("DNVP","SNVP")])
    ok <- EXP_VENT[[dd]] %in% labs
    msg(ok, sprintf("day %s ventilation contains '%s' (found: %s)",
                    dd, EXP_VENT[[dd]], paste(labs, collapse=",")))
    fail <- fail || !ok
  }

  # (e) flush WWeqR present on DNVP days (switching already in the column)
  for (dd in c("18","26")) {
    has_flush <- any(abs(d$WWeqR[d$day == dd] - EXP_FLUSH_WWEQR) < 1e-4, na.rm = TRUE)
    msg(has_flush, sprintf("day %s contains flush WWeqR=2.798%% (DNVP switching baked in)", dd))
    fail <- fail || !has_flush
  }

  # (f) C_ext constant 420 ppm
  cext_vals <- unique(d$CO2o[!is.na(d$CO2o)])
  ok <- length(cext_vals) == 1 && abs(cext_vals - EXP_CEXT) < TOL
  msg(ok, sprintf("C_ext constant = %s ppm (expected %g)",
                  paste(cext_vals, collapse=","), EXP_CEXT)); fail <- fail || !ok

  if (fail) stop("[ANCHOR FAILURE] Track B data integrity anchors not reproduced. Halting (anchor-first rule).")
  cat("  --> ALL TRACK B ANCHORS REPRODUCED. Safe to proceed.\n")
  invisible(TRUE)
}

# ---- build analysis objects -------------------------------------------------
build_trackB <- function(path = CO2_CSV, verbose = TRUE) {
  d <- load_co2(path)
  if (verbose) cat("Track B data loaded. Verifying anchors:\n")
  verify_anchors(d)

  anc <- d[d$sod >= ANCHOR_START_SEC & d$sod <= ANCHOR_END_SEC, ]
  days <- split(anc, anc$day)

  list(all = d, anchor = anc, days = days)
}

# ---- run when sourced directly ---------------------------------------------
if (sys.nframe() == 0) {
  tb <- build_trackB()
  cat(sprintf("\nData frames ready: full=%d rows, anchor=%d rows, days={%s}\n",
              nrow(tb$all), nrow(tb$anchor), paste(names(tb$days), collapse=", ")))
}

# =============================================================================
# load_calibrated_series.R
# -----------------------------------------------------------------------------
# Single source of truth for the deviation scripts. Builds, from the CALIBRATED
# Track A EnergyPlus runs (NOT the uncalibrated "Option A" series), the objects
#   meas   : measured 5-min series      (k, T_meas, T_out, heating)
#   df     : simulated 5-min series     (k, T_out, v3, v6)
#   D      : full merge over RunPeriod  (k, T_out, v3, v6, T_meas, win)
#   D_win  : weekend windows only       (Fri 22:00 -> Mon 07:00, >=80% coverage)
#
# Calibrated sources:
#   v3 -> runs/v3_trackA/eplusout.csv  (reproduces the canonical -5.41%/8.04%)
#   v6 -> runs/v6_trackA/eplusout.csv  (reproduces +6.05%/9.21% nocturnal)
#   measured -> meas_clean.csv         (openpyxl export, To mean 13.97 = EPW)
#
# ANCHOR PRINCIPLE: after loading, the script re-derives the nocturnal NMBE of
# each model and checks it against its reference metric. If a series does not
# reproduce its anchor, a warning is raised - no derived analysis should ever
# be run on a simulated series that has not first reproduced a validated metric
# (this is the safeguard that exposed the two artefacts earlier in the project).
# =============================================================================

suppressWarnings({library(lubridate)})

# --- EnergyPlus Date/Time parser (no year; 24:00:00 = midnight next day) -----
.parse_eplus_dt <- function(x, year = 2021){
  x <- trimws(x)
  p  <- regmatches(x, regexec("(\\d{2})/(\\d{2})\\s+(\\d{2}):(\\d{2}):(\\d{2})", x))
  g  <- function(i) as.integer(vapply(p, function(z) z[i], character(1)))
  m  <- g(2); d <- g(3); hh <- g(4); mi <- g(5)
  roll <- hh == 24
  base <- as.POSIXct(sprintf("%04d-%02d-%02d 00:00:00", year, m, d), tz = "UTC")
  t <- base + (hh %% 24)*3600 + mi*60
  t[roll] <- t[roll] + 24*3600
  t
}

# --- read the study-zone mean air temperature from an eplusout.csv -----------
.read_zone_temp <- function(path){
  stopifnot(file.exists(path))
  ec <- read.csv(path, check.names = FALSE)
  tcol <- grep("ZXCASESTUDYCLASSROOM:CLASSROOM.*Zone Mean Air Temperature",
               names(ec), ignore.case = TRUE, value = TRUE)[1]
  if(is.na(tcol)) stop("Study-zone Zone Mean Air Temperature column not found in ", path)
  data.frame(k = round_date(.parse_eplus_dt(ec[["Date/Time"]]), "5 mins"),
             Tz = as.numeric(ec[[tcol]]))
}

# --- (1) measured clean series ----------------------------------------------
stopifnot(file.exists("meas_clean.csv"))
meas <- read.csv("meas_clean.csv", stringsAsFactors = FALSE)
meas$k <- round_date(as.POSIXct(meas$k, tz = "UTC"), "5 mins")
meas <- meas[is.finite(meas$T_meas) & !is.na(meas$k), ]

# --- (2) calibrated simulated series (v3, v6) -------------------------------
s3 <- .read_zone_temp("runs/v3_trackA/eplusout.csv"); names(s3)[2] <- "v3"
s6 <- .read_zone_temp("runs/v6_trackA/eplusout.csv"); names(s6)[2] <- "v6"

# --- (3) df = (k, T_out, v3, v6); T_out from the measured series (= EPW) -----
df <- merge(s3, s6, by = "k")
df <- merge(df, meas[, c("k","T_out")], by = "k")
df <- df[order(df$k), c("k","T_out","v3","v6")]

# --- (4) D = full merge with measured indoor; add weekend window id ----------
D <- merge(df, meas[, c("k","T_meas")], by = "k"); D <- D[order(D$k), ]
wd <- wday(D$k); hh <- hour(D$k) + minute(D$k)/60
in_win <- (wd==6 & hh>=22) | (wd==7) | (wd==1) | (wd==2 & hh<7)
grp <- cumsum(c(TRUE, diff(as.integer(in_win))==1))
D$win <- ifelse(in_win, grp, NA)
D_win <- D[!is.na(D$win), ]
# keep windows covering >=80% of the nominal 57 h span
NOMINAL_H <- 57
span_h <- tapply(as.numeric(D_win$k), D_win$win, function(t) (max(t)-min(t))/3600)
good   <- names(span_h)[span_h >= 0.80 * NOMINAL_H]
D_win  <- D_win[as.character(D_win$win) %in% good, ]

# --- (5) ANCHOR VERIFICATION: nocturnal NMBE must reproduce reference --------
noct <- (hh >= 22) | (hh < 7)
Dn <- D[noct, ]
nmbe <- function(col) 100 * mean(Dn[[col]] - Dn$T_meas) / mean(Dn$T_meas)
cv   <- function(col) 100 * sqrt(mean((Dn[[col]] - Dn$T_meas)^2)) / mean(Dn$T_meas)
a3 <- nmbe("v3"); a6 <- nmbe("v6")
cat(sprintf("[load_calibrated] merged N=%d | nocturnal n=%d | weekend windows=%d\n",
            nrow(D), nrow(Dn), length(unique(D_win$win))))
cat(sprintf("  v3 nocturnal: NMBE %+.2f%% (ref -5.41%%)  CV-RMSE %.2f%% (ref 8.04%%)\n", a3, cv("v3")))
cat(sprintf("  v6 nocturnal: NMBE %+.2f%% (ref +6.05%%)  CV-RMSE %.2f%% (ref 9.21%%)\n", a6, cv("v6")))
tol <- 0.15
if(abs(a3 - (-5.41)) > tol) warning("v3 series does NOT reproduce its -5.41% anchor (got ",
                                     sprintf("%.2f%%", a3), ") - check runs/v3_trackA/eplusout.csv")
if(abs(a6 - ( 6.05)) > 0.25) warning("v6 series does NOT reproduce its +6.05% anchor (got ",
                                     sprintf("%.2f%%", a6), ") - check runs/v6_trackA/eplusout.csv")
if(abs(a3 - (-5.41)) <= tol && abs(a6 - 6.05) <= 0.25)
  cat("  [OK] both series reproduce their validated anchors -> safe to use.\n")

.CALIB_SERIES_LOADED <- TRUE

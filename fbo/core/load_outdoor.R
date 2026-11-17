# ============================================================================
# load_outdoor.R  --  T3 FBO framework, Phase 0: outdoor pollutant-series loader.
#
# Parses the Junta de Andalucia hourly CSVs (wide: one row per day, columns
# H01..H24) for station 41091009 (Torneo, Seville), year 2024, for PM2.5/NO2/O3,
# builds a continuous full-year HOURLY series, gap-fills it, and expands it to the
# 5-min grid consumed by g5.R (columns Cout_PM25, Cout_NO2, Cout_O3).
#
# Assumptions (documented):
#   * H01..H24 -> hours 0..23 of the day (H01 = [00:00,01:00)).
#   * timestamps handled in UTC (the network reports fixed 24-h days, no DST); a
#     constant offset vs the EPW local time is immaterial here because the pollutant
#     series is a decoupled representative boundary condition (co-timing not preserved).
#   * gap-fill: interior gaps <= maxgap_h are linearly interpolated; any remaining
#     gap (incl. whole-day and edge gaps) is filled with a month-by-hour climatology.
#   * PM2.5 zeros: kept as-is by default; set zeros_to_na=TRUE to treat them as gaps
#     (pending confirmation of the Junta missing-data convention).
#
# Base R only. Run `Rscript load_outdoor.R` for a report + self-checks on the files.
# ============================================================================

MAGNITUDE <- c(PM25 = 9L, NO2 = 8L, O3 = 14L)   # Junta magnitud codes


# ---------------------------------------------------------------------------
# 1. Parse one wide CSV -> long hourly data.frame(timestamp, value)
# ---------------------------------------------------------------------------
read_species_wide <- function(path, station = "41091009", year = 2024) {
  d <- read.csv2(path, stringsAsFactors = FALSE, fileEncoding = "latin1")
  d$stcode <- sprintf("%02d%03d%03d", d$PROVINCIA, d$MUNICIPIO, d$ESTACION)
  d <- d[d$stcode == station, , drop = FALSE]
  if (nrow(d) == 0) stop("No rows for station ", station, " in ", path)

  hcols <- sprintf("H%02d", 1:24)
  H <- as.matrix(d[, hcols]); suppressWarnings(storage.mode(H) <- "numeric")
  H[H < 0] <- NA                                   # negatives = invalid

  dts <- as.Date(sprintf("%d-%02d-%02d", year, d$MES, d$DIA))
  ts  <- as.POSIXct(rep(dts, each = 24), tz = "UTC") + rep(0:23, times = nrow(d)) * 3600
  data.frame(timestamp = ts, value = as.vector(t(H)))   # row-major: day1 H01..H24, ...
}

# full-year hourly index and left-join
align_full_year <- function(df, year = 2024) {
  idx <- seq(as.POSIXct(sprintf("%d-01-01 00:00", year), tz = "UTC"),
             as.POSIXct(sprintf("%d-12-31 23:00", year), tz = "UTC"), by = "hour")
  out <- data.frame(timestamp = idx)
  out$value <- df$value[match(out$timestamp, df$timestamp)]
  out
}


# ---------------------------------------------------------------------------
# 2. Gap-filling (linear for short gaps, month-hour climatology for the rest)
# ---------------------------------------------------------------------------
.na_runs <- function(is_na) {
  r <- rle(is_na); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  k <- which(r$values)
  data.frame(start = starts[k], end = ends[k], len = r$lengths[k])
}

fill_series <- function(v, ts, maxgap_h = 6) {
  n <- length(v); flag <- ifelse(is.na(v), "clim", "obs")   # provisional

  # (a) linear interpolation of interior gaps of length <= maxgap_h
  runs <- .na_runs(is.na(v))
  for (i in seq_len(nrow(runs))) {
    s <- runs$start[i]; e <- runs$end[i]
    if (runs$len[i] <= maxgap_h && s > 1 && e < n) {
      v[s:e] <- approx(c(s - 1, e + 1), c(v[s - 1], v[e + 1]), xout = s:e)$y
      flag[s:e] <- "interp"
    }
  }

  # (b) month-by-hour climatology for anything still missing
  need <- is.na(v)
  if (any(need)) {
    mo <- as.integer(format(ts, "%m")); hr <- as.integer(format(ts, "%H"))
    key <- mo * 100L + hr
    mh  <- tapply(v, key, mean, na.rm = TRUE)          # month-hour means
    hrm <- tapply(v, hr,  mean, na.rm = TRUE)          # hour means (fallback)
    fv  <- mh[as.character(key[need])]
    fb  <- hrm[as.character(hr[need])]
    fv[is.na(fv)] <- fb[is.na(fv)]
    v[need] <- fv                                       # flag stays "clim"
  }
  list(value = v, flag = flag)
}


# ---------------------------------------------------------------------------
# 3. Expand hourly -> 5-min grid (hold within the hour, or linear)
# ---------------------------------------------------------------------------
expand_5min <- function(ts_h, v_h, method = c("hold", "linear")) {
  method <- match.arg(method)
  n <- length(v_h)
  ts5 <- seq(ts_h[1], by = "5 min", length.out = n * 12L)
  if (method == "hold") {
    v5 <- rep(v_h, each = 12L)
  } else {
    v5 <- approx(x = as.numeric(ts_h), y = v_h, xout = as.numeric(ts5), rule = 2)$y
  }
  data.frame(timestamp = ts5, value = v5)
}


# ---------------------------------------------------------------------------
# 4. Main loader
# ---------------------------------------------------------------------------
load_outdoor <- function(paths, station = "41091009", year = 2024,
                         maxgap_h = 6, five_min = TRUE, method = "hold",
                         zeros_to_na = FALSE, verbose = TRUE) {
  # paths: named list/vector with PM25, NO2, O3 file paths
  species <- c("PM25", "NO2", "O3")
  hourly <- data.frame(timestamp = align_full_year(
    data.frame(timestamp = as.POSIXct(character()), value = numeric()), year)$timestamp)
  report <- list()

  for (sp in species) {
    raw <- read_species_wide(paths[[sp]], station, year)
    if (zeros_to_na && sp == "PM25") raw$value[raw$value == 0] <- NA
    al  <- align_full_year(raw, year)
    n_missing <- sum(is.na(al$value))
    longest   <- { r <- .na_runs(is.na(al$value)); if (nrow(r)) max(r$len) else 0L }
    fl  <- fill_series(al$value, al$timestamp, maxgap_h)
    hourly[[paste0("Cout_", sp)]] <- fl$value
    hourly[[paste0("flag_", sp)]] <- fl$flag
    report[[sp]] <- list(missing_hours = n_missing, longest_gap_h = longest,
                         pct_obs   = 100 * mean(fl$flag == "obs"),
                         pct_interp= 100 * mean(fl$flag == "interp"),
                         pct_clim  = 100 * mean(fl$flag == "clim"),
                         range = range(fl$value))
  }

  if (verbose) {
    cat("== load_outdoor report  (station", station, "| year", year, ") ==\n")
    for (sp in species) with(report[[sp]], cat(sprintf(
      "  %-5s: missing %4d h | longest gap %3d h | obs %.1f%% interp %.1f%% clim %.1f%% | range [%.1f, %.1f]\n",
      sp, missing_hours, longest_gap_h, pct_obs, pct_interp, pct_clim, range[1], range[2])))
  }

  if (!five_min) return(list(hourly = hourly, report = report))

  five <- data.frame(timestamp =
    seq(hourly$timestamp[1], by = "5 min", length.out = nrow(hourly) * 12L))
  for (sp in species)
    five[[paste0("Cout_", sp)]] <-
      expand_5min(hourly$timestamp, hourly[[paste0("Cout_", sp)]], method)$value

  list(hourly = hourly, fivemin = five, report = report)
}


# ---------------------------------------------------------------------------
# 5. Self-check on the real files
# ---------------------------------------------------------------------------
run_report_outdoor <- function(dir = ".") {
  paths <- list(
    PM25 = file.path(dir, "PM25_HH_2024_SEV_Torneo.csv"),
    NO2  = file.path(dir, "NO2_HH_2024_SEV_Torneo.csv"),
    O3   = file.path(dir, "O3_HH_2024_SEV_Torneo.csv"))
  res <- load_outdoor(paths, five_min = TRUE)

  e <- new.env(); e$ok <- 0L; e$fail <- 0L
  chk <- function(c, m) if (isTRUE(c)) e$ok <- e$ok+1L else { e$fail <- e$fail+1L; cat("  FAIL:", m, "\n") }

  cat("\n[checks]\n")
  chk(nrow(res$hourly) == 8784, paste("hourly length", nrow(res$hourly), "expected 8784 (leap year)"))
  chk(nrow(res$fivemin) == 8784 * 12, paste("5-min length", nrow(res$fivemin)))
  for (sp in c("PM25","NO2","O3")) {
    col <- paste0("Cout_", sp)
    chk(!any(is.na(res$hourly[[col]])), paste(sp, "has NA after fill"))
    chk(!any(is.na(res$fivemin[[col]])), paste(sp, "5-min has NA"))
    chk(all(res$fivemin[[col]] >= 0), paste(sp, "negative values"))
  }
  # hold-expansion: first 12 five-min values equal the first hourly value
  chk(all(abs(res$fivemin$Cout_NO2[1:12] - res$hourly$Cout_NO2[1]) < 1e-9), "5-min hold mismatch")
  # timestamps monotonic & 5-min spaced
  dts <- diff(as.numeric(res$fivemin$timestamp))
  chk(all(dts == 300), "5-min grid not uniformly 300 s")

  cat(sprintf("\nRESULT: %d passed, %d failed\n", e$ok, e$fail))
  if (e$fail > 0) stop("load_outdoor validation failed")
  invisible(res)
}

if (sys.nframe() == 0L) run_report_outdoor(dir = Sys.getenv("TORNEO_DIR", "."))

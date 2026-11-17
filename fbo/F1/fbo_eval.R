# ============================================================================
# fbo_eval.R  --  T3 FBO Phase 1: objective evaluator U(x) in PURE R (option B).
#
# The function the LHS / Bayesian optimizer calls per candidate. Given the 5
# operational set-points x and a scenario (orientation via Building North Axis),
# it: edits the IDF (text, make_v4.R idiom) -> runs EnergyPlus (system2) ->
# reads eplusout.csv (load_calibrated_series.R idiom) -> builds the df ->
# indicators.R / g5.R -> fis_scoring.R -> returns the scalar utility U.
#
# No Python / reticulate / eppy: EnergyPlus is the only external process.
# Consistent with Track A's own pure-R idioms (make_v4, load_calibrated_series).
#
# Run inside the FIS project (FIS.Rproj) so `here` finds core/ and data/.
# ============================================================================

library(here)
source(here::here("core", "indicators.R"))    # compute_indicators()
source(here::here("core", "g5.R"))            # compute_g5()
source(here::here("core", "fis_scoring.R"))   # default_variables(), build_scoring_fis(), scoring_utility()
source(here::here("core", "load_outdoor.R"))  # load_outdoor()

FBO <- new.env()

# ---- one-time setup: paths, template, scoring FIS, outdoor series -----------
fbo_setup <- function(base_idf,
                      eplus_exe = "C:/EnergyPlusV9-4-0/energyplus.exe",
                      epw       = "C:/EP_work/Seville_measured_NovDec2021.epw",
                      work_dir  = "C:/EP_work/FBO",
                      data_dir  = here::here("data"),
                      year = 2024, V = 126.7031) {
  stopifnot(file.exists(base_idf), file.exists(eplus_exe))
  FBO$eplus_exe <- eplus_exe; FBO$epw <- epw; FBO$work_dir <- work_dir; FBO$V <- V
  if (!dir.exists(work_dir)) dir.create(work_dir, recursive = TRUE)
  FBO$template <- readLines(base_idf, warn = FALSE)

  # scoring FIS, built once  (confirm signature against your fis_scoring.R)
  FBO$sf <- tryCatch(build_scoring_fis(default_variables()),
                     error = function(e) build_scoring_fis())

  # outdoor pollutant series (Torneo 2024), full-year 5-min, once
  CSV <- list(PM25 = file.path(data_dir, "PM25_HH_2024_SEV_Torneo.csv"),
              NO2  = file.path(data_dir, "NO2_HH_2024_SEV_Torneo.csv"),
              O3   = file.path(data_dir, "O3_HH_2024_SEV_Torneo.csv"))
  od <- load_outdoor(CSV, year = year, five_min = TRUE, verbose = FALSE)$fivemin
  od$key <- format(od$timestamp, "%m/%d %H:%M")   # month-day-HH:MM (year-agnostic BC)
  FBO$outdoor <- od
  cat("[fbo_setup] template", length(FBO$template), "lines | outdoor", nrow(od), "rows | ready\n")
  invisible(TRUE)
}

# ---- IDF text edits (make_v4.R idiom) --------------------------------------
.set_sched <- function(L, name, value) {
  pat <- sprintf("(Schedule:Constant,\\s*%s,\\s*Any Number,\\s*)[-0-9.eE+]+;", name)
  i <- grep(pat, L)
  if (length(i) != 1L) stop("Schedule:Constant not unique for ", name)
  L[i] <- sub(pat, sprintf("\\1%s;", format(value, scientific = FALSE, trim = TRUE)), L[i]); L
}
.set_north <- function(L, deg) {
  i <- grep("!- North Axis", L)
  if (length(i) != 1L) stop("North Axis field not unique")
  L[i] <- sub("^\\s*[-0-9.]+\\s*,", sprintf("   %s,", deg), L[i]); L
}

# ---- run EnergyPlus (system2) ----------------------------------------------
.run_eplus <- function(idf_path, out_dir) {
  if (dir.exists(out_dir)) unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE)
  code <- system2(FBO$eplus_exe,
                  c("-w", shQuote(FBO$epw), "-r", "-d", shQuote(out_dir), shQuote(idf_path)),
                  stdout = FALSE, stderr = FALSE)
  csv <- file.path(out_dir, "eplusout.csv")
  if (!file.exists(csv)) stop("EnergyPlus produced no eplusout.csv (exit ", code, ") for ", idf_path)
  csv
}

# ---- read the series indicators.R / g5.R need (load_calibrated idiom) -------
.read_series <- function(csv) {
  ec <- read.csv(csv, check.names = FALSE)
  g <- function(pat) { j <- grep(pat, names(ec)); j <- j[grepl("TimeStep", names(ec)[j])]
                       if (length(j) < 1L) stop("column not found: ", pat); ec[[j[1]]] }
  data.frame(
    dt   = trimws(ec[[1]]),
    top  = g("CLASSROOM:Zone Operative Temperature"),
    cin  = g("CLASSROOM:Zone Air CO2 Concentration"),
    tout = g("Site Outdoor Air Drybulb Temperature"),
    inf  = g("CLASSROOM:AFN Zone Infiltration Air Change Rate"),  # AFN total outdoor air (incl. windows)
    peo  = g("CLASSROOM:Zone People Sensible Heating Rate"),
    s    = as.integer(round(g("FBO Window State"))),
    stringsAsFactors = FALSE)
}

# ---- objective: U for one candidate ----------------------------------------
# x        : list(tau_CO2, delta_hys, d_slot, f_open, tau_AQ)
# north    : Building North Axis [deg] for the scenario orientation
evaluate_U <- function(x, north = 0, tag = "cand") {
  L <- FBO$template
  L <- .set_sched(L, "x_tau_CO2_sched",   x$tau_CO2)
  L <- .set_sched(L, "x_delta_hys_sched", x$delta_hys)
  L <- .set_sched(L, "x_d_slot_sched",    x$d_slot)
  L <- .set_sched(L, "x_f_open_sched",    x$f_open)
  L <- .set_sched(L, "x_tau_AQ_sched",    x$tau_AQ)
  L <- .set_north(L, north)
  idf <- file.path(FBO$work_dir, paste0(tag, ".idf"))
  writeLines(L, idf)

  d   <- .read_series(.run_eplus(idf, file.path(FBO$work_dir, tag)))

  # calendar day-of-year (compute_indicators groups by `day`: g2 theta_rm per day,
  # g3 daily EDM mean, g4 switches per occupied session) -- NOT a constant
  mmdd <- regmatches(d$dt, regexec("^(\\d{2})/(\\d{2})", d$dt))
  mm <- as.integer(vapply(mmdd, function(z) z[2], "")); dd <- as.integer(vapply(mmdd, function(z) z[3], ""))
  day <- as.integer(as.Date(sprintf("2021-%02d-%02d", mm, dd)) - as.Date("2020-12-31"))

  # T_occ = REAL school calendar: exclude non-lectivo weeks
  #   Christmas break 22 Dec - 7 Jan  |  summer break 22 Jun - 10 Sep (Jul/Aug already 0)
  notlect <- (mm == 12 & dd >= 22) | (mm == 1 & dd <= 7) |
             (mm == 6  & dd >= 22) | (mm %in% c(7, 8))   | (mm == 9 & dd <= 10)
  occ <- (d$peo > 0) & !notlect
  Q   <- d$inf * FBO$V                              # ach -> m3/h  (compute_g5 uses a = Q/V in h^-1)

  # g1..g4
  ind <- compute_indicators(data.frame(day = day, theta_op = d$top, C_in = d$cin,
                                        C_out = 420, theta_out = d$tout, s = d$s, occ = occ))
  # g5: align E+ timesteps to outdoor Cout_* by month-day-HH:MM (year-agnostic BC)
  key <- sub("^(\\d{2})/(\\d{2})\\s+(\\d{2}):(\\d{2}):\\d{2}$", "\\1/\\2 \\3:\\4", d$dt)
  oi  <- match(key, FBO$outdoor$key)
  if (anyNA(oi)) {                                   # fill unmatched (e.g. E+ '24:00' edge)
    last <- which(!is.na(oi))[1]
    for (i in seq_along(oi)) if (is.na(oi[i])) oi[i] <- last else last <- oi[i]
  }
  g5  <- compute_g5(data.frame(Q = Q,
                               Cout_PM25 = FBO$outdoor$Cout_PM25[oi],
                               Cout_NO2  = FBO$outdoor$Cout_NO2[oi],
                               Cout_O3   = FBO$outdoor$Cout_O3[oi],
                               occ = occ), V = FBO$V)$g5

  U <- scoring_utility(FBO$sf, ind$g1, ind$g2, ind$g3, ind$g4, g5)
  list(U = U, g1 = ind$g1, g2 = ind$g2, g3 = ind$g3, g4 = ind$g4, g5 = g5,
       x = x, north = north)
}

# ---- usage (integration test on the 17-Dec FBO IDF) ------------------------
# fbo_setup(base_idf = "F1/v4_G1_17dec.idf",
#           epw = "C:/EP_work/Seville_measured_NovDec2021.epw")
# res <- evaluate_U(list(tau_CO2=1000, delta_hys=100, d_slot=15, f_open=0.30, tau_AQ=3),
#                   north = 0, tag = "test01")
# str(res)

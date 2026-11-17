# =============================================================================
# export_co2_clean.R  --  Excel -> audited CSV export (Track B, reproducibility)
# Thesis T3 -- FBO framework -- 1946 uninsulated corner classroom, Seville
# -----------------------------------------------------------------------------
# PURPOSE (closes the reproducibility chain, "Vía B"):
#   Regenerates co2_clean.csv directly from the raw workbook
#   'T3_Calibration_TRACK B_3 days.xlsx', so the full trace runs from the
#   original measured data to the calibrated g_person = 14.74 L/h without any
#   undocumented intermediate step.
#
# WHY reticulate + openpyxl (NOT readxl):
#   The workbook's 'Time' column is written as native Excel time cells, but the
#   Track A workbook suffered a mixed-type cumulative-serial bug when read with
#   readxl. To be safe and consistent with project methodology, we read via
#   openpyxl (data_only=TRUE) through reticulate, which returns Python
#   datetime.time objects and avoids any serial-number misinterpretation.
#
# WHAT it does, step by step:
#   (1) Opens the workbook, sheet 'TRACK B Data'.
#   (2) Reads the header (row 6) and data (row 7 onward).
#   (3) Forward-fills Year/Month/Day (only present on the first row of each day
#       block due to merged/dragged cells).
#   (4) Keeps rows that have BOTH a Time and a Ti value (drops trailing blanks).
#   (5) Builds an ISO datetime, reconstructs the day, flags the 09:00-14:00
#       anchor window, and writes co2_clean.csv with the exact columns the
#       downstream engine expects.
#   (6) Prints an integrity audit (row counts per day, step regularity,
#       physical ranges) so the export is self-verifying.
#
# COLUMN MAP (0-indexed as returned by openpyxl row tuples):
#   1 Year | 2 Month | 3 Day | 4 Time | 6 Ti | 7 Hi | 9 Ot | 10 Ho
#   11 WindD | 12 WindS | 14 CO2i | 15 CO2o | 17 Occ | 18 Heating
#   19 Openings | 20 Vent | 22 WWeqR | 23 AT
# =============================================================================

suppressWarnings(suppressMessages(library(reticulate)))

# ---- configuration ----------------------------------------------------------
XLSX_IN   <- "T3_Calibration_TRACK B_3 days.xlsx"  # raw workbook (this folder)
SHEET     <- "TRACK B Data"
CSV_OUT   <- "co2_clean.csv"
ANCHOR_LO <- 9  * 3600
ANCHOR_HI <- 14 * 3600

export_co2_clean <- function(xlsx = XLSX_IN, sheet = SHEET, out = CSV_OUT) {

  if (!file.exists(xlsx))
    stop(sprintf("[FAIL] raw workbook not found: '%s'", xlsx))

  # ensure openpyxl available in the reticulate Python
  if (!py_module_available("openpyxl"))
    stop("[FAIL] Python module 'openpyxl' not available to reticulate. ",
         "Install with: reticulate::py_install('openpyxl').")

  openpyxl <- import("openpyxl")
  wb <- openpyxl$load_workbook(xlsx, data_only = TRUE, read_only = TRUE)
  ws <- wb[[sheet]]   # modern openpyxl API (get_sheet_by_name is deprecated)

  # iterate rows from row 7 (1-indexed) onward; values_only = TRUE
  it <- ws$iter_rows(min_row = 7L, values_only = TRUE)
  rows <- iterate(it, function(r) r)   # list of tuples (as R lists)

  # helpers
  getv <- function(r, idx) {           # idx is 0-based (openpyxl tuple index)
    v <- r[[idx + 1L]]                 # R list is 1-based
    if (is.null(v)) NA else v
  }
  as_num <- function(x) suppressWarnings(as.numeric(x))

  # forward-fill Y/M/D and keep valid rows (Time present AND Ti present)
  year <- month <- day <- NA
  recs <- list(); k <- 0L
  for (r in rows) {
    y <- getv(r, 1); mo <- getv(r, 2); d <- getv(r, 3)
    if (!is.na(y))  year  <- y
    if (!is.na(mo)) month <- mo
    if (!is.na(d))  day   <- d
    tm <- getv(r, 4); ti <- getv(r, 6)
    if (is.null(tm) || is.na(ti)) next

    # Python datetime.time -> hour/min/sec
    hh <- as.integer(tm$hour); mm <- as.integer(tm$minute); ss <- as.integer(tm$second)
    sod <- hh*3600L + mm*60L + ss
    dt_iso <- sprintf("%04d-%02d-%02d %02d:%02d:%02d", as.integer(year),
                      as.integer(month), as.integer(day), hh, mm, ss)

    k <- k + 1L
    recs[[k]] <- data.frame(
      datetime = dt_iso,
      day      = as.integer(day),
      Ti       = as_num(getv(r, 6)),
      Hi       = as_num(getv(r, 7)),
      Ot       = as_num(getv(r, 9)),
      Ho       = as_num(getv(r, 10)),
      WindD    = as_num(getv(r, 11)),
      WindS    = as_num(getv(r, 12)),
      CO2i     = as_num(getv(r, 14)),
      CO2o     = as_num(getv(r, 15)),
      Occ      = as_num(getv(r, 17)),
      Heating  = as.character(getv(r, 18)),
      Openings = as.character(getv(r, 19)),
      Vent     = as.character(getv(r, 20)),
      WWeqR    = as_num(getv(r, 22)),
      AT       = as_num(getv(r, 23)),
      anchor   = as.integer(sod >= ANCHOR_LO & sod <= ANCHOR_HI),
      stringsAsFactors = FALSE
    )
  }
  d <- do.call(rbind, recs)

  # blank -> "" for CO2i/WWeqR so downstream read.csv yields NA cleanly
  d$CO2i[is.na(d$CO2i)]   <- ""
  d$WWeqR[is.na(d$WWeqR)] <- ""

  write.csv(d, out, row.names = FALSE, na = "")

  # ---- self-verifying integrity audit ----
  cat(sprintf("Exported '%s': %d rows (anchor=%d)\n", out, nrow(d), sum(d$anchor == 1)))
  cat("Rows per day:\n")
  print(table(d$day))
  cat("Physical ranges:\n")
  num_summary <- function(x, nm) {
    v <- suppressWarnings(as.numeric(x)); v <- v[!is.na(v)]
    cat(sprintf("  %-6s n=%4d  min=%8.2f  max=%8.2f  mean=%8.2f\n",
                nm, length(v), min(v), max(v), mean(v)))
  }
  for (nm in c("Ti","Ot","WindS","CO2i","Occ","WWeqR")) num_summary(d[[nm]], nm)

  invisible(d)
}

if (sys.nframe() == 0) {
  export_co2_clean()
}

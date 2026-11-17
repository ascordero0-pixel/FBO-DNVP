# =============================================================================
# unpack_runs.R -- decompress the archived EnergyPlus series of Track A
# -----------------------------------------------------------------------------
# The calibrated v3 / v6 EnergyPlus outputs are stored gzip-compressed in the
# repository (runs/<model>/eplusout.csv.gz). The Track A scripts read the plain
# runs/<model>/eplusout.csv, so run this once from calibration/trackA:
#
#   source("unpack_runs.R")
#
# Base R only; works on Windows, macOS and Linux. Existing CSVs are kept.
# =============================================================================
for (gz in list.files("runs", pattern = "eplusout\\.csv\\.gz$", recursive = TRUE,
                      full.names = TRUE)) {
  csv <- sub("\\.gz$", "", gz)
  if (file.exists(csv)) { message("[skip] ", csv, " already present"); next }
  con_in <- gzfile(gz, "rb"); con_out <- file(csv, "wb")
  repeat {
    buf <- readBin(con_in, "raw", n = 1e7)
    if (!length(buf)) break
    writeBin(buf, con_out)
  }
  close(con_in); close(con_out)
  message("[ok] ", csv, " (", format(file.size(csv), big.mark = ","), " bytes)")
}

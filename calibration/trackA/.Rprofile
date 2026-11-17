# =====================================================================
# .Rprofile  —  T3 Calibration project
# Place in the ROOT of the project folder:
#   ...\04_Calibration protocols\Calibration R\.Rprofile
#
# Purpose: pin the Python interpreter for reticulate at session start,
# BEFORE any package or chunk can claim a different one. This is more
# robust than calling use_python() by hand each time.
#
# Keep this file MINIMAL. No library() calls, no data loading — it runs
# very early in startup when packages may not yet be available.
# =====================================================================

Sys.setenv(
  RETICULATE_PYTHON =
    "C:/Users/usuario/AppData/Local/Programs/Python/Python313/python.exe"
)

# Optional, harmless confirmation in interactive sessions only:
if (interactive()) {
  message("[.Rprofile] RETICULATE_PYTHON pinned to Python 3.13 (eppy env).")
}

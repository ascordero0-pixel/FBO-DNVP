# Run once in the environment used for the study (from the repository root):
#   Rscript environment/capture_session_info.R
# Writes environment/sessionInfo_author.txt (R, OS and package versions).
pk <- c("here", "ggplot2", "lubridate", "readxl", "patchwork", "scales",
        "DiceKriging", "DiceOptim", "mco", "reticulate", "rstudioapi")
suppressWarnings(invisible(lapply(pk, requireNamespace, quietly = TRUE)))
out <- file.path("environment", "sessionInfo_author.txt")
sink(out)
cat("Captured:", format(Sys.time(), "%Y-%m-%d %H:%M %Z"), "\n\n")
print(sessionInfo())
cat("\nPackage versions:\n")
for (p in pk) cat(sprintf("  %-12s %s\n", p,
  tryCatch(as.character(packageVersion(p)), error = function(e) "not installed")))
ep <- "C:/EnergyPlusV9-4-0/energyplus.exe"
if (file.exists(ep)) cat("\nEnergyPlus:", system2(ep, "--version", stdout = TRUE), "\n")
sink()
message("written: ", out)

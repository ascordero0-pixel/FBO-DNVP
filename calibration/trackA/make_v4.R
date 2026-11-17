# =============================================================================
# make_v4.R  --  Generate the corrected, publishable IDF (v4) from v3
# -----------------------------------------------------------------------------
# PURPOSE: persist in the IDF file the Wall_3 correction that the calibration
#   engine otherwise applies only in runtime (ce$normalize_wall3_to_outdoors).
#   The inherited T1 DesignBuilder export tagged the NORTH façade (Wall_3) with a
#   spurious ExteriorNaturalVentedCavity + OtherSideConditionsModel; it must be a
#   plain 'Outdoors' façade like its twin Wall_5 (South, same construction/area).
#
#   v4 is the RECOMMENDED IDF to publish: it is self-contained, so a reviewer
#   reproduces NMBE = -5.41 % by opening the IDF directly, with no hidden runtime
#   normalisation. Proven equivalent to v3+normalize (both -> -5.41 / 8.04).
#
# SAFETY: v3 is NEVER modified. Output is a new file v4.
#
# USAGE: source("make_v4.R")   # writes T1_aula_base_94_v4.idf and verifies it
# =============================================================================

make_v4 <- function(src = "T1_aula_base_94_v3.idf",
                    dst = "T1_aula_base_94_v4.idf") {
  stopifnot(file.exists(src))
  idf <- readLines(src, warn = FALSE)

  # --- locate the BuildingSurface:Detailed object of Wall_3 ---
  name_ln <- grep("CLASSROOM_Wall_3_0_0,", idf)[1]
  if (is.na(name_ln)) stop("Wall_3 surface not found in ", src)
  kw_ln  <- max(grep("BuildingSurface:Detailed", idf[1:name_ln], ignore.case = TRUE))
  end_ln <- kw_ln + which(grepl(";", idf[kw_ln:(kw_ln + 40)]))[1] - 1
  rng    <- kw_ln:end_ln

  # 1) OtherSideConditionsModel -> Outdoors, ,
  for (i in rng) if (grepl("OtherSideConditionsModel", idf[i]))
    idf[i] <- "      Outdoors, ,                                 !- Outside Face Environment"

  # 2) delete the dangling 'Outside Boundary Condition Object' line (the ...BoundaryObject)
  kill <- rng[grepl("BoundaryObject", idf[rng]) & grepl("Boundary Condition Object", idf[rng])]
  if (length(kill)) idf <- idf[-kill]

  # 3) exposure NoSun/NoWind -> Sun/Wind (defensive; Wall_3 already exposed)
  nl  <- grep("CLASSROOM_Wall_3_0_0,", idf)[1]
  kw2 <- max(grep("BuildingSurface:Detailed", idf[1:nl], ignore.case = TRUE))
  rng2 <- kw2:(kw2 + which(grepl(";", idf[kw2:(kw2 + 40)]))[1] - 1)
  for (i in rng2) { idf[i] <- sub("\\bNoSun\\b","SunExposed",idf[i])
                    idf[i] <- sub("\\bNoWind\\b","WindExposed",idf[i]) }

  # 4) remove the two spurious SurfaceProperty objects of Wall_3
  drop_object <- function(lines, keyword) {
    s <- grep(keyword, lines, ignore.case = TRUE)
    for (k in rev(s)) { e <- k + which(grepl(";", lines[k:(k + 15)]))[1] - 1
      if (any(grepl("Wall_3_0_0", lines[k:e]))) lines <- lines[-(k:e)] }
    lines
  }
  idf <- drop_object(idf, "SurfaceProperty:OtherSideConditionsModel")
  idf <- drop_object(idf, "SurfaceProperty:ExteriorNaturalVentedCavity")

  writeLines(idf, dst)

  # --- verification ---
  oscm <- length(grep("SurfaceProperty:OtherSideConditionsModel", idf, ignore.case=TRUE))
  natv <- length(grep("SurfaceProperty:ExteriorNaturalVentedCavity", idf, ignore.case=TRUE))
  bref <- length(grep("Wall_3_0_0 BoundaryObject", idf))
  nl3  <- grep("CLASSROOM_Wall_3_0_0,", idf)[1]
  kw3  <- max(grep("BuildingSurface:Detailed", idf[1:nl3], ignore.case=TRUE))
  en3  <- kw3 + which(grepl(";", idf[kw3:(kw3+40)]))[1] - 1
  blk3 <- idf[kw3:en3]
  cand <- blk3[grepl("Outside Face Environment", blk3) & !grepl("^\\s*!", blk3)]
  bc3  <- trimws(gsub("[,;].*", "", cand[1]))

  cat("=========  make_v4: Wall_3 (v4)  =========\n")
  cat(idf[kw3:en3], sep = "\n")
  cat(sprintf("\n[OK] v4 written: %s  (%d lines vs %d in v3)\n",
              dst, length(idf), length(readLines(src))))
  cat(sprintf("[chk] OSCM=%d | NatVent=%d | BoundaryObject refs=%d | Wall_3 BC=%s\n",
              oscm, natv, bref, bc3))
  pass <- (oscm==0 && natv==0 && bref==0 && grepl("Outdoors", bc3, ignore.case=TRUE))
  cat(sprintf("[chk] v4 integrity: %s\n", if (pass) "PASS" else "FAIL"))
  if (!pass) stop("make_v4: verification FAILED — inspect the block above")
  cat("==========================================\n")
  invisible(dst)
}

if (interactive()) make_v4()

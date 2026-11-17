# =============================================================================
# trackA_5min_metrics.R
# Re-run v3 (documented Track A) producing the values needed to fill §3.1:
#   * NMBE / CV-RMSE at BOTH 5-min and hourly resolution (nocturnal 22-07)
#   * identified thermal time constant tau (lumped-capacitance regression)
#   * envelope conductance UA (from EnergyPlus Envelope Summary)
#   * nocturnal / global infiltration (ren/h) for the record
#
# REQUIRES: source("run_calibration.R"); phase1_setup(). OneDrive PAUSED.
# Reuses the documented v3 config: crack=1, cap=1, STEP heat=20, real occupancy.
# =============================================================================
library(readxl); library(lubridate)
ce <- get("ce",envir=.STATE); cfg <- get("cfg",envir=.STATE)
ZUP <- "ZXCASESTUDYCLASSROOM:CLASSROOM"; IDF_V3 <- "T1_aula_base_94_v3.idf"
INT <- paste0(ZUP,":Zone Mean Air Temperature [C](TimeStep:ON)")
# NOTE: exterior column name is resolved robustly AFTER csv is loaded (see below),
# not hardcoded here -- EnergyPlus can emit "(TimeStep)" or "(TimeStep:ON)" depending
# on version/config, and a silent exact-match miss leaves T_out all-NA and tau=NA.

# ---- build v3 (documented) ---------------------------------------------------
idf <- ce$load_idf(cfg, IDF_V3)
ce$set_timestep(idf,12L); ce$scale_afn_crack_coeffs(idf,1.0); ce$set_capacitance_multiplier(idf,1.0)
ce$add_heating_gain(idf,20.0); ce$close_windows_for_trackA(idf)
tryCatch(ce$normalize_wall3_to_outdoors(idf),error=function(e)message("  wall3: ",conditionMessage(e)))
have<-vapply(seq_len(reticulate::py_len(idf$idfobjects[["ScheduleTypeLimits"]])),
  function(i) reticulate::py_to_r(idf$idfobjects[["ScheduleTypeLimits"]][[as.integer(i-1)]]$Name),character(1))
if(!("Fraction"%in%have)) idf$newidfobject("ScheduleTypeLimits",Name="Fraction",
  Lower_Limit_Value="0",Upper_Limit_Value="1",Numeric_Type="Continuous",Unit_Type="Dimensionless")
occ<-c("Through: 12/31","For: Weekdays","Until: 09:10","0","Until: 12:00","1",
  "Until: 12:30","0","Until: 14:00","1","Until: 24:00","0","For: AllOtherDays","Until: 24:00","0")
do.call(idf$newidfobject,c(list("Schedule:Compact",Name="T3_OCC_REAL",Schedule_Type_Limits_Name="Fraction"),
  setNames(as.list(occ),paste0("Field_",seq_along(occ)))))
ppl<-idf$idfobjects[["People"]]
for(i in seq_len(reticulate::py_len(ppl))){o<-ppl[[as.integer(i-1)]]
  reticulate::py_set_attr(o,"Number_of_People_Schedule_Name","T3_OCC_REAL")}
ce$add_output_variables(idf)
# exterior drybulb at timestep (duplicate Output:Variable is non-fatal)
idf$newidfobject("Output:Variable",Key_Value="*",Variable_Name="Site Outdoor Air Drybulb Temperature",Reporting_Frequency="Timestep")
# NOTE: the base IDF ALREADY has OutputControl:Table:Style + Output:Table:SummaryReports
# (DesignBuilder export). We must NOT touch them — EnergyPlus allows a single property
# each, and modifying/adding turns them invalid. The existing tabular output is reused.

out<-file.path(PATHS$runs_dir,"v3_5min.idf"); idf$saveas(out,encoding="utf-8")
unlink(file.path(PATHS$runs_dir,"v3_5min"),recursive=TRUE,force=TRUE)
od<-ce$run_energyplus(cfg,out,"v3_5min"); cat("[RUN] done.\n")

# ---- read sim 5-min interior + exterior --------------------------------------
csv<-read.csv(file.path(PATHS$runs_dir,"v3_5min","eplusout.csv"),check.names=FALSE)
parse_dt<-function(x){x<-trimws(x);mo<-as.integer(substr(x,1,2));da<-as.integer(substr(x,4,5))
  tt<-sub("^\\d{2}/\\d{2}\\s+","",x);hh<-as.integer(substr(tt,1,2));mi<-as.integer(substr(tt,4,5))
  nd<-which(hh==24);hh[nd]<-0;ts<-as.POSIXct(sprintf("2021-%02d-%02d %02d:%02d:00",mo,da,hh,mi),tz="UTC");ts[nd]<-ts[nd]+24*3600;ts}
ts<-parse_dt(csv[[1]])
EXT <- grep("Outdoor Air Drybulb.*TimeStep\\)$", names(csv), value=TRUE)[1]  # robust match
if (is.na(EXT)) stop("[FAIL] Exterior drybulb column not found in eplusout.csv -- ",
                      "check Output:Variable reporting frequency / naming.")
cat(sprintf("[EXT] using column: %s\n", EXT))
sim<-data.frame(timestamp=ts, T_sim=as.numeric(csv[[INT]]),
                T_out=if(EXT%in%names(csv)) as.numeric(csv[[EXT]]) else NA)

# ---- measured 5-min ----------------------------------------------------------
d<-read_excel(PATHS$xlsx_track_a,sheet=1,skip=6,col_names=FALSE);o<-if(ncol(d)>=22) 0L else -1L
ff<-function(z){for(i in seq_along(z)) if(is.na(z[i])&&i>1) z[i]<-z[i-1];z}
Y<-ff(d[[2+o]]);M<-ff(d[[3+o]]);hm<-format(d[[5+o]],"%H:%M");bad<-is.na(hm)|hm=="NA"
if(any(bad)) hm[bad]<-substr(as.character(d[[5+o]][bad]),12,16)
mts<-as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00",as.integer(Y),as.integer(M),
  as.integer(d[[4+o]]),as.integer(substr(hm,1,2)),as.integer(substr(hm,4,5))),tz="UTC")
meas<-data.frame(timestamp=mts,T_meas=as.numeric(d[[7+o]]));meas<-meas[!is.na(meas$T_meas),]

# align at 5-min (floor both) --------------------------------------------------
sim$k<-floor_date(sim$timestamp,"5 mins"); meas$k<-floor_date(meas$timestamp,"5 mins")
M5<-merge(sim, aggregate(T_meas~k,meas,mean), by="k"); M5<-M5[complete.cases(M5[,c("T_sim","T_meas")]),]
hh<-as.integer(format(M5$k,"%H")); M5$night<-(hh>=22|hh<7)

metric<-function(df){e<-df$T_sim-df$T_meas;c(NMBE=mean(e)/mean(df$T_meas)*100,
  CVRMSE=sqrt(mean(e^2))/mean(df$T_meas)*100,n=nrow(df))}
n5<-M5[M5$night,]; m5<-metric(n5)
# hourly
M5$hr<-floor_date(M5$k,"hour")
H<-merge(aggregate(T_sim~hr,M5,mean),aggregate(T_meas~hr,M5,mean),by="hr")
hh2<-as.integer(format(H$hr,"%H")); Hn<-H[hh2>=22|hh2<7,]; mh<-metric(Hn)

# ---- thermal time constant tau (lumped regression, nocturnal) ----------------
tau <- NA
tryCatch({
  nn <- n5[order(n5$k), ]; n <- nrow(nn)
  dt_h <- as.numeric(difftime(nn$k[-1], nn$k[-n], units="hours"))   # length n-1, hours
  dTdt <- diff(nn$T_sim) / dt_h                                     # degC/h
  dTd  <- nn$T_sim[-n] - nn$T_out[-n]                               # (Tin-Tout) at step start
  ok <- is.finite(dTdt) & is.finite(dTd) & abs(dt_h - 5/60) < 1e-3  # consecutive 5-min only
  cat(sprintf("[tau] usable steps: %d of %d nocturnal points\n", sum(ok), n))
  if (sum(ok) > 10) {
    fit <- lm(dTdt[ok] ~ dTd[ok]); slope <- unname(coef(fit)[2]); tau <- -1/slope
    cat(sprintf("[tau] dT/dt = %.4f*(Tin-Tout) %+.3f  (R2=%.3f)\n",
        slope, unname(coef(fit)[1]), summary(fit)$r.squared))
    cat(sprintf("[tau] identified thermal time constant tau = %.2f h\n", tau))
  } else cat("[tau] too few usable steps (check exterior column was read)\n")
}, error=function(e) cat("[tau] skipped:", conditionMessage(e), "\n"))

# ---- infiltration (from check function if available) -------------------------
inf<-tryCatch(ce$summarize_infiltration(od), error=function(e) NULL)
if(!is.null(inf)) print(reticulate::py_to_r(inf))

# ---- UA from Envelope Summary ------------------------------------------------
UA<-NA
tblf<-list.files(file.path(PATHS$runs_dir,"v3_5min"),
  pattern="tbl.*\\.(csv|htm|html)$|Table.*\\.(csv|htm|html)$",full.names=TRUE,ignore.case=TRUE)
if(length(tblf)>0){
  raw<-readLines(tblf[1]); 
  # opaque exterior: rows with 'U-Factor with Film' col and 'Net Area'; fenestration: 'Glass U-Factor' + area
  # crude parse: find lines under 'Opaque Exterior' / 'Exterior Fenestration'
  cat("\n[UA] Envelope Summary found:", basename(tblf[1]),"\n")
  cat("    (sumar manualmente U-Factor x Net Area de muros+cubierta+ventanas exteriores)\n")
  cat("    primeras filas de envolvente:\n")
  idx<-grep("Opaque Exterior|Exterior Fenestration|U-Factor|Net Area|Glass U-Factor",raw)
  for(i in head(idx,40)) cat("     ",substr(raw[i],1,160),"\n")
} else cat("\n[UA] no tabular file found — coger UA del survey T1 [11]\n")

# ---- report ------------------------------------------------------------------
cat("\n========  TRACK A v3 — values for §3.1  ========\n")
cat(sprintf("  Nights / period      : 50 free-running nights (1 Nov - 20 Dec 2021)\n"))
cat(sprintf("  5-min nocturnal      : NMBE = %+.2f %%   CV-RMSE = %.2f %%   (n = %d)\n", m5["NMBE"],m5["CVRMSE"],m5["n"]))
cat(sprintf("  Hourly nocturnal     : NMBE = %+.2f %%   CV-RMSE = %.2f %%   (n = %d)\n", mh["NMBE"],mh["CVRMSE"],mh["n"]))
cat(sprintf("  Thermal time const.  : tau = %.2f h\n", tau))
cat(sprintf("  Mean indoor T (night): %.2f C ; sensor +-0.3C floor = %.2f %% CV-RMSE\n",
  mean(n5$T_meas), 0.3/mean(n5$T_meas)*100))
cat("================================================\n")
write.csv(data.frame(res=c("5min_NMBE","5min_CVRMSE","hourly_NMBE","hourly_CVRMSE","tau_h","mean_T_night"),
  value=c(m5["NMBE"],m5["CVRMSE"],mh["NMBE"],mh["CVRMSE"],tau,mean(n5$T_meas))),
  file.path(PATHS$runs_dir,"trackA_5min_values.csv"),row.names=FALSE)
message("[OK] saved trackA_5min_values.csv")

# =============================================================================
# close_trackA_v6.R   (RAMP+TAIL heating, blinded upper-zone reading)
# 2-zone model (v6) for Track A. Model kept INTACT (both classrooms heated &
# occupied). Changes vs previous version:
#   * Heating profile: STEP -> RAMP+TAIL (water-radiator signature),
#     energy-conserved, applied to BOTH zones.
#   * Result reading: BLINDED to the exact case-study column at timestep res:
#       'ZXCASESTUDYCLASSROOM:CLASSROOM:Zone Mean Air Temperature [C](TimeStep:ON)'
#
# REQUIRES: source("run_calibration.R"); phase1_setup(); phase2_check_inputs().
#           OneDrive PAUSED.
# OUTPUTS (runs/): v6_trackA.idf · trackA_v6_metrics.csv · Fig_trackA_v6.(pdf|png)
#
# --- PARAMETERS ---------------------------------------------------------------
W_M2     <- 20      # daily heating energy equivalent to 20 W/m2 flat 07-18h
TAU_UP   <- 3       # ramp-up time constant (h)  -> water heating up
TAU_DOWN <- 3       # tail time constant (h)     -> residual water heat
HEAT_ON  <- 7;  HEAT_OFF <- 18
CRACK    <- 1.0; CAP <- 1.0
IDF_V6   <- "T1_aula_base_94_v6.idf"
ZUP      <- "ZXCASESTUDYCLASSROOM:CLASSROOM"          # case study (measured)
ZONES    <- c(ZUP, "ZXCSSTDYCLSSRMBLW:CLASSROOMBELOW")# both heated (intact)
COL_EXACT <- paste0(ZUP, ":Zone Mean Air Temperature [C](TimeStep:ON)")
# =============================================================================

library(ggplot2); library(readxl); library(lubridate)
ce <- get("ce", envir=.STATE); cfg <- get("cfg", envir=.STATE)
idf <- ce$load_idf(cfg, IDF_V6)

# ── 1) REAL occupancy schedule, assigned to BOTH People ───────────────────────
have <- vapply(seq_len(reticulate::py_len(idf$idfobjects[["ScheduleTypeLimits"]])),
   function(i) reticulate::py_to_r(idf$idfobjects[["ScheduleTypeLimits"]][[as.integer(i-1)]]$Name),
   character(1))
if (!("Fraction" %in% have))
  idf$newidfobject("ScheduleTypeLimits", Name="Fraction", Lower_Limit_Value="0",
                   Upper_Limit_Value="1", Numeric_Type="Continuous", Unit_Type="Dimensionless")
occ <- c("Through: 12/31","For: Weekdays","Until: 09:10","0","Until: 12:00","1",
         "Until: 12:30","0","Until: 14:00","1","Until: 24:00","0",
         "For: AllOtherDays","Until: 24:00","0")
do.call(idf$newidfobject, c(list("Schedule:Compact", Name="T3_OCC_REAL",
        Schedule_Type_Limits_Name="Fraction"),
        setNames(as.list(occ), paste0("Field_", seq_along(occ)))))
ppl <- idf$idfobjects[["People"]]
for (i in seq_len(reticulate::py_len(ppl))) {
  o <- ppl[[as.integer(i-1)]]
  reticulate::py_set_attr(o, "Number_of_People_Schedule_Name", "T3_OCC_REAL")
}
cat("[1] Occupancy 'T3_OCC_REAL' -> ", reticulate::py_len(ppl), " People object(s).\n", sep="")

# ── 2) Normalize BOTH Wall_3 (drop NatVentCavity -> Outdoors) ─────────────────
drop_all <- function(key){ lst<-idf$idfobjects[[key]]; n<-reticulate::py_len(lst)
  for(k in seq_len(n)) idf$removeidfobject(lst[[0L]]); n }
n_cav <- drop_all("SurfaceProperty:ExteriorNaturalVentedCavity")
n_osm <- tryCatch(drop_all("SurfaceProperty:OtherSideConditionsModel"), error=function(e) 0)
surfs <- idf$idfobjects[["BuildingSurface:Detailed"]]; n_fixed <- 0
for (i in seq_len(reticulate::py_len(surfs))) {
  s <- surfs[[as.integer(i-1)]]
  if (reticulate::py_to_r(s$Outside_Boundary_Condition)=="OtherSideConditionsModel") {
    reticulate::py_set_attr(s,"Outside_Boundary_Condition","Outdoors")
    reticulate::py_set_attr(s,"Outside_Boundary_Condition_Object","")
    reticulate::py_set_attr(s,"Sun_Exposure","SunExposed")
    reticulate::py_set_attr(s,"Wind_Exposure","WindExposed"); n_fixed<-n_fixed+1
  }
}
cat(sprintf("[2] Wall_3: dropped %d NatVentCavity + %d OSCM; %d wall(s)->Outdoors.\n",
            n_cav,n_osm,n_fixed))

# ── 3) RAMP+TAIL heating profile (energy conserved), BOTH zones ───────────────
prof <- function(h){
  if (h<=HEAT_ON) return(0)
  if (h< HEAT_ON+TAU_UP)  return((h-HEAT_ON)/TAU_UP)
  if (h< HEAT_OFF)        return(1)
  if (h< HEAT_OFF+TAU_DOWN) return(1-(h-HEAT_OFF)/TAU_DOWN)
  0 }
grid <- seq(0,24,by=1/120)
integral <- mean(vapply(grid, prof, numeric(1)))*24          # unit*hours
D <- W_M2 * (HEAT_OFF-HEAT_ON) / integral                    # conserve energy
# build step schedule (30-min steps on ramp/tail, midpoint values)
edges <- sort(unique(c(seq(HEAT_ON, HEAT_ON+TAU_UP, .5),
                       HEAT_OFF, seq(HEAT_OFF, HEAT_OFF+TAU_DOWN, .5), 24)))
edges <- edges[edges>HEAT_ON & edges<=24]
fld <- c("Through: 12/31","For: Weekdays","Until: 07:00","0"); prev<-HEAT_ON
for (e in edges){
  v <- prof((prev+e)/2)                                       # midpoint value
  hh<-floor(e); mm<-round((e-hh)*60); if(mm==60){hh<-hh+1;mm<-0}
  fld <- c(fld, sprintf("Until: %02d:%02d", hh, mm), sprintf("%.4f", max(v,0)))
  prev<-e
}
fld <- c(fld,"For: AllOtherDays","Until: 24:00","0")
do.call(idf$newidfobject, c(list("Schedule:Compact", Name="T3_HEAT_RAMP",
        Schedule_Type_Limits_Name="Fraction"),
        setNames(as.list(fld), paste0("Field_", seq_along(fld)))))
for (z in ZONES){
  sfx <- gsub("[^A-Za-z0-9]","_", z)
  idf$newidfobject("ElectricEquipment", Name=paste0("T3_HEAT_",sfx),
    Zone_or_ZoneList_Name=z, Schedule_Name="T3_HEAT_RAMP",
    Design_Level_Calculation_Method="Watts/Area",
    Watts_per_Zone_Floor_Area=sprintf("%.4f", D), Fraction_Latent="0",
    Fraction_Radiant="0.30", Fraction_Lost="0", EndUse_Subcategory="T3Heating")
}
cat(sprintf("[3] RAMP+TAIL heating (tau_up=%.1f, tau_down=%.1f) in %d zone(s).\n",
            TAU_UP,TAU_DOWN,length(ZONES)))
cat(sprintf("    energy-conserved peak power D = %.2f W/m2 (flat-equiv = %g W/m2).\n", D, W_M2))
cat("    profile (h -> fraction):")
for (h in c(7,8,9,10,12,15,17,18,19,20,21)) cat(sprintf(" %d:%.2f", h, prof(h)))
cat("\n")

# ── 4) Track A free-running: crack, cap, windows closed, CO2 EMS off ──────────
ce$set_timestep(idf, 12L)
ce$scale_afn_crack_coeffs(idf, CRACK)
ce$set_capacitance_multiplier(idf, CAP)
ce$close_windows_for_trackA(idf)
n_pcm <- tryCatch(drop_all("EnergyManagementSystem:ProgramCallingManager"), error=function(e) 0)
ce$add_output_variables(idf)
cat(sprintf("[4] Track A patches (crack=%g, cap=%g); removed %d EMS calling mgr(s).\n",
            CRACK,CAP,n_pcm))

# ── SAVE + verify before running ──────────────────────────────────────────────
out_idf <- file.path(PATHS$runs_dir, "v6_trackA.idf")
idf$saveas(out_idf, encoding="utf-8")
ln <- readLines(out_idf, warn=FALSE)
cav_left <- sum(grepl("ExteriorNaturalVentedCavity", ln))
osm_left <- sum(grepl("OtherSideConditionsModel", ln))
cat(sprintf("\n[SAVED] %s\n   NatVentCavity left: %d (=0) | OSCM left: %d (=0)\n",
            out_idf, cav_left, osm_left))
if (cav_left>0 || osm_left>0) stop("Normalization incomplete - STOP. Paste this to me.")

# ── RUN ───────────────────────────────────────────────────────────────────────
unlink(file.path(PATHS$runs_dir,"v6_trackA"), recursive=TRUE, force=TRUE)
out_dir <- ce$run_energyplus(cfg, out_idf, "v6_trackA")
cat("[RUN] done:", out_dir, "\n")

# ── BLINDED reading of the exact case-study column (timestep) ─────────────────
csvf <- file.path(PATHS$runs_dir, "v6_trackA", "eplusout.csv")
csv  <- read.csv(csvf, check.names=FALSE)
if (!(COL_EXACT %in% names(csv)))
  stop("Exact column not found:\n  ", COL_EXACT,
       "\nAvailable:\n  ", paste(grep("Zone Mean Air Temperature",names(csv),value=TRUE),collapse="\n  "))
cat("[READ] column ->", COL_EXACT, "\n")
x  <- trimws(csv[[1]])                                  # Date/Time, e.g. '11/01  00:05:00'
mo <- as.integer(substr(x,1,2)); da <- as.integer(substr(x,4,5))
tt <- sub("^\\d{2}/\\d{2}\\s+","",x); hh <- as.integer(substr(tt,1,2)); mi <- as.integer(substr(tt,4,5))
nd <- which(hh==24); hh[nd] <- 0
ts <- as.POSIXct(sprintf("2021-%02d-%02d %02d:%02d:00",mo,da,hh,mi),tz="UTC")
ts[nd] <- ts[nd] + 24*3600
sim <- data.frame(timestamp=ts, value=as.numeric(csv[[COL_EXACT]]))
sim <- sim[!is.na(sim$value),]
sim$hour <- floor_date(sim$timestamp,"hour")
sim_h <- aggregate(value~hour,sim,mean,na.rm=TRUE); names(sim_h)[2]<-"T_sim"

# ── measured ──────────────────────────────────────────────────────────────────
d<-read_excel(PATHS$xlsx_track_a,sheet=1,skip=6,col_names=FALSE); o<-if(ncol(d)>=22) 0L else -1L
ff<-function(z){for(i in seq_along(z)) if(is.na(z[i])&&i>1) z[i]<-z[i-1]; z}
Y<-ff(d[[2+o]]);M<-ff(d[[3+o]]);hm<-format(d[[5+o]],"%H:%M");bad<-is.na(hm)|hm=="NA"
if(any(bad)) hm[bad]<-substr(as.character(d[[5+o]][bad]),12,16)
mts<-as.POSIXct(sprintf("%04d-%02d-%02d %02d:%02d:00",as.integer(Y),as.integer(M),
   as.integer(d[[4+o]]),as.integer(substr(hm,1,2)),as.integer(substr(hm,4,5))),tz="UTC")
meas<-data.frame(timestamp=mts,T_meas=as.numeric(d[[7+o]])); meas<-meas[!is.na(meas$T_meas),]
meas$hour<-floor_date(meas$timestamp,"hour")
meas_h<-aggregate(T_meas~hour,meas,mean,na.rm=TRUE)

# ── metrics (night + 24h) + diurnal ───────────────────────────────────────────
H<-merge(sim_h,meas_h,by="hour"); H<-H[complete.cases(H),]
hh2<-as.integer(format(H$hour,"%H")); H$night<-(hh2>=22|hh2<7)
mt<-function(s){e<-s$T_sim-s$T_meas; c(NMBE=mean(e)/mean(s$T_meas)*100,
   CVRMSE=sqrt(mean(e^2))/mean(s$T_meas)*100,n=nrow(s))}
mn<-mt(H[H$night,]); m24<-mt(H)
res<-data.frame(window=c("NIGHT 22-07","FULL 24h"),NMBE=c(mn["NMBE"],m24["NMBE"]),
   CVRMSE=c(mn["CVRMSE"],m24["CVRMSE"]),n=c(mn["n"],m24["n"]))
cat("\n============  TRACK A on v6 (ramp+tail, blinded upper zone)  ============\n")
print(format(res,digits=4),row.names=FALSE)
cat(sprintf("\n  ASHRAE G14 (night): %s\n",
    ifelse(abs(mn["NMBE"])<=10 && mn["CVRMSE"]<=30,"PASS","FAIL")))
H$hod<-hh2; ps<-aggregate(T_sim~hod,H,mean)
pm<-aggregate(T_meas~hod,within(meas_h,{hod<-as.integer(format(hour,"%H"))}),mean)
cat(sprintf("  Diurnal peak: sim %d:00 vs measured %d:00\n",
    ps$hod[which.max(ps$T_sim)], pm$hod[which.max(pm$T_meas)]))
cat("========================================================================\n")
write.csv(res, file.path(PATHS$runs_dir,"trackA_v6_metrics.csv"), row.names=FALSE)

# ── figure (case-study zone only) ─────────────────────────────────────────────
f<-ggplot()+
  geom_line(data=ps,aes(hod,T_sim,color="Simulated (case study, ramp+tail)"),linewidth=0.9)+
  geom_line(data=pm,aes(hod,T_meas,color="Measured"),linewidth=1.2)+
  scale_color_manual(NULL,values=c("Measured"="#E0922F",
     "Simulated (case study, ramp+tail)"="#0E7C7B"))+
  scale_x_continuous("Hour of day",breaks=seq(0,23,3))+
  scale_y_continuous("Mean air temperature (\u00B0C)")+
  labs(title="Track A on v6 - diurnal cycle (case-study zone, ramp+tail heating)",
       subtitle=sprintf("crack=%g \u00B7 cap=%g \u00B7 ramp+tail tau_up=%.1f/tau_down=%.1f \u00B7 occ 09:10-12:00 / 12:30-14:00",
                        CRACK,CAP,TAU_UP,TAU_DOWN))+
  theme_minimal(base_size=11)+theme(legend.position="top",panel.grid.minor=element_blank())
ggsave(file.path(PATHS$runs_dir,"Fig_trackA_v6.png"),f,width=20,height=11,units="cm",dpi=300)
ggsave(file.path(PATHS$runs_dir,"Fig_trackA_v6.pdf"),f,width=20,height=11,units="cm")
message("[OK] Saved v6_trackA.idf + trackA_v6_metrics.csv + Fig_trackA_v6 in ", PATHS$runs_dir)

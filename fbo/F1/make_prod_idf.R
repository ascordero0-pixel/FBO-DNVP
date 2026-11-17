# =============================================================================
# make_prod_idf.R  --  Build the production annual IDF (v4_FBO_annual_prod) from
#                      the base annual IDF (v4_FBO_annual), self-documented.
# -----------------------------------------------------------------------------
# PURPOSE: persist, in ONE reproducible step, the heating/comfort model that the
#   evaluator-closure session (2026-07) added to make U(x) non-degenerate. The
#   base annual IDF ran FREE-RUNNING all year (no heating) -> winter collapsed g2
#   -> U was flat. This script derives the publishable production IDF and is the
#   canonical record of the derivation chain (the intermediate IDFs are retired).
#
# DERIVATION (base -> prod), each step justified in claude/T3_Fase1_inventario_scripts.md:
#   (1) Sizing OFF          -- deterministic; the baseboard uses a FIXED capacity.
#   (2) Thermostatic heating, RADIANT-convective baseboard:
#         - fixed 3000 W  = design load w/ windows closed at 21 C in / 1.9 C out
#           (Winter Design Day): UA 99.3*dT ~1896 W + infiltration 0.64 ach ~519 W
#           ~2.4 kW, +15% margin. Realistic -> does NOT mask CO2-purge overcooling.
#         - radiant fraction 0.4 to the 5 zone surfaces (+0.15 on people): a real
#           radiator warms surfaces -> raises MRT -> raises operative temp.
#         - SEASONAL setback: 20 C night / 21 C occupied (07:00-18:00), Nov-Apr;
#           heating OFF (setpoint 5 C) May-Oct so summer night cooling is unimpeded.
#           Availability 24/7; the setpoint schedule carries the seasonality.
#   (3) Zone Mean Air Temperature output (Timestep) for air-vs-operative diagnostics.
#
# NOTE: the school-calendar T_occ (exclude 22 Dec-7 Jan, 22 Jun-10 Sep) is applied
#   in fbo_eval.R::evaluate_U for the scoring. To bake it into this IDF's occupancy
#   Schedule:Year as well (for the definitive runs), see the TODO at the bottom.
#
# The EMS control law, the 5 x-schedules and the AQI_ext Schedule:File are INHERITED
# untouched from the base -- this script only ADDS the heating stack + one output.
#
# USAGE:  source(here::here("F1","make_prod_idf.R")); make_prod_idf()
# =============================================================================

make_prod_idf <- function(src = here::here("F1", "v4_FBO_annual.idf"),
                          dst = here::here("F1", "v4_FBO_annual_prod.idf")) {
  stopifnot(file.exists(src))
  L <- readLines(src, warn = FALSE)

  ZONE <- "ZXCASESTUDYCLASSROOM:CLASSROOM"

  # --- (1) SimulationControl: zone & system sizing -> No (fixed-capacity, deterministic)
  i <- grep("^\\s*SimulationControl", L)[1]
  if (is.na(i)) stop("SimulationControl not found in ", src)
  L[i + 1] <- sub("^\\s*[A-Za-z]+", "   No", L[i + 1])   # Do the zone sizing calculation
  L[i + 2] <- sub("^\\s*[A-Za-z]+", "   No", L[i + 2])   # Do the system sizing calculation

  # --- (2) heating block, inserted right after the 'ScheduleTypeLimits, On/Off' line
  blk <- c(
    "",
    "! ===================== FBO production heating (make_prod_idf.R) =====================",
    "!  Thermostatic RADIANT-convective baseboard. Fixed 3000 W (design load, windows",
    "!  closed, 21/1.9 C). Radiant fraction 0.4 -> raises MRT/operative. Seasonal setback:",
    "!  20 C night / 21 C occupied (Nov-Apr); heating OFF (5 C) May-Oct. Nocturnal free",
    "!  running (Track A) preserved in summer; CO2/AFN untouched (no outdoor air injected).",
    "! ===================================================================================",
    "",
    "Schedule:Compact, FBO_HeatAvail, On/Off,",
    "   Through: 12/31, For: AllDays, Until: 24:00, 1;",
    "",
    "Schedule:Compact, FBO_HeatCtrlType, Control Type,",
    "   Through: 12/31, For: AllDays, Until: 24:00, 1;   !- 1 = SingleHeating",
    "",
    "Schedule:Compact, FBO_HeatSP, Temperature,",
    "   Through: 4/30,",
    "   For: Weekdays, Until: 07:00, 20.0, Until: 18:00, 21.0, Until: 24:00, 20.0,",
    "   For: AllOtherDays, Until: 24:00, 20.0,",
    "   Through: 10/31,",
    "   For: AllDays, Until: 24:00, 5.0,               !- heating OFF (warm season)",
    "   Through: 12/31,",
    "   For: Weekdays, Until: 07:00, 20.0, Until: 18:00, 21.0, Until: 24:00, 20.0,",
    "   For: AllOtherDays, Until: 24:00, 20.0;",
    "",
    "ZoneControl:Thermostat, FBO_CLASSROOM_Thermostat,",
    sprintf("   %s,                        !- Zone Name", ZONE),
    "   FBO_HeatCtrlType,                              !- Control Type Schedule Name",
    "   ThermostatSetpoint:SingleHeating,              !- Control 1 Object Type",
    "   FBO_CLASSROOM_HeatSP;                           !- Control 1 Name",
    "",
    "ThermostatSetpoint:SingleHeating, FBO_CLASSROOM_HeatSP,",
    "   FBO_HeatSP;                                     !- Setpoint Temperature Schedule Name",
    "",
    "ZoneHVAC:EquipmentConnections,",
    sprintf("   %s,                        !- Zone Name", ZONE),
    "   FBO_CLASSROOM_EquipList,                        !- Zone Conditioning Equipment List Name",
    "   ,                                               !- Zone Air Inlet Node (none: convective)",
    "   ,                                               !- Zone Air Exhaust Node",
    "   FBO_CLASSROOM_ZoneAirNode,                      !- Zone Air Node Name",
    "   ;                                               !- Zone Return Air Node (none)",
    "",
    "ZoneHVAC:EquipmentList, FBO_CLASSROOM_EquipList,",
    "   SequentialLoad,",
    "   ZoneHVAC:Baseboard:RadiantConvective:Electric, !- Zone Equipment 1 Object Type",
    "   FBO_CLASSROOM_Baseboard,                        !- Zone Equipment 1 Name",
    "   1, 1, , ;                                       !- cooling seq, heating seq, frac scheds",
    "",
    "ZoneHVAC:Baseboard:RadiantConvective:Electric, FBO_CLASSROOM_Baseboard,",
    "   FBO_HeatAvail,                                  !- Availability Schedule Name",
    "   HeatingDesignCapacity,                          !- Heating Design Capacity Method",
    "   3000,                                           !- Heating Design Capacity {W}",
    "   ,                                               !- Capacity Per Floor Area",
    "   ,                                               !- Fraction of Autosized Capacity",
    "   1.0,                                            !- Efficiency",
    "   0.4,                                            !- Fraction Radiant",
    "   0.15,                                           !- Fraction Radiant Energy Incident on People",
    sprintf("   %s_Floor_0_0_0,   0.35,", ZONE),
    sprintf("   %s_Wall_2_0_0,    0.15,", ZONE),
    sprintf("   %s_Wall_3_0_0,    0.15,", ZONE),
    sprintf("   %s_Wall_5_0_0,    0.15,", ZONE),
    sprintf("   %s_Ceiling_1_0_0, 0.05;", ZONE)
  )
  a <- grep("ScheduleTypeLimits,\\s*On/Off", L)[1]
  if (is.na(a)) stop("anchor 'ScheduleTypeLimits, On/Off' not found in ", src)
  L <- append(L, blk, after = a)

  # --- (3) air-temperature output (diagnostics)
  L <- c(L, "",
         sprintf("Output:Variable, %s, Zone Mean Air Temperature, Timestep;", ZONE))

  writeLines(L, dst)

  # ---------------------------- verification ----------------------------------
  t <- readLines(dst, warn = FALSE); txt <- paste(t, collapse = "\n")
  chk <- function(cond, msg) cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "OK" else "FAIL", msg))
  cat("========= make_prod_idf: verification =========\n")
  # sizing off
  sc <- grep("^\\s*SimulationControl", t)[1]
  chk(grepl("No", t[sc + 1]) && grepl("No", t[sc + 2]), "zone & system sizing OFF")
  # heating stack present
  chk(grepl("ZoneHVAC:Baseboard:RadiantConvective:Electric, FBO_CLASSROOM_Baseboard", txt),
      "radiant baseboard object present")
  chk(sum(gregexpr("RadiantConvective:Electric", txt)[[1]] > 0) == 2, "baseboard: object + list ref (=2)")
  chk(grepl("Until: 24:00, 5.0", txt), "seasonal OFF (May-Oct, setpoint 5)")
  chk(grepl("Until: 07:00, 20.0", txt), "setback 20 present")
  chk(grepl("Zone Mean Air Temperature", txt), "air-temperature output present")
  # radiant fractions sum (0.35+0.15+0.15+0.15+0.05) + people 0.15 = 1.0
  chk(abs((0.35 + 0.15 + 0.15 + 0.15 + 0.05) + 0.15 - 1.0) < 1e-9, "radiant fractions + people = 1.0")
  # control INHERITED intact
  chk(length(grep("EnergyManagementSystem:Actuator", t)) == 9, "9 EMS actuators intact")
  for (s in c("x_tau_CO2_sched","x_delta_hys_sched","x_d_slot_sched","x_f_open_sched","x_tau_AQ_sched"))
    chk(length(grep(s, t)) >= 2, paste("x-schedule intact:", s))
  chk(length(grep("AQI_ext", t)) >= 2, "AQI_ext Schedule:File intact")
  cat(sprintf("[OK] prod written: %s  (%d lines vs %d in base)\n", dst, length(t), length(L) - length(blk) - 2))
  cat("===============================================\n")
  invisible(dst)
}

# DESIGN DECISION (2026-07): the school-calendar T_occ (exclude 22 Dec-7 Jan &
#   22 Jun-10 Sep, on top of the Jul/Aug closure already in the occupancy schedule)
#   is applied in the SCORING layer (fbo_eval.R::evaluate_U), NOT baked into this
#   IDF's occupancy Schedule:Year. Rationale: T_occ is a post-processing occupied-
#   hours filter; the E+ run uses the full lectivo-year occupancy and the calendar
#   mask is applied when the indicators are computed. Clean, defensible, and avoids
#   Schedule:Year surgery (which also drives lighting/equipment). Documented in the
#   manuscript methods (occupancy / T_occ).

if (sys.nframe() == 0L && interactive()) make_prod_idf()

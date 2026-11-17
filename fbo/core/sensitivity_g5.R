# ============================================================================
# sensitivity_g5.R  --  Sensitivity of the outdoor-import indicator g5 to the
# assumed per-species indoor loss rates k_s, over their primary-literature ranges.
#
# Purpose: show that g5 is ROBUST to the (assumed, non-calibrated) k_s within the
# literature ranges -- turning the "assumed parameters" flank into a strength.
#
# Scenario: real outdoor series (Torneo 2024, load_outdoor.R) + a representative
# DNVP ventilation regime (windows open during occupied lective hours, infiltration
# otherwise) + lective occupancy mask. Q(t) is illustrative here; the definitive
# analysis will use the EnergyPlus Q of the 9 scenarios.
#
# Literature ranges (h^-1):
#   PM2.5 xi   : 0.10 - 0.49  (central 0.20)     Chen & Zhao 2011; He et al. 2005
#   NO2   k    : 0.05 - 0.28  (central ~0.094)   Grontoft & Raychaudhuri 2004
#   O3    k    : 0.94 - 2.83  (central ~1.887)   Weschler 2000
# ============================================================================
source(here::here("core", "g5.R"))
source(here::here("core", "load_outdoor.R"))
suppressMessages(library(ggplot2))

out_dir <- "runs"; if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
V <- 126.70

# --- 1. real outdoor series + representative scenario -----------------------
paths <- list(
  PM25 = file.path(Sys.getenv("TORNEO_DIR", "."), "PM25_HH_2024_SEV_Torneo.csv"),
  NO2  = file.path(Sys.getenv("TORNEO_DIR", "."), "NO2_HH_2024_SEV_Torneo.csv"),
  O3   = file.path(Sys.getenv("TORNEO_DIR", "."), "O3_HH_2024_SEV_Torneo.csv"))
od <- load_outdoor(paths, five_min = TRUE, verbose = FALSE)
five <- od$fivemin
ts <- five$timestamp

dow   <- as.integer(format(ts, "%u"))                 # 1=Mon .. 7=Sun
hourd <- as.integer(format(ts, "%H")) + as.integer(format(ts, "%M")) / 60
mo    <- as.integer(format(ts, "%m")); dy <- as.integer(format(ts, "%d"))
summer_closed <- (mo == 7) | (mo == 8) | (mo == 6 & dy >= 22) | (mo == 9 & dy <= 10)
occ <- (dow <= 5) & !summer_closed & (hourd >= 9) & (hourd < 14)   # lective 09:00-14:00

# representative DNVP: a = 4 h^-1 when occupied (windows open), 0.64 h^-1 otherwise (infiltration)
a  <- ifelse(occ, 4.0, 0.64)
df <- data.frame(Q = a * V, Cout_PM25 = five$Cout_PM25, Cout_NO2 = five$Cout_NO2,
                 Cout_O3 = five$Cout_O3, occ = occ)

g5_of <- function(k_pm, k_no2, k_o3) {
  pr <- default_species_params()
  pr$PM25$k <- k_pm; pr$NO2$k <- k_no2; pr$O3$k <- k_o3
  compute_g5(df, V = V, params = pr)$g5
}

# --- 2. ranges and central values -----------------------------------------
central <- with(default_species_params(), c(PM25 = PM25$k, NO2 = NO2$k, O3 = O3$k))
rng <- list(PM25 = c(0.10, 0.49), NO2 = c(0.05, 0.28), O3 = c(0.94, 2.83))
g5_central <- g5_of(central["PM25"], central["NO2"], central["O3"])

# --- 3. one-at-a-time sweep -----------------------------------------------
np <- 9
sweep <- do.call(rbind, lapply(names(rng), function(sp) {
  ks <- seq(rng[[sp]][1], rng[[sp]][2], length.out = np)
  g5 <- sapply(ks, function(k) {
    args <- as.list(central); args[[sp]] <- k
    g5_of(args$PM25, args$NO2, args$O3)
  })
  data.frame(species = sp, k = ks,
             frac = (ks - rng[[sp]][1]) / (rng[[sp]][2] - rng[[sp]][1]), g5 = g5)
}))

# --- 4. tornado (g5 span induced by each parameter, others central) --------
tornado <- do.call(rbind, lapply(names(rng), function(sp) {
  s <- sweep[sweep$species == sp, ]
  data.frame(species = sp, g5_min = min(s$g5), g5_max = max(s$g5),
             span = max(s$g5) - min(s$g5))
}))

# --- 5. envelope (all-low k = least removal = worst; all-high = best) ------
g5_worst <- g5_of(rng$PM25[1], rng$NO2[1], rng$O3[1])   # low k -> high C_in -> worst
g5_best  <- g5_of(rng$PM25[2], rng$NO2[2], rng$O3[2])   # high k -> low C_in -> best

# --- 6. report -------------------------------------------------------------
cat("== g5 sensitivity to k_s (Torneo 2024, representative DNVP) ==\n")
cat(sprintf("central g5 = %.2f  (band units, 1-6)\n\n", g5_central))
cat("one-at-a-time span (others at central):\n")
for (i in seq_len(nrow(tornado))) with(tornado[i, ],
  cat(sprintf("  %-5s: g5 in [%.2f, %.2f]  -> span %.2f  (k range %.3f-%.3f)\n",
              species, g5_min, g5_max, span, rng[[species]][1], rng[[species]][2])))
cat(sprintf("\nfull envelope over all-low / all-high k:  g5 in [%.2f, %.2f]  -> span %.2f\n",
            g5_best, g5_worst, g5_worst - g5_best))

# --- 7. figure -------------------------------------------------------------
sweep$species <- factor(sweep$species, levels = c("O3", "PM25", "NO2"),
                        labels = c("O3", "PM2.5", "NO2"))
cent_pts <- data.frame(
  species = factor(c("O3","PM2.5","NO2"), levels = c("O3","PM2.5","NO2")),
  frac = c((central["O3"]-rng$O3[1])/(rng$O3[2]-rng$O3[1]),
           (central["PM25"]-rng$PM25[1])/(rng$PM25[2]-rng$PM25[1]),
           (central["NO2"]-rng$NO2[1])/(rng$NO2[2]-rng$NO2[1])),
  g5 = g5_central)

p <- ggplot(sweep, aes(frac, g5, color = species)) +
  geom_step(linewidth = 0.9, direction = "hv") +
  geom_point(data = cent_pts, aes(frac, g5), shape = 21, fill = "white", size = 2.6, stroke = 0.9) +
  scale_color_manual(values = c("O3"="#d7301f","PM2.5"="#7fbc41","NO2"="#2c7fb8"), name = NULL) +
  scale_y_continuous(breaks = 1:6, limits = c(1, 6)) +
  scale_x_continuous(labels = scales::percent) +
  labs(x = "position of ks within its literature range",
       y = "g5  (P95 worst ICA sub-index, 1-6)",
       title = "Sensitivity of g5 to the assumed loss rates ks  (Torneo 2024)",
       subtitle = "open circle = central value; flat lines = robustness to ks") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(), legend.position = c(0.9, 0.85),
        plot.title = element_text(size = 11), plot.subtitle = element_text(size = 9))

ggsave(file.path(out_dir, "T3_g5_sensitivity.png"), p, width = 6.5, height = 3.6, dpi = 200)
ggsave(file.path(out_dir, "T3_g5_sensitivity.pdf"), p, width = 6.5, height = 3.6, device = cairo_pdf)
cat(sprintf("\nfigure saved to %s\n", normalizePath(out_dir)))

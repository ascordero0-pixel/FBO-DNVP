# =============================================================================
# F1/fig_bo_utility.R -- Figure 6: optimal utility U* of the nine scenarios
#   Bars   : U* = max U over the LHS seed + BO evaluations (Table 4), full precision
#   Marks  : best U of the LHS seed alone, re-scored with the current FIS
#   Bands  : utility categories of the scoring FIS (Poor < 0.33, Acceptable
#            0.33-0.67, Good >= 0.67)
# Reads  : F1/seed_9scen_results.csv, F1/bo/bo_*.csv, core/fis_scoring.R
# Writes : runs/Fig_BO_utility_by_scenario.{pdf,png}
# USAGE  : source(here::here("F1", "fig_bo_utility.R"))     (inside FIS.Rproj)
# Port to R of the original matplotlib figure (2026-07-31); same data and layout. No title or
# footnote (the manuscript caption carries them); aspect = the figure frame in the manuscript.
# =============================================================================
library(here); library(ggplot2)
source(here("core", "fis_scoring.R"))

seed <- read.csv(here("F1", "seed_9scen_results.csv"))
bo   <- do.call(rbind, lapply(list.files(here("F1", "bo"), "^bo_.*\\.csv$", full.names = TRUE), read.csv))
cols <- c("scn", "g1", "g2_pII", "g2_pIII", "g3", "g4", "g5")
sf   <- build_scoring_fis(default_variables())
rescore <- function(x) vapply(seq_len(nrow(x)), function(i) scoring_utility(sf,
  g1 = x$g1[i], g2 = c(x$g2_pII[i], x$g2_pIII[i]), g3 = x$g3[i], g4 = x$g4[i], g5 = x$g5[i]), numeric(1))
seed$U_now <- rescore(seed)
all  <- rbind(seed[cols], bo[cols]); all$U <- rescore(all)
# U* at full precision (labels rounded once, to 3 decimals); seed-only best for the marks
d <- merge(aggregate(U ~ scn, all, max), aggregate(U_now ~ scn, seed, max), by = "scn")
d$orient <- factor(sub("_.*", "", d$scn), levels = c("NE", "N", "E"))
d$ap     <- factor(sub(".*_", "", d$scn), levels = c("w154", "w216", "w410"))
d$x      <- as.integer(d$ap) + (as.integer(d$orient) - 2) * 0.25
cat("U* range / max orientation spread =",
    round(diff(range(d$U)) / max(tapply(d$U, d$ap, function(u) diff(range(u)))), 2), "\n")

ORI <- c(NE = "#0072B2", N = "#E69F00", E = "#009E73")
p <- ggplot(d) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0,    ymax = 0.33, fill = "#F7E1E1", alpha = 0.6) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.33, ymax = 0.67, fill = "#FBF1D9", alpha = 0.7) +
  annotate("rect", xmin = -Inf, xmax = Inf, ymin = 0.67, ymax = 0.72, fill = "#DDEFE3", alpha = 0.8) +
  geom_hline(yintercept = c(0.33, 0.67), linetype = "dashed", colour = "grey55", linewidth = 0.4) +
  geom_col(aes(x = x, y = U, fill = orient), width = 0.245, colour = "white", linewidth = 0.6) +
  geom_segment(aes(x = x - 0.11, xend = x + 0.11, y = U_now, yend = U_now, linetype = "LHS seed best"),
               colour = "grey15", linewidth = 0.7) +
  geom_text(aes(x = x, y = U + 0.012, label = sprintf("%.3f", U)), size = 3.6, colour = "grey15", vjust = 0) +
  annotate("text", x = 3.5, y = c(0.695, 0.495, 0.165), label = c("Good", "Acceptable", "Poor"),
           hjust = 0, fontface = "italic", size = 4.2, colour = c("#3D6B4A", "#8A6D1F", "#A44A4A")) +
  scale_fill_manual(values = ORI, labels = paste("North axis:", names(ORI)), name = NULL) +
  scale_linetype_manual(values = c("LHS seed best" = "solid"), name = NULL) +
  scale_x_continuous(breaks = 1:3, labels = c("1.54 %\n(small)", "2.16 %\n(medium)", "4.10 %\n(large)"),
                     limits = c(0.55, 4.05), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.72), breaks = seq(0, 0.7, 0.1), expand = c(0, 0)) +
  labs(x = expression("Window aperture size  (WWeqR"[max]*")"), y = expression("Optimal utility  "*italic(U)^"*")) +
  guides(fill = guide_legend(order = 1, ncol = 2, byrow = FALSE), linetype = guide_legend(order = 2)) +
  theme_classic(base_size = 13) +
  theme(axis.text = element_text(colour = "black"),
        panel.grid.major.y = element_line(colour = "#D5DEE8", linewidth = 0.3),
        legend.position = c(0.015, 0.985), legend.justification = c(0, 1),
        legend.box = "horizontal", legend.box.just = "top", legend.spacing.x = unit(0.2, "cm"),
        legend.background = element_blank(), legend.key = element_blank(),
        legend.key.size = unit(0.45, "cm"), legend.text = element_text(size = 11))

dir.create(here("runs"), showWarnings = FALSE)
ggsave(here("runs", "Fig_BO_utility_by_scenario.pdf"), p, width = 8.35, height = 4.74, device = cairo_pdf)
ggsave(here("runs", "Fig_BO_utility_by_scenario.png"), p, width = 8.35, height = 4.74, dpi = 300)
print(d[order(d$ap, d$orient), c("scn", "U", "U_now")], row.names = FALSE)

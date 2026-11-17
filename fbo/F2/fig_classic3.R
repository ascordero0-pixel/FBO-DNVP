# =============================================================================
# F2/fig_classic3.R -- Figure 10 (v2): optimized FBO protocol vs TWO fixed-threshold
#   baselines (1000 ppm, field default; 800 ppm, threshold of T1 [11]).
#   Two-panel dumbbell: (a) composite utility U, (b) control effort g4.
#   Reads: F2/out/classic_vs_opt.csv (opt + cls1000) and F2/out/classic800_vs_opt.csv.
#   Writes: F2/fig/Fig_classic3.pdf / .png  (does not touch Fig_classic.*)
# USAGE: source(here::here("F2", "fig_classic3.R"))
# =============================================================================
library(here); library(ggplot2)

a <- read.csv(here("F2", "out", "classic_vs_opt.csv"))
b <- read.csv(here("F2", "out", "classic800_vs_opt.csv"))
d <- merge(a[, c("scn", "U_opt", "U_cls", "g4_opt", "g4_cls")],
           b[, c("scn", "U_cls800", "g4_cls800")], by = "scn")

ord <- c("E_w410","N_w410","NE_w410","E_w216","N_w216","NE_w216","E_w154","N_w154","NE_w154")  # bottom -> top
d$scn <- factor(d$scn, levels = ord)

long <- rbind(
  data.frame(scn = d$scn, panel = "(a) Composite utility U",  proto = "Optimized (FBO)",      x = d$U_opt),
  data.frame(scn = d$scn, panel = "(a) Composite utility U",  proto = "Baseline 1000 ppm",    x = d$U_cls),
  data.frame(scn = d$scn, panel = "(a) Composite utility U",  proto = "Baseline 800 ppm [11]",x = d$U_cls800),
  data.frame(scn = d$scn, panel = "(b) Control effort g4",    proto = "Optimized (FBO)",      x = d$g4_opt),
  data.frame(scn = d$scn, panel = "(b) Control effort g4",    proto = "Baseline 1000 ppm",    x = d$g4_cls),
  data.frame(scn = d$scn, panel = "(b) Control effort g4",    proto = "Baseline 800 ppm [11]",x = d$g4_cls800))
long$proto <- factor(long$proto, levels = c("Baseline 800 ppm [11]", "Baseline 1000 ppm", "Optimized (FBO)"))
seg <- aggregate(x ~ scn + panel, long, function(v) c(lo = min(v), hi = max(v)))
seg <- data.frame(scn = seg$scn, panel = seg$panel, lo = seg$x[, "lo"], hi = seg$x[, "hi"])
band <- data.frame(scn = factor(ord, levels = ord), y = seq_along(ord))
vline <- data.frame(panel = "(a) Composite utility U", x = 0.33)   # Poor | Acceptable boundary of the FIS

p <- ggplot() +
  geom_blank(data = long, aes(x = x, y = scn)) +   # [repo patch P6] discrete y first
  geom_rect(data = band, aes(xmin = -Inf, xmax = Inf, ymin = y - 0.38, ymax = y + 0.38), fill = "#EEF1F3") +
  geom_vline(data = vline, aes(xintercept = x), colour = "#B23A48", linetype = "dashed", linewidth = 0.4) +
  geom_segment(data = seg, aes(x = lo, xend = hi, y = scn, yend = scn), colour = "#BFC5CA", linewidth = 0.9) +
  geom_point(data = long, aes(x = x, y = scn, fill = proto, colour = proto, shape = proto), size = 2.8, stroke = 0.9) +
  scale_shape_manual(values = c("Baseline 800 ppm [11]" = 21, "Baseline 1000 ppm" = 19, "Optimized (FBO)" = 19), name = NULL) +
  scale_colour_manual(values = c("Baseline 800 ppm [11]" = "#8C8C8C", "Baseline 1000 ppm" = "#8C8C8C", "Optimized (FBO)" = "#E08A1E"), name = NULL) +
  scale_fill_manual(values = c("Baseline 800 ppm [11]" = "white", "Baseline 1000 ppm" = "#8C8C8C", "Optimized (FBO)" = "#E08A1E"), name = NULL) +
  facet_wrap(~ panel, scales = "free_x") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(), strip.text = element_text(face = "bold", hjust = 0),
        legend.position = "bottom", axis.ticks.x = element_line(colour = "grey60"),
        panel.spacing = unit(1.2, "lines"))
# Panel-specific x labels (ggplot facets share one xlab; add as strip-independent annotation)
p <- p + geom_text(data = data.frame(panel = c("(a) Composite utility U", "(b) Control effort g4"),
                                     lab = c("Utility U", "Control effort (actions / session)"),
                                     x = c(0.365, 16.5), y = c(0.1, 0.1)),
                   aes(x = x, y = y, label = lab), size = 3.1, colour = "grey30", inherit.aes = FALSE)

dir.create(here("F2", "fig"), showWarnings = FALSE)
ggsave(here("F2", "fig", "Fig_classic3.pdf"), p, width = 10.5, height = 4.8, units = "in", device = cairo_pdf)
ggsave(here("F2", "fig", "Fig_classic3.png"), p, width = 10.5, height = 4.8, units = "in", dpi = 300)
cat("[fig_classic3] written F2/fig/Fig_classic3.{pdf,png}\n")
print(d)

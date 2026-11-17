# ============================================================================
# plot_mf_panel.R  --  Section 2.7 figure: membership functions of the three
# PRIMARY indicators (g1, g2, g5) as a single 3-panel figure. Subscripts are
# rendered with plotmath (locale/font-independent). Knots read from fis_scoring.R.
# ============================================================================
source(here::here("core", "fis_scoring.R"))
suppressMessages(library(ggplot2))
out_dir <- "runs"; if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

V <- default_variables(); RHO <- 0.5
cols <- c(Low = "#2c7fb8", Medium = "#7fbc41", High = "#d7301f")
tlab <- c(L = "Low", M = "Medium", H = "High")

prim_df <- function(var, panel, disp_hi) {
  grid <- seq(var$lo, disp_hi, length.out = 800); cur <- var_term_curves(var, grid)
  do.call(rbind, lapply(names(cur), function(t)
    data.frame(panel = panel, x = grid, mu = cur[[t]], term = tlab[[t]], row.names = NULL)))
}
df1 <- prim_df(V$g1, "g1", 1000)
df5 <- prim_df(V$g5, "g5", 5)
p   <- seq(0, 25, length.out = 800)   # g2 axis: covers the achievable p_II range (~12-24%) with the re-anchored [5,10] ramp
mem <- sapply(p, function(pi) var_memberships(V$g2, c(pi, RHO * pi)))
df2 <- do.call(rbind, lapply(c("L","M","H"), function(t)
  data.frame(panel = "g2", x = p, mu = mem[t, ], term = tlab[[t]], row.names = NULL)))

# --- plotmath strip labels (real subscripts, no Unicode dependency) --------
key <- c("g1","g2","g5")
plot_lab <- c(
  g1 = "(a)~~italic(g)[1]*\":  IAQ, CO\"[2]*\" excess (ppm)\"",
  g2 = "(b)~~italic(g)[2]*\":  thermal comfort (% out of Cat II)\"",
  g5 = "(c)~~italic(g)[5]*\":  outdoor import (ICA / EAQI index)\"")
relabel <- function(df) { df$panel <- factor(df$panel, levels = key, labels = plot_lab[key]); df }

curves <- relabel(rbind(df1, df2, df5))
curves$term <- factor(curves$term, levels = c("Low","Medium","High"))

bounds <- relabel(rbind(
  data.frame(panel="g1", x=c(550,800), lab=c("550","800")),
  data.frame(panel="g2", x=c(5,10),    lab=c("5%","10%")),   # RE-ANCHORED: 5% = EN Cat III limit, 10% = worst achievable
  data.frame(panel="g5", x=c(2.5,3.5), lab=c("2.5","3.5"))))
cats <- relabel(rbind(
  data.frame(panel="g1", x=c(275,675,900), lab=c("Cat I","Cat II","Cat III+")),
  data.frame(panel="g5", x=c(1.75,3,4.25), lab=c("Good-Fair","Moderate","Poor +"))))

p_out <- ggplot(curves, aes(x, mu, color = term, fill = term)) +
  geom_area(aes(group = term), position = "identity", alpha = 0.08, color = NA) +
  geom_line(linewidth = 0.8) +
  geom_vline(data = bounds, aes(xintercept = x), linetype = "dashed",
             color = "grey35", linewidth = 0.35) +
  geom_text(data = bounds, aes(x = x, y = 1.07, label = lab), inherit.aes = FALSE,
            size = 2.6, color = "grey25") +
  geom_text(data = cats, aes(x = x, y = 0.06, label = lab), inherit.aes = FALSE,
            size = 2.3, color = "grey45", fontface = "italic") +
  facet_wrap(~panel, scales = "free_x", labeller = label_parsed) +
  scale_color_manual(values = cols, name = NULL) +
  scale_fill_manual(values = cols, guide = "none") +
  scale_y_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1.16)) +
  labs(x = NULL, y = expression("membership  " * mu),
       caption = "Membership functions of the primary indicators (g1, g5 normatively anchored; g2 fail-ramp anchored at the EN 16798-1 Cat III limit, 5%, and graded to 10%). g2 shown at p_III = 0.5 x p_II.") +
  theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "bottom",
        legend.key.size = unit(0.9, "lines"),
        strip.text = element_text(size = 8.6),
        plot.caption = element_text(size = 7.6, hjust = 0, color = "grey35"))

ggsave(file.path(out_dir, "T3_Fig_MF_primaries_3panel.png"), p_out, width = 11, height = 3.1, dpi = 200)
ggsave(file.path(out_dir, "T3_Fig_MF_primaries_3panel.pdf"), p_out, width = 11, height = 3.1, device = cairo_pdf)
cat("saved 3-panel figure (plotmath subscripts) to", normalizePath(out_dir), "\n")

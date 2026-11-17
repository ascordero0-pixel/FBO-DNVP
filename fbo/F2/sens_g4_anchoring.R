# =============================================================================
# F2/sens_g4_anchoring.R -- Sensitivity of the published optima to the anchoring
#   of the effort indicator g4 (Table B7). Re-scores (NO simulation) every
#   evaluated point of the seed (450) and the BO logs (270) plus the 27 protocols
#   of Table 10 under three g4 anchorings, and reports for each scenario:
#   U(x*), best U among evaluated points, dU, and dominance over both baselines.
#   Pre-registered reading: dU <= 0.01 (eps of Section 3.3) and opt >= both
#   baselines in all 9 scenarios => ranking robust to the anchoring.
# USAGE: source(here::here("F2", "sens_g4_anchoring.R"))
# =============================================================================
library(here)
source(here("core", "fis_scoring.R"))

SCN <- c("NE_w154","N_w154","E_w154","NE_w216","N_w216","E_w216","NE_w410","N_w410","E_w410")
seed <- read.csv(here("F1", "seed_9scen_results.csv")); seed$src <- "seed"
bo   <- do.call(rbind, lapply(SCN, function(s) { d <- read.csv(here("F1", "bo", paste0("bo_", s, ".csv"))); d$src <- "bo"; d }))
cols <- c("scn","tag","U","g1","g2_pII","g2_pIII","g3","g4","g5","src")
E <- rbind(seed[, cols], bo[, cols])
c1 <- read.csv(here("F2", "out", "classic_vs_opt.csv"))
c8 <- read.csv(here("F2", "out", "classic800_vs_opt.csv"))

ANCH <- list("A0 (6.44/7.44/8.44, published)" = c(6.44, 8.44),
             "A1 (5/10/15, 1-3 actions per hour)" = c(5, 15),
             "A2 (10/15/20, stress)"              = c(10, 20))

fis_with_g4 <- function(b) { v <- default_variables(); v$g4 <- tercile_partition("g4", 0, 60, b[1], b[2], "actions/session"); build_scoring_fis(v) }
U_of <- function(sf, g1, pII, pIII, g3, g4, g5 = 3) scoring_utility(sf, g1, c(pII, pIII), g3, g4, g5)

res <- list()
for (an in names(ANCH)) {
  sf <- fis_with_g4(ANCH[[an]])
  E$U_new <- mapply(function(a,b,c,d,e,f) U_of(sf,a,b,c,d,e,f), E$g1, E$g2_pII, E$g2_pIII, E$g3, E$g4, E$g5)
  for (s in SCN) {
    d <- E[E$scn == s, ]; o <- c1[c1$scn == s, ]; o8 <- c8[c8$scn == s, ]
    Uopt <- U_of(sf, o$g1_opt, o$pII_opt, o$pIII_opt, o$g3_opt, o$g4_opt)
    Uc1  <- U_of(sf, o$g1_cls, o$pII_cls, o$pIII_cls, o$g3_cls, o$g4_cls)
    Uc8  <- U_of(sf, o8$g1_cls800, o8$pII_cls800, o8$pIII_cls800, o8$g3_cls800, o8$g4_cls800)
    b <- d[which.max(d$U_new), ]
    res[[length(res)+1]] <- data.frame(anchor = an, scn = s, U_opt = round(Uopt,4), U_best = round(b$U_new,4),
      dU = round(b$U_new - Uopt, 4), best_tag = b$tag, best_g4 = round(b$g4,1), best_g1 = round(b$g1),
      U_cls1000 = round(Uc1,4), U_cls800 = round(Uc8,4),
      opt_ge_both = (Uopt >= Uc1 - 1e-4) && (Uopt >= Uc8 - 1e-4))
  }
}
R <- do.call(rbind, res)
dir.create(here("F2", "out"), showWarnings = FALSE)
write.csv(R, here("F2", "out", "sens_g4_anchoring.csv"), row.names = FALSE)
for (an in names(ANCH)) { cat("\n===", an, "\n"); print(R[R$anchor == an, -1], row.names = FALSE) }
cat(sprintf("\nA1: max dU = %.4f | dominance kept in %d/9\n",
            max(R$dU[grepl("A1", R$anchor)]), sum(R$opt_ge_both[grepl("A1", R$anchor)])))

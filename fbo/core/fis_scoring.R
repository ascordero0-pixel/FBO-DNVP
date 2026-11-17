# ============================================================================
# fis_scoring.R  --  T3 FBO framework, Block A: the SCORING Fuzzy Inference System
#
# Native-R port (base R only, no external packages) of fis_scoring.py.
# Hierarchical Mamdani FIS (max-min inference, centroid defuzzification) mapping
# the five annual performance indicators g1..g5 to a scalar utility U in [0,1].
#
#   Level 1  sub-FIS Cost : (g3 overcooling, g4 effort) -> Cost {Low,Med,High}   (9 rules)
#   Level 2  main FIS     : (g1 CO2, g2 comfort, g5 outdoor-AQ, Cost) -> U  (23 rules)
#
# Anchoring (Table 6):
#   * Primaries g1, g5 : NORMATIVE single-axis (EN 16798-1 CO2 categories; ICA/EAQI).
#   * Primary   g2     : NORMATIVE nested-category ladder (Cat II +/-3, Cat III +/-4),
#                        input pair (p_II, p_III), fail-ramp RE-ANCHORED 5%..10%
#                        (5% = EN Cat III acceptability limit; see comfort_ladder_var).
#   * Cost pair g3, g4 : TERCILE placeholders until the post-LHS distribution exists
#                        (g4 analytical cap N_max = 300/5 = 60 actions/session).
#   * Output U vertices: Poor 0.20, Acceptable 0.50, Good 0.85.
#
# Usage from RStudio:
#   source("fis_scoring.R")
#   sf <- build_scoring_fis()
#   scoring_utility(sf, g1=200, g2=c(2,1), g3=50, g4=10, g5=1.5)   # -> U in [0,1]
# Running `Rscript fis_scoring.R` executes the validation harness (105 checks).
# ============================================================================


# ---------------------------------------------------------------------------
# 1. Membership-function primitives (vectorised over x)
# ---------------------------------------------------------------------------
trapmf <- function(x, a, b, c, d) {
  y <- numeric(length(x))
  if (b > a) { m <- x > a & x < b; y[m] <- (x[m] - a) / (b - a) }
  y[x >= b & x <= c] <- 1
  if (d > c) { m <- x > c & x < d; y[m] <- (d - x[m]) / (d - c) }
  if (a == b) y[x <= a] <- 1
  if (c == d) y[x >= d] <- 1
  y
}

trimf <- function(x, a, b, c) {
  y <- numeric(length(x))
  if (b > a) { m <- x > a & x < b; y[m] <- (x[m] - a) / (b - a) }
  y[x == b] <- 1
  if (c > b) { m <- x > b & x < c; y[m] <- (c - x[m]) / (c - b) }
  y
}


# ---------------------------------------------------------------------------
# 2. Linguistic variables (lists of term closures) and dispatch
# ---------------------------------------------------------------------------
new_var <- function(name, lo, hi, terms, unit = "") {
  list(kind = "var", name = name, lo = lo, hi = hi, unit = unit,
       term_names = names(terms), terms = terms,
       reps = sapply(terms, function(t) t$rep))
}

normative_primary <- function(name, lo, hi, b1, b2, unit = "") {
  terms <- list(
    L = list(fn = function(x) trapmf(x, lo, lo, b1, b2),      rep = 0.5 * (lo + b1)),
    M = list(fn = function(x) trimf(x, b1, 0.5 * (b1 + b2), b2), rep = 0.5 * (b1 + b2)),
    H = list(fn = function(x) trapmf(x, b1, b2, hi, hi),      rep = 0.5 * (b2 + hi))
  )
  new_var(name, lo, hi, terms, unit)
}

# tercile placeholder = same geometry with cut points t1<t2 (injected post-LHS)
tercile_partition <- function(name, lo, hi, t1, t2, unit = "") {
  v <- normative_primary(name, lo, hi, t1, t2, unit)
  v$placeholder <- TRUE
  v
}

output_utility <- function() {
  terms <- list(
    Poor = list(fn = function(x) trapmf(x, 0, 0, 0.20, 0.50), rep = 0.10),
    Acc  = list(fn = function(x) trimf(x, 0.20, 0.50, 0.85),   rep = 0.50),
    Good = list(fn = function(x) trapmf(x, 0.50, 0.85, 1, 1),  rep = 0.925)
  )
  new_var("U", 0, 1, terms)
}

# g2 comfort indicator: EN 16798-1 nested-category ladder (Option 2).
# Crisp input is the PAIR (p_II, p_III) = % occupied hours outside Cat II / Cat III.
# fail-ramp RE-ANCHORED [5% -> 10%] (T3, case-study specific): the a_strict end
# is the EN 16798-1 Cat III long-term acceptability limit (5%), taken as the
# "comfortable" reference, and comfort degrades over 5%..10% up to the worst
# control achievable in this passive 1946 Csa envelope (which cannot meet Cat III
# <5% by ventilation alone). This keeps a NORMATIVE anchor (the 5% EN limit) while
# letting g2 rank controls for the optimizer -- consistent with the data-adaptive
# (tercile) anchoring already used for the g3/g4 cost pair.
comfort_ladder_var <- function(name = "g2", a_strict = 5, a_std = 10,
                               axis_hi = 30, unit = "% occ. hours") {
  fail <- function(p) max(0, min(1, (p - a_strict) / (a_std - a_strict)))
  list(kind = "ladder", name = name, lo = 0, hi = axis_hi, unit = unit,
       a_strict = a_strict, a_std = a_std, term_names = c("L", "M", "H"),
       reps = list(L = c(3, 1), M = c(12, 4), H = c(20, 12)), fail = fail)
}

# memberships(): dispatch on variable kind. val is a scalar (var) or c(p_II,p_III) (ladder)
var_memberships <- function(var, val) {
  if (!is.null(var$kind) && var$kind == "ladder") {
    fII <- var$fail(val[1]); fIII <- var$fail(val[2])
    return(c(L = 1 - fII, M = fII * (1 - fIII), H = fIII))
  }
  vapply(var$terms, function(t) t$fn(val), numeric(1))
}

var_term_curves <- function(var, grid) lapply(var$terms, function(t) t$fn(grid))


# ---------------------------------------------------------------------------
# 3. Generic single-level Mamdani engine
# ---------------------------------------------------------------------------
new_rule <- function(antecedents, consequent, tag = "")
  list(antecedents = antecedents, consequent = consequent, tag = tag)

new_fis <- function(inputs, output, rules, n_grid = 1001) {
  grid <- seq(output$lo, output$hi, length.out = n_grid)
  list(inputs = inputs, output = output, rules = rules,
       grid = grid, curves = var_term_curves(output, grid))
}

.fire <- function(fis, crisp) {
  # returns named list: mu (per input) already computed
  mu <- lapply(names(fis$inputs), function(n) var_memberships(fis$inputs[[n]], crisp[[n]]))
  names(mu) <- names(fis$inputs)
  mu
}

fis_infer <- function(fis, crisp) {
  mu <- .fire(fis, crisp)
  agg <- numeric(length(fis$grid))
  for (r in fis$rules) {
    strengths <- numeric(0)
    for (vn in names(r$antecedents))
      strengths <- c(strengths, max(mu[[vn]][r$antecedents[[vn]]]))
    w <- if (length(strengths)) min(strengths) else 1
    if (w <= 0) next
    agg <- pmax(agg, pmin(w, fis$curves[[r$consequent]]))
  }
  denom <- sum(agg)
  if (denom <= 1e-12) return(0.5 * (fis$output$lo + fis$output$hi))
  sum(fis$grid * agg) / denom
}

fis_rule_activations <- function(fis, crisp) {
  mu <- .fire(fis, crisp)
  acts <- list()
  for (r in fis$rules) {
    strengths <- numeric(0)
    for (vn in names(r$antecedents))
      strengths <- c(strengths, max(mu[[vn]][r$antecedents[[vn]]]))
    w <- if (length(strengths)) min(strengths) else 1
    prev <- acts[[r$consequent]]; if (is.null(prev)) prev <- 0
    acts[[r$consequent]] <- max(prev, w)
  }
  acts
}

fis_dominant <- function(fis, crisp) {
  acts <- fis_rule_activations(fis, crisp)
  names(acts)[which.max(unlist(acts))]
}


# ---------------------------------------------------------------------------
# 4. Rule bases (verbatim from T3_FIS_reglas.xlsx)
# ---------------------------------------------------------------------------
build_cost_subfis <- function(g3, g4, cost_out) {
  R <- new_rule
  rules <- list(
    R(list(g3 = "L", g4 = "L"), "L", "SC1"),
    R(list(g3 = "L", g4 = "M"), "L", "SC2"),
    R(list(g3 = "L", g4 = "H"), "M", "SC3"),
    R(list(g3 = "M", g4 = "L"), "L", "SC4"),
    R(list(g3 = "M", g4 = "M"), "M", "SC5"),
    R(list(g3 = "M", g4 = "H"), "H", "SC6"),
    R(list(g3 = "H", g4 = "L"), "H", "SC7"),   # CORREGIDO 2026-08: g3 Alto -> Coste Alto (la energia no se enmascara con g4 bajo)
    R(list(g3 = "H", g4 = "M"), "H", "SC8"),
    R(list(g3 = "H", g4 = "H"), "H", "SC9")
  )
  new_fis(list(g3 = g3, g4 = g4), cost_out, rules)
}

build_main_fis <- function(g1, g2, g5, cost_in, U) {
  R <- new_rule
  rules <- list(
    R(list(g2 = "H"), "Poor", "R1"),
    R(list(g1 = "H", g5 = "H"), "Poor", "R2"),
    R(list(g1 = "H", g2 = "L", g5 = "L"), "Acc", "R3"),
    R(list(g1 = "H", g2 = "L", g5 = "M"), "Acc", "R4"),
    R(list(g1 = "H", g2 = "M", g5 = "L"), "Acc", "R5"),
    R(list(g1 = "H", g2 = "M", g5 = "M"), "Acc", "R6"),
    R(list(g1 = "L", g2 = "L", g5 = "H"), "Acc", "R7"),
    R(list(g1 = "M", g2 = "L", g5 = "H"), "Acc", "R8"),
    R(list(g1 = "L", g2 = "M", g5 = "H"), "Acc", "R9"),
    R(list(g1 = "M", g2 = "M", g5 = "H"), "Acc", "R10"),
    R(list(g1 = "L", g2 = "L", g5 = "L", Cost = "L"), "Good", "R11"),
    R(list(g1 = "L", g2 = "L", g5 = "L", Cost = "M"), "Good", "R12"),
    R(list(g1 = "L", g2 = "L", g5 = "L", Cost = "H"), "Acc", "R13"),
    R(list(g1 = "M", g2 = "L", g5 = "L", Cost = "L"), "Good", "R14"),
    R(list(g1 = "M", g2 = "L", g5 = "L", Cost = c("M", "H")), "Acc", "R15"),
    R(list(g1 = "L", g2 = "M", g5 = "L", Cost = "L"), "Good", "R16"),
    R(list(g1 = "L", g2 = "M", g5 = "L", Cost = c("M", "H")), "Acc", "R17"),
    R(list(g1 = "L", g2 = "L", g5 = "M", Cost = "L"), "Good", "R18"),
    R(list(g1 = "L", g2 = "L", g5 = "M", Cost = c("M", "H")), "Acc", "R19"),
    R(list(g1 = "M", g2 = "M", g5 = "L"), "Acc", "R20"),
    R(list(g1 = "M", g2 = "L", g5 = "M"), "Acc", "R21"),
    R(list(g1 = "L", g2 = "M", g5 = "M"), "Acc", "R22"),
    R(list(g1 = "M", g2 = "M", g5 = "M"), "Acc", "R23")
  )
  new_fis(list(g1 = g1, g2 = g2, g5 = g5, Cost = cost_in), U, rules)
}


# ---------------------------------------------------------------------------
# 5. Default variable set + hierarchical scoring FIS
# ---------------------------------------------------------------------------
default_variables <- function() {
  list(
    g1   = normative_primary("g1", 0, 1350, 550, 800, "ppm CO2 excess"),
    g2   = comfort_ladder_var("g2"),                                  # nested Cat II/III ladder
    g5   = normative_primary("g5", 1, 6, 2.5, 3.5, "ICA/EAQI index"),
    g3   = tercile_partition("g3", 0, 300, 5.40, 12.02, "degC.min/day"),   # terciles post-LHS (seed, global)
    g4   = tercile_partition("g4", 0,  60, 6.44,  8.44, "actions/session"), # terciles post-LHS (seed, global)
    Cost = normative_primary("Cost", 0, 1, 1/3, 2/3, "adim."),
    U    = output_utility()
  )
}

build_scoring_fis <- function(vars = default_variables()) {
  list(vars = vars,
       sub  = build_cost_subfis(vars$g3, vars$g4, vars$Cost),
       main = build_main_fis(vars$g1, vars$g2, vars$g5, vars$Cost, vars$U))
}

scoring_cost <- function(sf, g3, g4) fis_infer(sf$sub, list(g3 = g3, g4 = g4))

scoring_utility <- function(sf, g1, g2, g3, g4, g5, return_cost = FALSE) {
  cst <- scoring_cost(sf, g3, g4)
  u <- fis_infer(sf$main, list(g1 = g1, g2 = g2, g5 = g5, Cost = cst))
  if (return_cost) list(U = u, Cost = cst) else u
}

scoring_dominant_block <- function(sf, g1, g2, g3, g4, g5) {
  cst <- scoring_cost(sf, g3, g4)
  fis_dominant(sf$main, list(g1 = g1, g2 = g2, g5 = g5, Cost = cst))
}


# ---------------------------------------------------------------------------
# 6. Validation harness
# ---------------------------------------------------------------------------
.classify_U <- function(u) if (u < 0.32) "Poor" else if (u > 0.70) "Good" else "Acc"

run_tests <- function() {
  sf <- build_scoring_fis(); v <- sf$vars
  env <- new.env(); env$ok <- 0L; env$fail <- 0L
  check <- function(cond, msg) {
    if (isTRUE(cond)) env$ok <- env$ok + 1L
    else { env$fail <- env$fail + 1L; cat("  FAIL:", msg, "\n") }
  }
  rep <- lapply(c("g1", "g2", "g3", "g4", "g5"), function(n) v[[n]]$reps)
  names(rep) <- c("g1", "g2", "g3", "g4", "g5")

  cat("[N] normative anchors sanity (g1,g5 single-axis; g2 ladder)\n")
  for (nm in c("g1", "g5")) {
    b <- if (nm == "g1") c(550, 800) else c(2.5, 3.5)
    mb1 <- var_memberships(v[[nm]], b[1]); mb2 <- var_memberships(v[[nm]], b[2])
    mc  <- var_memberships(v[[nm]], mean(b))
    check(abs(mb1["L"] - 1) < 1e-9 && mb1["M"] < 1e-9 && mb1["H"] < 1e-9,
          paste(nm, "b1 not pure Low"))
    check(abs(mb2["H"] - 1) < 1e-9 && mb2["M"] < 1e-9 && mb2["L"] < 1e-9,
          paste(nm, "b2 not pure High"))
    check(abs(mc["M"] - 1) < 1e-9 && abs(mc["L"] - 0.5) < 1e-9 && abs(mc["H"] - 0.5) < 1e-9,
          paste(nm, "centre not (0.5,1,0.5)"))
  }
  g2v <- v$g2
  for (nm in c("L", "M", "H")) {
    pair <- g2v$reps[[nm]]; m <- var_memberships(g2v, pair)
    check(abs(m[nm] - 1) < 1e-9 && (sum(m) - m[nm]) < 1e-9,
          paste("g2 ladder", nm, "not pure"))
  }
  check(var_memberships(g2v, c(8, 4))["M"] >= 0.49, "g2 (8,4) not Medium-leaning")

  cat("[T1] sub-FIS Cost matrix\n")
  mat <- list("L,L"="L","L,M"="L","L,H"="M","M,L"="M","M,M"="M",
              "M,H"="H","H,L"="H","H,M"="H","H,H"="H")   # (H,L) -> H tras la correccion de SC7
  for (k in names(mat)) {
    ab <- strsplit(k, ",")[[1]]
    cc <- scoring_cost(sf, rep$g3[[ab[1]]], rep$g4[[ab[2]]])
    lab <- if (cc < 0.32) "L" else if (cc > 0.68) "H" else "M"
    check(lab == mat[[k]], paste0("Cost(", k, ")=", round(cc, 3), " -> ", lab))
  }
  dgL <- scoring_cost(sf, rep$g3[["L"]], rep$g4[["L"]])
  dgM <- scoring_cost(sf, rep$g3[["M"]], rep$g4[["M"]])
  dgH <- scoring_cost(sf, rep$g3[["H"]], rep$g4[["H"]])
  check(dgL < dgM && dgM < dgH, "Cost not monotone on diagonal")

  cat("[T2] trap R2 -> Poor\n")
  check(.classify_U(scoring_utility(sf, rep$g1[["H"]], rep$g2[["L"]], rep$g3[["L"]],
                                    rep$g4[["L"]], rep$g5[["H"]])) == "Poor", "trap not Poor")
  cat("[T3] optimum R11 -> Good\n")
  uopt <- scoring_utility(sf, rep$g1[["L"]], rep$g2[["L"]], rep$g3[["L"]], rep$g4[["L"]], rep$g5[["L"]])
  check(.classify_U(uopt) == "Good", paste("optimum U =", round(uopt, 3)))
  cat("[T4] dominance R1 (g2 High) -> Poor\n")
  check(.classify_U(scoring_utility(sf, rep$g1[["L"]], rep$g2[["H"]], rep$g3[["L"]],
                                    rep$g4[["L"]], rep$g5[["L"]])) == "Poor", "g2-High not Poor")
  cat("[T5] Cost tie-break (#M=0)\n")
  uch <- scoring_utility(sf, rep$g1[["L"]], rep$g2[["L"]], rep$g3[["L"]], rep$g4[["L"]], rep$g5[["L"]])
  upr <- scoring_utility(sf, rep$g1[["L"]], rep$g2[["L"]], rep$g3[["H"]], rep$g4[["H"]], rep$g5[["L"]])
  check(uch > upr + 0.10 && .classify_U(uch) == "Good" && .classify_U(upr) == "Acc",
        paste0("tie-break cheap=", round(uch,3), " pricey=", round(upr,3)))

  cat("[T6] completeness over 27 x 3 corners (dominant-consequent check)\n")
  expected <- function(g1t, g2t, g5t, ct) {
    if (g2t == "H") return("Poor")
    if (g1t == "H" && g5t == "H") return("Poor")
    if (g1t == "H" || g5t == "H") return("Acc")
    nM <- sum(c(g1t, g2t, g5t) == "M")
    if (nM == 0) return(if (ct %in% c("L","M")) "Good" else "Acc")
    if (nM == 1) return(if (ct == "L") "Good" else "Acc")
    "Acc"
  }
  Tt <- c("L","M","H"); ub <- list(Poor=numeric(0), Acc=numeric(0), Good=numeric(0))
  for (g1t in Tt) for (g2t in Tt) for (g5t in Tt) for (ct in Tt) {
    dom <- scoring_dominant_block(sf, rep$g1[[g1t]], rep$g2[[g2t]],
                                  rep$g3[[ct]], rep$g4[[ct]], rep$g5[[g5t]])
    exp <- expected(g1t, g2t, g5t, ct)
    check(dom == exp, paste0("(", g1t, ",", g2t, ",", g5t, ",C~", ct, ") ", dom, " exp ", exp))
    u <- scoring_utility(sf, rep$g1[[g1t]], rep$g2[[g2t]], rep$g3[[ct]], rep$g4[[ct]], rep$g5[[g5t]])
    ub[[exp]] <- c(ub[[exp]], u)
  }
  cat("\n  centroid U ranges per intended block:\n")
  for (b in c("Poor","Acc","Good"))
    cat(sprintf("    %-5s: [%.3f, %.3f]  (n=%d)\n", b, min(ub[[b]]), max(ub[[b]]), length(ub[[b]])))
  cat(sprintf("\nRESULT: %d passed, %d failed  (pure-corner Poor~%.3f, Good~%.3f)\n",
              env$ok, env$fail,
              scoring_utility(sf, rep$g1[["H"]], rep$g2[["L"]], rep$g3[["L"]], rep$g4[["L"]], rep$g5[["H"]]),
              uopt))
  if (env$fail > 0) stop("validation failed")
  invisible(env$ok)
}

if (sys.nframe() == 0L) run_tests()

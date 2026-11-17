# =============================================================================
# figA_calibration_v3v6.R   (Annex A - calibration validation figure)
# -----------------------------------------------------------------------------
# Compares the evolution of the measured indoor temperature with the two
# candidate models (calibrated single-zone v3, two-zone v6) to support the
# validation of the retained model. Two panels:
#   (A) a readable ~10-day representative window (weekdays + a weekend) with the
#       three 5-min series, so v3 tracking the measurement and v6 overshooting
#       (+2.1 C) are directly visible; weekend nights shaded.
#   (B) the per-night mean bias of each model across ALL 50 nights, with the
#       pooled-mean reference lines, proving the pattern is systematic over the
#       whole campaign rather than an artefact of the shown window.
#
# Per-night mean bias:  bias_D^(m) = (1/n_D) * sum_{t in D} ( T^(m)_t - T_meas_t )
# Pooled mean bias (reference lines): v3 -> NMBE -5.41% (~ -1.0 C),
#                                     v6 -> NMBE +6.05% (~ +1.1 C).
#
# Window override:  ANNEX_WIN_START <- "YYYY-MM-DD"; ANNEX_WIN_END <- "YYYY-MM-DD"
# Uses the CALIBRATED series via load_calibrated_series.R (anchored).
# Output: PNG + PDF in OUT_DIR/deviation_figs.
# =============================================================================

suppressWarnings({library(lubridate)})
if(!exists(".CALIB_SERIES_LOADED") || !isTRUE(.CALIB_SERIES_LOADED))
  source("load_calibrated_series.R")
stopifnot(exists("D"))

# locale-independent English date labels (system locale may be Spanish)
.en_mon <- c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
fmt_en  <- function(d){ d <- as.Date(d); sprintf("%d %s", mday(d), .en_mon[month(d)]) }

# --- per-night mean bias over the whole campaign (Panel B) -------------------
hh   <- hour(D$k) + minute(D$k)/60
noct <- hh >= 22 | hh < 7
Dn   <- D[noct, ]
hhn  <- hour(Dn$k) + minute(Dn$k)/60
ngt  <- as.Date(Dn$k, tz = "UTC"); ngt[hhn < 7] <- ngt[hhn < 7] - 1
Dn$night <- ngt
nights <- sort(unique(Dn$night))
pn <- data.frame(
  night = nights,
  b3 = tapply(Dn$v3 - Dn$T_meas, Dn$night, mean),
  b6 = tapply(Dn$v6 - Dn$T_meas, Dn$night, mean))
pn$wend   <- wday(pn$night) %in% c(6, 7, 1)
pooled3 <- mean(Dn$v3 - Dn$T_meas); pooled6 <- mean(Dn$v6 - Dn$T_meas)

# --- representative window for Panel A ---------------------------------------
# Explicit override if provided; otherwise auto-select the 10-day window whose
# v3 NOCTURNAL NMBE is closest to the campaign headline (-5.41 %), so Panel (a)
# illustrates typical performance and cannot be (mis)read as cherry-picked.
if(exists("ANNEX_WIN_START") && exists("ANNEX_WIN_END")){
  ws <- as.POSIXct(ANNEX_WIN_START, tz="UTC"); we <- as.POSIXct(ANNEX_WIN_END, tz="UTC")
} else {
  WLEN   <- 10  # days
  starts <- seq(as.POSIXct("2021-11-01", tz="UTC"),
                as.POSIXct("2021-12-20", tz="UTC") - WLEN*86400, by = "1 day")
  score <- vapply(seq_along(starts), function(i){
    s  <- starts[i]; w <- D[D$k >= s & D$k < s + WLEN*86400, ]
    hw <- hour(w$k) + minute(w$k)/60; wn <- w[hw >= 22 | hw < 7, ]
    if(nrow(wn) < 800) return(NA_real_)
    100*mean(wn$v3 - wn$T_meas)/mean(wn$T_meas)
  }, numeric(1))
  best <- which.min(abs(score - (-5.41)))
  ws <- starts[best]; we <- ws + WLEN*86400
  cat(sprintf("Auto-selected representative window: v3 nocturnal NMBE = %+.2f%% (target -5.41%%)\n",
              score[best]))
}
W <- D[D$k >= ws & D$k < we, ]; W <- W[order(W$k), ]
# nocturnal mask within the window (calibration regime); daytime -> NA (gaps)
hhW   <- hour(W$k) + minute(W$k)/60
noctW <- hhW >= 22 | hhW < 7
Wn    <- W[noctW, ]
mbn   <- mean(Wn$T_meas)
nmbe_wn <- function(col) 100*mean(Wn[[col]] - Wn$T_meas)/mbn
cat(sprintf("Annex window %s..%s (NOCTURNAL): v3 NMBE %+.2f%%, v6 NMBE %+.2f%%\n",
            fmt_en(ws), fmt_en(we), nmbe_wn("v3"), nmbe_wn("v6")))
cat(sprintf("Whole-campaign pooled nocturnal bias: v3 %+.2f C, v6 %+.2f C\n", pooled3, pooled6))

# --- weekend shading spans within the window --------------------------------
wd <- wday(W$k); hhw <- hour(W$k) + minute(W$k)/60
inwk <- (wd==6 & hhw>=22) | (wd==7) | (wd==1) | (wd==2 & hhw<7)
rr <- rle(inwk); ce <- cumsum(rr$lengths); cs <- c(1, head(ce,-1)+1)
wk_spans <- lapply(which(rr$values), function(i) c(W$k[cs[i]], W$k[ce[i]]))

col_meas<-"black"; col_v3<-"#1f6feb"; col_v6<-"#d1242f"
faint <- function(c,a=0.14) adjustcolor(c, alpha.f=a)
lab_T <- expression("Indoor T ("*degree*"C)")
lab_b <- expression("Per-night mean bias  (sim "-" meas)  ("*degree*"C)")
OUTFIG <- if(exists("OUT_DIR")) file.path(OUT_DIR,"deviation_figs") else file.path(getwd(),"deviation_figs")
dir.create(OUTFIG, recursive=TRUE, showWarnings=FALSE)

draw <- function(){
  layout(matrix(1:2, ncol=1), heights=c(1.55, 1.25))
  # (A) representative window, NOCTURNAL segments only (daytime masked -> gaps)
  par(mar=c(3.0, 4.8, 2.6, 1.0))
  mask <- function(x){ x[!noctW] <- NA; x }
  yl <- range(c(Wn$T_meas, Wn$v3, Wn$v6), na.rm=TRUE)
  plot(W$k, mask(W$T_meas), type="n", ylim=yl, xaxt="n", xlab="", ylab=lab_T,
       main=sprintf("(a)  Nocturnal temperature: measured vs models  (%s - %s 2021)",
                    fmt_en(ws), fmt_en(we-86400)))
  for(s in wk_spans) rect(s[1], yl[1], s[2], yl[2], col=faint("grey40"), border=NA)
  lines(W$k, mask(W$T_meas), col=col_meas, lwd=1.7)
  lines(W$k, mask(W$v3),     col=col_v3,   lwd=1.5)
  lines(W$k, mask(W$v6),     col=col_v6,   lwd=1.5)
  axis.POSIXct(1, at=seq(ws, we, by="1 day"), format="%d %b")
  legend("topright", bty="n", horiz=TRUE, lwd=2, col=c(col_meas,col_v3,col_v6),
         legend=c("Measured","v3 (1-zone)","v6 (2-zone)"), cex=0.9)
  legend("bottomleft", bty="n", cex=0.85,
         legend=c(sprintf("window NMBE:  v3 %+.1f%%,  v6 %+.1f%%", nmbe_wn("v3"), nmbe_wn("v6")),
                  "shaded = weekend"))
  # (B) per-night mean bias across the whole campaign
  par(mar=c(3.4, 4.8, 2.6, 1.0))
  yb <- range(c(pn$b3, pn$b6), na.rm=TRUE)
  plot(pn$night, pn$b3, type="n", ylim=yb, xlab="", ylab=lab_b, xaxt="n",
       main="(b)  Per-night mean bias across the campaign (50 nights)")
  .bt <- as.Date(c("2021-11-01","2021-11-15","2021-12-01","2021-12-15"))
  axis(1, at=as.numeric(.bt), labels=fmt_en(.bt))
  abline(h=0, col="grey20", lty=3)
  abline(h=pooled3, col=col_v3, lty=2); abline(h=pooled6, col=col_v6, lty=2)
  lines(pn$night, pn$b3, col=col_v3, lwd=1.4, type="b", pch=16, cex=0.6)
  lines(pn$night, pn$b6, col=col_v6, lwd=1.4, type="b", pch=17, cex=0.6)
  mtext("Night (evening date, 2021)", side=1, line=2.3, cex=0.95)
  legend("topleft", bty="n", cex=0.85, lwd=c(1.4,1.4,1,1), pch=c(16,17,NA,NA),
         lty=c(1,1,2,2), col=c(col_v3,col_v6,col_v3,col_v6),
         legend=c("v3 per-night bias","v6 per-night bias",
                  sprintf("v3 pooled = %+.2f C (NMBE -5.41%%)", pooled3),
                  sprintf("v6 pooled = %+.2f C (NMBE +6.05%%)", pooled6)))
}

png(file.path(OUTFIG,"FigA_calibration_v3v6.png"), width=11, height=8.2, units="in", res=200); draw(); dev.off()
pdf(file.path(OUTFIG,"FigA_calibration_v3v6.pdf"), width=11, height=8.2, encoding="ISOLatin1", useDingbats=FALSE); draw(); dev.off()
cat("Saved FigA_calibration_v3v6.{png,pdf} to:\n  ", OUTFIG, "\n")

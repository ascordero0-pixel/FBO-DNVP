# =============================================================================
# deviation_nocturnal.R
# -----------------------------------------------------------------------------
# Deviation (residual) of the v3 (single-zone) and v6 (two-zone) models with
# respect to the MEASURED indoor temperature, restricted to the NOCTURNAL
# regime 22:00 -> 07:00, for every 5-min step of the whole RunPeriod
# (1 Nov -> 20 Dec 2021). This is the regime Track A calibrates, so these
# metrics are directly comparable with Track A (NMBE ~ -5.4 %).
#
# Residuals (sign: negative = model colder than measured):
#     r_v3(t) = T_v3(t) - T_meas(t) ;  r_v6(t) = T_v6(t) - T_meas(t)
#
# A point is nocturnal iff  hour>=22  OR  hour<7 . Each night (22:00 day D ->
# 07:00 day D+1) is labelled by its evening date D: early-morning points
# (hour<7) are reassigned to the previous day. The bold overlay is the
# PER-NIGHT MEAN residual (not a moving average, which would wrongly bridge the
# ~15 h daytime gap between nights). Weekend nights (Fri/Sat/Sun evenings, when
# the staircase heat input appears) are shaded; global metrics are split into
# weekday vs weekend nights.
#
# REQUIRES in memory: df (k, T_out, v3, v6) and meas (k, T_meas, ...).
# Output: PNG + PDF in OUT_DIR/deviation_figs.
# =============================================================================

suppressWarnings({library(lubridate)})
# Ensure df/meas come from the CALIBRATED Track A series (not the Option A run).
if(!exists(".CALIB_SERIES_LOADED") || !isTRUE(.CALIB_SERIES_LOADED))
  source("load_calibrated_series.R")
stopifnot(exists("df"), exists("meas"))

# --- full merged 5-min series over the RunPeriod ----------------------------
DF <- merge(df, meas[, c("k","T_meas")], by = "k")
DF <- DF[order(DF$k), ]
t0 <- as.POSIXct("2021-11-01 00:00:00", tz = "UTC")
t1 <- as.POSIXct("2021-12-21 00:00:00", tz = "UTC")
DF <- DF[DF$k >= t0 & DF$k < t1, ]
DF$r_v3 <- DF$v3 - DF$T_meas
DF$r_v6 <- DF$v6 - DF$T_meas

# --- nocturnal mask (22:00 -> 07:00) ----------------------------------------
hh   <- hour(DF$k) + minute(DF$k)/60
noct <- (hh >= 22) | (hh < 7)

# --- night label: evening date D (early-morning -> previous day) ------------
night <- as.Date(DF$k, tz = "UTC")
night[hh < 7] <- night[hh < 7] - 1
DF$night <- night

# --- nocturnal subset and per-night means -----------------------------------
DFn <- DF[noct, ]
mbar_n <- mean(DFn$T_meas)
nmbe   <- function(r) 100 * mean(r) / mbar_n
cvrmse <- function(r) 100 * sqrt(mean(r^2)) / mbar_n

nights <- sort(unique(DFn$night))
pnight <- data.frame(
  night = nights,
  tc = as.POSIXct(tapply(as.numeric(DFn$k), DFn$night, mean),
                  origin = "1970-01-01", tz = "UTC"),
  m3 = tapply(DFn$r_v3, DFn$night, mean),
  m6 = tapply(DFn$r_v6, DFn$night, mean))
# CANONICAL weekday/weekend partition: a night (22:00->07:00) is an indivisible
# physical unit, labelled by its EVENING date D (early-morning points, hour<7,
# were reassigned to D-1 above). A night is "weekend" iff D is Fri(6)/Sat(7)/
# Sun(1) (lubridate wday: 1=Sun). This is the reference criterion for the
# manuscript. Do NOT classify by wday(timestamp), which would split a single
# night at midnight and can even invert which regime a model fits better.
pnight$wend <- wday(pnight$night) %in% c(6, 7, 1)
DFn$wend    <- wday(DFn$night)   %in% c(6, 7, 1)

# --- console metrics: overall + weekday/weekend split -----------------------
cat(sprintf("NOCTURNAL 22:00-07:00  (1 Nov - 20 Dec)  N=%d  nights=%d  mean meas=%.2f C\n",
            nrow(DFn), length(nights), mbar_n))
cat(sprintf("  ALL nights      v3: NMBE %+.2f%%  CV-RMSE %.2f%%   |  v6: NMBE %+.2f%%  CV-RMSE %.2f%%\n",
            nmbe(DFn$r_v3), cvrmse(DFn$r_v3), nmbe(DFn$r_v6), cvrmse(DFn$r_v6)))
for(lbl in c("weekday","weekend")){
  s <- if(lbl=="weekday") DFn[!DFn$wend, ] else DFn[DFn$wend, ]
  mb <- mean(s$T_meas)
  f <- function(r) c(100*mean(r)/mb, 100*sqrt(mean(r^2))/mb)
  a3 <- f(s$v3 - s$T_meas); a6 <- f(s$v6 - s$T_meas)
  cat(sprintf("  %-7s nights  v3: NMBE %+.2f%%  CV-RMSE %.2f%%   |  v6: NMBE %+.2f%%  CV-RMSE %.2f%%  (n=%d)\n",
              lbl, a3[1], a3[2], a6[1], a6[2], nrow(s)))
}

# --- masked series (daytime -> NA so lines break at the gaps) ---------------
mask <- function(x){ x[!noct] <- NA; x }
r3m <- mask(DF$r_v3); r6m <- mask(DF$r_v6)
Tm  <- mask(DF$T_meas); v3m <- mask(DF$v3); v6m <- mask(DF$v6)

# --- weekend spans for shading ----------------------------------------------
inwk <- (wday(DF$k)==6 & hh>=22) | (wday(DF$k)==7) | (wday(DF$k)==1) | (wday(DF$k)==2 & hh<7)
rr <- rle(inwk); ce <- cumsum(rr$lengths); cs <- c(1, head(ce,-1)+1)
wk_spans <- lapply(which(rr$values), function(i) c(DF$k[cs[i]], DF$k[ce[i]]))

col_meas<-"black"; col_v3<-"#1f6feb"; col_v6<-"#d1242f"
faint <- function(c,a=0.20) adjustcolor(c, alpha.f=a)
lab_Tin <- expression("Indoor T ("*degree*"C)")
lab_res <- expression("Deviation  (model "-" measured)  ("*degree*"C)")
OUTFIG <- if(exists("OUT_DIR")) file.path(OUT_DIR,"deviation_figs") else file.path(getwd(),"deviation_figs")
dir.create(OUTFIG, recursive=TRUE, showWarnings=FALSE)
tick_axis <- function() axis.POSIXct(1, at=seq(t0,t1,by="5 days"), format="%d-%b")
shade <- function(yb,yt) for(s in wk_spans) rect(s[1],yb,s[2],yt,col=faint("grey40",0.12),border=NA)

draw <- function(){
  layout(matrix(1:2, ncol=1), heights=c(1.4, 2.1))
  # (A) nocturnal absolute temperatures
  par(mar=c(2.4,4.8,2.6,1.2)); yl<-range(c(DFn$T_meas,DFn$v3,DFn$v6),na.rm=TRUE)
  plot(DF$k, Tm, type="n", ylim=yl, xaxt="n", xlab="", ylab=lab_Tin,
       main="Nocturnal indoor T (22:00-07:00): measured vs models - 1 Nov to 20 Dec 2021")
  shade(yl[1],yl[2])
  lines(DF$k, Tm, col=col_meas, lwd=0.7)
  lines(DF$k, v3m, col=faint(col_v3,0.8), lwd=0.6)
  lines(DF$k, v6m, col=faint(col_v6,0.8), lwd=0.6)
  tick_axis()
  legend("topright", bty="n", horiz=TRUE, lwd=2, col=c(col_meas,col_v3,col_v6),
         legend=c("Measured","v3","v6"), cex=0.9)
  # (B) nocturnal residuals + per-night means
  par(mar=c(3.6,4.8,1.2,1.2)); yb<-range(c(DFn$r_v3,DFn$r_v6),na.rm=TRUE)
  plot(DF$k, r3m, type="n", ylim=yb, xaxt="n", xlab="", ylab=lab_res)
  shade(yb[1],yb[2]); abline(h=0, col="grey20", lty=3)
  lines(DF$k, r3m, col=faint(col_v3,0.22), lwd=0.5)     # raw nocturnal 5-min
  lines(DF$k, r6m, col=faint(col_v6,0.22), lwd=0.5)
  lines(pnight$tc, pnight$m3, col=col_v3, lwd=2.2, type="b", pch=16, cex=0.6)  # per-night mean
  lines(pnight$tc, pnight$m6, col=col_v6, lwd=2.2, type="b", pch=16, cex=0.6)
  tick_axis(); mtext("Date (2021)", side=1, line=2.3, cex=0.95)
  legend("bottomleft", bty="n", ncol=2, lwd=c(0.8,2.2,0.8,2.2),
         col=c(faint(col_v3,0.5),col_v3,faint(col_v6,0.5),col_v6),
         legend=c("v3 raw 5-min","v3 per-night mean","v6 raw 5-min","v6 per-night mean"))
  legend("topright", bty="n", cex=0.9,
         legend=c(sprintf("v3: NMBE %+.1f%%  CV-RMSE %.1f%%", nmbe(DFn$r_v3), cvrmse(DFn$r_v3)),
                  sprintf("v6: NMBE %+.1f%%  CV-RMSE %.1f%%", nmbe(DFn$r_v6), cvrmse(DFn$r_v6)),
                  "shaded = weekend nights"))
}
png(file.path(OUTFIG,"deviation_nocturnal.png"), width=15, height=8.5, units="in", res=150); draw(); dev.off()
pdf(file.path(OUTFIG,"deviation_nocturnal.pdf"), width=15, height=8.5, encoding="ISOLatin1", useDingbats=FALSE); draw(); dev.off()
cat("Saved deviation_nocturnal.{png,pdf} to:\n  ", OUTFIG, "\n")

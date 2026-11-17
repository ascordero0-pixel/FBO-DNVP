# tau by per-night exponential fit: T(t) = Tinf + (T0 - Tinf) * exp(-t/tau)
# Robust to distributed mass: fits each night's cooling separately, averages.
# Uses n5 (nocturnal 5-min sim+meas) already in memory + exterior column from CSV.

csv <- read.csv(file.path(PATHS$runs_dir,"v3_5min","eplusout.csv"), check.names=FALSE)
EXT <- grep("Outdoor Air Drybulb.*TimeStep", names(csv), value=TRUE)[1]
parse_dt<-function(x){x<-trimws(x);mo<-as.integer(substr(x,1,2));da<-as.integer(substr(x,4,5))
  tt<-sub("^\\d{2}/\\d{2}\\s+","",x);hh<-as.integer(substr(tt,1,2));mi<-as.integer(substr(tt,4,5))
  nd<-which(hh==24);hh[nd]<-0;ts<-as.POSIXct(sprintf("2021-%02d-%02d %02d:%02d:00",mo,da,hh,mi),tz="UTC");ts[nd]<-ts[nd]+24*3600;ts}
ext <- data.frame(k=lubridate::floor_date(parse_dt(csv[[1]]),"5 mins"), T_out=as.numeric(csv[[EXT]]))

dat <- merge(n5[,c("k","T_sim")], ext, by="k"); dat <- dat[order(dat$k),]
# assign a "night id": a night is 22:00 -> 07:00; label by the date at 22:00
hh <- as.integer(format(dat$k,"%H"))
# night index increments at each >2h gap
gap <- c(TRUE, as.numeric(difftime(dat$k[-1],dat$k[-nrow(dat)],units="hours")) > 2)
dat$night <- cumsum(gap)

fit_night <- function(g){
  g <- g[order(g$k),]
  # cooling phase only: from the start until temperature stops falling materially
  if(nrow(g) < 24) return(NA)                       # need >=2 h of 5-min data
  t_h <- as.numeric(difftime(g$k, g$k[1], units="hours"))
  Tout <- mean(g$T_out, na.rm=TRUE)                 # use night-mean outdoor as asymptote ref
  y <- g$T_sim - Tout                               # excess temperature above outdoor
  y <- y[y > 0.3]                                   # while clearly above outdoor
  if(length(y) < 24) return(NA)
  tt <- t_h[seq_along(y)]
  # linearize: ln(y) = ln(y0) - t/tau  -> slope = -1/tau
  ok <- is.finite(log(y))
  if(sum(ok) < 24) return(NA)
  f <- lm(log(y[ok]) ~ tt[ok])
  sl <- unname(coef(f)[2]); r2 <- summary(f)$r.squared
  if(sl >= 0) return(NA)                            # must decay
  c(tau = -1/sl, r2 = r2, npts = sum(ok))
}

res <- lapply(split(dat, dat$night), fit_night)
res <- do.call(rbind, res[!vapply(res, function(z) all(is.na(z)), logical(1))])
res <- as.data.frame(res)
res <- res[is.finite(res$tau) & res$tau > 0 & res$tau < 48 & res$r2 > 0.7, ]  # keep good fits

cat(sprintf("\nPer-night exponential fits kept: %d nights (R2 > 0.7)\n", nrow(res)))
cat(sprintf("  tau  : median = %.2f h | mean = %.2f h | SD = %.2f h | IQR = [%.2f, %.2f]\n",
    median(res$tau), mean(res$tau), sd(res$tau),
    quantile(res$tau,0.25), quantile(res$tau,0.75)))
cat(sprintf("  fit  : median R2 = %.3f | typical npts = %d\n", median(res$r2), round(median(res$npts))))
cat("\n  (use the MEDIAN tau as the identified value; the spread reflects\n",
    "  night-to-night outdoor variability, not model error)\n")

write.csv(res, file.path(PATHS$runs_dir,"tau_per_night.csv"), row.names=FALSE)

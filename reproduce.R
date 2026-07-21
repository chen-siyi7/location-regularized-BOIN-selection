#!/usr/bin/env Rscript
# ============================================================================
# Location-regularized BOIN selection -- full reproducible analysis (base R).
# Run from a terminal (e.g. macOS):
#     Rscript reproduce.R          # fast defaults (a few minutes)
#     Rscript reproduce.R full     # manuscript-scale replicates (slow)
# Prints Tables 2-7 to the console and writes Figures 1-5 as PDF files.
# No external packages required (base R only).
# ============================================================================

source("boin_core.R")     # BOIN conduct + location-regularized selection
source("designs.R")       # 3+3, CRM, biased-coin, k-in-a-row
source("comparators.R")   # logistic-regression / Bayesian model-based selection

## ---- settings -------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
TEST  <- "test" %in% args           # tiny smoke-test scale
QUICK <- !("full" %in% args)
PHI <- 0.25; N <- 30; J <- 6
n_fixed <- if (TEST) 800 else if (QUICK) 4000 else 10000
n_rand  <- if (TEST) 400 else if (QUICK) 1500 else 10000
n_train <- if (TEST)  12 else if (QUICK)   50 else 150
n_test  <- if (TEST)  15 else if (QUICK)   80 else 200
k_pos   <- if (TEST)   8 else if (QUICK)   30 else 55
k_sweep <- if (TEST)  10 else if (QUICK)   40 else 70
n_cross <- if (TEST)  10 else if (QUICK)   40 else 110
n_comp  <- if (TEST) 300 else if (QUICK)  800 else 1500
tab  <- elim_table(PHI, N)
SKEL <- c(0.05, 0.12, 0.25, 0.40, 0.55, 0.70)
cat(sprintf("Mode: %s\n\n", if (TEST) "TEST (tiny smoke test)" else if (QUICK)
            "QUICK (fast; pass 'full' for manuscript-scale)" else "FULL (manuscript-scale)"))
hdr <- function(s) cat("\n", strrep("=", 70), "\n", s, "\n", strrep("=", 70), "\n", sep = "")

## ---- cache helpers --------------------------------------------------------
cache_design <- function(p, design, n_sim, seed) {
  set.seed(seed); bnd <- boin_boundaries(PHI); Jl <- length(p)
  nm <- matrix(0L, n_sim, Jl); xm <- matrix(0L, n_sim, Jl)
  crv <- if (design == "crm") matrix(0, n_sim, Jl) else NULL; own <- integer(n_sim)
  for (t in seq_len(n_sim)) {
    tr <- switch(design,
      boin = simulate_one_trial(p, PHI, N, bnd = bnd, tab = tab),
      "3p3" = sim_3p3(p), bcd = sim_bcd(p, PHI, N),
      krow = sim_krow(p, PHI, N), crm = sim_crm(p, PHI, N, SKEL))
    nm[t, ] <- tr$n; xm[t, ] <- tr$x
    if (design == "crm") crv[t, ] <- tr$curve
    if (design == "3p3") own[t] <- tr$mtd
  }
  list(n = nm, x = xm, curve = crv, own = own, mtd = which.min(abs(p - PHI)))
}
pc_hit <- function(s, m) mean(s == m) * 100
pc_std <- function(c) pc_hit(if (is.null(c$curve))
  vapply(seq_len(nrow(c$n)), function(t) select_boin(c$n[t,], c$x[t,], PHI, tab = tab), 0L)
  else vapply(seq_len(nrow(c$n)), function(t) select_curve(c$curve[t,], c$n[t,], PHI), 0L), c$mtd)
pc_prior <- function(c, rho, alpha) pc_hit(if (is.null(c$curve))
  vapply(seq_len(nrow(c$n)), function(t) select_lr(c$n[t,], c$x[t,], PHI, rho, alpha, tab = tab), 0L)
  else vapply(seq_len(nrow(c$n)), function(t) select_curve(c$curve[t,], c$n[t,], PHI, rho, alpha), 0L), c$mtd)

## ---- validation -----------------------------------------------------------
hdr("Validation")
tab60 <- elim_table(0.30, 60)
cn <- run_paired(c(0.05,0.15,0.30,0.45,0.60), 0.30, 60, 4000, 0, 2.5, seed = 6, tab = tab60)
cat("rho = 0 recovers BOIN exactly:",
    isTRUE(all.equal(cn["boin"], cn["lr"], check.names = FALSE)), "\n")

## ---- Table 2: fixed case-study scenarios ---------------------------------
hdr("Table 2 -- fixed case-study scenarios, (rho,alpha)=(3.0,2.5)")
FIXED <- list(
 "1 Moderate ramp"=c(0.05,0.10,0.25,0.40,0.55,0.70), "2 Steep cliff"=c(0.02,0.08,0.25,0.55,0.80,0.95),
 "3 Gradual rise"=c(0.08,0.14,0.20,0.26,0.34,0.45),   "5 Sparse DLT"=c(0.04,0.06,0.12,0.25,0.40,0.58),
 "6 Flat near tgt"=c(0.15,0.20,0.24,0.28,0.33,0.40),  "7 Very steep"=c(0.01,0.03,0.25,0.60,0.85,0.95),
 "4 Low-mid"=c(0.12,0.25,0.40,0.55,0.70,0.85), "S1 High MTD"=c(0.02,0.05,0.08,0.15,0.25,0.45),
 "S2 All overdose"=c(0.30,0.45,0.55,0.70,0.80,0.90), "S3 Non-monotone"=c(0.10,0.25,0.15,0.30,0.45,0.60))
t2 <- data.frame()
for (nm in names(FIXED)) {
  r <- run_paired(FIXED[[nm]], PHI, N, n_fixed, 3.0, 2.5, seed = 42, tab = tab)
  t2 <- rbind(t2, data.frame(Scenario=nm, BOIN=round(r["boin"],1),
                             prior=round(r["lr"],1), dPC=round(r["lr"]-r["boin"],1)))
}
rownames(t2) <- NULL; print(t2, row.names = FALSE)
cat("Mean recommended-domain (rows 1-6):", round(mean(t2$dPC[1:6]),1), "pp\n")

## ---- random scenarios (shared) -------------------------------------------
pool <- gen_random_scenarios(2500, J, PHI, seed = 2024)
interior <- Filter(function(s) s$mtd %in% c(3,4), pool)
set.seed(7); interior <- interior[sample(length(interior))]
train <- interior[seq_len(n_train)]; test <- interior[n_train + seq_len(n_test)]

## ---- Table 3: out-of-sample tuning & validation --------------------------
hdr("Table 3 -- out-of-sample tuning and validation")
trc <- lapply(seq_along(train), function(i) cache_design(train[[i]]$p, "boin", n_rand, 5000+i))
grid <- expand.grid(rho=c(2,3,4), alpha=c(2.0,2.5,3.0))
grid$dPC <- apply(grid, 1, function(g) mean(sapply(trc, function(c) pc_prior(c,g[1],g[2]) - pc_std(c))))
best <- grid[which.max(grid$dPC), ]
tec <- lapply(seq_along(test), function(i) cache_design(test[[i]]$p, "boin", n_rand, 8000+i))
oos <- function(rho, alpha) { d <- sapply(tec, function(c) pc_prior(c,rho,alpha) - pc_std(c))
  c(dPC = round(mean(d),2), SE = round(sd(d)/sqrt(length(d)),2)) }
cat(sprintf("tuned (rho,alpha) = (%.1f, %.1f), in-sample dPC = %.2f\n", best$rho, best$alpha, best$dPC))
t3 <- rbind(`tuned` = oos(best$rho,best$alpha),
            `recommended (3.0,2.5)` = oos(3.0,2.5),
            `conservative (2.0,2.0)` = oos(2.0,2.0))
print(t3)

## ---- Table 4: parameter sensitivity by MTD position ----------------------
hdr("Table 4 -- parameter sensitivity by true-MTD position")
posp <- gen_random_scenarios(5000, J, PHI, seed = 555)
bypos <- lapply(1:J, function(j) Filter(function(s) s$mtd==j, posp))
cache_pos <- lapply(1:J, function(j) {
  sc <- bypos[[j]]; if (!length(sc)) return(list()); sc <- sc[seq_len(min(k_pos, length(sc)))]
  lapply(seq_along(sc), function(i) cache_design(sc[[i]]$p, "boin", n_rand, 20000+100*j+i)) })
grp_dpc <- function(js, rho, alpha) { cs <- do.call(c, cache_pos[js])
  mean(sapply(cs, function(c) pc_prior(c,rho,alpha) - pc_std(c))) }
settings <- list(c(1,2.0), c(2,2.0), c(3,2.5), c(3,3.0), c(4,3.0))
t4 <- do.call(rbind, lapply(settings, function(s) data.frame(rho=s[1], alpha=s[2],
  interior=round(grp_dpc(c(3,4),s[1],s[2]),1), boundary_adj=round(grp_dpc(c(2,5),s[1],s[2]),1),
  boundary=round(grp_dpc(c(1,6),s[1],s[2]),1))))
print(t4, row.names = FALSE)

## ---- Table 5: sweeps ------------------------------------------------------
hdr("Table 5 -- dose-count, target-rate, sample-size sweeps")
interior_scen <- function(Jl, phi, seed) { lo <- 2; hi <- Jl-3
  pl <- gen_random_scenarios(3000, Jl, phi, seed)
  sc <- Filter(function(s) s$mtd >= (lo+1) & s$mtd <= (hi+1), pl); sc[seq_len(min(k_sweep, length(sc)))] }
sweep_dpc <- function(scen, phi, Nn, Jl, seed0) { tb <- elim_table(phi, Nn)
  d <- sapply(seq_along(scen), function(i) { set.seed(seed0+i); bnd <- boin_boundaries(phi)
    mtd <- scen[[i]]$mtd; p <- scen[[i]]$p; hs <- 0; hl <- 0
    for (t in 1:n_rand) { tr <- simulate_one_trial(p, phi, Nn, bnd = bnd, tab = tb)
      iso <- iso_estimate(tr$n, tr$x, phi)
      hs <- hs + (select_boin(tr$n,tr$x,phi,iso,tab=tb)==mtd)
      hl <- hl + (select_lr(tr$n,tr$x,phi,3,2.5,iso,tab=tb)==mtd) }
    100*(hl-hs)/n_rand })
  c(mean = mean(d), se = sd(d)/sqrt(length(d))) }
Js <- c(5,6,7,8); dJ <- sapply(Js, function(Jl) sweep_dpc(interior_scen(Jl,0.25,1000+Jl),0.25,30,Jl,1e6+Jl))
phis <- c(0.10,0.15,0.20,0.25,0.30,0.35,0.40)
dphi <- sapply(phis, function(ph) sweep_dpc(interior_scen(6,ph,round(2000+ph*100)),ph,30,6,2e6+round(ph*1000)))
Ns <- c(18,24,30,36,42,48); scN <- interior_scen(6,0.25,3030)
dN <- sapply(Ns, function(nn) sweep_dpc(scN,0.25,nn,6,3e6+nn))
cat("J-sweep   dPC:", paste0(Js,":",round(dJ["mean",],1), collapse="  "), "\n")
cat("phi-sweep dPC:", paste0(phis,":",round(dphi["mean",],1), collapse="  "), "\n")
cat("N-sweep   dPC:", paste0(Ns,":",round(dN["mean",],1), collapse="  "), "\n")

## ---- Table 6: cross-design transfer --------------------------------------
hdr("Table 6 -- transfer across allocation designs")
cross <- interior[seq_len(n_cross)]; designs <- c("3p3","bcd","krow","boin","crm")
t6 <- data.frame(); dpc_cross <- numeric(); se_cross <- numeric()
for (dz in designs) {
  cs <- lapply(seq_along(cross), function(i) cache_design(cross[[i]]$p, dz, if (QUICK) 1000 else 3000, 9000+i))
  sstd <- sapply(cs, pc_std); spri <- sapply(cs, function(c) pc_prior(c,3.0,2.5)); dvec <- spri - sstd
  own <- if (dz=="3p3") round(mean(sapply(cs, function(c) pc_hit(c$own, c$mtd))),1) else NA
  t6 <- rbind(t6, data.frame(design=dz, own_rule=own, std=round(mean(sstd),1),
                             prior=round(mean(spri),1), dPC=round(mean(dvec),1)))
  dpc_cross[dz] <- mean(dvec); se_cross[dz] <- sd(dvec)/sqrt(length(dvec))
}
print(t6, row.names = FALSE)

## ---- Table 7: model-based comparators ------------------------------------
hdr("Table 7 -- model-based selection comparators")
methods <- c("BOIN","logreg (Liao)","bayes-logit (Wakayama)","bayes-cloglog (Wakayama)","BOIN+prior","bayes-logit+prior")
acc <- matrix(0, length(cross), length(methods))
for (i in seq_along(cross)) {
  p <- cross[[i]]$p; mtd <- cross[[i]]$mtd; set.seed(12000+i); h <- numeric(length(methods))
  for (t in 1:n_comp) { tr <- simulate_one_trial(p, PHI, N, tab = tab)
    h[1] <- h[1] + (select_boin(tr$n,tr$x,PHI,tab=tab)==mtd)
    h[2] <- h[2] + (select_logreg(tr$n,tr$x,PHI)==mtd)
    h[3] <- h[3] + (select_bayes(tr$n,tr$x,PHI,"logit")==mtd)
    h[4] <- h[4] + (select_bayes(tr$n,tr$x,PHI,"cloglog")==mtd)
    h[5] <- h[5] + (select_lr(tr$n,tr$x,PHI,3,2.5,tab=tab)==mtd)
    h[6] <- h[6] + (select_bayes(tr$n,tr$x,PHI,"logit",3,2.5)==mtd) }
  acc[i, ] <- 100*h/n_comp
}
pcm <- colMeans(acc)
print(data.frame(method=methods, mean_PC=round(pcm,1), dPC_vs_BOIN=round(pcm-pcm[1],1)), row.names = FALSE)

## ---- Figures (all five, publication style) --------------------------------
hdr("Writing figures 1-5")
OI <- c("#0072B2","#E69F00","#009E73","#CC79A7","#56B4E9","#D55E00","#000000")  # colorblind-safe
ebar <- function(x,m,se) suppressWarnings(arrows(x, m-se, x, m+se, angle=90, code=3,
                                                 length=0.03, col="grey30"))

## Figure 1 -- location penalty shape
pdf("fig1_prior_shape.pdf", 6.5, 4.2)
z <- seq(0.001,0.999,length=400)
plot(z,(3/J)*1.5*(log(z)+log(1-z)), type="l", lwd=2, col=OI[6], ylim=c(-6,0.5),
     xlab="normalized dose position  z = (j-0.5)/J", ylab="location-prior contribution")
lines(z,(2/J)*(log(z)+log(1-z)), lty=2, lwd=1.5, col="grey40")
lines(z,(4/J)*2*(log(z)+log(1-z)), lty=3, lwd=1.5, col=OI[1])
zj <- (1:J-0.5)/J; pj <- (3/J)*1.5*(log(zj)+log(1-zj))
abline(v=zj, col="grey92"); points(zj, pj, pch=19, col=OI[6], cex=1.2)
text(zj, -6, paste0("d",1:J), cex=0.7, col="grey40", pos=3)
legend("bottom", c("conservative (2,2)","recommended (3,2.5)","aggressive (4,3)"),
       lty=c(2,1,3), lwd=c(1.5,2,1.5), col=c("grey40",OI[6],OI[1]), bty="n", cex=0.8)
dev.off()

## Figure 2 -- dose-toxicity scenario curves
pdf("fig2_scenarios.pdf", 9, 4); par(mfrow=c(1,2), mar=c(4,4,2.2,1))
doses <- 1:6
rec_nm <- c("1 Moderate ramp","2 Steep cliff","3 Gradual rise","5 Sparse DLT","6 Flat near tgt","7 Very steep")
plot(NA, xlim=c(1,6), ylim=c(0,1), xlab="dose level", ylab="true DLT probability",
     main="(a) Recommended-domain & Low-mid", cex.main=0.95)
for (k in seq_along(rec_nm)) { p <- FIXED[[rec_nm[k]]]; lines(doses,p,type="o",pch=19,cex=0.7,col=OI[k],lwd=1.3)
  m <- which.min(abs(p-PHI)); points(m,p[m],pch=1,cex=2,lwd=1.5,col=OI[k]) }
p <- FIXED[["4 Low-mid"]]; lines(doses,p,type="o",pch=19,cex=0.7,lty=2,col=OI[7],lwd=1.3)
m <- which.min(abs(p-PHI)); points(m,p[m],pch=1,cex=2,lwd=1.5,col=OI[7])
abline(h=PHI,lty=3,col="grey50")
legend("topleft", c(rec_nm,"4 Low-mid"), col=c(OI[1:6],OI[7]), lty=c(rep(1,6),2), pch=19, cex=0.6, bty="n")
str_nm <- c("S1 High MTD","S2 All overdose","S3 Non-monotone"); sc <- OI[c(6,2,4)]
plot(NA, xlim=c(1,6), ylim=c(0,1), xlab="dose level", ylab="true DLT probability",
     main="(b) External stress scenarios", cex.main=0.95)
for (k in seq_along(str_nm)) { p <- FIXED[[str_nm[k]]]; lines(doses,p,type="o",pch=19,cex=0.7,col=sc[k],lwd=1.3)
  m <- which.min(abs(p-PHI)); points(m,p[m],pch=1,cex=2,lwd=1.5,col=sc[k]) }
abline(h=PHI,lty=3,col="grey50")
legend("topleft", str_nm, col=sc, lty=1, pch=19, cex=0.7, bty="n")
dev.off()

## Figure 3 -- gain by true-MTD position (recommended vs conservative)
pdf("fig3_position.pdf", 6.5, 4.2)
dpos_rec <- sapply(1:J, function(j) grp_dpc(j,3.0,2.5))
dpos_con <- sapply(1:J, function(j) grp_dpc(j,2.0,2.0))
yl <- range(c(dpos_rec,dpos_con,0)); yl <- yl + c(-1,1)*diff(yl)*0.12
bp <- barplot(rbind(dpos_rec,dpos_con), beside=TRUE, names.arg=paste0("d",1:J),
              col=c(OI[6],OI[5]), ylab=expression(Delta*P[C]*" (pp)"), xlab="true MTD position", ylim=yl)
rect(bp[1,3]-0.5, yl[1], bp[2,4]+0.5, yl[2], col=rgb(0,0.62,0.45,0.08), border=NA)
barplot(rbind(dpos_rec,dpos_con), beside=TRUE, col=c(OI[6],OI[5]), add=TRUE, axes=FALSE, names.arg=rep("",J))
abline(h=0)
legend("bottomleft", c("recommended (3.0,2.5)","conservative (2.0,2.0)"),
       fill=c(OI[6],OI[5]), bty="n", cex=0.8)
dev.off()

## Figure 4 -- dose-count, target-rate, sample-size sweeps
pdf("fig4_sweeps.pdf", 9.5, 3.2); par(mfrow=c(1,3), mar=c(4,4,2,1))
plot(Js, dJ["mean",], type="b", pch=19, col=OI[1], ylim=c(0,max(dJ["mean",]+dJ["se",])*1.15),
     xlab="number of doses J", ylab=expression(Delta*P[C]*" (pp)"), main="(a)"); ebar(Js,dJ["mean",],dJ["se",]); abline(h=0,col="grey80")
plot(phis, dphi["mean",], type="b", pch=15, col=OI[3], ylim=c(0,max(dphi["mean",]+dphi["se",])*1.15),
     xlab=expression("target rate "*varphi), ylab="", main="(b)"); ebar(phis,dphi["mean",],dphi["se",]); abline(h=0,col="grey80")
plot(Ns, dN["mean",], type="b", pch=17, col=OI[6], ylim=c(0,max(dN["mean",]+dN["se",])*1.15),
     xlab="sample size N", ylab="", main="(c)"); ebar(Ns,dN["mean",],dN["se",]); abline(h=0,col="grey80")
dev.off()

## Figure 5 -- transfer across allocation designs
pdf("fig5_crossdesign.pdf", 6.5, 4.2)
ord <- c("3p3","bcd","krow","boin","crm"); mm <- dpc_cross[ord]; ss <- se_cross[ord]
bp <- barplot(mm, names.arg=c("3+3","BCD","k-in-a-row","BOIN","CRM"), col=OI[3],
              ylab=expression(Delta*P[C]*" from location prior (pp)"), ylim=c(0,max(mm+ss)*1.18))
suppressWarnings(arrows(bp, mm-ss, bp, mm+ss, angle=90, code=3, length=0.04, col="grey30")); abline(h=0)
text(bp, mm+ss, sprintf("+%.1f", mm), pos=3, cex=0.8)
dev.off()

cat("Wrote fig1_prior_shape.pdf, fig2_scenarios.pdf, fig3_position.pdf,",
    "fig4_sweeps.pdf, fig5_crossdesign.pdf\n")
cat("\nDone.\n")

# Model-based selection comparators (base R), applied to BOIN-conduct (n,x) tables:
#   select_logreg  - ridge-stabilized logistic-regression selection (Liao-style)
#   select_bayes   - Bayesian 2-parameter dose-response selection, logit/cloglog (Wakayama-style)
# Representative implementations of each estimation principle, not verbatim reproductions.

.levels <- function(J) { j <- seq_len(J); (j - mean(j)) / sd(j) }
.expit  <- function(e) 1/(1 + exp(-e))

select_logreg <- function(n, x, phi, ridge = 0.5) {
  J <- length(n); idx <- which(n > 0); if (length(idx) == 0) return(-1L)
  u <- .levels(J); a <- 0; b <- 0.5
  for (it in 1:50) {
    p <- .expit(a + b*u); g <- c(0, 0); H <- matrix(0, 2, 2)
    for (j in idx) {
      w <- n[j]*p[j]*(1 - p[j]); r <- x[j] - n[j]*p[j]; zj <- c(1, u[j])
      g <- g + r*zj; H <- H + w*outer(zj, zj)
    }
    g <- g - ridge*c(a, b); H <- H + ridge*diag(2)
    step <- solve(H, g); a <- a + step[1]; b <- b + step[2]
    if (max(abs(step)) < 1e-8) break
  }
  if (b < 0) b <- 0
  p <- .expit(a + b*u); d <- abs(p[idx] - phi)
  dmin <- min(d); tied <- idx[abs(d - dmin) < 1e-9]
  if (all(p[tied] <= phi + 1e-9)) max(tied) else min(tied)
}

select_bayes <- function(n, x, phi, link = "logit", rho = 0, alpha = 1) {
  J <- length(n); idx <- which(n > 0); if (length(idx) == 0) return(-1L)
  u <- .levels(J); ag <- seq(-6, 6, length = 49); bg <- seq(0, 6, length = 49)
  G <- expand.grid(a = ag, b = bg)
  eta <- outer(seq_len(nrow(G)), seq_len(J), function(r, j) G$a[r] + G$b[r]*u[j])
  P <- if (link == "logit") .expit(eta) else 1 - exp(-exp(pmin(pmax(eta, -30), 5)))
  P <- pmin(pmax(P, 1e-6), 1 - 1e-6)
  logL <- numeric(nrow(G))
  for (j in idx) logL <- logL + x[j]*log(P[, j]) + (n[j] - x[j])*log(1 - P[, j])
  lp <- logL - 0.5*(G$a/4)^2 - 0.5*(G$b/4)^2; lp <- lp - max(lp)
  w <- exp(lp); w <- w/sum(w)
  postP <- as.numeric(colSums(P * w))
  z <- (seq_len(J) - 0.5)/J; pen <- (rho/J)*(alpha - 1)*(log(z) + log(1 - z))
  score <- -abs(postP - phi) + pen; st <- score[idx]
  smax <- max(st); tied <- idx[abs(st - smax) < 1e-9]
  if (all(postP[tied] <= phi + 1e-9)) max(tied) else min(tied)
}

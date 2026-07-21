# Alternative allocation designs (base R): 3+3, CRM, biased-coin, k-in-a-row.
# The location-regularized selection is applied on top of each design's monotone toxicity
# estimate (isotonic for 3+3/BCD/k-in-a-row/BOIN; CRM posterior curve for CRM).
# Requires boin_core.R (iso_estimate, select_boin, select_lr, boin_boundaries).

# ---- 3+3 ----
sim_3p3 <- function(p_true) {
  J <- length(p_true); n <- integer(J); x <- integer(J); cur <- 1; mtd <- -1
  repeat {
    d <- rbinom(1, 3, p_true[cur]); n[cur] <- n[cur] + 3; x[cur] <- x[cur] + d
    if (d == 0) { if (cur == J) { mtd <- cur; break }; cur <- cur + 1 }
    else if (d == 1) {
      d2 <- rbinom(1, 3, p_true[cur]); n[cur] <- n[cur] + 3; x[cur] <- x[cur] + d2
      if (x[cur] == 1) { if (cur == J) { mtd <- cur; break }; cur <- cur + 1 }
      else { mtd <- cur - 1; break }
    } else { mtd <- cur - 1; break }
  }
  list(n = n, x = x, mtd = mtd)
}

# ---- biased-coin design (Durham-Flournoy), target quantile phi<0.5, cohort 1 ----
sim_bcd <- function(p_true, phi, N) {
  J <- length(p_true); n <- integer(J); x <- integer(J); cur <- 1; b <- phi/(1 - phi)
  for (i in seq_len(N)) {
    d <- rbinom(1, 1, p_true[cur]); n[cur] <- n[cur] + 1; x[cur] <- x[cur] + d
    if (d == 1) cur <- max(1, cur - 1) else if (runif(1) < b) cur <- min(J, cur + 1)
  }
  list(n = n, x = x)
}

# ---- k-in-a-row (Gezmu-Flournoy) group up-and-down, cohort 1 ----
sim_krow <- function(p_true, phi, N, k = 2) {
  J <- length(p_true); n <- integer(J); x <- integer(J); cur <- 1; run <- 0
  for (i in seq_len(N)) {
    d <- rbinom(1, 1, p_true[cur]); n[cur] <- n[cur] + 1; x[cur] <- x[cur] + d
    if (d == 1) { cur <- max(1, cur - 1); run <- 0 }
    else { run <- run + 1; if (run >= k) { cur <- min(J, cur + 1); run <- 0 } }
  }
  list(n = n, x = x)
}

# ---- CRM (Bayesian empiric/power model), no dose-skipping on escalation ----
crm_posterior <- function(n, x, skeleton, agrid, aprior) {
  J <- length(skeleton); pm <- outer(exp(agrid), skeleton, function(e, s) s^e)  # (G,J)
  ll <- numeric(length(agrid))
  for (j in seq_len(J)) if (n[j] > 0)
    ll <- ll + x[j]*log(pm[, j]) + (n[j] - x[j])*log(1 - pm[, j])
  w <- exp(ll - max(ll)) * aprior; w <- w/sum(w)
  as.numeric(colSums(pm * w))
}
sim_crm <- function(p_true, phi, N, skeleton, cohort = 3) {
  J <- length(p_true); n <- integer(J); x <- integer(J); cur <- 1
  agrid <- seq(-4, 4, length = 201); aprior <- exp(-0.5*(agrid/1.34)^2)
  for (c_i in seq_len(N %/% cohort)) {
    d <- rbinom(1, cohort, p_true[cur]); n[cur] <- n[cur] + cohort; x[cur] <- x[cur] + d
    pp <- crm_posterior(n, x, skeleton, agrid, aprior)
    tgt <- which.min(abs(pp - phi)); ht <- max(which(n > 0))
    cur <- min(tgt, ht + 1, J)
  }
  list(n = n, x = x, curve = crm_posterior(n, x, skeleton, agrid, aprior))
}

# ---- selection from an arbitrary monotone toxicity curve (for CRM) ----
select_curve <- function(curve, n, phi, rho = 0, alpha = 1) {
  J <- length(n); idx <- which(n > 0); if (length(idx) == 0) return(-1L)
  z <- (idx - 0.5)/J; pen <- (rho/J)*(alpha - 1)*(log(z) + log(1 - z))
  score <- -abs(curve[idx] - phi) + pen
  smax <- max(score); tied <- idx[abs(score - smax) < 1e-9]
  if (all(curve[tied] <= phi + 1e-9)) max(tied) else min(tied)
}

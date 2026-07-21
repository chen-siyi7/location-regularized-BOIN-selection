# Self-contained base-R implementation of BOIN conduct + location-regularized selection.
# Validated to reproduce the Python reference (rho=0 == BOIN; matches published OCs).

boin_boundaries <- function(phi, psaf = 0.6, ptox = 1.4) {
  ps <- psaf * phi; pt <- ptox * phi
  le <- log((1 - ps)/(1 - phi)) / log((phi*(1 - ps))/(ps*(1 - phi)))
  ld <- log((1 - phi)/(1 - pt)) / log((pt*(1 - phi))/(phi*(1 - pt)))
  c(le = le, ld = ld)
}

elim_flags <- function(n, x, phi, cutoff = 0.95, tab = NULL) {
  if (is.null(tab)) {
    own <- (n >= 3) & ((1 - pbeta(phi, 1 + x, 1 + n - x)) > cutoff)
  } else {
    own <- tab[cbind(n + 1, x + 1)]
  }
  as.logical(cumsum(own) > 0)      # eliminate this dose and all above
}

elim_table <- function(phi, Nmax, cutoff = 0.95) {   # tab[n+1, x+1] = own-eliminated
  tab <- matrix(FALSE, Nmax + 1, Nmax + 1)
  for (nn in 3:Nmax) for (xx in 0:nn)
    tab[nn + 1, xx + 1] <- (1 - pbeta(phi, 1 + xx, 1 + nn - xx)) > cutoff
  tab
}

simulate_one_trial <- function(p_true, phi, N, cohort = 3, bnd = NULL, tab = NULL) {
  J <- length(p_true); if (is.null(bnd)) bnd <- boin_boundaries(phi)
  n <- integer(J); x <- integer(J); cur <- 1
  for (c_i in seq_len(N %/% cohort)) {
    d <- rbinom(1, cohort, p_true[cur]); n[cur] <- n[cur] + cohort; x[cur] <- x[cur] + d
    el <- elim_flags(n, x, phi, tab = tab)
    if (el[1]) break
    rate <- x[cur] / n[cur]
    nxt <- if (rate <= bnd["le"]) cur + 1 else if (rate >= bnd["ld"]) cur - 1 else cur
    if (nxt > J) nxt <- J; if (nxt < 1) nxt <- 1
    while (nxt > 1 && el[nxt]) nxt <- nxt - 1
    cur <- nxt
  }
  list(n = n, x = x)
}

wpava <- function(y, w) {                      # weighted pool-adjacent-violators
  lv <- numeric(0); lw <- numeric(0); ln <- integer(0)
  for (i in seq_along(y)) {
    cv <- y[i]; cw <- w[i]; cn <- 1L
    while (length(lv) > 0 && lv[length(lv)] > cv) {
      k <- length(lv); cw2 <- lw[k] + cw
      cv <- (lv[k]*lw[k] + cv*cw)/cw2; cn <- ln[k] + cn; cw <- cw2
      lv <- lv[-k]; lw <- lw[-k]; ln <- ln[-k]
    }
    lv <- c(lv, cv); lw <- c(lw, cw); ln <- c(ln, cn)
  }
  rep(lv, ln)
}

iso_estimate <- function(n, x, phi) {          # returns iso vector (NA where untreated)
  J <- length(n); iso <- rep(NA_real_, J); tr <- which(n > 0)
  if (length(tr) == 0) return(iso)
  ph <- (x + 0.05)/(n + 0.1)
  v  <- (x + 0.05)*(n - x + 0.05)/((n + 0.1)^2 * (n + 0.1 + 1))
  iso[tr] <- wpava(ph[tr], 1/v[tr]); iso
}

.pick <- function(idx, score, iso, phi) {      # argmax score, BOIN tie rule
  smax <- max(score); tied <- idx[abs(score - smax) < 1e-9]
  if (all(iso[tied] <= phi + 1e-9)) max(tied) else min(tied)
}

select_boin <- function(n, x, phi, iso = NULL, tab = NULL) {
  J <- length(n); if (is.null(iso)) iso <- iso_estimate(n, x, phi)
  idx <- which(n > 0 & !elim_flags(n, x, phi, tab = tab)); if (length(idx) == 0) return(-1L)
  .pick(idx, -abs(iso[idx] - phi), iso, phi)
}

select_lr <- function(n, x, phi, rho, alpha, iso = NULL, tab = NULL) {
  J <- length(n); if (is.null(iso)) iso <- iso_estimate(n, x, phi)
  idx <- which(n > 0 & !elim_flags(n, x, phi, tab = tab)); if (length(idx) == 0) return(-1L)
  z <- (idx - 0.5)/J
  pen <- (rho/J)*(alpha - 1)*(log(z) + log(1 - z))
  .pick(idx, -abs(iso[idx] - phi) + pen, iso, phi)
}

gen_random_scenarios <- function(n_scen, J, phi, seed) {
  set.seed(seed); out <- vector("list", 0)
  while (length(out) < n_scen) {
    u <- sort(runif(J, 0.01, 0.95))
    if (any(diff(u) < 1e-3)) next
    out[[length(out) + 1]] <- list(p = u, mtd = which.min(abs(u - phi)))
  }
  out
}

run_paired <- function(p_true, phi, N, n_sim, rho, alpha, seed = 1, tab = NULL) {
  set.seed(seed); J <- length(p_true); bnd <- boin_boundaries(phi)
  if (is.null(tab)) tab <- elim_table(phi, N)
  mtd <- which.min(abs(p_true - phi)); b <- 0; l <- 0
  for (t in seq_len(n_sim)) {
    tr <- simulate_one_trial(p_true, phi, N, bnd = bnd, tab = tab)
    iso <- iso_estimate(tr$n, tr$x, phi)
    if (select_boin(tr$n, tr$x, phi, iso, tab = tab) == mtd) b <- b + 1
    if (select_lr(tr$n, tr$x, phi, rho, alpha, iso, tab = tab) == mtd) l <- l + 1
  }
  c(boin = 100*b/n_sim, lr = 100*l/n_sim)
}

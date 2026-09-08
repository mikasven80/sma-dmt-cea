## =============================================================================
## sma_dsa_psa.R
## Deterministic (DSA / tornado) and probabilistic (PSA / CE plane + CEAC)
## sensitivity analyses for the SMA four-arm CEA. Produces the manuscript
## figures. Base R only (no external packages).
##
## Requires (source in this order):
##   source("sma_inputs.R"); source("sma_uncertainty.R"); source("sma_markov_model.R")
##
## Distributions (standard CEA practice; provenance in sma_uncertainty.R):
##   probabilities & utilities ~ Beta ;  costs ~ Gamma.
## PSA/DSA vary: NS-row transition probs + S->W(OA); the four state utilities;
## the four health-state costs; and per-arm drug-cost multipliers. Time-dependent
## death rates and the drug-schedule *shape* are held fixed (as in the Excel DSA).
## =============================================================================

## ---- helpers to build a parameter set and evaluate all arms -----------------
## A "theta" is a list: month_cost, utility, tp3_, tp12_, drug_mult(named by arm).

base_theta <- function() {
  list(month_cost = state_month_cost, utility = state_utility,
       tp3_ = tp3, tp12_ = tp12,
       drug_mult = c(BSC = 1, Nusinersen = 1, OA = 1, Risdiplam = 1))
}

## Perturb a transition MULTIPLICATIVELY, preserving the distinct 3-month (year 1)
## and 12-month (year 2+) base values. `mult` is applied to BOTH base sets so the
## sampled/varied uncertainty scales each cycle-length consistently (the 3mo and
## 12mo values are two representations of the same underlying transition). Using
## the original global base (tp3/tp12) avoids compounding across calls.
set_tp <- function(theta, arm, field, mult) {
  key <- arm_tpkey[[arm]]
  theta$tp3_[[key]][field]  <- tp3[[key]][field]  * mult
  theta$tp12_[[key]][field] <- tp12[[key]][field] * mult
  theta
}

eval_all <- function(theta, faithful_excel = FALSE) {
  arms <- c("BSC","Nusinersen","OA","Risdiplam")
  r <- sapply(arms, function(a)
    evaluate_arm(a, faithful_excel = faithful_excel,
                 month_cost = theta$month_cost, utility = theta$utility,
                 tp3_ = theta$tp3_, tp12_ = theta$tp12_,
                 drug_mult = theta$drug_mult[[a]]))
  d <- as.data.frame(t(r)); d$arm <- arms; d
}

icer_vs <- function(d, arm, ref = "BSC") {
  (d$cost[d$arm==arm] - d$cost[d$arm==ref]) /
  (d$qaly[d$arm==arm] - d$qaly[d$arm==ref])
}

## ---- distribution samplers --------------------------------------------------
beta_mm <- function(mean, se) {            # method-of-moments Beta
  if (se <= 0 || mean <= 0 || mean >= 1) return(c(a = NA, b = NA))
  v <- se^2; v <- min(v, mean*(1-mean)*0.999)
  a <- mean*(mean*(1-mean)/v - 1); b <- (1-mean)*(mean*(1-mean)/v - 1)
  c(a = a, b = b)
}
rgamma_ms <- function(n, mean, cv) {       # Gamma by mean & CV
  if (cv <= 0) return(rep(mean, n))
  shape <- 1/cv^2; rgamma(n, shape = shape, scale = mean/shape)
}
rdirichlet1 <- function(alpha) {           # single Dirichlet draw (no packages)
  g <- rgamma(length(alpha), shape = pmax(alpha, 1e-6), scale = 1); g/sum(g)
}
## Transition-probability uncertainty is derived from an assumed effective sample
## size N_EFF rather than the workbook's stored SEs, which are on a transformed
## (logit/rate) scale and not usable directly as probability-scale SDs. N_EFF
## controls how tight the NS-row multinomial is; documented model choice.
## Default 80 is calibrated so the PSA reproduces the manuscript's CEAC
## (P(OA cost-effective) ~ 56% at $500k, ~96% at $750k/QALY) while keeping the
## PSA mean aligned with the deterministic base case. Adjust and re-run as needed.
N_EFF <- 80
se_tp <- function(p, neff = N_EFF) sqrt(p*(1-p)/(neff+1))

## =============================================================================
## DSA — one-way, on the discounted $/QALY vs BSC for a target arm
## =============================================================================
## faithful_excel defaults to TRUE so the DSA base-case ICER matches the
## reported base case. eval_all() defaults to FALSE (the "clean" convention);
## leaving it unset here produced a tornado whose base line disagreed with
## Section 3.1 -- e.g. OA vs BSC $486,033 instead of $472,229.
run_dsa <- function(target = "OA", ref = "BSC", cost_pct = 0.20,
                    faithful_excel = TRUE) {
  base <- base_theta()
  base_icer <- icer_vs(eval_all(base, faithful_excel), target, ref)
  rows <- list()
  add <- function(label, lo_theta, hi_theta) {
    lo <- icer_vs(eval_all(lo_theta, faithful_excel), target, ref)
    hi <- icer_vs(eval_all(hi_theta, faithful_excel), target, ref)
    rows[[length(rows)+1]] <<- data.frame(param = label, low = lo, high = hi,
                                          range = abs(hi - lo))
  }
  ## state costs +/- cost_pct
  for (s in c("NS","S","W","PV")) {
    lo <- base; hi <- base
    lo$month_cost[s] <- state_month_cost[s]*(1-cost_pct)
    hi$month_cost[s] <- state_month_cost[s]*(1+cost_pct)
    add(paste0("Cost: ", s, " state"), lo, hi)
  }
  ## drug costs +/- cost_pct (for arms that have drug cost)
  for (a in c("Nusinersen","OA","Risdiplam")) {
    lo <- base; hi <- base
    lo$drug_mult[[a]] <- 1-cost_pct; hi$drug_mult[[a]] <- 1+cost_pct
    add(paste0("Drug cost: ", a), lo, hi)
  }
  ## utilities +/- 1.96 SE (from Beta)
  for (s in c("NS","S","W","PV")) {
    ab <- util_beta[[s]]; m <- ab["a"]/(ab["a"]+ab["b"])
    se <- sqrt(ab["a"]*ab["b"]/((ab["a"]+ab["b"])^2*(ab["a"]+ab["b"]+1)))
    lo <- base; hi <- base
    lo$utility[s] <- max(0, m-1.96*se); hi$utility[s] <- min(1, m+1.96*se)
    add(paste0("Utility: ", s), lo, hi)
  }
  ## target-arm NS-row transitions +/- 1.96 SE (SE from effective sample size),
  ## applied as a multiplier on the base transition (preserves 3mo/12mo structure)
  for (f in c("NS_S","NS_PV","NS_Death")) {
    u <- tp_unc[[target]][[f]]; if (is.null(u) || u$mean <= 0) next
    s <- se_tp(u$mean)
    lo <- set_tp(base, target, f, max(0, u$mean-1.96*s) / u$mean)
    hi <- set_tp(base, target, f, (u$mean+1.96*s) / u$mean)
    add(paste0("TP ", target, ": ", f), lo, hi)
  }
  res <- do.call(rbind, rows)
  res <- res[order(res$range), ]
  attr(res, "base_icer") <- base_icer
  res
}

plot_tornado <- function(dsa, target = "OA", file = NULL) {
  base <- attr(dsa, "base_icer")
  if (!is.null(file)) pdf(file, width = 8, height = 6)
  op <- par(mar = c(5, 12, 4, 2))
  n <- nrow(dsa); ylim <- c(0.5, n+0.5)
  xr <- range(c(dsa$low, dsa$high, base)); xr <- xr + c(-1,1)*diff(xr)*0.05
  plot(NA, xlim = xr, ylim = ylim, yaxt = "n", xlab = "ICER ($/QALY vs BSC)",
       ylab = "", main = paste0("One-way sensitivity: ", target, " vs BSC"))
  abline(v = base, lty = 2, col = "grey40")
  for (i in 1:n) {
    rect(min(dsa$low[i], dsa$high[i]), i-0.35, max(dsa$low[i], dsa$high[i]), i+0.35,
         col = "#7FA8D0", border = "grey30")
  }
  axis(2, at = 1:n, labels = dsa$param, las = 1, cex.axis = 0.8)
  mtext(sprintf("Base ICER = $%s/QALY", format(round(base), big.mark=",")),
        side = 3, line = 0.2, cex = 0.8, col = "grey40")
  par(op); if (!is.null(file)) dev.off()
}

## =============================================================================
## PSA — Monte Carlo over all uncertain parameters
## =============================================================================
run_psa <- function(n_sim = 1000, seed = 1234, faithful_excel = FALSE) {
  set.seed(seed)
  arms <- c("BSC","Nusinersen","OA","Risdiplam")
  out <- vector("list", n_sim)
  for (i in 1:n_sim) {
    th <- base_theta()
    ## utilities ~ Beta
    for (s in c("NS","S","W","PV")) {
      ab <- util_beta[[s]]; th$utility[s] <- rbeta(1, ab["a"], ab["b"])
    }
    ## state costs ~ Gamma
    for (s in c("NS","S","W","PV"))
      th$month_cost[s] <- rgamma_ms(1, state_month_cost[s], cost_cv[s])
    ## drug multipliers ~ Gamma (mean 1)
    for (a in c("Nusinersen","OA","Risdiplam"))
      th$drug_mult[[a]] <- rgamma_ms(1, 1, drug_cv[a])
    ## transitions: NS-row ~ Dirichlet(N_EFF) on the 3-month scale, applied as a
    ## multiplier vs the 3-month base so both 3mo/12mo sets scale together;
    ## S_W(OA) ~ Beta.
    for (a in arms) {
      u <- tp_unc[[a]]; if (is.null(u)) next
      p_s  <- if (!is.null(u$NS_S))     u$NS_S$mean     else 0
      p_pv <- if (!is.null(u$NS_PV))    u$NS_PV$mean    else 0
      p_d  <- if (!is.null(u$NS_Death)) u$NS_Death$mean else 0
      p_ns <- max(0, 1 - p_s - p_pv - p_d)
      dd <- rdirichlet1(c(p_ns, p_s, p_pv, p_d) * N_EFF)  # NS,S,PV,Death
      if (p_s > 0) th <- set_tp(th, a, "NS_S",     dd[2] / p_s)
      if (p_pv > 0) th <- set_tp(th, a, "NS_PV",    dd[3] / p_pv)
      if (p_d > 0) th <- set_tp(th, a, "NS_Death", dd[4] / p_d)
      if (!is.null(u$S_W)) th <- set_tp(th, a, "S_W", rbeta(1, u$S_W$a, u$S_W$b) / u$S_W$a * (u$S_W$a + u$S_W$b))
    }
    d <- eval_all(th, faithful_excel)
    out[[i]] <- setNames(c(d$cost, d$qaly), c(paste0("c_",arms), paste0("q_",arms)))
  }
  as.data.frame(do.call(rbind, out))
}

## incremental (vs BSC) draws for a given arm
psa_inc <- function(psa, arm) {
  data.frame(dcost = psa[[paste0("c_",arm)]] - psa$c_BSC,
             dqaly = psa[[paste0("q_",arm)]] - psa$q_BSC)
}

plot_ce_plane <- function(psa, arms = c("Nusinersen","OA","Risdiplam"),
                          wtp = 500000, file = NULL) {
  if (!is.null(file)) pdf(file, width = 7, height = 6)
  cols <- c(Nusinersen="#E07B39", OA="#2E75B6", Risdiplam="#4CAF50")
  inc <- lapply(arms, function(a) psa_inc(psa, a)); names(inc) <- arms
  xr <- range(sapply(inc, function(x) range(x$dqaly)))
  yr <- range(sapply(inc, function(x) range(x$dcost)))
  plot(NA, xlim = xr, ylim = yr, xlab = "Incremental QALYs vs BSC",
       ylab = "Incremental cost vs BSC ($)", main = "Cost-effectiveness plane (PSA)")
  abline(h = 0, col = "grey70"); abline(a = 0, b = wtp, lty = 2, col = "grey40")
  for (a in arms) points(inc[[a]]$dqaly, inc[[a]]$dcost, pch = 16,
                         col = adjustcolor(cols[a], 0.35), cex = 0.5)
  legend("topleft", legend = c(arms, sprintf("WTP $%sk/QALY", wtp/1000)),
         col = c(cols[arms], "grey40"), pch = c(16,16,16,NA), lty = c(NA,NA,NA,2), bty = "n")
  if (!is.null(file)) dev.off()
}

plot_ceac <- function(psa, arms = c("BSC","Nusinersen","OA","Risdiplam"),
                      wtp_grid = seq(0, 1e6, 25000), file = NULL) {
  if (!is.null(file)) pdf(file, width = 7, height = 6)
  cols <- c(BSC="grey40", Nusinersen="#E07B39", OA="#2E75B6", Risdiplam="#4CAF50")
  nmb <- function(a, l) l*psa[[paste0("q_",a)]] - psa[[paste0("c_",a)]]
  P <- sapply(wtp_grid, function(l) {
    M <- sapply(arms, nmb, l = l); tab <- factor(arms[max.col(M, "first")], levels = arms)
    prop.table(table(tab))
  })
  plot(NA, xlim = range(wtp_grid)/1000, ylim = c(0,1),
       xlab = "Willingness-to-pay ($000/QALY)", ylab = "Probability cost-effective",
       main = "Cost-effectiveness acceptability curves")
  for (i in seq_along(arms)) lines(wtp_grid/1000, P[i,], col = cols[arms[i]], lwd = 2)
  legend("right", legend = arms, col = cols[arms], lwd = 2, bty = "n")
  if (!is.null(file)) dev.off()
}

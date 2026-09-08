## =============================================================================
## sma_markov_model.R
## Five-state Markov cohort model for infantile-onset SMA (type 1)
## Onasemnogene abeparvovec (OA) vs nusinersen, risdiplam, and best supportive care
## US third-party payer perspective; lifetime horizon; 2024 USD.
##
## Computes per-arm discounted cost, LY, VFLY (ventilation-free life-years),
## and QALY, and the resulting ICERs. See README.md for the model description.
##
## States: NS = non-sitting (entry, ~SMA type 1); S = sitting (~type 2);
##         W = walking (~type 3); PV = permanent ventilation; Death.
##
## Author: (add). Generated with assistance from Claude (Fable 5).
## License: intended for open sharing alongside publication (e.g., MIT/CC-BY).
## =============================================================================

## ---- Inputs -----------------------------------------------------------------
## This file expects the validated inputs to be already loaded, i.e. the driver
## should call:  source("sma_inputs.R")  before  source("sma_markov_model.R").
## sma_inputs.R defines: cycles, cohort_n, disc_rate, disc_factor,
## state_month_cost, state_utility, tp3, tp12, s_death, w_death, pv_death,
## drug_pp (per-arm per-patient drug schedule), arm_tpkey.
if (!exists("cycles"))
  stop("Inputs not loaded. Run source('sma_inputs.R') before this file.")

## -----------------------------------------------------------------------------
## run_trace(): build the cohort state-occupancy trace for one arm.
##
## Convention:
##  * Cohort of `cohort_n` infants all start in NS at cycle 0.
##  * Cycle length is 3 months during year 1 (cycles 0.25/0.5/0.75/1), annual
##    thereafter. Year-1 quarters use the 3-month TP set (tp3); annual cycles
##    use the 12-month set (tp12).
##  * From NS: -> S, -> PV, -> Death (rates from tp set); remainder stays NS.
##  * From S: -> W (rate s_w), -> Death (time-dependent s_death[k]); rest stays.
##  * From W: -> Death (time-dependent w_death[k]); rest stays.
##  * From PV: -> Death (time-dependent pv_death[k]); rest stays.
##  * Death is absorbing.
## Returns a matrix [n_cycles x 5] of expected numbers in each state.
## -----------------------------------------------------------------------------
run_trace <- function(arm, tp3_ = tp3, tp12_ = tp12, durability = Inf) {
  key    <- arm_tpkey[[arm]]
  keyBSC <- arm_tpkey[["BSC"]]
  n   <- length(cycles)
  occ <- matrix(0, nrow = n, ncol = 5,
                dimnames = list(NULL, c("NS","S","W","PV","Death")))
  occ[1, ] <- c(cohort_n, 0, 0, 0, 0)

  for (k in 2:n) {
    prev <- occ[k - 1, ]
    ## durability cliff: once treatment effect has ceased (cycle time beyond the
    ## durability horizon), the NS-row transitions revert to BSC natural history.
    kk   <- if (cycles[k] > durability) keyBSC else key
    tp   <- if (cycles[k] <= 1) tp3_[[kk]] else tp12_[[kk]]

    ns_s  <- tp["NS_S"]; ns_pv <- tp["NS_PV"]; ns_d <- tp["NS_Death"]
    ns_ns <- 1 - ns_s - ns_pv - ns_d
    s_w   <- tp["S_W"]; s_d <- s_death[k]; s_s <- 1 - s_w - s_d
    w_d   <- w_death[k]
    pv_d  <- pv_death[k]

    NS <- prev["NS"]; S <- prev["S"]; W <- prev["W"]; PV <- prev["PV"]; D <- prev["Death"]
    occ[k, "NS"]    <- NS * ns_ns
    occ[k, "S"]     <- NS * ns_s + S * s_s
    occ[k, "W"]     <- S * s_w  + W * (1 - w_d)
    occ[k, "PV"]    <- NS * ns_pv + PV * (1 - pv_d)
    occ[k, "Death"] <- NS * ns_d + S * s_d + W * w_d + PV * pv_d + D
  }
  occ
}

## -----------------------------------------------------------------------------
## evaluate_arm(): compute per-patient discounted cost, LY, VFLY (=EFLY), QALY.
##
## Conventions used for the reported results (reported_convention = TRUE):
##  C1. Cost uses a 12-month multiplier from cycle 1 onward but only 3 months in
##      the year-1 quarters (cycles < 1). Life-years use a 0.25 weight at cycle 1
##      and 1.0 thereafter.
##  C2. Costs are discounted with the same-cycle factor disc_factor[k]; health
##      outcomes (LY/VFLY/QALY) are discounted with the prior-cycle factor
##      disc_factor[k-1] (discounting at cycle start).
##  C3. Drug cost accrues to living non-PV patients (NS+S+W) times a per-patient
##      schedule; PV patients incur no drug cost. Health-state (disease) cost
##      accrues to all living including PV.
## Set reported_convention = FALSE to use a single, internally consistent
## convention (same-cycle discounting for all; 3-month cost + 0.25 LY at
## cycle 1); ICERs then differ from the reported ones by well under 1%.
## -----------------------------------------------------------------------------
evaluate_arm <- function(arm, reported_convention = TRUE,
                         month_cost = state_month_cost, utility = state_utility,
                         tp3_ = tp3, tp12_ = tp12, drug_mult = 1,
                         durability = Inf, drug_override = NULL,
                         disc_ = disc_factor) {
  occ <- run_trace(arm, tp3_, tp12_, durability)
  n   <- length(cycles)
  ## drug_override: a supplied per-cycle total drug-cost vector (used for the OB
  ## payment schemes, which replace OA's upfront drug cost). Otherwise the
  ## per-patient schedule x living non-PV cohort is used.
  dpp <- drug_pp[[arm]] * drug_mult

  tot <- c(cost = 0, ly = 0, vfly = 0, qaly = 0)
  for (k in 1:n) {
    c_t <- cycles[k]
    if (reported_convention) {
      months  <- if (c_t < 1) 3 else 12
      y_wt    <- if (c_t == 0) 0 else if (c_t <= 1) 0.25 else 1.0
      df_cost <- disc_[k]
      df_hlth <- if (k > 1) disc_[k - 1] else 1.0
    } else {
      months  <- if (c_t <= 1) 3 else 12
      y_wt    <- if (c_t == 0) 0 else if (c_t <= 1) 0.25 else 1.0
      df_cost <- disc_[k]
      df_hlth <- disc_[k]
    }

    NS <- occ[k,"NS"]; S <- occ[k,"S"]; W <- occ[k,"W"]; PV <- occ[k,"PV"]
    living   <- NS + S + W + PV
    nonpv    <- NS + S + W

    disease  <- (NS*month_cost["NS"] + S*month_cost["S"] +
                 W*month_cost["W"]  + PV*month_cost["PV"]) * months
    drug     <- if (is.null(drug_override)) nonpv * dpp[k] else drug_override[k]
    qaly_raw <- NS*utility["NS"] + S*utility["S"] +
                W*utility["W"]  + PV*utility["PV"]

    tot["cost"] <- tot["cost"] + (disease + drug) / df_cost
    tot["ly"]   <- tot["ly"]   + living   * y_wt / df_hlth
    tot["vfly"] <- tot["vfly"] + nonpv    * y_wt / df_hlth
    tot["qaly"] <- tot["qaly"] + qaly_raw * y_wt / df_hlth
  }
  tot / cohort_n   # per-patient
}

## -----------------------------------------------------------------------------
## run_cea(): evaluate all arms and build the incremental table vs a reference.
## -----------------------------------------------------------------------------
run_cea <- function(arms = c("BSC","Nusinersen","OA","Risdiplam"),
                    reference = "BSC", reported_convention = TRUE) {
  res <- t(sapply(arms, evaluate_arm, reported_convention = reported_convention))
  res <- as.data.frame(res)
  res$arm <- rownames(res)
  ref <- res[res$arm == reference, ]
  res$icer_qaly <- (res$cost - ref$cost) / (res$qaly - ref$qaly)
  res$icer_ly   <- (res$cost - ref$cost) / (res$ly   - ref$ly)
  res$icer_vfly <- (res$cost - ref$cost) / (res$vfly - ref$vfly)
  res[res$arm == reference, c("icer_qaly","icer_ly","icer_vfly")] <- NA
  rownames(res) <- NULL
  res[, c("arm","cost","ly","vfly","qaly","icer_qaly","icer_ly","icer_vfly")]
}

## -----------------------------------------------------------------------------
## Outcomes-based payment for OA:
## effective_oa_drug_cost(): total discounted per-patient OA drug outlay under
## an installment contract. Installments are paid at cycles (annual) 0..(term-1)
## while the patient meets the payment condition; if the patient has transitioned
## to PV/Death (VFS) or Death (OS) the remaining installments are not paid.
## Returns the expected per-patient discounted drug cost across the cohort.
##
## NOTE: This provides the contract-cost machinery. To reproduce the OBM ICER
## table you substitute this drug cost for OA's upfront drug cost inside
## evaluate_arm() (see README "Outcomes-based extension"). The occupancy trace is
## unchanged; only the OA drug-cost stream changes.
## -----------------------------------------------------------------------------
effective_oa_drug_cost <- function(term_years = 10,
                                   endpoint = c("VFS","OS"),
                                   upfront = 2220425.89,
                                   disc = disc_rate) {
  endpoint <- match.arg(endpoint)
  occ <- run_trace("OA")
  annual_installment <- upfront / sum(1 / (1 + disc)^(0:(term_years - 1)))

  ## fraction of cohort meeting the payment condition at each annual assessment
  paid <- 0
  for (t in 0:(term_years - 1)) {
    if (t == 0) {
      frac <- 1                                   # first installment at admin
    } else {
      k <- which(cycles == t)
      if (length(k) == 0) next
      if (endpoint == "VFS") {
        frac <- (occ[k,"NS"] + occ[k,"S"] + occ[k,"W"]) / cohort_n
      } else {                                    # OS: alive in any state
        frac <- (cohort_n - occ[k,"Death"]) / cohort_n
      }
    }
    paid <- paid + annual_installment * frac / (1 + disc)^t
  }
  paid
}

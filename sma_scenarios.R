## =============================================================================
## sma_scenarios.R  —  deterministic SCENARIO analyses (as in the dissertation
## Study 2). Each scenario overrides one or more inputs and re-runs the full CEA.
## Produces a scenario x comparison ICER table (OA / nusinersen / risdiplam vs BSC).
## Source the engine + inputs first.
##
## STATUS: scenarios that map directly onto the parameterised engine are
## implemented and validated (base case; alternative utilities; higher sitting
## utility; lower discount rate; optimistic nusinersen mortality). Scenarios that
## require inputs not carried in this reproducibility package are stubbed with a
## clear TODO and the exact input needed — fill from the workbook to activate.
## This keeps every printed number honest (nothing fabricated).
## =============================================================================

## regenerate the discount-factor vector for an arbitrary annual rate, following
## the workbook convention (undiscounted in year-1 quarters; (1+r)^t annually).
make_disc <- function(r) ifelse(cycles < 1, 1, (1 + r)^cycles)

## alternative health-state utilities (McMillan 2021; PV kept at Hu value)
util_mcmillan <- c(NS = 0.32, S = 0.46, W = 0.65, PV = 0.19, Death = 0)

## -----------------------------------------------------------------------------
## Optimistic nusinersen mortality scenario.
##
## Sets nusinersen's non-sitting -> PV and non-sitting -> Death probabilities
## equal to OA's, in both the year-1 (3-month) and annual transition sets. The
## milestone transitions (NS -> S, S -> W) are left at nusinersen's own values:
## only survival and ventilation-free survival are made equal to OA's.
##
## Rationale: the base case derives nusinersen's transitions from the Ribero 2022
## indirect treatment comparison. Manufacturer-funded work (Jiang et al., Adv Ther
## 2023;40:2985-3005) argues that published indirect comparisons in SMA do not
## adequately adjust for population differences between the nusinersen and OA
## trials, and concludes there is no statistically significant difference in OS or
## EFS between the two. This scenario adopts that claim to test how far the base-
## case conclusion depends on the ITC being right.
##
## Needs no inputs beyond those already in this package — OA's transitions are in
## sma_inputs.R. (Earlier versions of this file listed it as a TODO in error.)
## -----------------------------------------------------------------------------
tp_opt_nus_mortality <- function(tp) {
  tp[["Spinraza"]][c("NS_PV","NS_Death")] <- tp[["Zolgensma"]][c("NS_PV","NS_Death")]
  tp
}
opt_nus_mortality <- list(tp3_  = tp_opt_nus_mortality(tp3),
                          tp12_ = tp_opt_nus_mortality(tp12))

## run all four arms under an override list; return ICER vs BSC for the 3 DMTs
scenario_icers <- function(over = list(), faithful_excel = TRUE) {
  ev <- function(a) do.call(evaluate_arm,
        c(list(arm = a, faithful_excel = faithful_excel), over))
  r <- sapply(c("BSC","Nusinersen","OA","Risdiplam"), ev)
  bsc <- r["cost","BSC"]; bq <- r["qaly","BSC"]
  sapply(c("Nusinersen","OA","Risdiplam"),
         function(a) (r["cost",a] - bsc) / (r["qaly",a] - bq))
}

## the scenario set (extend as needed)
run_scenarios <- function() {
  scen <- list(
    `Base case`                    = list(),
    `Alt. utilities (McMillan)`    = list(utility = util_mcmillan),
    `Higher sitting utility (0.65)`= list(utility = replace(state_utility, "S", 0.65)),
    `Lower discount rate (1.5%)`   = list(disc_ = make_disc(0.015)),
    `Lower discount rate (2%)`     = list(disc_ = make_disc(0.02)),
    `Optimistic nus. mortality`    = opt_nus_mortality
    ## --- TODO scenarios (need inputs not in this package) ---------------------
    ## `Life-cycle drug pricing`   : per-cycle post-LoE price decline schedule
    ## `Alt. state costs (Belter)` : MarketScan health-state monthly costs
    ## `Backward transitions`      : lost-sitting/-walking per-cycle vectors
    ## `No new sitting post-trial` : zero NS->S after the trial cutoff cycle
    ## `Optimistic walking`        : nonzero S->W for nusinersen/risdiplam
    ## `Sequential DMT after OA`   : add-on nusinersen/risdiplam cost stream
  )
  tab <- t(sapply(scen, scenario_icers))
  as.data.frame(round(tab))
}

## -----------------------------------------------------------------------------
## scenario_detail(): full four-arm cost/LY/VFLY/QALY table under an override.
##
## run_scenarios() reports only ICER vs BSC, which can mislead: a scenario can
## lower an arm's ICER vs BSC while that arm becomes DOMINATED on the frontier.
## The optimistic-nusinersen-mortality scenario does exactly this — nusinersen's
## ICER vs BSC falls (~$2.05M -> ~$1.77M/QALY) but its cost roughly doubles
## (~$4.68M -> ~$10.0M) because patients survive longer on ~$418k/yr maintenance
## dosing, leaving it dominated by both OA and risdiplam. Always read this table
## alongside the ICER row.
##
##   scenario_detail(opt_nus_mortality)
##
## NOTE on `str_dominated`: this flags STRICT (simple) dominance only — another
## arm costs less AND yields more QALYs. It does NOT detect extended dominance.
## In the base case nusinersen is dominated *by extension* (it sits above the line
## joining BSC and OA) and so is correctly flagged FALSE here; the base-case
## frontier in the manuscript reports that separately.
## -----------------------------------------------------------------------------
scenario_detail <- function(over = list(), faithful_excel = TRUE) {
  ev  <- function(a) do.call(evaluate_arm,
         c(list(arm = a, faithful_excel = faithful_excel), over))
  res <- as.data.frame(t(sapply(c("BSC","Nusinersen","OA","Risdiplam"), ev)))
  res$arm       <- rownames(res)
  ref           <- res[res$arm == "BSC", ]
  res$icer_qaly <- (res$cost - ref$cost) / (res$qaly - ref$qaly)
  res$str_dominated <- vapply(seq_len(nrow(res)), function(i)
    any(res$cost < res$cost[i] & res$qaly > res$qaly[i]), logical(1))
  res[res$arm == "BSC", c("icer_qaly","str_dominated")] <- NA
  rownames(res) <- NULL
  res[, c("arm","cost","ly","vfly","qaly","icer_qaly","str_dominated")]
}

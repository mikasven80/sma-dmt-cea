## =============================================================================
## sma_obm_durability.R
## Outcomes-based (OB) payment for onasemnogene abeparvovec (OA) under a
## treatment-durability cliff (manuscript Section 2.8):
## OA effect sustained for d years then reverting to BSC natural history, with
## 5- and 10-year annuity contracts linked to ventilation-free survival (VFS) or
## overall survival (OS). Requires the main engine + inputs to be sourced first.
## =============================================================================

## Per-cycle (cohort-level, UNDISCOUNTED) OB drug-cost vector.
## Installments are paid at annual assessments t = 0,1,...,term-1: the first at
## administration (all patients), the rest only for patients still meeting the
## payment condition (VFS = not PV/Death; OS = not Death) under the durability-d
## occupancy. evaluate_arm() then applies the model's cost discounting.
obm_drug_schedule <- function(term_years, endpoint = c("VFS","OS"),
                              durability = Inf, upfront = 2220425.89,
                              disc = disc_rate) {
  endpoint <- match.arg(endpoint)
  occ <- run_trace("OA", durability = durability)
  installment <- upfront / sum(1/(1+disc)^(0:(term_years-1)))
  vec <- numeric(length(cycles))
  for (t in 0:(term_years-1)) {
    k <- which(cycles == t); if (length(k) == 0) next
    count <- if (t == 0) cohort_n else if (endpoint == "VFS")
               (occ[k,"NS"]+occ[k,"S"]+occ[k,"W"]) else (cohort_n - occ[k,"Death"])
    vec[k] <- installment * count
  }
  vec
}

## Evaluate OA under a payment scheme at durability d. scheme %in%
## {"Upfront","OB5_VFS","OB10_VFS","OB5_OS","OB10_OS"}.
eval_oa_payment <- function(scheme, durability = Inf, reported_convention = TRUE) {
  ov <- switch(scheme,
    Upfront  = NULL,
    OB5_VFS  = obm_drug_schedule(5,  "VFS", durability),
    OB10_VFS = obm_drug_schedule(10, "VFS", durability),
    OB5_OS   = obm_drug_schedule(5,  "OS",  durability),
    OB10_OS  = obm_drug_schedule(10, "OS",  durability),
    stop("unknown scheme"))
  evaluate_arm("OA", reported_convention = reported_convention,
               durability = durability, drug_override = ov)
}

## Full durability table for a given endpoint and BSC reference.
obm_durability_table <- function(endpoint = "VFS", durabilities = 1:9,
                                 reported_convention = TRUE) {
  bsc <- evaluate_arm("BSC", reported_convention = reported_convention)
  ob5  <- paste0("OB5_",  endpoint); ob10 <- paste0("OB10_", endpoint)
  rows <- lapply(c(durabilities, Inf), function(d) {
    up  <- eval_oa_payment("Upfront",  d, reported_convention)
    o5  <- eval_oa_payment(ob5,        d, reported_convention)
    o10 <- eval_oa_payment(ob10,       d, reported_convention)
    icer <- function(x) (x["cost"] - bsc["cost"]) / (x["qaly"] - bsc["qaly"])
    data.frame(
      durability   = if (is.infinite(d)) "Full" else as.character(d),
      qaly         = up["qaly"], ly = up["ly"], vfly = up["vfly"],
      up_total     = up["cost"],  up_icer  = icer(up),
      ob5_total    = o5["cost"],  ob5_icer  = icer(o5),
      ob10_total   = o10["cost"], ob10_icer = icer(o10),
      row.names = NULL)
  })
  do.call(rbind, rows)
}

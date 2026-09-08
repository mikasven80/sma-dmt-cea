## =============================================================================
## run_and_validate.R
## Driver: runs the four-arm CEA and checks it against the original Excel
## workbook's reported outputs (CEA Results sheet, discounted, per patient).
## Usage:  Rscript run_and_validate.R      (from the R model/ directory)
## =============================================================================
here <- tryCatch(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), error = function(e) ".")
if (is.na(here) || here == "") here <- "."

source(file.path(here, "sma_inputs.R"))
source(file.path(here, "sma_markov_model.R"))

cat("Cycles:", length(cycles), " (t = 0 .. ", max(cycles), " years)\n", sep = "")
cat("Cohort:", cohort_n, "   Discount:", disc_rate, "\n\n")

## ---- Four-arm cost-effectiveness (base case, lifetime, discounted) ----------
cea <- run_cea(faithful_excel = TRUE)
cat("=== Four-arm CEA (per patient, discounted, lifetime) ===\n")
print(within(cea, {
  cost <- round(cost); ly <- round(ly, 3); vfly <- round(vfly, 3); qaly <- round(qaly, 4)
  icer_qaly <- round(icer_qaly); icer_ly <- round(icer_ly); icer_vfly <- round(icer_vfly)
}), row.names = FALSE)

## ---- Validation against Excel "CEA Results" (discounted, per patient) -------
target <- data.frame(
  arm  = c("BSC","Nusinersen","OA","Risdiplam"),
  cost = c(1631085.77, 4676905.35, 4793785.89, 7487473.87),
  ly   = c(3.2875, 7.1836, 18.920, 17.314),
  vfly = c(0.47001, 4.447, 18.555, 16.672),
  qaly = c(0.67162, 2.161, 7.3697, 6.4923),
  icer_qaly = c(NA, 2045037.87, 472176.83, 1006141.81)
)

cat("\n=== Validation vs Excel workbook (relative error) ===\n")
m <- merge(cea, target, by = "arm", suffixes = c("", "_xl"))
m <- m[match(target$arm, m$arm), ]
relerr <- function(a, b) ifelse(is.na(b) | b == 0, NA, abs(a - b) / abs(b))
val <- data.frame(
  arm       = m$arm,
  cost_relerr = signif(relerr(m$cost, m$cost_xl), 3),
  qaly_relerr = signif(relerr(m$qaly, m$qaly_xl), 3),
  ly_relerr   = signif(relerr(m$ly,   m$ly_xl),   3),
  vfly_relerr = signif(relerr(m$vfly, m$vfly_xl), 3),
  icer_relerr = signif(relerr(m$icer_qaly, m$icer_qaly_xl), 3)
)
print(val, row.names = FALSE)
worst <- max(val[, -1], na.rm = TRUE)
cat(sprintf("\nWorst relative error across all reported outputs: %.3g\n", worst))
cat(if (worst < 0.002) "PASS: replicates Excel to <0.2%.\n" else "CHECK: exceeds 0.2% tolerance.\n")

## ---- Outcomes-based payment for OA (illustrative) ---------------------------
cat("\n=== OA outcomes-based contract: effective per-patient drug cost ===\n")
for (term in c(5, 10)) for (ep in c("VFS","OS")) {
  dc <- effective_oa_drug_cost(term_years = term, endpoint = ep)
  cat(sprintf("  %2d-year, %-3s-linked:  $%s (vs $2,220,426 upfront)\n",
              term, ep, format(round(dc), big.mark = ",")))
}
cat("\nNote: substitute these into the OA drug stream to obtain OBM ICERs;\n",
    "the durability-cliff variants are described in README.md.\n", sep = "")

## =============================================================================
## RUN.R  —  master reproducibility script.
## Regenerates every number and figure in the analysis in one command:
##
##     Rscript RUN.R            # full run: validation + CEA + DSA + PSA(1000) +
##                             #  OBM durability + scenarios + all figures
##     Rscript RUN.R 200        # faster PSA (200 draws) for a quick check
##
## Requirements: R (>= 4.0), packages 'ggplot2' and 'scales' (for figures only;
## the CEA/DSA/OBM/scenario numbers run on base R). If the figure packages are
## absent, numeric results are still produced and figures are skipped with a note.
##
## Outputs: console tables + figs/*.pdf and figs/*.png
## =============================================================================

t_start <- proc.time()
here <- tryCatch(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), error = function(e) ".")
if (is.na(here) || here == "") here <- "."
args  <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[1]) else 1000L

rule <- function(x) cat("\n", strrep("=", 78), "\n", x, "\n", strrep("=", 78), "\n", sep = "")

## ---- 1. Load engine + inputs ------------------------------------------------
for (f in c("sma_inputs.R", "sma_uncertainty.R", "sma_markov_model.R",
            "sma_dsa_psa.R", "sma_obm_durability.R", "sma_scenarios.R"))
  source(file.path(here, f))

## ---- 2. Base-case CEA + validation vs the Excel workbook --------------------
rule("1. Base-case four-arm CEA (per patient, discounted) + validation vs Excel")
cea <- run_cea(faithful_excel = TRUE)
print(within(cea, { cost <- round(cost); ly <- round(ly, 3); vfly <- round(vfly, 3)
  qaly <- round(qaly, 4); icer_qaly <- round(icer_qaly)
  icer_ly <- round(icer_ly); icer_vfly <- round(icer_vfly) }), row.names = FALSE)

target <- data.frame(
  arm = c("BSC","Nusinersen","OA","Risdiplam"),
  cost = c(1631085.77, 4676905.35, 4793785.89, 7487473.87),
  qaly = c(0.67162, 2.161, 7.3697, 6.4923),
  icer = c(NA, 2045037.87, 472176.83, 1006141.81))
m <- merge(cea, target, by = "arm", suffixes = c("", "_xl"))
worst <- max(abs(m$cost - m$cost_xl)/m$cost_xl,
             abs(m$qaly - m$qaly_xl)/m$qaly_xl, na.rm = TRUE)
cat(sprintf("\nValidation: worst relative error vs Excel = %.3g  -> %s\n",
            worst, if (worst < 0.002) "PASS (<0.2%)" else "CHECK"))

## ---- 3. Scenario analyses ---------------------------------------------------
rule("2. Scenario analyses (ICER $/QALY vs BSC)")
print(run_scenarios())

## The optimistic-nusinersen-mortality scenario lowers nusinersen's ICER vs BSC
## while making it dominated on the frontier; print the full table so the ICER
## row above is not read on its own.
cat("\n  Detail — optimistic nusinersen mortality (NS->PV/Death set to OA's):\n")
print(within(scenario_detail(opt_nus_mortality), {
  cost <- round(cost); ly <- round(ly, 2); vfly <- round(vfly, 2)
  qaly <- round(qaly, 4); icer_qaly <- round(icer_qaly) }), row.names = FALSE)

## ---- 4. Outcomes-based payment: durability tables --------------------------
rule("3. Outcomes-based payment for OA — ICER by durability ($/QALY vs BSC)")
tab_vfs <- obm_durability_table("VFS")
tab_os  <- obm_durability_table("OS")
print(within(tab_vfs, { up_icer <- round(up_icer); ob5_icer <- round(ob5_icer)
  ob10_icer <- round(ob10_icer); qaly <- round(qaly, 3) })[
  , c("durability","qaly","up_icer","ob5_icer","ob10_icer")], row.names = FALSE)

## ---- 5. DSA + PSA + figures -------------------------------------------------
rule("4. Sensitivity analyses and figures")
have_gg <- requireNamespace("ggplot2", quietly = TRUE) &&
           requireNamespace("scales",  quietly = TRUE)
if (have_gg) {
  source(file.path(here, "sma_figures_gg.R"))
  source(file.path(here, "sma_figures_structure_frontier.R"))
  source(file.path(here, "sma_figures_mono.R"))   # monochrome overrides
  source(file.path(here, "sma_dsa_supplement.R"))
  figdir <- file.path(here, "figs"); dir.create(figdir, showWarnings = FALSE)
  dsa <- run_dsa("OA", "BSC")
  set.seed(1234); psa <- run_psa(n_sim, faithful_excel = TRUE)
  ## main-body figures 1-4
  save_fig(fig_model_structure(),         "fig1_model_structure", 7.6, 4.2, figdir)
  save_fig(fig_ceac_manuscript(),         "fig2_ceac_manuscript", 9.2, 4.4, figdir)
  save_fig(fig_frontier(cea),             "fig2_ce_frontier",     7.2, 5.0, figdir)  # supplement
  ## supplement / sensitivity figures
  save_fig(fig_tornado(dsa, "OA"),        "tornado_OA",          8, 5.2, figdir)
  dsa_panels <- dsa_all()
  save_fig(fig_tornado_panels(dsa_panels), "figS1_tornado_panels", 8.2, 9.4, figdir)
  dsa_table(dsa_panels, file = file.path(here, "dsa_results.csv"))
  save_fig(fig_ce_plane(psa, wtp = 5e5),  "ce_plane",            7, 6,   figdir)
  save_fig(fig_ceac(psa),                 "ceac",                7, 5,   figdir)
  save_fig(fig_obm_icer(tab_vfs, tab_os), "obm_icer_durability", 8, 5.2, figdir)
  save_fig(fig_obm_saving(tab_vfs),       "obm_saving",          8, 5.2, figdir)
  cat("Figures written to", figdir, "(.pdf + .png each)\n",
      "  main body : fig1_model_structure, fig2_ceac_manuscript,",
      "obm_icer_durability (Fig 3)\n",
      "  supplement: figS1_tornado_panels, fig2_ce_frontier, obm_saving, ce_plane, ceac\n")
} else {
  cat("ggplot2/scales not installed - numeric results above are complete;",
      "figures skipped.\n  Install with: install.packages(c('ggplot2','scales'))\n")
}

rule(sprintf("Done in %.1f s. See README.md for interpretation and provenance.",
             (proc.time() - t_start)["elapsed"]))

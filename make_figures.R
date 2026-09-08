## =============================================================================
## make_figures.R  —  produce all manuscript figures (ggplot2) + print the
## sensitivity/OBM tables. Usage:  Rscript make_figures.R [n_sim]  (default 1000)
## Outputs (figs/): tornado_OA, ce_plane, ceac, obm_icer_durability, obm_saving
##                  (each as .pdf and .png)
## =============================================================================
here <- tryCatch(dirname(sub("^--file=", "",
          grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), error = function(e) ".")
if (is.na(here) || here == "") here <- "."
for (f in c("sma_inputs.R","sma_uncertainty.R","sma_markov_model.R",
            "sma_dsa_psa.R","sma_obm_durability.R","sma_figures_gg.R"))
  source(file.path(here, f))

args  <- commandArgs(trailingOnly = TRUE)
n_sim <- if (length(args) >= 1) as.integer(args[1]) else 1000L
figdir <- file.path(here, "figs"); dir.create(figdir, showWarnings = FALSE)

## ---- DSA / tornado ----------------------------------------------------------
cat("DSA (OA vs BSC)...\n")
dsa_oa <- run_dsa("OA", "BSC")
save_fig(fig_tornado(dsa_oa, "OA"), "tornado_OA", w = 8, h = 5.2, dir = figdir)

## ---- PSA / CE plane + CEAC --------------------------------------------------
cat("PSA (", n_sim, " sims)...\n", sep = "")
set.seed(1234); psa <- run_psa(n_sim, faithful_excel = TRUE)
save_fig(fig_ce_plane(psa, wtp = 5e5), "ce_plane", w = 7, h = 6, dir = figdir)
save_fig(fig_ceac(psa),                "ceac",     w = 7, h = 5, dir = figdir)

## ---- OBM durability ---------------------------------------------------------
cat("OBM durability tables...\n")
tab_vfs <- obm_durability_table("VFS")
tab_os  <- obm_durability_table("OS")
save_fig(fig_obm_icer(tab_vfs, tab_os),   "obm_icer_durability", w = 8, h = 5.2, dir = figdir)
save_fig(fig_obm_saving(tab_vfs),         "obm_saving",          w = 8, h = 5.2, dir = figdir)

cat("\nOBM durability (VFS) — key rows (per patient, discounted):\n")
show <- within(tab_vfs, {
  up_icer <- round(up_icer); ob10_icer <- round(ob10_icer); ob5_icer <- round(ob5_icer)
  qaly <- round(qaly, 3)
})[, c("durability","qaly","up_icer","ob5_icer","ob10_icer")]
print(show, row.names = FALSE)

cat("\nFigures written to", figdir, "\n")
cat("Files: tornado_OA, ce_plane, ceac, obm_icer_durability, obm_saving (.pdf/.png)\n")

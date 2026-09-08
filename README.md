# Cost-effectiveness of disease-modifying therapies for infantile-onset SMA

R code for the analysis reported in:

> Abuloha S, Goodin A, Jiao T, Corti M, Svensson M. *Cost-Effectiveness of
> Disease-Modifying Therapies for Infantile-Onset Spinal Muscular Atrophy in the
> United States under an Outcomes-Based Payment Design.* (manuscript under review)

The model compares onasemnogene abeparvovec (OA), nusinersen, risdiplam, and best
supportive care (BSC) for infantile-onset (type 1) spinal muscular atrophy from a
US third-party payer perspective over a lifetime horizon, and evaluates 5- and
10-year outcomes-based annuity contracts for OA across treatment-durability
assumptions.

## How to run
```bash
Rscript RUN.R          # base-case CEA, scenarios, outcomes-based payment, DSA, PSA (1,000 draws), figures
Rscript RUN.R 200      # same, with 200 PSA draws for a quick check
```
`RUN.R` is the single entry point and runs in a few seconds. The numeric results
need only base R (≥ 4.0). The figures additionally need `ggplot2` and `scales`
(`install.packages(c("ggplot2", "scales"))`); if they are absent, the numbers are
still produced and the figures are skipped.

## Files
| File | Purpose |
|---|---|
| `RUN.R` | Master script: runs everything below and writes figures to `figs/`. |
| `sma_inputs.R` | All base-case inputs: transition probabilities, time-dependent per-cycle death rates, per-patient drug-cost schedules, health-state costs, utilities, discount vector, cycle timeline. |
| `sma_uncertainty.R` | Uncertainty inputs for the sensitivity analyses (Beta parameters for utilities, cost coefficients of variation, transition uncertainty). |
| `sma_markov_model.R` | Model engine: cohort trace, economic evaluation, four-arm CEA table, and the outcomes-based payment machinery for OA. |
| `sma_obm_durability.R` | Outcomes-based payment analysis: OA effect sustained for *d* years then reverting to natural history, with 5-/10-year annuities linked to ventilation-free survival (VFS) or overall survival (OS). |
| `sma_scenarios.R` | Deterministic scenario analyses. |
| `sma_dsa_psa.R` | Deterministic (tornado) and probabilistic sensitivity analyses. |
| `sma_dsa_supplement.R` | Three-panel tornado figure and the underlying table (`dsa_results.csv`) for the supplement. |
| `sma_figures_gg.R`, `sma_figures_structure_frontier.R`, `sma_figures_mono.R` | Figure functions (ggplot2). `sma_figures_mono.R` provides the greyscale versions used in the manuscript. |
| `run_and_validate.R`, `make_figures.R` | Stand-alone drivers for the base-case CEA and the figures, respectively (both are also called by `RUN.R`). |
| `psa_excel_draws.csv` | The 1,000 probabilistic draws (cost and QALYs per arm) used for the manuscript's acceptability curves and cost-effectiveness plane. |
| `dsa_results.csv` | One-way sensitivity analysis results (all three DMTs versus BSC). |
| `figs/` | Generated figures (PDF and PNG). |

## Model structure
Five health states: non-sitting (entry state, equivalent to SMA type 1), sitting
(type 2), walking (type 3), permanent ventilation, and death. Cycle length is
three months in year 1 and annual thereafter, with a lifetime horizon (79 years)
and 3% discounting of costs and outcomes. The treatment-specific transitions
(non-sitting → sitting / permanent ventilation / death, and sitting → walking)
come from the published indirect treatment comparison (Ribero et al. 2022) and
STR1VE-US for OA; transitions to death from the sitting, walking, and
permanent-ventilation states are state-specific, time-dependent, and identical
across arms (parametric survival models fitted to digitised Kaplan–Meier data
from Zerres et al. 1997 and Gregoretti et al. 2013, and US life tables for the
walking state). Health-state costs are from the HealthCare Integrated Research
Database (2024 USD) and utilities from Hu et al. 2022. Full sources and
distributions are given in Table 1 of the manuscript; see the header comments in
`sma_markov_model.R` for the transition equations.

## Sensitivity and scenario analyses
- **Deterministic (one-way)**: health-state and drug costs varied by ±20% (no
  usable variance in the source), utilities and treatment-specific transitions by
  ±1.96 standard errors. Output: `figs/figS1_tornado_panels`, `dsa_results.csv`.
- **Probabilistic**: the manuscript's acceptability curves (Figure 2) and
  cost-effectiveness plane (Supplementary Figure S2) are drawn from the 1,000
  draws in `psa_excel_draws.csv` (`fig_ceac_manuscript()`,
  `fig_ce_plane(read_excel_psa())`). `run_psa()` regenerates a PSA from the
  distributions in `sma_uncertainty.R`: utilities ~ Beta, costs ~ Gamma,
  non-sitting-row transition probabilities ~ Dirichlet, sitting → walking (OA) ~
  Beta, with transition uncertainty governed by an assumed effective sample size
  (`N_EFF`, default 80).
- **Scenarios** (`run_scenarios()`): base case, alternative utilities
  (McMillan et al. 2021), a higher sitting-state utility, lower discount rates
  (1.5% and 2%), and optimistic nusinersen mortality (nusinersen assigned OA's
  survival and ventilation-free survival, leaving its milestone transitions
  unchanged). `scenario_detail()` returns the full four-arm table for any
  override. The remaining manuscript scenarios (life-cycle drug pricing,
  alternative health-state costs, backward transitions, no new sitting after
  trial cut-off, optimistic walking) are not included in this package; the
  corresponding functions document the input each would require.
- **Outcomes-based payment** (`obm_durability_table("VFS")`, `("OS")`): expected
  per-patient OA drug outlay and ICER versus BSC under upfront payment and under
  5- and 10-year annuities whose installments are paid only while the patient
  meets the contract endpoint at each annual assessment, for durability of 1–9
  years and a fully sustained effect. Output: `figs/obm_icer_durability`
  (Figure 3) and `figs/obm_saving`.

## Figures produced by `RUN.R`
| File in `figs/` | Manuscript |
|---|---|
| `fig1_model_structure` | Figure 1 |
| `fig2_ceac_manuscript` | Figure 2 |
| `obm_icer_durability` | Figure 3 |
| `figS1_tornado_panels` | Supplementary Figure S1 |
| `ce_plane` | Supplementary Figure S2 |
| `fig2_ce_frontier`, `tornado_OA`, `ceac`, `obm_saving` | Additional outputs, not shown in the manuscript |

## Citation and licence
Licensed under the MIT License (see `LICENSE`). If you use this code, please
cite the manuscript above and the archived release of this repository
(Zenodo DOI to be added on release).

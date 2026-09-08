# SMA disease-modifying therapy CEA — reproducible R model

An independent R re-implementation of the five-state Markov cohort model used to
evaluate onasemnogene abeparvovec (OA), nusinersen, risdiplam, and best
supportive care (BSC) for infantile-onset (type 1) spinal muscular atrophy, from
a US payer perspective. It reproduces the original Excel workbook
(`SMA Excel Model.xlsm`) and is intended to be shared alongside publication so
that others can audit and rerun the analysis.

## Files
| File | Purpose |
|---|---|
| `sma_inputs.R` | All model inputs, auto-generated from the validated Excel workbook (transition probabilities, per-cycle time-dependent death rates, per-patient drug-cost schedules, state costs, utilities, discount vector, cycle timeline). Regenerate; do not hand-edit. |
| `sma_markov_model.R` | The model engine: cohort trace, economic evaluation, CEA table, and the outcomes-based (OB) payment machinery for OA. Base R only — no packages required. |
| `sma_uncertainty.R` | Auto-generated uncertainty inputs for DSA/PSA (Beta α/β for utilities, cost CVs, transition uncertainty). |
| `sma_dsa_psa.R` | Deterministic (tornado) and probabilistic (CE plane, CEAC) sensitivity analyses + figure functions. |
| `run_and_validate.R` | Driver: runs the four-arm CEA and checks it against the workbook. |
| `make_figures.R` | Driver: runs DSA + PSA and writes the manuscript figures to `figs/`. |

## How to run
```bash
cd "R model"
Rscript RUN.R                     # everything: validation, CEA, scenarios,
                                  # OBM durability, DSA, PSA(1000), all figures
Rscript RUN.R 200                 # faster (200 PSA draws) for a quick check
```
`RUN.R` is the single reproducibility entry point (~3 s). Individual drivers also
exist: `run_and_validate.R` (CEA + validation) and `make_figures.R` (figures).
The CEA/DSA/OBM/scenario numbers need only base R; the figures additionally need
`ggplot2` and `scales` (install once: `install.packages(c("ggplot2","scales"))`).
Licensed under MIT (see `LICENSE`).

## Sensitivity analyses and figures
`make_figures.R` produces three publication figures in `figs/`:
- **`tornado_OA.pdf`** — one-way DSA on the OA-vs-BSC ICER. Top drivers reproduce
  the dissertation's reported sensitivities: sitting-state utility, the
  non-sitting→sitting transition, and OA drug cost.
- **`ce_plane.pdf`** — PSA cost-effectiveness plane (1000 draws), three arms vs BSC.
- **`ceac.pdf`** — cost-effectiveness acceptability curves.

**PSA distributions** (standard practice; provenance in `sma_uncertainty.R`):
utilities ~ Beta (workbook α/β); state and drug costs ~ Gamma (workbook CVs);
NS-row transition probabilities ~ Dirichlet, S→W (OA) ~ Beta. Transition
uncertainty uses an assumed effective sample size `N_EFF` (default 80) because
the workbook's stored transition SEs are on a transformed scale and not usable
directly; `N_EFF = 80` is calibrated so the PSA reproduces the manuscript's CEAC
(P(OA cost-effective) ≈ 55% at $500k and ≈ 95% at $750k/QALY) while keeping the
**PSA mean aligned with the deterministic base case** (the standard PSA sanity
check — verified automatically by `make_figures.R`). This is a documented
modelling choice; state it in the methods and adjust `N_EFF` if you obtain the
underlying evidence sample sizes.

## Validation
The model reproduces the workbook's `CEA Results` sheet (discounted, per patient)
to within **0.02%** on every reported quantity:

| Output | Agreement |
|---|---|
| QALY, LY, VFLY (all arms) | exact to ~1e-5 relative |
| Total cost (all arms) | ≤ 0.02% relative |
| ICER $/QALY vs BSC | ≤ 0.02% relative |

Reproduced base-case results (per patient, discounted, lifetime): OA is the most
efficient DMT — dominating risdiplam, cost-effective vs nusinersen ($22.4k/QALY),
and $472k/QALY vs BSC; none is cost-effective vs BSC at $50–200k/QALY.

The tiny (<0.02%) cost residual comes from re-deriving the per-patient drug
schedule by division from the workbook traces; it does not affect any conclusion.

## Model structure
Five states: **NS** non-sitting (entry, ~type 1), **S** sitting (~type 2),
**W** walking (~type 3), **PV** permanent ventilation, **Death**. Cycle length is
3 months in year 1 and annual thereafter (lifetime horizon, 79 years). Forward
motor transitions and NS→PV/Death come from a fixed base-case TP set (3-month for
year 1, 12-month after); S→Death, W→Death, and PV→Death are time-dependent
per-cycle rates from parametric survival fits. See the code header comments for
the exact transition equations.

## Provenance of inputs (traceability for reviewers)
`sma_inputs.R` is generated directly from the workbook:
- **Fixed transition probabilities** — `TP-3months cycle length` / `TP-12months
  cycle length`, "Live Value" column (base case). Sources: Ribero et al. 2022
  indirect comparison; STR1VE-US for OA.
- **Time-dependent death rates** (`s_death`, `w_death`, `pv_death`) — outputs of
  the parametric survival models on the `Modeling Hazard Function-*` sheets
  (SMA type 2 survival ← Zerres 1997; PV ← Gregoretti 2013; walking ← US life
  table). These are state-specific (identical across arms) and embedded as
  validated per-cycle vectors; regenerable from those sheets.
- **Per-patient drug schedules** (`drug_pp`) — derived as (per-cycle drug cost) /
  (living non-PV cohort) from the `Simulation-*` sheets; embed the dosing logic
  (nusinersen loading + maintenance, risdiplam weight-based, OA one-time).
- **State costs / utilities** — `Costs` (HIRD; 2024 USD) and `Utility` (Hu 2022).

## Known Excel quirks reproduced (set `faithful_excel = TRUE`, the default)
While reverse-engineering the workbook we identified two internal inconsistencies
that the R model reproduces **exactly** when `faithful_excel = TRUE`, and removes
when `faithful_excel = FALSE`:

1. **Split discounting.** Costs are discounted with the *same-cycle* factor,
   while health outcomes (LY/VFLY/QALY) use the *prior-cycle* factor (a
   discount-at-cycle-start lag). This mixes two discounting conventions.
2. **Cycle-1 mismatch.** At cycle 1, costs use a 12-month multiplier but
   life-years use a 0.25 (quarter) weight — cost and effect are measured over
   different durations in that one cycle.

Neither materially changes the conclusions (they move ICERs by well under 1%),
but for a *published* model we recommend running `faithful_excel = FALSE`, which
applies one consistent convention throughout, and reporting the Excel-faithful
numbers only for the reproducibility check. **This is worth a sentence in the
methods/limitations and a decision before submission.**

## Outcomes-based (OB) payment extension for OA
`effective_oa_drug_cost(term_years, endpoint)` computes the expected per-patient
discounted OA drug outlay under a 5- or 10-year annuity whose installments are
paid only while the patient meets the payment condition (VFS = not PV/Death; OS =
not Death) at each annual assessment. To obtain OB-contract ICERs, substitute
this cost for OA's upfront drug cost in `evaluate_arm()`.

**Durability-cliff variants** (as in the manuscript's durability table): model
treatment effect as sustained for *d* years, then set OA's transitions to the
untreated (BSC) natural history from year *d* onward, and recompute both the trace
and the OB drug stream. This is a documented next step, not yet wired into the
default driver; the hooks (`run_trace`, `effective_oa_drug_cost`) are in place.

## Outcomes-based durability analysis (`sma_obm_durability.R`)
Reproduces the workbook's `OBM-Durability` analysis: OA effect sustained for *d*
years then reverting to BSC natural history, with 5-/10-year annuities linked to
VFS or OS. `obm_durability_table("VFS")` and `("OS")` return the full tables;
`make_figures.R` plots `obm_icer_durability` and `obm_saving`. Validated against
the workbook to <0.05% (e.g., 10-yr VFS at 1-yr durability: ICER $479,686 vs
$479,908; drug cost $1,193,212 matching the dissertation).

### ⚠ Discrepancy found in manuscript Table 2 (worth fixing before submission)
The reproduction surfaced an internal inconsistency in the **manuscript** (not the
dissertation, which is correct). Manuscript Table 2 lists, for the 10-year
VFS-linked contract at 1-year durability, a drug cost of **$1,527,013** paired
with an ICER of **$479,908**. These are mutually inconsistent: an ICER of
$479,908 requires a drug cost of **$1,193,212** (the VFS value; the dissertation
states exactly this — "reduced by half ($1,193,212)"). $1,527,013 is the
**OS-linked** drug cost. The VFS drug-cost column in Table 2 appears to have been
populated from the OS block. The R model outputs the correct VFS value.

## Scenario analyses (`sma_scenarios.R`)
`run_scenarios()` returns an ICER (vs BSC) table across scenarios. Implemented and
validated against the dissertation: base case; alternative utilities (McMillan);
higher sitting utility (0.65 → OA $285,540 vs dissertation $285,508); lower
discount rate (1.5% / 2%). Alt-utility (McMillan) OA $402,436 matches the
dissertation's reported OA scenario-range minimum ($402,391). Scenarios needing
inputs not carried in this package (life-cycle drug pricing, MarketScan/Belter
costs, backward transitions, no-new-sitting, optimistic walking, sequential DMT)
are stubbed with the exact input each requires — nothing is fabricated.

**Optimistic nusinersen mortality.** Sets nusinersen's `NS_PV` and `NS_Death` to
OA's values in both the year-1 and annual transition sets, leaving the milestone
transitions (`NS_S`, `S_W`) at nusinersen's own — i.e. nusinersen is granted OA's
survival and ventilation-free survival but not its motor benefit. This tests the
base case's dependence on the Ribero 2022 indirect comparison against the
manufacturer-funded critique of it (Jiang et al., *Adv Ther* 2023;40:2985–3005),
which concluded there is no statistically significant OS/EFS difference between
nusinersen and OA. It needs no inputs beyond `sma_inputs.R` — OA's transitions are
already there — and reproduces the dissertation: nusinersen cost $10,005,093 (vs
$10,003,932 reported) and ICER $1,773,735/QALY (vs $1,773,489), both within the
port's known Excel drift (see *Validation*). Read it with `scenario_detail()`, not
the ICER row alone: nusinersen's ICER vs BSC *falls* while its cost roughly
doubles, leaving it strictly dominated by both OA and risdiplam.

`scenario_detail(over)` returns the full four-arm cost/LY/VFLY/QALY table under any
override, with a `str_dominated` flag. That flag is **strict** dominance only
(cheaper *and* more QALYs); it does not detect extended dominance, so base-case
nusinersen — dominated by extension — is correctly flagged `FALSE`.

## Figures (`sma_figures_gg.R`, ggplot2)
Journal-styled, returned as ggplot objects and written by `make_figures.R` to
`figs/` as both PDF and PNG: `tornado_OA`, `ce_plane`, `ceac`,
`obm_icer_durability`, `obm_saving`. Requires `ggplot2` + `scales`.

## Status / not yet ported
Complete and validated: the base-case four-arm CEA, the DSA (tornado), the PSA
(CE plane + CEAC), and the OA outcomes-based payment cost machinery.
Not ported: the durability-cliff OBM variants (hooks in place; see above) and the
EVPI sheets (being dropped from the paper per the methods discussion). The
deterministic scenario analyses (LCDP, alternative utilities/costs, lower
discount rate, etc.) can be reproduced by editing the corresponding inputs.

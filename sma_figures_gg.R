## =============================================================================
## sma_figures_gg.R  —  journal-ready ggplot2 figures for the SMA CEA.
## Requires: ggplot2, scales. Source the model + DSA/PSA + OBM files first.
## Each function returns a ggplot object; save_fig() writes PDF + PNG.
## =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

## ---- shared theme + palette -------------------------------------------------
arm_cols <- c(BSC = "#6C757D", Nusinersen = "#E07B39",
              OA = "#2E75B6", Risdiplam = "#4CAF50")

theme_sma <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = "grey92"),
      panel.border     = element_rect(colour = "grey70"),
      plot.title       = element_text(face = "bold", size = rel(1.05)),
      plot.subtitle    = element_text(colour = "grey35", size = rel(0.9)),
      axis.title       = element_text(colour = "grey20"),
      legend.key       = element_blank(),
      legend.background= element_blank(),
      strip.background = element_rect(fill = "grey95", colour = "grey70"),
      plot.margin      = margin(10, 14, 8, 8)
    )
}

save_fig <- function(p, file, w = 7, h = 5, dir = "figs") {
  dir.create(dir, showWarnings = FALSE)
  ## base pdf device (portable; cairo is unavailable on some installs)
  suppressWarnings(ggsave(file.path(dir, paste0(file, ".pdf")), p,
                          width = w, height = h, device = pdf))
  suppressWarnings(ggsave(file.path(dir, paste0(file, ".png")), p,
                          width = w, height = h, dpi = 200))
  invisible(p)
}

dollar_k <- label_number(scale_cut = cut_short_scale(), prefix = "$")

## ---- 1. Tornado (DSA) -------------------------------------------------------
fig_tornado <- function(dsa, target = "OA", ref = "BSC") {
  base <- attr(dsa, "base_icer")
  d <- dsa
  d$param <- factor(d$param, levels = d$param[order(d$range)])
  d$lo <- pmin(d$low, d$high); d$hi <- pmax(d$low, d$high)
  ggplot(d, aes(y = param)) +
    geom_vline(xintercept = base, linetype = 2, colour = "grey45") +
    geom_segment(aes(x = lo, xend = hi, yend = param),
                 linewidth = 7, colour = arm_cols[[target]], alpha = 0.85,
                 lineend = "butt") +
    scale_x_continuous(labels = dollar_k) +
    labs(title = paste0("One-way sensitivity analysis: ", target, " vs ", ref),
         subtitle = paste0("Base-case ICER = ", dollar(round(base)), "/QALY"),
         x = "ICER ($/QALY gained)", y = NULL) +
    theme_sma()
}

## ---- 2. Cost-effectiveness plane (PSA) --------------------------------------
fig_ce_plane <- function(psa, arms = c("Nusinersen","OA","Risdiplam"), wtp = 5e5) {
  df <- do.call(rbind, lapply(arms, function(a)
    data.frame(arm = a,
               dqaly = psa[[paste0("q_",a)]] - psa$q_BSC,
               dcost = psa[[paste0("c_",a)]] - psa$c_BSC)))
  df$arm <- factor(df$arm, levels = arms)
  ggplot(df, aes(dqaly, dcost, colour = arm)) +
    geom_hline(yintercept = 0, colour = "grey80") +
    geom_abline(slope = wtp, intercept = 0, linetype = 2, colour = "grey45") +
    geom_point(alpha = 0.28, size = 0.7) +
    annotate("text", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5, size = 3,
             colour = "grey45", label = paste0("WTP = ", dollar(wtp), "/QALY")) +
    scale_colour_manual(values = arm_cols, name = NULL) +
    scale_y_continuous(labels = dollar_k) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2))) +
    labs(title = "Cost-effectiveness plane",
         subtitle = paste0("Probabilistic sensitivity analysis (", nrow(psa), " simulations), vs BSC"),
         x = "Incremental QALYs vs BSC", y = "Incremental cost vs BSC") +
    theme_sma() + theme(legend.position = c(0.14, 0.85))
}

## ---- 3. CEAC ----------------------------------------------------------------
fig_ceac <- function(psa, arms = c("BSC","Nusinersen","OA","Risdiplam"),
                     wtp_grid = seq(0, 1e6, 20000)) {
  nmb <- function(a, l) l*psa[[paste0("q_",a)]] - psa[[paste0("c_",a)]]
  P <- sapply(wtp_grid, function(l) {
    M <- sapply(arms, nmb, l = l)
    prop.table(table(factor(arms[max.col(M, "first")], levels = arms)))
  })
  df <- do.call(rbind, lapply(seq_along(arms), function(i)
    data.frame(arm = arms[i], wtp = wtp_grid, p = P[i, ])))
  df$arm <- factor(df$arm, levels = arms)
  ggplot(df, aes(wtp, p, colour = arm)) +
    geom_line(linewidth = 0.9) +
    scale_colour_manual(values = arm_cols, name = NULL) +
    scale_x_continuous(labels = dollar_k) +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(title = "Cost-effectiveness acceptability curves",
         subtitle = "Probability each strategy is optimal, by willingness-to-pay",
         x = "Willingness-to-pay ($/QALY)", y = "Probability cost-effective") +
    theme_sma() + theme(legend.position = c(0.85, 0.55))
}

## ---- 4. OBM: ICER vs durability (standard vs OB contracts) ------------------
fig_obm_icer <- function(tab_vfs, tab_os) {
  d <- rbind(
    data.frame(dur = tab_vfs$durability, icer = tab_vfs$up_icer,  scheme = "Standard upfront"),
    data.frame(dur = tab_vfs$durability, icer = tab_vfs$ob10_icer, scheme = "10-year OB (VFS)"),
    data.frame(dur = tab_vfs$durability, icer = tab_vfs$ob5_icer,  scheme = "5-year OB (VFS)"),
    data.frame(dur = tab_os$durability,  icer = tab_os$ob10_icer,  scheme = "10-year OB (OS)"))
  ord <- c(as.character(1:9), "Full")
  d$dur <- factor(d$dur, levels = ord)
  d$scheme <- factor(d$scheme, levels = c("Standard upfront","10-year OB (OS)","5-year OB (VFS)","10-year OB (VFS)"))
  cols <- c("Standard upfront"="grey35","10-year OB (OS)"="#E07B39",
            "5-year OB (VFS)"="#8E7CC3","10-year OB (VFS)"="#2E75B6")
  ggplot(d, aes(dur, icer, colour = scheme, group = scheme)) +
    geom_line(linewidth = 0.9) + geom_point(size = 1.8) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_y_continuous(labels = dollar_k) +
    labs(title = "Effect of outcomes-based payment on the OA vs BSC ICER",
         subtitle = "By assumed treatment durability (years of sustained effect before cliff)",
         x = "Treatment durability (years)", y = "ICER ($/QALY gained)") +
    theme_sma() + theme(legend.position = c(0.75, 0.8))
}

## ---- 5. OBM: financial risk mitigation (drug-cost saving) vs durability -----
fig_obm_saving <- function(tab_vfs, upfront_drug = 2220425.89) {
  ## OB drug cost = total - disease; disease = upfront_total - upfront_drug
  disease <- tab_vfs$up_total - upfront_drug
  d <- rbind(
    data.frame(dur = tab_vfs$durability, saving = upfront_drug - (tab_vfs$ob10_total - disease), scheme = "10-year OB (VFS)"),
    data.frame(dur = tab_vfs$durability, saving = upfront_drug - (tab_vfs$ob5_total  - disease), scheme = "5-year OB (VFS)"))
  ord <- c(as.character(1:9), "Full"); d$dur <- factor(d$dur, levels = ord)
  ggplot(d, aes(dur, saving, fill = scheme)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65) +
    scale_fill_manual(values = c("5-year OB (VFS)"="#8E7CC3","10-year OB (VFS)"="#2E75B6"), name = NULL) +
    scale_y_continuous(labels = dollar_k) +
    labs(title = "Financial risk mitigation from outcomes-based payment (OA)",
         subtitle = "Expected reduction in per-patient drug outlay vs standard upfront, by durability",
         x = "Treatment durability (years)", y = "Drug-cost saving vs upfront") +
    theme_sma() + theme(legend.position = c(0.8, 0.85))
}

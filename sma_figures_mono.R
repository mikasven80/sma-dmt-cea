## =============================================================================
## sma_figures_mono.R  —  monochrome (greyscale) overrides for all figures.
##
## Source AFTER sma_figures_gg.R and sma_figures_structure_frontier.R. Every
## function defined here shadows the colour version of the same name, so the
## rest of the pipeline (RUN.R, make_figures.R) needs no change.
##
## Design rule: no series is distinguished by shade alone. Line series carry a
## distinct LINETYPE, point series a distinct SHAPE, and bars a distinct FILL
## plus border. The figures therefore survive greyscale printing and are legible
## to readers with colour-vision deficiency.
## =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

## ---- shared monochrome scales ----------------------------------------------
arm_cols <- c(BSC = "grey62", Nusinersen = "grey40",
              OA = "black",   Risdiplam = "grey20")
arm_ltys <- c(BSC = "dotted", Nusinersen = "dashed",
              OA = "solid",   Risdiplam = "dotdash")
arm_shps <- c(BSC = 15, Nusinersen = 17, OA = 16, Risdiplam = 4)

scheme_cols <- c("Standard upfront" = "grey55", "10-year OB (OS)" = "grey30",
                 "5-year OB (VFS)"  = "grey45", "10-year OB (VFS)" = "black")
scheme_ltys <- c("Standard upfront" = "solid",  "10-year OB (OS)" = "dotdash",
                 "5-year OB (VFS)"  = "dashed", "10-year OB (VFS)" = "solid")
scheme_shps <- c("Standard upfront" = 1, "10-year OB (OS)" = 17,
                 "5-year OB (VFS)"  = 15, "10-year OB (VFS)" = 16)

## ---- Figure 1: model structure ---------------------------------------------
fig_model_structure <- function() {
  bw <- 1.5; bh <- 0.62
  nodes <- data.frame(
    key   = c("NS", "S", "W", "PV", "D"),
    label = c("Non-sitting", "Sitting", "Walking",
              "Permanent\nventilation (PV)", "Death"),
    x     = c(0, 4.4, 8.8, 0, 8.8),
    y     = c(3, 3,   3,   0,  0),
    fill  = c("grey97", "grey97", "grey97", "grey88", "grey78"),
    stringsAsFactors = FALSE)
  MG <- "Motor-milestone gain"; PVk <- "Progression to ventilation"; DK <- "Death"

  edges <- data.frame(rbind(
    c(     bw,            3,        4.4 - bw,        3     ),
    c(4.4 + bw,           3,        8.8 - bw,        3     ),
    c(      0,       3 - bh,               0,     0 + bh   ),
    c(    8.8,       3 - bh,             8.8,     0 + bh   ),
    c(     bw,            0,        8.8 - bw,        0     ),
    c(bw*0.55,       3 - bh,        8.8 - bw,   0 + bh*0.55),
    c(4.4 + bw*0.55, 3 - bh,   8.8 - bw*0.55,   0 + bh     )), stringsAsFactors = FALSE)
  names(edges) <- c("x", "y", "xend", "yend")
  edges$kind <- factor(c(MG, MG, PVk, DK, DK, DK, DK), levels = c(MG, PVk, DK))

  mkloop <- function(k, up) {
    n <- nodes[nodes$key == k, ]
    data.frame(x = n$x - bw*0.30, y = n$y + if (up) bh else -bh,
               xend = n$x + bw*0.30, yend = n$y + if (up) bh else -bh,
               curv = if (up) -1.5 else 1.5)
  }
  loops <- rbind(mkloop("NS", TRUE), mkloop("S", TRUE), mkloop("W", TRUE),
                 mkloop("PV", FALSE), mkloop("D", FALSE))
  ar <- arrow(length = unit(0.15, "cm"), type = "closed")

  ggplot() +
    geom_curve(data = loops[loops$curv < 0, ], aes(x, y, xend = xend, yend = yend),
               curvature = -1.5, linewidth = 0.3, colour = "grey55", arrow = ar) +
    geom_curve(data = loops[loops$curv > 0, ], aes(x, y, xend = xend, yend = yend),
               curvature = 1.5, linewidth = 0.3, colour = "grey55", arrow = ar) +
    geom_segment(data = edges,
                 aes(x, y, xend = xend, yend = yend, linetype = kind),
                 linewidth = 0.55, colour = "black", arrow = ar) +
    geom_rect(data = nodes,
              aes(xmin = x - bw, xmax = x + bw, ymin = y - bh, ymax = y + bh),
              fill = nodes$fill, colour = "black", linewidth = 0.55) +
    geom_text(data = nodes, aes(x, y, label = label),
              size = 3.4, fontface = "bold", colour = "black", lineheight = 0.95) +
    annotate("text", x = -bw - 0.22, y = 3, label = "entry\nstate", hjust = 1,
             size = 2.9, fontface = "italic", colour = "grey35", lineheight = 0.95) +
    scale_linetype_manual(values = setNames(c("solid", "dashed", "dotted"),
                                            c(MG, PVk, DK)), name = NULL) +
    coord_equal(xlim = c(-3.4, 10.7), ylim = c(-1.5, 4.5), clip = "off") +
    guides(linetype = guide_legend(override.aes = list(linewidth = 0.7))) +
    theme_void(base_size = 11) +
    theme(legend.position  = "bottom",
          legend.text      = element_text(size = 9, colour = "grey20"),
          legend.key.width = unit(1.1, "cm"),
          plot.margin      = margin(4, 6, 2, 6))
}

## ---- Figure 2: two-panel CEAC ----------------------------------------------
fig_ceac_manuscript <- function(psa = read_psa_draws(),
                                wtp_grid = seq(0, 1e6, 5000), ref_wtp = 2e5) {
  all_arms <- c("BSC", "Nusinersen", "OA", "Risdiplam")
  dmts     <- c("Nusinersen", "OA", "Risdiplam")
  df <- rbind(
    ceac_df(psa, all_arms, "A. All strategies (vs best supportive care)", wtp_grid),
    ceac_df(psa, dmts,     "B. Disease-modifying therapies only",         wtp_grid))
  df$arm   <- factor(df$arm, levels = all_arms)
  df$panel <- factor(df$panel, levels = unique(df$panel))

  ggplot(df, aes(wtp, p, colour = arm, linetype = arm)) +
    geom_vline(xintercept = ref_wtp, linetype = "dotted",
               colour = "grey45", linewidth = 0.4) +
    geom_line(linewidth = 0.85) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = arm_cols, name = NULL, drop = FALSE) +
    scale_linetype_manual(values = arm_ltys, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = dollar_k, breaks = seq(0, 1e6, 250000),
                       expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_continuous(labels = percent, limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(x = "Cost-effectiveness threshold ($/QALY)",
         y = "Probability cost-effective") +
    theme_sma() +
    theme(legend.position = "bottom",
          panel.spacing   = unit(1.0, "lines"),
          strip.text      = element_text(face = "bold", size = rel(0.92), hjust = 0,
                                         margin = margin(4, 4, 4, 4)))
}

## ---- Figure 3: OBM ICER by durability --------------------------------------
fig_obm_icer <- function(tab_vfs, tab_os) {
  d <- rbind(
    data.frame(dur = tab_vfs$durability, icer = tab_vfs$up_icer,   scheme = "Standard upfront"),
    data.frame(dur = tab_vfs$durability, icer = tab_vfs$ob10_icer, scheme = "10-year OB (VFS)"),
    data.frame(dur = tab_vfs$durability, icer = tab_vfs$ob5_icer,  scheme = "5-year OB (VFS)"),
    data.frame(dur = tab_os$durability,  icer = tab_os$ob10_icer,  scheme = "10-year OB (OS)"))
  ord <- c(as.character(1:9), "Full")
  d$dur <- factor(d$dur, levels = ord)
  d$scheme <- factor(d$scheme, levels = c("Standard upfront", "10-year OB (OS)",
                                          "5-year OB (VFS)", "10-year OB (VFS)"))
  ggplot(d, aes(dur, icer, colour = scheme, linetype = scheme,
                shape = scheme, group = scheme)) +
    geom_line(linewidth = 0.85) + geom_point(size = 2.1, fill = "white") +
    scale_colour_manual(values = scheme_cols, name = NULL) +
    scale_linetype_manual(values = scheme_ltys, name = NULL) +
    scale_shape_manual(values = scheme_shps, name = NULL) +
    scale_y_continuous(labels = dollar_k) +
    labs(x = "Treatment durability (years)",
         y = "ICER ($/QALY gained)") +
    theme_sma() + theme(legend.position = c(0.75, 0.8))
}

## ---- supplement: OBM drug-cost saving --------------------------------------
fig_obm_saving <- function(tab_vfs, upfront_drug = 2220425.89) {
  disease <- tab_vfs$up_total - upfront_drug
  d <- rbind(
    data.frame(dur = tab_vfs$durability,
               saving = upfront_drug - (tab_vfs$ob10_total - disease), scheme = "10-year OB (VFS)"),
    data.frame(dur = tab_vfs$durability,
               saving = upfront_drug - (tab_vfs$ob5_total  - disease), scheme = "5-year OB (VFS)"))
  ord <- c(as.character(1:9), "Full"); d$dur <- factor(d$dur, levels = ord)
  ggplot(d, aes(dur, saving, fill = scheme)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.65,
             colour = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("5-year OB (VFS)" = "grey78",
                                 "10-year OB (VFS)" = "grey30"), name = NULL) +
    scale_y_continuous(labels = dollar_k) +
    labs(x = "Treatment durability (years)",
         y = "Drug-cost saving vs upfront") +
    theme_sma() + theme(legend.position = c(0.8, 0.85))
}

## ---- supplement: tornado ----------------------------------------------------
fig_tornado <- function(dsa, target = "OA", ref = "BSC") {
  base <- attr(dsa, "base_icer")
  d <- dsa
  d$param <- factor(d$param, levels = d$param[order(d$range)])
  d$lo <- pmin(d$low, d$high); d$hi <- pmax(d$low, d$high)
  ggplot(d, aes(y = param)) +
    geom_vline(xintercept = base, linetype = 2, colour = "grey45") +
    geom_segment(aes(x = lo, xend = hi, yend = param), linewidth = 7,
                 colour = "grey55", lineend = "butt") +
    scale_x_continuous(labels = dollar_k) +
    labs(x = "ICER ($/QALY gained)",
         y = NULL) +
    theme_sma()
}

## ---- supplement: CE plane ---------------------------------------------------
fig_ce_plane <- function(psa, arms = c("Nusinersen", "OA", "Risdiplam"), wtp = 5e5) {
  d <- do.call(rbind, lapply(arms, function(a) data.frame(
    arm = a,
    dq  = psa[[paste0("q_", a)]] - psa$q_BSC,
    dc  = psa[[paste0("c_", a)]] - psa$c_BSC)))
  d$arm <- factor(d$arm, levels = arms)
  ggplot(d, aes(dq, dc, colour = arm, shape = arm)) +
    geom_abline(slope = wtp, intercept = 0, linetype = 2, colour = "grey45") +
    geom_point(alpha = 0.55, size = 1.1, stroke = 0.5) +
    scale_colour_manual(values = arm_cols, name = NULL) +
    scale_shape_manual(values = arm_shps, name = NULL) +
    scale_y_continuous(labels = dollar_k) +
    annotate("text", x = Inf, y = Inf, hjust = 1.05, vjust = 1.6, size = 3,
             colour = "grey30", label = paste0("Threshold = ", dollar(wtp), "/QALY")) +
    labs(x = "Incremental QALYs vs BSC",
         y = "Incremental cost vs BSC") +
    theme_sma() + theme(legend.position = c(0.16, 0.86)) +
    guides(colour = guide_legend(override.aes = list(alpha = 1, size = 2.2)))
}

## ---- supplement: single-panel CEAC (regenerated PSA) ----------------------------
fig_ceac <- function(psa, arms = c("BSC", "Nusinersen", "OA", "Risdiplam"),
                     wtp_grid = seq(0, 1e6, 20000)) {
  df <- ceac_df(psa, arms, "x", wtp_grid)
  df$arm <- factor(df$arm, levels = arms)
  ggplot(df, aes(wtp, p, colour = arm, linetype = arm)) +
    geom_line(linewidth = 0.85) +
    scale_colour_manual(values = arm_cols, name = NULL) +
    scale_linetype_manual(values = arm_ltys, name = NULL) +
    scale_x_continuous(labels = dollar_k) +
    scale_y_continuous(labels = percent, limits = c(0, 1)) +
    labs(x = "Cost-effectiveness threshold ($/QALY)",
         y = "Probability cost-effective") +
    theme_sma() + theme(legend.position = c(0.85, 0.55))
}

## ---- supplement: CE frontier ------------------------------------------------
fig_frontier <- function(res) {
  fr <- ce_frontier(res)
  res$on_frontier <- res$arm %in% fr$arm
  res$status <- ifelse(res$on_frontier, "On frontier", NA)
  for (i in which(!res$on_frontier))
    res$status[i] <- if (any(res$cost < res$cost[i] & res$qaly > res$qaly[i]))
      "Dominated" else "Extendedly dominated"
  seg <- data.frame(x = head(fr$qaly, -1), y = head(fr$cost, -1),
                    xend = tail(fr$qaly, -1), yend = tail(fr$cost, -1))
  seg$icer <- (seg$yend - seg$y) / (seg$xend - seg$x)
  seg$mx <- (seg$x + seg$xend) / 2; seg$my <- (seg$y + seg$yend) / 2
  res$nx <- ifelse(res$arm %in% c("Risdiplam", "Nusinersen"), -0.35, 0.35)
  res$ny <- ifelse(res$arm == "OA", -0.42e6, 0.34e6)
  res$hj <- ifelse(res$nx < 0, 1, 0)

  ggplot(res, aes(qaly, cost)) +
    geom_segment(data = seg, aes(x, y, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = "black", linewidth = 0.7) +
    geom_text(data = seg, aes(mx, my, label = paste0(dollar_k(round(icer)), "/QALY")),
              inherit.aes = FALSE, vjust = 2.1, size = 3.1,
              colour = "black", fontface = "bold") +
    geom_point(aes(shape = status), size = 3.2, stroke = 1.1, colour = "black",
               fill = "white") +
    geom_text(aes(label = arm, x = qaly + nx, y = cost + ny, hjust = hj),
              size = 3.3, fontface = "bold", colour = "black") +
    scale_shape_manual(values = c("On frontier" = 16, "Dominated" = 4,
                                  "Extendedly dominated" = 21), name = NULL) +
    scale_y_continuous(labels = dollar_k, limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.10))) +
    scale_x_continuous(expand = expansion(mult = c(0.07, 0.13))) +
    labs(x = "Quality-adjusted life-years (QALY)",
         y = "Total cost") +
    theme_sma() +
    theme(legend.position = c(0.99, 0.03), legend.justification = c(1, 0),
          legend.text = element_text(size = 8.5))
}

## =============================================================================
## sma_figures_structure_frontier.R
##   fig_model_structure() — Figure 1, five-state Markov model structure
##   fig_ceac_panels()     — Figure 2, two-panel CEAC (vs BSC / excluding BSC)
##   fig_frontier()        — cost-effectiveness frontier (supplement)
##
## Replaces the two dissertation-era raster figures with ggplot2 versions in the
## same theme as the rest of the figure set (sma_figures_gg.R). Source AFTER
## sma_figures_gg.R (uses theme_sma(), arm_cols, save_fig(), dollar_k).
##
## Terminology note: the state is labelled "PV" (permanent ventilation) to match
## the manuscript. The dissertation figure used "PAV"; that term is not used in
## the manuscript and is deliberately not reproduced here.
## =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

## ---- Figure 1: model structure ---------------------------------------------
## Transition set mirrors Section 2.2 exactly:
##   forward milestones NS -> S -> W; NS -> PV; all live states -> Death;
##   NO direct S/W -> PV; NO backward transitions in the base case.

fig_model_structure <- function() {

  bw <- 1.5; bh <- 0.62                       # box half-width / half-height
  ## Layout: milestone states left-to-right along the top; PV bottom-left
  ## (reachable only from non-sitting); Death bottom-right (absorbing).
  nodes <- data.frame(
    key   = c("NS", "S", "W", "PV", "D"),
    label = c("Non-sitting", "Sitting", "Walking",
              "Permanent\nventilation (PV)", "Death"),
    x     = c(0, 4.4, 8.8, 0, 8.8),
    y     = c(3, 3,   3,   0,  0),
    fill  = c("#EAF2FA", "#EAF2FA", "#EAF2FA", "#FBEEE3", "#EFEFEF"),
    line  = c("#2E75B6", "#2E75B6", "#2E75B6", "#E07B39", "#6C757D"),
    stringsAsFactors = FALSE
  )
  MG <- "Motor-milestone gain"; PVk <- "Progression to ventilation"; DK <- "Death"

  ## Explicit endpoints, all trimmed to box borders. Chosen so that no arrow
  ## crosses another and none passes through a box.
  edges <- data.frame(rbind(
    #        x,              y,        xend,            yend
    c(     bw,             3,        4.4 - bw,        3     ),  # NS -> S
    c(4.4 + bw,            3,        8.8 - bw,        3     ),  # S  -> W
    c(      0,        3 - bh,               0,     0 + bh   ),  # NS -> PV
    c(    8.8,        3 - bh,             8.8,     0 + bh   ),  # W  -> Death
    c(     bw,             0,        8.8 - bw,        0     ),  # PV -> Death
    c(bw*0.55,        3 - bh,        8.8 - bw,   0 + bh*0.55),  # NS -> Death
    c(4.4 + bw*0.55,  3 - bh,   8.8 - bw*0.55,   0 + bh     )), # S  -> Death
    stringsAsFactors = FALSE)
  names(edges) <- c("x", "y", "xend", "yend")
  edges$kind <- factor(c(MG, MG, PVk, DK, DK, DK, DK),
                       levels = c(MG, PVk, DK))

  ## self-loops: above the top-row states, below PV and Death
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
    geom_curve(data = loops[loops$curv < 0, ],
               aes(x, y, xend = xend, yend = yend),
               curvature = -1.5, linewidth = 0.32, colour = "grey60", arrow = ar) +
    geom_curve(data = loops[loops$curv > 0, ],
               aes(x, y, xend = xend, yend = yend),
               curvature = 1.5, linewidth = 0.32, colour = "grey60", arrow = ar) +
    geom_segment(data = edges,
                 aes(x, y, xend = xend, yend = yend, colour = kind),
                 linewidth = 0.6, arrow = ar) +
    geom_rect(data = nodes,
              aes(xmin = x - bw, xmax = x + bw, ymin = y - bh, ymax = y + bh),
              fill = nodes$fill, colour = nodes$line, linewidth = 0.6) +
    geom_text(data = nodes, aes(x, y, label = label),
              size = 3.4, fontface = "bold", colour = "grey15", lineheight = 0.95) +
    annotate("text", x = -bw - 0.22, y = 3, label = "entry\nstate", hjust = 1,
             size = 2.9, fontface = "italic", colour = "grey45", lineheight = 0.95) +
    scale_colour_manual(values = setNames(c("#2E75B6", "#E07B39", "#6C757D"),
                                          c(MG, PVk, DK)), name = NULL) +
    coord_equal(xlim = c(-3.4, 10.7), ylim = c(-1.5, 4.5), clip = "off") +
    guides(colour = guide_legend(override.aes = list(linewidth = 1.1))) +
    theme_void(base_size = 11) +
    theme(legend.position   = "bottom",
          legend.text       = element_text(size = 9, colour = "grey25"),
          legend.key.width  = unit(1.0, "cm"),
          plot.margin       = margin(4, 6, 2, 6))
}

## ---- Figure 2: cost-effectiveness frontier ---------------------------------
## Frontier computed from the model output, not hard-coded: strictly dominated
## options are removed first, then extendedly dominated options iteratively.

ce_frontier <- function(res) {
  d <- res[order(res$cost), c("arm", "cost", "qaly")]
  keep <- rep(TRUE, nrow(d))
  for (i in seq_len(nrow(d)))                       # strict dominance
    if (any(d$cost < d$cost[i] & d$qaly > d$qaly[i])) keep[i] <- FALSE
  d <- d[keep, ]
  repeat {                                          # extended dominance
    if (nrow(d) < 3) break
    ## icer[k] = ICER of option k vs option k-1 (icer[1] is NA for the anchor).
    ## If icer[k] > icer[k+1], option k is extendedly dominated -> drop row k.
    icer <- c(NA, diff(d$cost) / diff(d$qaly))
    bad  <- which(diff(icer[-1]) < 0) + 1            # row index in d, not +1 again
    if (!length(bad)) break
    d <- d[-bad[1], ]
  }
  d
}

fig_frontier <- function(res) {
  fr  <- ce_frontier(res)
  res$on_frontier <- res$arm %in% fr$arm

  ## status label for each off-frontier option
  res$status <- ifelse(res$on_frontier, "On frontier", NA)
  for (i in which(!res$on_frontier))
    res$status[i] <- if (any(res$cost < res$cost[i] & res$qaly > res$qaly[i]))
      "Dominated" else "Extendedly dominated"

  seg <- data.frame(x = head(fr$qaly, -1), y = head(fr$cost, -1),
                    xend = tail(fr$qaly, -1), yend = tail(fr$cost, -1))
  seg$icer <- (seg$yend - seg$y) / (seg$xend - seg$x)
  seg$mx <- (seg$x + seg$xend) / 2; seg$my <- (seg$y + seg$yend) / 2

  ## nudge labels so they don't sit on the markers
  res$nx <- ifelse(res$arm %in% c("Risdiplam", "Nusinersen"), -0.35, 0.35)
  res$ny <- ifelse(res$arm == "OA", -0.42e6, 0.34e6)
  res$hj <- ifelse(res$nx < 0, 1, 0)

  ggplot(res, aes(qaly, cost)) +
    geom_segment(data = seg, aes(x, y, xend = xend, yend = yend),
                 inherit.aes = FALSE, colour = "#2E75B6", linewidth = 0.8) +
    geom_text(data = seg,
              aes(mx, my, label = paste0(dollar_k(round(icer)), "/QALY")),
              inherit.aes = FALSE, angle = 0, vjust = 2.1, size = 3.1,
              colour = "#2E75B6", fontface = "bold") +
    geom_point(aes(shape = status, colour = arm), size = 3.4, stroke = 1.1) +
    geom_text(aes(label = arm, x = qaly + nx, y = cost + ny, hjust = hj),
              size = 3.3, fontface = "bold", colour = "grey15") +
    scale_colour_manual(values = arm_cols, guide = "none") +
    scale_shape_manual(values = c("On frontier" = 16, "Dominated" = 4,
                                  "Extendedly dominated" = 1), name = NULL) +
    scale_y_continuous(labels = dollar_k, limits = c(0, NA),
                       expand = expansion(mult = c(0.02, 0.10))) +
    scale_x_continuous(expand = expansion(mult = c(0.07, 0.13))) +
    labs(title = "Lifetime cost-effectiveness frontier",
         subtitle = "Discounted cost and QALYs per patient; US payer perspective, 2024 USD",
         x = "Quality-adjusted life-years (QALY)", y = "Total cost") +
    theme_sma() +
    theme(legend.position = c(0.99, 0.03),
          legend.justification = c(1, 0),
          legend.text = element_text(size = 8.5))
}

## ---- Figure 2: two-panel cost-effectiveness acceptability curves -----------
## Panel A answers "is any DMT cost-effective versus best supportive care?"
## Panel B answers "given that a DMT will be used, which one?" — the BSC-excluded
## comparison reported in the second half of Section 3.2. Both panels plot the
## probability each strategy is optimal (highest net monetary benefit) at each
## willingness-to-pay value, over the same PSA draws.

ceac_df <- function(psa, arms, panel, wtp_grid) {
  nmb <- function(a, l) l * psa[[paste0("q_", a)]] - psa[[paste0("c_", a)]]
  P <- sapply(wtp_grid, function(l) {
    M <- sapply(arms, nmb, l = l)
    prop.table(table(factor(arms[max.col(M, "first")], levels = arms)))
  })
  do.call(rbind, lapply(seq_along(arms), function(i)
    data.frame(arm = arms[i], wtp = wtp_grid, p = P[i, ], panel = panel)))
}

fig_ceac_panels <- function(psa, wtp_grid = seq(0, 1e6, 10000),
                            ref_wtp = 2e5) {
  all_arms <- c("BSC", "Nusinersen", "OA", "Risdiplam")
  dmts     <- c("Nusinersen", "OA", "Risdiplam")

  df <- rbind(
    ceac_df(psa, all_arms, "A. All strategies (vs best supportive care)", wtp_grid),
    ceac_df(psa, dmts,     "B. Disease-modifying therapies only",         wtp_grid))
  df$arm   <- factor(df$arm, levels = all_arms)
  df$panel <- factor(df$panel, levels = unique(df$panel))

  ggplot(df, aes(wtp, p, colour = arm)) +
    geom_vline(xintercept = ref_wtp, linetype = 2,
               colour = "grey55", linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = arm_cols, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = dollar_k,
                       breaks = seq(0, 1e6, 250000),
                       expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_continuous(labels = percent, limits = c(0, 1),
                       breaks = seq(0, 1, 0.25)) +
    labs(title = "Cost-effectiveness acceptability curves",
         subtitle = paste0("Probability each strategy is optimal, by willingness-to-pay; ",
                           "dashed line = $200,000/QALY"),
         x = "Willingness-to-pay ($/QALY)", y = "Probability cost-effective") +
    theme_sma() +
    theme(legend.position  = "bottom",
          panel.spacing    = unit(1.0, "lines"),
          strip.text       = element_text(face = "bold", size = rel(0.92),
                                          hjust = 0, margin = margin(4, 4, 4, 4)))
}

## ---- Figure 2 (manuscript): CEAC built from the WORKBOOK's PSA draws --------
## IMPORTANT: the R port's own run_psa() does NOT reproduce the workbook PSA.
## Per README, transition uncertainty in the port uses an assumed effective
## sample size (N_EFF) because the workbook's stored transition SEs are on a
## transformed scale; N_EFF is a calibrated approximation. Section 3.2 of the
## manuscript reports the WORKBOOK's PSA, so the manuscript figure is built from
## the workbook's own 1,000 draws, exported to psa_excel_draws.csv
## (sheet "PSA", cols BE/BF, BI/BJ, BM/BN, BQ/BR = QALY/Cost per arm).
## Use fig_ceac_panels() for the R-port PSA; use this for the manuscript.

read_excel_psa <- function(path = "psa_excel_draws.csv") read.csv(path)

fig_ceac_manuscript <- function(psa = read_excel_psa(),
                                wtp_grid = seq(0, 1e6, 5000), ref_wtp = 2e5) {
  all_arms <- c("BSC", "Nusinersen", "OA", "Risdiplam")
  dmts     <- c("Nusinersen", "OA", "Risdiplam")
  df <- rbind(
    ceac_df(psa, all_arms, "A. All strategies (vs best supportive care)", wtp_grid),
    ceac_df(psa, dmts,     "B. Disease-modifying therapies only",         wtp_grid))
  df$arm   <- factor(df$arm, levels = all_arms)
  df$panel <- factor(df$panel, levels = unique(df$panel))

  ggplot(df, aes(wtp, p, colour = arm)) +
    geom_vline(xintercept = ref_wtp, linetype = 2,
               colour = "grey55", linewidth = 0.4) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = arm_cols, name = NULL, drop = FALSE) +
    scale_x_continuous(labels = dollar_k, breaks = seq(0, 1e6, 250000),
                       expand = expansion(mult = c(0.02, 0.03))) +
    scale_y_continuous(labels = percent, limits = c(0, 1),
                       breaks = seq(0, 1, 0.25)) +
    labs(title = "Cost-effectiveness acceptability curves",
         subtitle = paste0("Probability each strategy is optimal, by willingness-to-pay; ",
                           "dashed line = $200,000/QALY"),
         x = "Willingness-to-pay ($/QALY)", y = "Probability cost-effective") +
    theme_sma() +
    theme(legend.position = "bottom",
          panel.spacing   = unit(1.0, "lines"),
          strip.text      = element_text(face = "bold", size = rel(0.92),
                                         hjust = 0, margin = margin(4, 4, 4, 4)))
}

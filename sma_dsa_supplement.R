## =============================================================================
## sma_dsa_supplement.R  —  deterministic sensitivity analysis outputs for the
## supplementary material: a three-panel monochrome tornado figure and the
## underlying numeric table.
##
## Source AFTER sma_figures_gg.R, sma_figures_structure_frontier.R and
## sma_figures_mono.R.
##
## Parameters are varied as described in Section 2.9: health-state and drug
## costs by +/-20% (no usable variance in the source), and utilities and the
## treatment-specific non-sitting transitions by +/-1.96 standard errors.
## =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales) })

## ---- readable parameter labels ---------------------------------------------
STATE <- c(NS = "non-sitting", S = "sitting", W = "walking",
           PV = "permanent ventilation")
ARM   <- c(Nusinersen = "nusinersen", OA = "onasemnogene abeparvovec",
           Risdiplam = "risdiplam")

pretty_param <- function(x) {
  out <- x
  for (k in names(STATE)) {
    out <- sub(paste0("^Cost: ", k, " state$"),
               paste0("Health-state cost: ", STATE[[k]]), out)
    out <- sub(paste0("^Utility: ", k, "$"),
               paste0("Utility: ", STATE[[k]]), out)
  }
  for (k in names(ARM))
    out <- sub(paste0("^Drug cost: ", k, "$"),
               paste0("Drug cost: ", ARM[[k]]), out)
  out <- sub("^TP [A-Za-z]+: NS_S$",      "Transition: non-sitting to sitting", out)
  out <- sub("^TP [A-Za-z]+: NS_PV$",     "Transition: non-sitting to PV", out)
  out <- sub("^TP [A-Za-z]+: NS_Death$",  "Transition: non-sitting to death", out)
  out
}

## ---- collect DSA for several targets ---------------------------------------
dsa_all <- function(targets = c("Nusinersen", "OA", "Risdiplam"), ref = "BSC") {
  panels <- c(Nusinersen = "A. Nusinersen vs BSC",
              OA         = "B. Onasemnogene abeparvovec vs BSC",
              Risdiplam  = "C. Risdiplam vs BSC")
  out <- list(); bases <- list()
  for (tg in targets) {
    d <- run_dsa(tg, ref)
    d$param  <- pretty_param(as.character(d$param))
    d$target <- tg
    d$panel  <- panels[[tg]]
    d$base   <- attr(d, "base_icer")
    out[[tg]] <- d
    bases[[tg]] <- data.frame(panel = panels[[tg]], base = attr(d, "base_icer"))
  }
  res <- do.call(rbind, out)
  res$panel <- factor(res$panel, levels = unname(panels[targets]))
  attr(res, "bases") <- do.call(rbind, bases)
  rownames(res) <- NULL
  res
}

## ---- Supplementary Figure: three-panel tornado ------------------------------
fig_tornado_panels <- function(d = dsa_all()) {
  bases <- attr(d, "bases")
  bases$panel <- factor(bases$panel, levels = levels(d$panel))
  d$lo <- pmin(d$low, d$high); d$hi <- pmax(d$low, d$high)
  ## order bars within each panel by influence
  d <- d[order(d$panel, d$range), ]
  d$key <- factor(paste(d$panel, d$param, sep = "|"),
                  levels = paste(d$panel, d$param, sep = "|"))

  ggplot(d) +
    geom_vline(data = bases, aes(xintercept = base),
               linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    geom_segment(aes(y = key, yend = key, x = lo, xend = hi),
                 linewidth = 5.5, colour = "grey58", lineend = "butt") +
    facet_wrap(~ panel, ncol = 1, scales = "free") +
    scale_y_discrete(labels = function(z) sub("^.*\\|", "", z)) +
    scale_x_continuous(labels = dollar_k) +
    labs(x = "ICER ($/QALY gained)",
         y = NULL) +
    theme_sma() +
    theme(strip.text   = element_text(face = "bold", size = rel(0.9), hjust = 0,
                                      margin = margin(4, 4, 4, 4)),
          panel.spacing = unit(0.9, "lines"),
          axis.text.y   = element_text(size = rel(0.82)))
}

## ---- Supplementary Table: DSA numbers --------------------------------------
dsa_table <- function(d = dsa_all(), file = "dsa_results.csv") {
  t <- data.frame(
    Comparison   = sub("^[A-C]\\. ", "", as.character(d$panel)),
    Parameter    = d$param,
    Base_ICER    = round(d$base),
    ICER_low     = round(pmin(d$low, d$high)),
    ICER_high    = round(pmax(d$low, d$high)),
    Range        = round(d$range))
  t <- t[order(t$Comparison, -t$Range), ]
  write.csv(t, file, row.names = FALSE)
  invisible(t)
}

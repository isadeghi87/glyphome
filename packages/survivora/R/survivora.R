## survivora: Prognostic Landscape Cartography
##
## Treats hazard ratio as a topographic surface on a 2-D embedding.
## Bootstrap stability is rendered as contours, not as error bars.
##
## Visual grammar:
##   - Background terrain      = smoothed local log-hazard surface
##   - Stability contours      = bootstrap reproducibility of hazard sign
##   - Locality-resolved KM    = small-multiple Kaplan–Meier per region
##   - Forest of cluster HR    = effect size + bootstrap CI + reproducibility
##
## Figures generated:
##   A. Prognostic landscape   (terrain + stability + region anchors)
##   B. Locality-resolved KM   (per-region survival small multiples)
##   C. Reproducibility forest (cluster HR with bootstrap CI + agreement bar)

# ── Palettes ─────────────────────────────────────────────────────────────────
survivora_hazard_pal <- function() {
  # diverging: protective (blue) ↔ neutral (cream) ↔ adverse (red)
  c("#1d3557", "#457b9d", "#a8dadc", "#f1faee", "#fcbf49", "#e76f51", "#7a1f1f")
}

survivora_region_pal <- function() {
  c(R1 = "#264653", R2 = "#2a9d8f", R3 = "#e9c46a",
    R4 = "#f4a261", R5 = "#e76f51", R6 = "#6a994e")
}

# ── Data generators ──────────────────────────────────────────────────────────
survivora_landscape_data <- function(seed = 11) {
  set.seed(seed)
  # Six prognostic regions (clusters) on a 2-D embedding
  centers <- data.frame(
    region = c("R1","R2","R3","R4","R5","R6"),
    cx = c(-3.2, -1.0,  1.5,  3.5, -2.0,  2.0),
    cy = c( 2.5,  3.2,  2.8, -0.5, -2.5, -2.8),
    # log hazard ratio per region: protective (negative) to adverse (positive)
    log_hr = c(-0.85, -0.40,  0.15,  0.95, -0.10,  0.65),
    n      = c(60L, 80L, 70L, 90L, 55L, 75L),
    # bootstrap reproducibility of HR sign across resamples (0..1)
    repro  = c(0.92, 0.74, 0.45, 0.88, 0.60, 0.81),
    stringsAsFactors = FALSE
  )
  do.call(rbind, lapply(seq_len(nrow(centers)), function(i) {
    r <- centers[i, ]; nn <- r$n
    dx <- stats::rnorm(nn, 0, 0.55); dy <- stats::rnorm(nn, 0, 0.55)
    # event times follow region log-HR
    base_haz <- 0.08
    rate <- base_haz * exp(r$log_hr)
    t <- stats::rexp(nn, rate)
    cens <- stats::rexp(nn, 0.04)
    time <- pmin(t, cens, 60)
    event <- as.integer(t <= cens & t <= 60)
    data.frame(
      x = r$cx + dx, y = r$cy + dy, region = r$region,
      log_hr_local = r$log_hr + stats::rnorm(nn, 0, 0.10),
      reproducibility = pmin(1, pmax(0, r$repro + stats::rnorm(nn, 0, 0.05))),
      time = time, event = event,
      stringsAsFactors = FALSE
    )
  }))
}

# Smoothed hazard terrain on a regular grid
survivora_hazard_terrain <- function(df, gridn = 90L) {
  rx <- range(df$x) + c(-0.6, 0.6); ry <- range(df$y) + c(-0.6, 0.6)
  gx <- seq(rx[1], rx[2], length.out = gridn)
  gy <- seq(ry[1], ry[2], length.out = gridn)
  grid <- expand.grid(x = gx, y = gy)
  # Gaussian-kernel weighted local mean of log_hr_local
  bw <- 0.9
  sapply <- base::sapply
  hr <- sapply(seq_len(nrow(grid)), function(i) {
    d2 <- (df$x - grid$x[i])^2 + (df$y - grid$y[i])^2
    w  <- exp(-d2 / (2 * bw^2))
    sum(w * df$log_hr_local) / (sum(w) + 1e-9)
  })
  repro <- sapply(seq_len(nrow(grid)), function(i) {
    d2 <- (df$x - grid$x[i])^2 + (df$y - grid$y[i])^2
    w  <- exp(-d2 / (2 * bw^2))
    sum(w * df$reproducibility) / (sum(w) + 1e-9)
  })
  grid$log_hr <- hr
  grid$reproducibility <- repro
  grid
}

# Simple Kaplan–Meier estimator (no external deps)
survivora_km <- function(time, event) {
  ord <- order(time); time <- time[ord]; event <- event[ord]
  ut <- sort(unique(time[event == 1]))
  if (length(ut) == 0) return(data.frame(t = c(0, max(time)), s = c(1, 1)))
  at_risk <- sapply(ut, function(tt) sum(time >= tt))
  d <- sapply(ut, function(tt) sum(time == tt & event == 1))
  s <- cumprod(1 - d / at_risk)
  data.frame(t = c(0, ut), s = c(1, s))
}

# Bootstrap HR (reference = pooled hazard) per region
survivora_region_hr <- function(df, B = 200L) {
  base_rate <- sum(df$event) / sum(df$time)
  out <- do.call(rbind, lapply(unique(df$region), function(rg) {
    sub <- df[df$region == rg, ]
    obs_rate <- sum(sub$event) / sum(sub$time)
    obs_hr <- obs_rate / base_rate
    boots <- replicate(B, {
      idx <- sample(seq_len(nrow(sub)), nrow(sub), replace = TRUE)
      sb <- sub[idx, ]
      (sum(sb$event) / sum(sb$time)) / base_rate
    })
    data.frame(region = rg, hr = obs_hr,
               lo = stats::quantile(boots, 0.025, names = FALSE),
               hi = stats::quantile(boots, 0.975, names = FALSE),
               agree = mean(sign(log(boots)) == sign(log(obs_hr))),
               n = nrow(sub),
               stringsAsFactors = FALSE)
  }))
  out
}

# ── Plot A: Prognostic landscape ────────────────────────────────────────────
survivora_plot_landscape <- function() {
  df <- survivora_landscape_data()
  grid <- survivora_hazard_terrain(df)

  centers <- aggregate(cbind(x, y) ~ region, data = df, FUN = mean)

  pal <- survivora_hazard_pal()

  # Stability contour: only show ridge where reproducibility > 0.7
  grid_stable <- grid[grid$reproducibility > 0.70, ]

  p_terrain <- ggplot() +
    geom_raster(data = grid, aes(x = x, y = y, fill = log_hr), interpolate = TRUE) +
    geom_contour(data = grid, aes(x = x, y = y, z = log_hr),
                 colour = "#ffffff", alpha = 0.35, linewidth = 0.25, bins = 8) +
    geom_contour(data = grid_stable, aes(x = x, y = y, z = reproducibility),
                 breaks = c(0.7, 0.85), colour = "grey25",
                 linewidth = 0.55, linetype = "dashed") +
    geom_point(data = df, aes(x = x, y = y),
               colour = "white", fill = "grey15",
               shape = 21, size = 0.9, stroke = 0.15, alpha = 0.55) +
    geom_label(data = centers, aes(x = x, y = y, label = region),
               size = 3.2, label.size = 0.2, fill = "white", alpha = 0.85) +
    scale_fill_gradientn(colours = pal, name = "log HR",
                         limits = c(-1, 1), oob = scales::squish) +
    coord_fixed() +
    theme_void(base_size = 11) +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10)) +
    labs(title = "A. Prognostic landscape",
         subtitle = "Background terrain = local log hazard ratio; dashed contour = bootstrap-stable region (>=0.7 sign-agreement)")
  p_terrain
}

# ── Plot B: Locality-resolved KM small multiples ────────────────────────────
survivora_plot_locality_km <- function() {
  df <- survivora_landscape_data()
  km_all <- do.call(rbind, lapply(unique(df$region), function(rg) {
    sub <- df[df$region == rg, ]
    k <- survivora_km(sub$time, sub$event)
    k$region <- rg
    k
  }))
  km_pool <- survivora_km(df$time, df$event)
  km_pool_rep <- do.call(rbind, lapply(unique(df$region), function(rg) {
    k <- km_pool; k$region <- rg; k
  }))

  reg_pal <- survivora_region_pal()

  ggplot(km_all, aes(x = t, y = s)) +
    geom_step(data = km_pool_rep, colour = "grey70", linewidth = 0.5) +
    geom_step(aes(colour = region), linewidth = 1.0) +
    facet_wrap(~ region, nrow = 2) +
    scale_colour_manual(values = reg_pal, guide = "none") +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 10) +
    theme(strip.text = element_text(face = "bold"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 13),
          plot.subtitle = element_text(colour = "grey35", size = 10)) +
    labs(title = "B. Locality-resolved Kaplan–Meier",
         subtitle = "Coloured curve = region survival; grey = pooled cohort reference",
         x = "Time", y = "Survival")
}

# ── Plot C: Reproducibility forest ──────────────────────────────────────────
survivora_plot_repro_forest <- function() {
  df <- survivora_landscape_data()
  hr <- survivora_region_hr(df, B = 300L)
  hr <- hr[order(hr$hr), ]
  hr$region <- factor(hr$region, levels = hr$region)
  hr$direction <- ifelse(hr$hr < 1, "Protective", "Adverse")

  p_forest <- ggplot(hr, aes(y = region, x = hr)) +
    geom_vline(xintercept = 1, colour = "grey50", linetype = "dashed", linewidth = 0.4) +
    geom_segment(aes(x = lo, xend = hi, yend = region, colour = direction),
                 linewidth = 1.3) +
    geom_point(aes(fill = direction, size = n), shape = 21, colour = "grey20", stroke = 0.4) +
    scale_x_log10() +
    scale_colour_manual(values = c(Protective = "#1d3557", Adverse = "#7a1f1f"),
                        guide = "none") +
    scale_fill_manual(values = c(Protective = "#457b9d", Adverse = "#e76f51"),
                      name = "Direction") +
    scale_size_continuous(range = c(3, 7), name = "Cohort n") +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold", size = 12)) +
    labs(title = "C1. Hazard ratio per landscape region (bootstrap 95% CI)",
         x = "Hazard ratio (log scale)", y = NULL)

  p_agree <- ggplot(hr, aes(y = region, x = agree, fill = agree)) +
    geom_col(width = 0.6) +
    geom_vline(xintercept = 0.7, colour = "grey30", linetype = "dashed", linewidth = 0.4) +
    scale_fill_gradientn(colours = c("#e9c46a", "#a8dadc", "#1d3557"),
                         limits = c(0, 1), name = "Repro.") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.minor = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold", size = 12)) +
    labs(title = "C2. Bootstrap sign-agreement",
         x = "Fraction of bootstraps agreeing on HR sign", y = NULL)

  patchwork::wrap_plots(p_forest, p_agree, nrow = 1, widths = c(1.4, 1)) +
    patchwork::plot_annotation(
      title = "C. Reproducibility forest",
      subtitle = "Effect size and stability shown together — fragile regions visible at a glance"
    )
}

# ── Combined demo composer ──────────────────────────────────────────────────
survivora_plot_demo <- function() {
  pA <- survivora_plot_landscape()
  pB <- survivora_plot_locality_km()
  pC <- survivora_plot_repro_forest()
  patchwork::wrap_plots(pA, pB, pC, ncol = 1, heights = c(1.4, 1.0, 1.0)) +
    patchwork::plot_annotation(
      title = "survivora: Prognostic landscape cartography",
      subtitle = "Hazard as terrain · stability as contour · effect size as forest",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 15))
    )
}

# ── Save demo ────────────────────────────────────────────────────────────────
survivora_save_demo <- function(output_dir) {
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  ggplot2::ggsave(file.path(output_dir, "survivora_01_landscape.png"),
                  plot = survivora_plot_landscape(),
                  width = 11, height = 8, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "survivora_02_locality_km.png"),
                  plot = survivora_plot_locality_km(),
                  width = 12, height = 7, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "survivora_03_repro_forest.png"),
                  plot = survivora_plot_repro_forest(),
                  width = 11, height = 6, dpi = 300, bg = "white")
}

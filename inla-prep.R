fit_inla_models <- function(survey_short = c("QCS", "HS", "WCVI", "WCHG"),
  include_depth = FALSE, n_knots = 50, posterior_samples = 500,
  bootstrap_reps = 500) {

  survey_short <- match.arg(survey_short)
  survey <- switch(survey_short,
    QCS = "Queen Charlotte Sound Synoptic Survey",
    HS = "Hecate Strait Synoptic Survey",
    WCVI = "West Coast Vancouver Island Synoptic Survey",
    WCHG = "West Coast Haida Gwaii Synoptic Survey")

  d <- tidy_survey_sets(d_survey_sets,
    survey = survey,
    years = 2000:2017
  )
  sum(is.na(d$depth))
  d <- d[!is.na(d$depth), , drop = FALSE]
  d <- scale_survey_predictors(d)
  head(d)

  fake_yrs <- data.frame(year = unique(d$year),
    fake_year = seq_along(unique(d$year)))
  d <- inner_join(d, fake_yrs, by = "year")

  pg_out <- make_prediction_grid(d, region = survey_short)

  d$X <- d$X / 10 # for computational purposes
  d$Y <- d$Y / 10 # for computational purposes
  dat <- rename(d, orig_year = year, year = fake_year)

  cell_area <- pg_out$cell_area
  pg <- pg_out$grid
  pg$X <- pg$X / 10
  pg$Y <- pg$Y / 10

  # Fit models:

  binary_fitted <- FALSE
  m_inla_bin <- NA
  m_inla_pos <- NA

  giveup_i <- 0

  while (is.na(m_inla_bin[[1]]) && giveup_i <= 5) {
    m_inla_bin <- tryCatch({fit_inla(dat,
      response = "present", n_knots = n_knots, fit_model = TRUE,
      family = "binomial", plot = FALSE, kmeans = FALSE,
      include_depth = include_depth
    )}, error = function(e) NA)
    giveup_i <- giveup_i + 1
  }
  if (length(m_inla_bin) == 1) return(NA)

  dpos <- filter(dat, present == 1)
  dpos$density <- dpos$density * 1e3 # for computational purposes; too small otherwise

  giveup_i <- 0
  while (is.na(m_inla_pos[[1]]) && giveup_i <= 5) {
    m_inla_pos <- tryCatch({fit_inla(dpos,
      n_knots = n_knots, fit_model = TRUE,
      family = "gamma",
      response = "density", kmeans = FALSE,
      include_depth = include_depth
    )}, error = function(e) NA)
    giveup_i <- giveup_i + 1
  }
  if (length(m_inla_pos) == 1) return(NA)

  # Draw from posterior and project onto grid:

  pg_one_year <- filter(pg, year == min(pg$year))
  p_bin <- predict_inla(m_inla_bin,
    pred_grid = pg_one_year, samples = posterior_samples,
    include_depth = include_depth
  )
  p_pos <- predict_inla(m_inla_pos,
    pred_grid = pg_one_year, samples = posterior_samples,
    include_depth = include_depth
  )

  # Calculate probability * Gamma median values on the grid:

  out <- stats::plogis(p_bin) * exp(p_pos)
  out_median <- apply(out, c(1, 3), median)

  pg_all <- bind_rows(replicate(dim(p_bin)[3], pg_one_year, simplify = FALSE))
  pg_all$year <- rep(unique(pg$year), each = nrow(pg_one_year))
  out_long <- reshape2::melt(out_median)
  pg_all$pred <- out_long$value

  # Spatiotemporal plot:

  coast <- gfplot:::load_coastline(range(d$lon), range(d$lat), utm_zone = 9)
  isobath <- gfplot:::load_isobath(range(d$lon), range(d$lat),
    bath = c(100, 200, 500), utm_zone = 9
  )
  map_padding <- c(-5, 5)

  g <- ggplot(pg_all, aes(X * 10, Y * 10)) +
    geom_raster(aes(fill = sqrt(pred))) +
    scale_fill_viridis_c(option = "C") +
    facet_wrap(~ year) +
    geom_polygon(
      data = coast, aes_string(x = "X", y = "Y", group = "PID"),
      fill = "grey80"
    ) +
    geom_point(
      data = filter(d, present == 1),
      aes(x = X * 10, y = Y * 10, size = density), inherit.aes = FALSE,
      pch = 21, col = "white", alpha = 0.3
    ) +
    geom_point(
      data = filter(d, present == 0),
      aes(x = X * 10, y = Y * 10), inherit.aes = FALSE,
      pch = 4, col = "white", alpha = 0.2
    )

  g <- g + geom_path(
    data = isobath, aes_string(
      x = "X", y = "Y",
      group = "paste(PID, SID)"
    ),
    inherit.aes = FALSE, lwd = 0.4, col = "grey70", alpha = 0.4
  )

  g <- g + xlab("UTM 9N Easting (km)") + ylab("UTM 9N Northing (km)") +
    scale_size(range = c(0, 9)) +
    guides(size = FALSE, fill = FALSE) +
    theme_pbs() +
    coord_equal(
      expand = FALSE, xlim = range(pg_all$X * 10) + map_padding,
      ylim = range(pg_all$Y * 10) + map_padding
    )

  # Annual index plot:

  yr_est <- apply(out, c(2, 3), sum)
  yr_est <- reshape2::melt(yr_est) %>%
    rename(i = Var1, fake_year = Var2) %>%
    inner_join(fake_yrs, by = "fake_year")

  # Calculate bootstrapped design-based index:

  surv <- d_survey_sets %>%
    filter(survey_series_desc == survey) %>%
    select(-sample_id) %>%
    unique()

  calc_bio <- function(dat, i = seq_len(nrow(dat))) {
    dat[i, , drop = FALSE] %>%
      group_by(year, survey_id, area_km2, grouping_code) %>%
      summarise(density = mean(density_kgpm2 * 1e6)) %>%
      group_by(year) %>%
      summarise(biomass = sum(density * area_km2)) %>%
      pull(biomass)
  }

  boot_biomass <- function(dat, reps = bootstrap_reps) {
    out <- dat %>%
      group_by(year, species_common_name, survey_series_desc) %>%
      do({
        b <- boot::boot(., statistic = calc_bio, strata = .$grouping_code, R = reps)
        suppressWarnings(bci <- boot::boot.ci(b, type = "perc"))
        dplyr::tibble(
          mean_boot = mean(b$t),
          median_boot = median(b$t),
          lwr = bci$percent[[4]],
          upr = bci$percent[[5]],
          cv = sd(b$t) / mean(b$t),
          biomass = calc_bio(.)
        )
      })
  }

  out_boot <- boot_biomass(surv)

  # Compare:

  yr_est_summ <- group_by(yr_est, year) %>%
    summarise(
      lwr = quantile(value * 1e3 * cell_area, probs = 0.025),
      upr = quantile(value * 1e3 * cell_area, probs = 0.975),
      biomass = quantile(value * 1e3 * cell_area, probs = 0.5)
    ) %>%
    mutate(type = "spatiotemporal")

  both <- bind_rows(
    yr_est_summ,
    mutate(out_boot, type = "design-based bootstrap")
  )

  g_index <- ggplot(both, aes(as.numeric(year), biomass,
    fill = type, ymin = lwr, ymax = upr
  )) +
    geom_ribbon(alpha = 0.2, colour = NA) +
    geom_line(lwd = 1, aes(colour = type)) +
    theme_pbs() +
    ylim(0, NA) +
    scale_fill_manual(values = c("#474747", "#2188ff")) +
    scale_colour_manual(values = c("#474747", "#2188ff")) +
    ylab("Biomass") + xlab("Year") +
    labs(fill = "Type", colour = "Type")

  list(
    index_plot = g_index, map_plot = g, index_data = both,
    model_bin = m_inla_bin, model_pos = m_inla_pos,
    prediction_grid = pg_all, data = d,
    isobath = isobath, coast = coast,
    survey = survey, survey_short = survey_short
  )
}

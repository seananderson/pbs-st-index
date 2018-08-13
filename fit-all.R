library(dplyr)
library(ggplot2)
library(gfplot)
library(dplyr)
library(INLA)
source("inla-prep.R")
source("fit-inla.R")

# -------------------------------
# setup:
species <- "pacific cod"
species_short <- "pcod"
include_depth <- TRUE
depth_text <- if (include_depth) "depth" else "no-depth"

# -------------------------------
# data:
# faked right now
# gfplot::get_survey_sets(species)

d_survey_sets <- readRDS("../gfsynopsis/report/data-cache/pbs-survey-sets.rds")
d_survey_sets <- dplyr::filter(d_survey_sets, species_common_name == species)

# -------------------------------
# fit models:

out <- lapply(c("QCS", "HS", "WCVI", "WCHG"), function(surv) {
  message("\nFitting model for the survey ", surv)
  x <- fit_inla_models(surv,
    bootstrap_reps = 800, include_depth = include_depth,
    n_knots = 50, posterior_samples = 1000
  )
  x
})

# -------------------------------
# cache output:
# readr::save_rds(out, file = "2018-03-28-models.rds")
# out <- readr::read_rds("2018-03-28-models.rds")

# -------------------------------
# index plots:

index_dat <- purrr::map_df(out, function(x) {
  data.frame(x$index_data, survey = x$survey_short, stringsAsFactors = FALSE)
})

g_index <- ggplot(index_dat, aes(as.numeric(year), biomass,
  fill = type, ymin = lwr, ymax = upr
)) +
  geom_ribbon(alpha = 0.2, colour = NA) +
  geom_line(lwd = 1, aes(colour = type)) +
  theme_pbs() +
  facet_wrap(~survey) +
  ylim(0, NA) +
  scale_fill_manual(values = c("#474747", "#2188ff")) +
  scale_colour_manual(values = c("#474747", "#2188ff")) +
  ylab("Biomass") + xlab("Year") +
  labs(fill = "Type", colour = "Type") +
  coord_cartesian(expand = FALSE)
g_index
ggsave(paste0("figs/", species_short, "-index-", depth_text, ".pdf"),
  width = 9, height = 5.5)

# -------------------------------
# maps:

d <- d_survey_sets %>% dplyr::filter(survey_series_desc %in%
  c(
    "Queen Charlotte Sound Synoptic Survey",
    "Hecate Strait Synoptic Survey",
    "West Coast Vancouver Island Synoptic Survey",
    "West Coast Haida Gwaii Synoptic Survey"
  )) %>%
  dplyr::mutate(present = ifelse(density_kgpm2 > 1, 1, 0))

coast <- gfplot:::load_coastline(range(d$longitude),
  range(d$latitude),
  utm_zone = 9
)
isobath <- gfplot:::load_isobath(range(d$longitude),
  range(d$latitude),
  bath = c(100, 200, 500), utm_zone = 9
)
map_padding <- c(-5, 5)

prediction_data <- lapply(out, function(x) {
  mutate(x$prediction_grid, survey_short = x$survey_short)
}) %>% dplyr::bind_rows()

# grid:
g <- ggplot(prediction_data, aes(X * 10, Y * 10, group = survey_short)) +
  geom_raster(aes(fill = sqrt(pred))) +
  scale_fill_viridis_c(option = "C") +
  facet_wrap(~ year) +
  geom_polygon(
    data = coast, aes_string(x = "X", y = "Y", group = "PID"),
    fill = "grey80"
  )

## include_tow dots?

# g <- g + geom_point(
#   data = filter(out[[1]]$data, present == 1),
#   aes(x = X * 10, y = Y * 10, size = density), inherit.aes = FALSE,
#   pch = 21, col = "white", alpha = 0.3
# ) +
#   geom_point(
#     data = filter(out[[1]]$data, present == 0),
#     aes(x = X * 10, y = Y * 10), inherit.aes = FALSE,
#     pch = 4, col = "white", alpha = 0.2
#   )

# isobath:
g <- g + geom_path(
  data = isobath, aes_string(
    x = "X", y = "Y",
    group = "paste(PID, SID)"
  ),
  inherit.aes = FALSE, lwd = 0.4, col = "grey70", alpha = 0.4
)

# aesthetics:
g <- g + xlab("UTM 9N Easting (km)") + ylab("UTM 9N Northing (km)") +
  scale_size(range = c(0, 9)) +
  guides(size = FALSE, fill = FALSE) +
  theme_light() +
  coord_equal(
    expand = FALSE, xlim = range(prediction_data$X * 10) + map_padding,
    ylim = range(prediction_data$Y * 10) + map_padding
  )

ggsave(plot = g,
  filename = paste0("figs/", species_short, "-map-", depth_text, ".pdf"),
  width = 12, height = 14)

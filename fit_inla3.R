fit_inla3 <- function(d, response = "present", n_knots = 50,
                      family = "binomial", plot = FALSE) {
  library(dplyr)
  library(INLA)

  coords <- select(d, X, Y) %>% unique()

  # Use kmeans() to calculate centers
  km <- stats::kmeans(x = coords, centers = n_knots)
  mesh <- INLA::inla.mesh.create(km$centers)
  spde <- INLA::inla.spde2.matern(mesh, alpha = 3 / 2)

  # bnd <- inla.nonconvex.hull(km$centers, 20, 60)
  # mesh <- inla.mesh.2d(loc = km$centers, boundary = bnd,
  #   max.edge = c(25, 200), cutoff = c(25, 25))
  # spde <- inla.spde2.matern(mesh, alpha = 3 / 2)

  if (plot) {
    plot(mesh)
    points(coords)
  }

  d$year <- d$year - min(d$year) + 1
  k <- max(d$year)
  dat <- data.frame(
    y = d[, response, drop = TRUE],
    time = d$year,
    xcoo = d$X,
    ycoo = d$Y
  )

  # Make a design matrix where the first year is the intercept
  YEARS.lab <- paste0("Y", seq(1, k))
  dat[YEARS.lab] <- 0
  dat[, YEARS.lab[1]] <- 1
  for (j in seq_along(YEARS.lab)) {
    dat[dat$time == j, YEARS.lab[j]] <- 1
  }

  dat$depth <- d$depth_scaled
  dat$depth2 <- d$depth_scaled2

  # Construct index for iid / ar1 model
  iset <- INLA::inla.spde.make.index(
    "i2D",
    n.spde = mesh$n,
    n.group = k
  )

  # Make the covariates
  X.1 <- dat[, -c(1:4)]
  Covar.names <- colnames(X.1)
  XX.list <- as.list(X.1)
  effect.list <- list()
  effect.list[[1]] <- c(iset)
  for (Z in seq_len(ncol(X.1)))
    effect.list[[Z + 1]] <- XX.list[[Z]]
  names(effect.list) <- c("1", Covar.names)

  # Make projection points stack
  A <- INLA::inla.spde.make.A(
    mesh = mesh,
    loc = cbind(dat$xcoo, dat$ycoo),
    group = dat$time
  )
  A.list <- list()
  A.list[[1]] <- A
  for (Z in seq_len(ncol(X.1)))
    A.list[[Z + 1]] <- 1

  # Make projection points stack
  sdat <- INLA::inla.stack(
    tag = "stdata",
    data = list(y = dat$y),
    A = A.list,
    effects = effect.list
  )

  formula <- as.formula(paste0(
    "y ~ -1 + ",
    paste(Covar.names, collapse = "+"),
    "+ f(i2D, model=spde, group = i2D.group, control.group = list(model='ar1'))"
  ))

  model <- try(INLA::inla(formula,
    family = family,
    data = INLA::inla.stack.data(sdat),
    control.predictor = list(compute = TRUE, A = INLA::inla.stack.A(sdat)),
    control.compute = list(config = TRUE),
    verbose = TRUE,
    debug = FALSE,
    keep = FALSE
  ), silent = TRUE)

  list(
    model = model, mesh = mesh, spde = spde, data = dat, formula = formula,
    iset = iset
  )
}

predict_inla3 <- function(obj, pred_grid, samples = 100L) {
  library(INLA)

  mesh <- obj$mesh
  model <- obj$model
  iset <- obj$iset

  inla.mcmc <- INLA::inla.posterior.sample(n = samples, model)

  na <- rownames(inla.mcmc[[1]]$latent)
  re.indx <- grep("^i2D", na)
  depth.indx1 <- grep("^depth$", na)
  depth.indx2 <- grep("^depth2$", na)
  fe.indx <- grep("^Y[0-9]+$", na)

  # Read in locations and knots, form projection matrix
  knot_locs <- mesh$loc[, c(1, 2)]

  gridLocs <- pred_grid[, c("X", "Y")]
  # Projections will be stored as an array
  projectedLatentGrid <- array(0, dim = c(nrow(gridLocs), samples, length(fe.indx)))

  projMatrix <- INLA::inla.spde.make.A(mesh, loc = as.matrix(gridLocs))

  # YEARS.lab <- paste0("Y", seq_along(fe.indx))
  # mm <- matrix(nrow = nrow(pred_grid), ncol = length(fe.indx) + 2, data = 0)
  # colnames(mm) <- c(YEARS.lab, "depth", "depth2")
  #
  # mm[,1] <- 1
  #
  # for (j in YEARS.lab) {
  #   mm[colnames(mm) == YEARS.lab, YEARS.lab[j]] <- 1
  # }

  # dat$depth <- d$depth_scaled
  # dat$depth2 <- d$depth_scaled2

  for (yr in 1) {
    # Grab random effects from this year
    indx <- which(iset$i2D.group == yr)
    for (i in seq_len(samples)) {
      projectedLatentGrid[, i, yr] <-
        as.numeric(projMatrix %*% inla.mcmc[[i]]$latent[re.indx][indx]) +
        pred_grid[, "depth_scaled"] * inla.mcmc[[i]]$latent[depth.indx1] +
        pred_grid[, "depth_scaled2"] * inla.mcmc[[i]]$latent[depth.indx2] +
        inla.mcmc[[i]]$latent[fe.indx[[yr]]]
    }
  }

  # Loop over all subsequent years, and do MCMC projections for that year
  for (yr in seq(2, length(fe.indx))) {

    # Grab random effects from this year
    indx <- which(iset$i2D.group == yr)
    for (i in seq_len(samples)) {
      projectedLatentGrid[, i, yr] <-
        as.numeric(projMatrix %*% inla.mcmc[[i]]$latent[re.indx][indx]) +
        pred_grid[, "depth_scaled"] * inla.mcmc[[i]]$latent[depth.indx1] +
        pred_grid[, "depth_scaled2"] * inla.mcmc[[i]]$latent[depth.indx2] +
        inla.mcmc[[i]]$latent[fe.indx[[1]]] +
        inla.mcmc[[i]]$latent[fe.indx[[yr]]]
    }
  }

  projectedLatentGrid
}













#
#
#
#
#
#
#
#
#
#
#
#   b1 <- plyr::laply(seq_len(samples), function(i) inla.mcmc.pos[[i]]$latent[b1_])
#   b2 <- plyr::laply(seq_len(samples), function(i) inla.mcmc.pos[[i]]$latent[b2_])
#
#   b1_bin <- plyr::laply(seq_len(samples), function(i) inla.mcmc.bin[[i]]$latent[b1_])
#   b2_bin <- plyr::laply(seq_len(samples), function(i) inla.mcmc.bin[[i]]$latent[b2_])
#
#   projMatrix_bin <- INLA::inla.spde.make.A(mesh_bin,
#     loc = as.matrix(pred_grid[, c("X", "Y")]))
#   projMatrix_pos <- INLA::inla.spde.make.A(mesh_pos,
#     loc = as.matrix(pred_grid[, c("X", "Y")]))
#
#   pp <- plyr::laply(seq_len(samples), function(i) {
#     as.numeric(projMatrix %*% inla.mcmc.pos[[i]]$latent[sf]) +
#       b1[i] * pred_grid$depth_scaled +
#       b2[i] * pred_grid$depth_scaled2
#   })
#   pb <- plyr::laply(seq_len(samples), function(i) {
#     as.numeric(projMatrix %*% inla.mcmc.bin[[i]]$latent[sf]) +
#       b1_bin[i] * pred_grid$depth_scaled +
#       b2_bin[i] * pred_grid$depth_scaled2
#   })
#   pc <- stats::plogis(pb) * exp(pp)
#
#   # spatial field only:
#   spp <- plyr::laply(seq_len(samples), function(i) {
#     as.numeric(projMatrix %*% inla.mcmc.pos[[i]]$latent[sf])
#   })
#   spb <- plyr::laply(seq_len(samples), function(i) {
#     as.numeric(projMatrix %*% inla.mcmc.bin[[i]]$latent[sf])
#   })
#
#   pred_grid$pred_delta <- apply(pc, 2L, stats::median)
#   pred_grid$pred_binary <- apply(stats::plogis(pb), 2L, stats::median)
#   pred_grid$pred_positive <- apply(exp(pp), 2L, stats::median)
#   pred_grid$spatial_field <- apply(pc, 2L, stats::median)
#   pred_grid$spatial_field_binary <- apply(spb, 2L, stats::median)
#   pred_grid$spatial_field_positive <- apply(spp, 2L, stats::median)
#
#   list(
#     pred = pred_grid, prediction_posterior = pc,
#     params = list(
#       b0 = b0, b1 = b1, b2 = b2, b0_bin = b0_bin, b1_bin = b1_bin,
#       b2_bin = b2_bin
#     )
#   )
# }

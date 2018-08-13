fit_inla <- function(dat, plot = FALSE, max.edge = c(1, 3), convex = 1.5) {
  library(INLA)

  pos_dat <- filter(dat, present == 1)
  coords <- cbind(dat$X10, dat$Y10)
  coords_pos <- cbind(pos_dat$X10, pos_dat$Y10)

  bnd <- inla.nonconvex.hull(coords, convex = convex)
  mesh6 <- inla.mesh.2d(
    boundary = bnd,
    max.edge = max.edge,
    cutoff = 0.01,
    offset = 1.5
  )

  if (plot) {
    plot(mesh6)
    points(dat$X, dat$Y, col = "red")
  }

  A.est6 <- inla.spde.make.A(mesh = mesh6, loc = coords)
  A.est6.pos <- inla.spde.make.A(mesh = mesh6, loc = coords_pos)
  spde <- inla.spde2.matern(mesh = mesh6, alpha = 1.5)
  formula <- y ~ -1 + intercept + depth_scaled + depth_scaled2 +
    f(spatial.field, model = spde)
  s.index <- inla.spde.make.index(name = "spatial.field", n.spde = spde$n.spde)

  stack.est.pos <- inla.stack(
    data = list(y = log(pos_dat$density)), A = list(A.est6.pos, 1, 1),
    effects = list(
      c(s.index, list(intercept = 1)),
      list(depth_scaled = pos_dat$depth_scaled),
      list(depth_scaled2 = pos_dat$depth_scaled2)
    ), tag = "est_pos"
  )
  output6.stack.pos <- inla(formula,
    data = inla.stack.data(stack.est.pos, spde = spde),
    family = "gaussian", control.predictor = list(
      A = inla.stack.A(stack.est.pos),
      compute = TRUE
    ), verbose = FALSE, control.compute = list(config = TRUE)
  )

  stack.est.bin <- inla.stack(
    data = list(y = dat$present), A = list(A.est6, 1, 1),
    effects = list(
      c(s.index, list(intercept = 1)),
      list(depth_scaled = dat$depth_scaled),
      list(depth_scaled2 = dat$depth_scaled2)
    ), tag = "est_bin"
  )
  output6.stack.bin <- inla(formula,
    data = inla.stack.data(stack.est.bin, spde = spde),
    family = "binomial", control.predictor = list(
      A = inla.stack.A(stack.est.bin),
      compute = TRUE
    ), verbose = FALSE, control.compute = list(config = TRUE)
  )

  list(pos = output6.stack.pos, bin = output6.stack.bin, mesh = mesh6)
}

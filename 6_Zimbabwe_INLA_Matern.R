
library(geojsonio)
library(sp)
library(INLA)

# Matern examples: 
# https://www.paulamoraga.com/book-geospatial/sec-geostatisticaldataexamplespatial.html?q=matern#building-the-spde-model-on-the-mesh-1
# https://becarioprecario.bitbucket.io/spde-gitbook/ch-intro.html#spde-model-definition
# https://groups.google.com/g/r-inla-discussion-group/c/ho1tPj7MtF0 - specify priors

#---------------------------------------------------------------------
# load data
#---------------------------------------------------------------------

spdf <- geojson_read("Zimbabwe_HIV/zwe2016phia.geojson",  what = "sp")
names(spdf)
class(spdf)

# plot borders
plot(spdf)

# plot centroids
points(spdf$center_x, spdf$center_y, col="red", pch=16)

length(spdf$center_x)

x_coords = spdf$center_x
y_coords = spdf$center_y

# x_coords = (x_coords - min(x_coords))/(max(x_coords) - min(x_coords))
# y_coords = (y_coords - min(y_coords))/(max(y_coords) - min(y_coords))

#plot(x_coords, y_coords, col="red", pch=16, main="Centroids - normalised")

y = round(spdf$y)
n_obs = round(spdf$n_obs)
estimate = spdf$estimate

# plot(x = c(0,0.3), y = c(0.,0.3), type="l")
# points(estimate, y/n_obs, pch=16, col="blue")

ids <- str_sub(spdf$area_id, start = -2)
ids <- str_replace(ids, "_", "")
ids <- as.numeric(ids)

df_model <- data.frame(ids=ids, x_coords, y_coords, y, n_obs, estimate)
head(df_model)

#---------------------------------------------------------------------
# build model - Matern
#---------------------------------------------------------------------

coo <- cbind(x_coords, y_coords)
mesh <- inla.mesh.2d(loc = coo, max.edge = c(0.1, 5),cutoff = 0.2)
mesh$n
plot(mesh)
points(coo, col = "red")

# alpha = nu + d/2 -> [0,2]
# alpha = 2
# nu = alpha - 1 = 1

# spde <- inla.spde2.matern(
#   # Mesh and smoothness parameter
#   mesh = mesh, alpha = 2, 
#   prior.variance.nominal = 0.1,
#   prior.range.nominal = 1,
#   constr = TRUE)

spde <- inla.spde2.matern(mesh = mesh, alpha = 2,
                          # hyper = list(theta1 = list(param = c(a1, a2), prior = "loggamma"),
                          #              theta2 = list(param = c(b1, b2), prior = "loggamma"))
                          #hyper = list(theta1 = list(param = c(a1, a2), prior = "loggamma"),
                          #             theta2 = list(param = c(b1, b2), prior = "loggamma"))
                          )


indexs <- inla.spde.make.index("s", spde$n.spde)
lengths(indexs)

A <- inla.spde.make.A(mesh = mesh, loc = coo)
dim(A)

# prediction points are the same as observation points
coop <- coo

Ap <- inla.spde.make.A(mesh = mesh, loc = coop)
dim(Ap)

# stack for estimation stk.e
stk.e <- inla.stack(
  tag = "est",
  data = list(y = df_model$y, numtrials = df_model$n_obs),
  A = list(1, A),
  effects = list(data.frame(b0 = rep(1, nrow(coo))), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA, n_obs = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = rep(1, nrow(coop))), s = indexs)
)

# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ b0 + f(s, model = spde)

res <- inla(formula,
            family = "binomial", 
            Ntrials = numtrials,
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(compute = TRUE, link = 1, A = inla.stack.A(stk.full))
)

idx <- inla.stack.index(stk.p, "pred")$data
fitted <- res$summary.fitted.values[idx, c("0.5quant")]


#---------------------------------------------------------------------
# Plot results
#---------------------------------------------------------------------

plot(x = c(0,0.3), y = c(0.,0.3), type="l", xlab="raw estimate (y/n_obs)", ylab="INLA estimate, Matern1")
points(estimate, fitted, pch = 16, col = "red")
legend(0.01, 0.24, legend = c("INLA estimates, Matern1"), 
       col = c("red"), pch = c(16, 16), 
       cex = 0.8, title = "", text.font = 4, bg = "lightblue")

spdf@data['fitted'] <- fitted
spplot(spdf, zcol="fitted")


#---------------------------------------------------------------------
# Hyperparameters
#---------------------------------------------------------------------
(hyper <- res$summary.hyperpar)
dim(hyper)

interpretable_res <- inla.spde2.result(inla = res, name = "s", spde = spde)
str(interpretable_res)

variance <- interpretable_res$marginals.variance.nominal$variance.nominal.1
range <- interpretable_res$marginals.range.nominal$range.nominal.1
#theta1 <- res$marginals.hyperpar$`Theta1 for s` # theta1 = log(tau)
theta2 <- res$marginals.hyperpar$`Theta2 for s` # theta2 = log(kappa)

lengthscale <- inla.tmarginal(function(x, nu = 1) {
  kappa <- exp(x)
  (2 * nu) ** (1 / 2) / kappa
}, theta2)


par(mfrow = c(2, 1), mar = c(3, 3, 1, 1), mgp = c(2, 1, 0))
plot(variance, type = "l", xlab = expression(sigma^2))
plot(lengthscale, type = "l",
     xlab = 'lengthscale', ylab = 'Posterior density',
     xlim = c(0,10))


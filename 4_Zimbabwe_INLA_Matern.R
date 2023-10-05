
library(geojsonio)
library(sp)
library(INLA)

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

x_coords = (x_coords - min(x_coords))/(max(x_coords) - min(x_coords))
y_coords = (y_coords - min(y_coords))/(max(y_coords) - min(y_coords))

plot(x_coords, y_coords, col="red", pch=16, main="Centroids")

y = round(spdf$y)
n_obs = round(spdf$n_obs)
estimate = spdf$estimate

plot(x = c(0,0.3), y = c(0.,0.3), type="l")
points(estimate, y/n_obs, pch=16, col="blue")

df_model <- data.frame(area_id=spdf$area_id, x_coords, y_coords, y, n_obs, estimate)
head(df_model)

#---------------------------------------------------------------------
# build model - Matern
#---------------------------------------------------------------------

coo <- cbind(x_coords, y_coords)
mesh <- inla.mesh.2d(
  loc = coo, max.edge = c(0.1, 5),
  cutoff = 0.01
)
mesh$n
plot(mesh)
points(coo, col = "red")

# alpha = nu + d/2
# nu = alpha - d/2 = alpha - 1
# nu = 5/2 = alpha -1, hence alpha = 5/2+1 ?
spde <- inla.spde2.matern(mesh = mesh, alpha = 2, constr = TRUE)

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
  effects = list(data.frame(b0 = 1), s = indexs)
)

# stack for prediction stk.p
stk.p <- inla.stack(
  tag = "pred",
  data = list(y = NA, n_obs = NA),
  A = list(1, Ap),
  effects = list(data.frame(b0 = 1), s = indexs)
)
# Error in (function (data, A, effects, tag = "", compress = TRUE, remove.unused = TRUE)  : 
#             Row count mismatch for A: 1,63



# stk.full has stk.e and stk.p
stk.full <- inla.stack(stk.e, stk.p)

formula <- y ~ b0 + f(s, model = spde)

res <- inla(formula,
            family = "binomial", 
            Ntrials = numtrials,
            control.family = list(link = "logit"),
            data = inla.stack.data(stk.full),
            control.predictor = list(
              compute = TRUE, link = 1,
              A = inla.stack.A(stk.full)
            )
)






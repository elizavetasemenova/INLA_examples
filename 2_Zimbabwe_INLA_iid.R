
library(geojsonio)
library(sp)
library(INLA)

#---------------------------------------------------------------------
# load data
#---------------------------------------------------------------------

spdf <- geojson_read("Zimbabwe_HIV/zwe2016phia.geojson",  what = "sp")
names(spdf)
class(spdf)

# # plot borders
# plot(spdf)
# 
# # plot centroids
# points(spdf$center_x, spdf$center_y, col="red", pch=16)
# 
# length(spdf$center_x)
# 
# x_coords = spdf$center_x
# y_coords = spdf$center_y
# 
# x_coords = (x_coords - min(x_coords))/(max(x_coords) - min(x_coords))
# y_coords = (y_coords - min(y_coords))/(max(y_coords) - min(y_coords))
# 
# plot(x_coords, y_coords, col="red", pch=16, main="Centroids")

y = round(spdf$y)
n_obs = round(spdf$n_obs)
estimate = spdf$estimate

# plot(x = c(0,0.3), y = c(0.,0.3), type="l")
# points(estimate, y/n_obs, pch=16, col="blue")

df_model <- data.frame(area_id=spdf$area_id, y, n_obs, estimate)
head(df_model)

#---------------------------------------------------------------------
# build model - iid
#---------------------------------------------------------------------

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.01)))
formula <- y ~  f(area_id, model = "iid", hyper = prior.prec)

res <- inla(formula,
            data = df_model,
            family = "binomial", 
            Ntrials = n_obs,
            control.predictor = list(compute = TRUE),
            control.compute = list(config = TRUE)
)

summary(res)

# Results summary
# res$summary.fixed
# res$summary.random
# res$summary.hyperpar
# res$marginals.fixed

fitted <- res$summary.fitted.values
fitted_0.5quant <- fitted$`0.5quant`


plot(x = c(0,0.3), y = c(0.,0.3), type="l", xlab="raw estimate (y/n_obs)", ylab="INLA estimate, iid")
points(estimate, fitted_0.5quant, pch=16, col="red")
legend(0.01, 0.27, 
       legend=c("Raw estimates", "INLA estimates, iid"),
       col=c("blue", "red"), 
       pch= c(16, 16),
       cex=0.8,
       title="", text.font=4, bg='lightblue')



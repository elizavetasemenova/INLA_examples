
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

y = round(spdf$y)
n_obs = round(spdf$n_obs)
estimate = spdf$estimate

plot(x = c(0,0.3), y = c(0.,0.3), type="l")
points(estimate, y/n_obs, pch=16, col="blue")

#df_model <- data.frame(area_id=spdf$area_id, x_coords, y_coords, y, n_obs, estimate)
df_model <- data.frame(area_id=spdf$area_id, y, n_obs, estimate)
head(df_model)

#---------------------------------------------------------------------
# build model - iid
#---------------------------------------------------------------------

# prior.prec <- list(prec = list(prior = "pc.prec",
#                                param = c(1, 0.01)))

plot(spdf)
library(spdep)
nb <- poly2nb(spdf)
head(nb)

nb2INLA("spdf.adj", nb)
g <- inla.read.graph(filename = "spdf.adj")

#formula <- y ~ f(area_id, model = "besag", graph = g, scale.model = TRUE) 
formula <- y ~ f(area_id, model = "besag", graph = g) 

res <- inla(formula,
            data = df_model,
            family = "binomial", 
            Ntrials = n_obs,
            control.predictor = list(compute = TRUE),
            control.compute = list(config = TRUE)
)

# Error in inla.core(formula = formula, family = family, contrasts = contrasts,  : 
#                      In f(area_id): 'covariate' must match 'values',  and both must either be 'numeric', or 'factor'/'character'.
#                    
#                    *** inla.core.safe:  inla.program has crashed: rerun to get better initial values. try=1/1 
#                    Error in inla.core(formula = formula, family = family, contrasts = contrasts,  : 
#                                         In f(area_id): 'covariate' must match 'values',  and both must either be 'numeric', or 'factor'/'character'.
#                                       Error in inla.core.safe(formula = formula, family = family, contrasts = contrasts,  : 
#                                                                 *** Failed to get good enough initial values. Maybe it is due to something else.

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


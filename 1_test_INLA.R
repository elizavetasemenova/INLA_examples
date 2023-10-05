#R.Version()
#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

sessionInfo()

library(geojsonio)
library(sp)
library(INLA)
#inla.upgrade()
inla.version()

#inla.upgrade(testing=TRUE)

#---------------------------------------------------------------------
# build model
#---------------------------------------------------------------------

# follow example here: https://www.paulamoraga.com/book-geospatial/sec-inla.html

prior.prec <- list(prec = list(prior = "pc.prec",
                               param = c(1, 0.01)))
formula <- r ~ f(hospital, model = "iid", hyper = prior.prec)

res <- inla(formula,
            data = Surg,
            family = "binomial", 
            Ntrials = n,
            control.predictor = list(compute = TRUE),
            control.compute = list(dic = TRUE)
)

summary(res)
          


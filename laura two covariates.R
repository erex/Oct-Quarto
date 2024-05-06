library(dsims)

# Multi-strata example (make sf shape)
s1 = matrix(c(0,0,0,2,1,2,1,0,0,0),ncol=2, byrow=TRUE)
s2 = matrix(c(1,0,1,2,2,2,2,0,1,0),ncol=2, byrow=TRUE)
pol1 = sf::st_polygon(list(s1))
pol2 = sf::st_polygon(list(s2))
sfc <- sf::st_sfc(pol1,pol2)
strata.names <- c("low", "high")

mytrunc <- 0.2
sf.pol <- sf::st_sf(strata = strata.names, geom = sfc)

region <- make.region(region.name = "Multi-strata Eg",
                      strata.name = strata.names,
                      shape = sf.pol)


mydesign <- make.design(region=region, design="systematic", 
                        truncation=mytrunc, samplers = 30)
somelines <- generate.transects(mydesign)
plot(region, somelines, covered.area=TRUE)


density <- make.density(region = region,
                        x.space = 0.22,
                        constant = c(20,50))

covs <- list()
covs$size <- list(list(distribution = "poisson", lambda = 25),
                  list(distribution = "poisson", lambda = 15))
covs$sex <- data.frame(level = rep(c("male", "female"),2),
                       prob = c(0.5, 0.5, 0.6, 0.4),
                       strata = c(rep("low",2),rep("high",2)))

# Define the population description (this time using the density to determine
# the population size)
popdesc <- make.population.description(region = region,
                                       density = density,
                                       covariates = covs,
                                       fixed.N = FALSE)

cov.param <- list()
cov.param$size <- c(log(1.02),log(1.005))
cov.param$sex <- data.frame(level = c("male", "female", "male", "female"),
                            param = c(log(1.5), 0, log(1.7), log(1.2)),
                            strata = c("low","low","high","high"))

# define the detecability
detect <- make.detectability(key.function = "hn",
                             scale.param = 0.08,
                             cov.param = cov.param,
                             truncation = mytrunc)

plot(detect, popdesc)


myds <- make.ds.analysis(truncation=mytrunc)

simmake <- make.simulation(reps=10, design=mydesign,
                           population.description=popdesc,
                           detectability=detect,
                           ds.analysis=myds)
seeit <- run.survey(simmake)
plot(seeit, region)
out <- run.simulation(simmake)

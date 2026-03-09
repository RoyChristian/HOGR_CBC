root <- getwd()
################################################################################
#                                 packages
################################################################################
library(cmdstanr)
library(sp)
library(sf)
library(reshape2)
library(tidyverse)
library(scales)
library(bayesplot)
library(bayestestR)
library(posterior)

################################################################################
#                                 functions
################################################################################
dist.func <- function(Coords=Coords){
  dist.list <- lapply(1:nrow(Coords), function(i){
    a <- Coords[i,1] - Coords[,1]
    b <- Coords[i,2] - Coords[,2]
    return(rbind(a,b))
  })
  dist.array  <- array(unlist(dist.list), dim = c(nrow(dist.list[[1]]), ncol(dist.list[[1]]), length(dist.list)))
  dist.perm <- aperm(dist.array, c(2,3,1))
  return(dist.perm)
}

# General inits function for standalone testing
inits.func = function(){
  mu_param = rnorm(8,0,0.025)
  sigma = runif(ncol(mean.ordered),0.01,0.1)
  phi = runif(1,1,3)
  L_rho = t(chol((runif(1,0.10,0.20) ^ outer(1:8, 1:8, function(x, y) abs(x - y)))))
  L_xi = runif(8,0.01,0.05)
  param_epsilon_raw = rnorm(8*ncol(mean.ordered),0,0.01)
  epsilon_N_raw = matrix(rnorm(nrow(mean.ordered)*ncol(mean.ordered),0,0.1),nrow=nrow(mean.ordered))
  return(list(mu_param=mu_param,
              sigma=sigma,
              phi_x=phi,
              phi_y=phi,
              L_rho=L_rho,
              L_xi=L_xi,
              param_epsilon_raw=param_epsilon_raw,
              epsilon_N_raw=epsilon_N_raw))  
}

prob90 = function (x) {
  test <- ci(x, ci = 0.90, method = "ETI")
  out <- c(median(x),test$CI_low, test$CI_high)
  names(out) <- c("y","ymin", "ymax")
  return(out)
}

sig90 = function (x) {
  test <- ci(x, ci = 0.90, method = "ETI")
  out <- between(0,test$CI_low, test$CI_high)
  return(!out)
}

extract_samples <- function(fit_obj) {
  vars <- fit_obj$metadata()$stan_variables
  draws <- posterior::as_draws_rvars(fit_obj$draws())
  
  lapply(vars, \(var_name){  
    posterior::draws_of(draws[[var_name]], with_chains = FALSE)
  }) |> setNames(vars)
}

################################################################################
#                               import covariates
################################################################################
ponds.df <- read.csv(file.path(root,"Data","new_ponds.csv"), header=T) 
breed.df <- read.csv(file.path(root, "Data/breed_condtions.csv"))
winter.df = read.csv(file.path(root, "Data/climate_table/winter_bystrata.csv")) %>% 
  rename(stratum2=Strata,year=Year)

### Keep only 1970 to 2017
ponds.df <- left_join(ponds.df %>% filter(Year %in% seq(1970,2017)),
                      ponds.df %>% mutate(Year=Year+1) %>% rename(Pondslag=Ponds))

################################################################################
#                      import data and calculate log sd
################################################################################
hogr.df = read.csv(file.path(root,"Data","new_extractions","annual_indices_per_analytical_stratum_recoded.csv"))
hogr.df$index_sdlog = sqrt((log(hogr.df$index_mean) - log(hogr.df$index_median))*2)
hogr.df = droplevels(hogr.df[hogr.df$year>=1970,])

### Add winter climate data
hogr.df = left_join(hogr.df, winter.df, by=c("stratum2","year"))

#### Join ponds and summer precipitation
ponds.df = left_join(ponds.df, breed.df)

################################################################################
#                               import shapefiles
################################################################################
CBC.shp <- st_read(file.path(root, "data", "shapefiles", "CBC_Strata_Lambert.shp"))

# Rename the 5th column to "Strata" 
names(CBC.shp)[5] <- "Strata"

################################################################################
#                               check projection
################################################################################
# st_crs(CBC.shp)

################################################################################
#                            create data matrices
################################################################################
median.mat = acast(hogr.df, year~stratum2, value.var = "index_median")
mean.mat = acast(hogr.df, year~stratum2, value.var = "index_mean")
sd.mat = acast(hogr.df, year~stratum2, value.var = "index_sdlog")
area.mat = acast(hogr.df, year~stratum2, value.var = "stratum_area")
ratio.mat = acast(hogr.df, year~stratum2, value.var = "circle_detection_ratio")

### Winter climate matrices
Win_prec.mat = acast(hogr.df, year~stratum2, value.var = "Winter_prec")
Win_mint.mat = acast(hogr.df, year~stratum2, value.var = "Winter_mint")

################################################################################
#                         extract distance information
################################################################################
# Filter to keep only relevant strata found in the data
CBC_kept.shp <- CBC.shp[as.character(CBC.shp$Strata) %in% levels(as.factor(hogr.df$stratum2)), ]
CBC_kept.shp <- CBC_kept.shp %>%
  mutate(across(where(is.factor), fct_drop))

### Extract centroids from polygons
centroids <- st_centroid(CBC_kept.shp)
coords_matrix <- st_coordinates(centroids)
coords.selected <- data.frame(strata = as.character(CBC_kept.shp$Strata),
                              x = coords_matrix[,1],
                              y = coords_matrix[,2],
                              stringsAsFactors = FALSE)

#### Reorder coordinates to match alphabetical matrix columns
coords.selected <- coords.selected[order(coords.selected$strata), ]

#### Estimate distance matrix in meters
distance.mat <- as.matrix(dist(coords.selected[, c("x", "y")]))

##### Create new order based on distance from reference stratum (Row 7: Alaska)
new.order <- order(distance.mat[7, ])

#### Reorder all matrices and coordinates to follow migration path
mean.ordered = mean.mat[,new.order]
sd.ordered = sd.mat[,new.order]
median.ordered = median.mat[,new.order]
area.ordered = area.mat[,new.order]
ratio.ordered = ratio.mat[,new.order]
coords.ordered = coords.selected[new.order,]

# Scale distances for the spatial kernel (km/45)
distance.ordered =  as.matrix(dist(coords.ordered[,2:3]/(1000*45))) 
coords.scaled = coords.ordered[,2:3]/(50000)

#### Reorder winter climate matrices 
Win_prec.ordered <- Win_prec.mat[,new.order]
Win_mint.ordered <- Win_mint.mat[,new.order]

### Center climate values for better convergence
scale_Win_prec.ordered <- apply(Win_prec.ordered/1e7,2,scale, scale=FALSE)
scale_Win_mint.ordered <- apply(Win_mint.ordered,2,scale, scale=FALSE)

################################################################################
#                                 setup Stan data
################################################################################
data.stan <- list(NYear = nrow(mean.ordered), 
                  NRegions = ncol(mean.ordered),
                  NParams = 8,
                  distance_vec = dist.func(coords.scaled),
                  Y = mean.ordered, 
                  log_median_Y1  = log(median.ordered[1,]),
                  sigma_obs = sd.ordered,
                  var_obs = sd.ordered^2,
                  Trend = scale(seq(1,nrow(mean.ordered))/20, scale=F)[,1],
                  PPRprec=scale(ponds.df$ppr_prec/1e7, scale = FALSE)[,1],
                  Borealprec=scale(ponds.df$boreal_prec/1e7, scale = FALSE)[,1],
                  PPRmaxt=scale(ponds.df$ppr_maxt, scale = FALSE)[,1],
                  Borealmaxt =scale(ponds.df$boreal_maxt, scale = FALSE)[,1],
                  Wprec = scale_Win_prec.ordered,
                  Wmint = scale_Win_mint.ordered)

################################################################################
#                                 run Stan model
################################################################################
full.model <- cmdstan_model(stan_file=file.path(root,"stan","sst_climate_linear_allyears.stan"))

fit <- full.model$sample(data = data.stan,
                         seed = 541, 
                         chains = 6, 
                         parallel_chains = 6,
                         iter_warmup = 1000,
                         iter_sampling = 2500, 
                         refresh = 300,
                         adapt_delta = 0.99,
                         max_treedepth = 20,
                         init = function(){
                           mu_param = rnorm(data.stan$NParams, 0, 0.025)
                           sigma = runif(data.stan$NRegions, 0.01, 0.1)
                           phi_x = runif(1, 1, 3)
                           phi_y = runif(1, 1, 3)
                           L_rho = t(chol((runif(1, 0.10, 0.20) ^ outer(1:data.stan$NParams, 1:data.stan$NParams, function(x, y) abs(x - y)))))
                           L_xi = runif(data.stan$NParams, 0.01, 0.05)
                           param_epsilon_raw = rnorm(data.stan$NParams * data.stan$NRegions, 0, 0.01)
                           epsilon_N_raw = matrix(rnorm(data.stan$NYear * data.stan$NRegions, 0, 0.1), nrow=data.stan$NYear)
                           return(list(mu_param=mu_param,
                                       sigma=sigma,
                                       phi_x=phi_x,
                                       phi_y=phi_y,
                                       L_rho=L_rho,
                                       L_xi=L_xi,
                                       param_epsilon_raw=param_epsilon_raw,
                                       epsilon_N_raw=epsilon_N_raw))  
                         })

output = extract_samples(fit)

################################################################################
#                                 save output
################################################################################
save.image(file.path(root, paste("output/climate_linear","final",".RData",sep="")))
fit$save_object(file = file.path(root, paste("output/fit_climate_linear","final","fit.RDS",sep="")))

################################################################################
#                                  diagnose
################################################################################
fit_diagnostic = fit$cmdstan_diagnose()

################################################################################
#                                  traceplots
################################################################################
fit_array = fit$draws()
mcmc_trace(fit_array, pars = "phi_x")
mcmc_trace(fit_array, pars = "phi_y")
mcmc_trace(fit_array, pars = paste("mu_param[", 1:data.stan$NParams, "]", sep=""))
mcmc_trace(fit_array, pars = paste("L_xi[", 1:data.stan$NParams, "]", sep=""))

################################################################################
#                             extract fit summary
################################################################################
fit_summary <- fit$summary()
fit_summary %>% filter(rhat > 1.05) # Modern Stan threshold is often 1.05
fit_summary %>% filter(ess_bulk < 400) # Minimum recommended by Stan team

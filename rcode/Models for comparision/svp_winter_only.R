root <- getwd()
################################################################################
#                                packages
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
#                               functions
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

extract_samples <- function(fit_obj) {
  vars <- fit_obj$metadata()$stan_variables
  draws <- posterior::as_draws_rvars(fit_obj$draws())
  lapply(vars, \(var_name){  
    posterior::draws_of(draws[[var_name]], with_chains = FALSE)
  }) |> setNames(vars)
}

################################################################################
#                                import covariates
################################################################################
ponds.df <- read.csv(file.path(root,"Data","new_ponds.csv"), header=T) 
breed.df <- read.csv(file.path(root, "Data/breed_condtions.csv"))
winter.df = read.csv(file.path(root, "Data/climate_table/winter_bystrata.csv")) %>% 
  rename(stratum2=Strata,year=Year)

### Keep only 1970 to 2017
ponds.df <- left_join(ponds.df %>% filter(Year %in% seq(1970,2017)),
                      ponds.df %>% mutate(Year=Year+1) %>% rename(Pondslag=Ponds))

################################################################################
#                                import data
################################################################################
hogr.df = read.csv(file.path(root,"Data","new_extractions","annual_indices_per_analytical_stratum_recoded.csv"))
hogr.df$index_sdlog = sqrt((log(hogr.df$index_mean) - log(hogr.df$index_median))*2)
hogr.df = droplevels(hogr.df[hogr.df$year>=1970,])

### Add winter climate
hogr.df = left_join(hogr.df, winter.df, by=c("stratum2","year"))
ponds.df = left_join(ponds.df, breed.df)

################################################################################
#                           import shapefiles (sf version)
################################################################################
CBC.shp <- st_read(file.path(root, "data", "shapefiles", "CBC_Strata_Lambert.shp"))
names(CBC.shp)[5] <- "Strata"

################################################################################
#                           create data matrices
################################################################################
median.mat = acast(hogr.df, year~stratum2, value.var = "index_median")
mean.mat   = acast(hogr.df, year~stratum2, value.var = "index_mean")
sd.mat     = acast(hogr.df, year~stratum2, value.var = "index_sdlog")

Win_prec.mat = acast(hogr.df, year~stratum2, value.var = "Winter_prec")
Win_mint.mat = acast(hogr.df, year~stratum2, value.var = "Winter_mint")

################################################################################
#                           extract distance information
################################################################################
CBC_kept.shp <- CBC.shp[as.character(CBC.shp$Strata) %in% levels(as.factor(hogr.df$stratum2)), ]
CBC_kept.shp <- CBC_kept.shp %>% mutate(across(where(is.factor), fct_drop))

### Extract centroids for distance matrix
centroids <- st_centroid(CBC_kept.shp)
coords_matrix <- st_coordinates(centroids)

coords.selected <- data.frame(strata = as.character(CBC_kept.shp$Strata),
                              x = coords_matrix[,1],
                              y = coords_matrix[,2],
                              stringsAsFactors = FALSE) %>% 
  arrange(strata)

distance.mat <- as.matrix(dist(coords.selected[, c("x", "y")]))
new.order    <- order(distance.mat[7, ]) # Distance from Alaska

#### Reorder everything to follow migration path
mean.ordered   = mean.mat[,new.order]
sd.ordered     = sd.mat[,new.order]
median.ordered = median.mat[,new.order]
coords.ordered = coords.selected[new.order,]

# Scaling for Stan
coords.scaled = coords.ordered[,2:3]/(50000)

#### Reorder and center winter climate 
Win_prec.ordered <- Win_prec.mat[,new.order]
Win_mint.ordered <- Win_mint.mat[,new.order]
scale_Win_prec.ordered <- apply(Win_prec.ordered/1e7, 2, scale, scale=FALSE)
scale_Win_mint.ordered <- apply(Win_mint.ordered, 2, scale, scale=FALSE)

################################################################################
#                         setup Stan data
################################################################################
Kept_out = 2 # Number of years held back for predictive check

data.stan <- list(NYear = nrow(mean.ordered), 
                  MAxYear = nrow(mean.ordered) - Kept_out, 
                  NRegions = ncol(mean.ordered),
                  NParams = 3,
                  distance_vec = dist.func(coords.scaled),
                  Y = mean.ordered, 
                  log_median_Y1  = log(median.ordered[1,]),
                  sigma_obs = sd.ordered,
                  var_obs = sd.ordered^2,
                  Wprec = scale_Win_prec.ordered,
                  Wmint = scale_Win_mint.ordered)

################################################################################
#                         run Stan model
################################################################################
full.model <- cmdstan_model(stan_file=file.path(root,"stan","sst_winter_only.stan"))

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
                           L_xi  = runif(data.stan$NParams, 0.01, 0.05)
                           param_epsilon_raw = rnorm(data.stan$NParams * data.stan$NRegions, 0, 0.01)
                           epsilon_N_raw = matrix(rnorm(data.stan$MAxYear * data.stan$NRegions, 0, 0.1), nrow=data.stan$MAxYear)
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
#                      save output
################################################################################
save.image(file.path(root, paste0("output/winter_only_N", (nrow(mean.ordered)-Kept_out), ".RData")))
fit$save_object(file = file.path(root, paste0("output/fit_winter_only_N", (nrow(mean.ordered)-Kept_out), "fit.RDS")))

################################################################################
#                       diagnose and summaries
################################################################################
fit_diagnostic = fit$cmdstan_diagnose()
fit_summary    = fit$summary()
fit_summary %>% filter(rhat > 1.05) 
fit_summary %>% filter(ess_bulk < 400)

################################################################################
#                      extract lppd for observed values
################################################################################
log_mean_exp <- function(x) {
  max_x <- max(x)
  return(max_x + log(mean(exp(x - max_x))))
}

# Calculate lppd per year and strata
lppd_matrix <- apply(output$log_lik, c(2, 3), log_mean_exp)
total_lppd  <- sum(lppd_matrix)

################################################################################
#                      extract lppd for predicted t + 1 
################################################################################
lppd_pred <- apply(output$lik_pred, 2, log_mean_exp)

pred.df <- data.frame(Strata = colnames(mean.ordered),
                      Year   = (nrow(mean.ordered) - Kept_out) + 1, 
                      lppd_pred = lppd_pred)

write.csv(pred.df, file = file.path(root, paste0("output/fit_winter_only_N", (nrow(mean.ordered)-Kept_out), "lppd.csv")), row.names = FALSE)
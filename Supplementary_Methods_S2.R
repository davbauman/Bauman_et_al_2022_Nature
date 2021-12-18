###########################################################################
###########################################################################
### PART I: Calculation of latent logit-survivorship per sp, plot, year ###
###########################################################################
###########################################################################

# Survival model outputs expressed in terms of mortality (1 - survival).

##########################################
### SURVIVAL WRAPPER - LOADS STAN MCMC ###
##########################################

rm(list = ls())
set.seed(1)

#################
### LIBRARIES ###
#################

library(rstan)
library(parallel)
library(doParallel)
library(tidyverse)
library(tidybayes)
library(brms)
library(data.table)
library(viridis)

rstan_options(auto_write = TRUE)

# Themes for figures:
# -------------------
my_theme_big <- theme_classic() +
  theme(panel.border = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        plot.subtitle = element_text(size = 16, colour = "black")) 

my_theme_inset <- theme_classic() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 24, colour = "black"),
        axis.title.x = element_blank())

## DATA:
########
# Import data

load('DATA/datafinal_all_sp.RData')   # All species

# Data preparation:
# *****************

# We use the datasets with mean trait values per species per plot:

data <- datafinal_all_sp %>%
  filter(is.na(dbh_0) == FALSE)

data2 <- data %>% 
  # Transform column 'dead' (dead = 1, surv = 0) into 'surv' (dead = 0, surv = 1):
  rowwise() %>% 
  mutate(dead = ifelse(dead == 0, 1, 0)) %>% 
  ungroup %>% 
  rename(surv = dead) %>% 
  # Rename species code 'code' into 'Sp' and dbh_0 into dbh:
  rename(Sp = code, dbh = dbh_0) %>% 
  # Express DBH in mm:
  mutate(dbh = dbh * 10)

# Time between two censuses, expressed as fraction of a year:
data2 <- data2 %>% 
  dplyr::select(plot:nbdays, everything()) %>% 
  # Create two index variables corresponding to 'year_0' and 'year_1' of the intervals:
  mutate(yrstart = year_0 - min(data$year_0) + 1,
         yrend = year_1 - min(data$year_0) + 1)

# Create Stan survival model:
# ***************************
write("
data {
  int<lower=0> n;                   // no. of observations
  int<lower=0> n_plot;              // no. of plots
  int<lower=0> T;                   // no. of Years
  real dbh[n];                      // scaled diameter in mm
  int surv[n];                      // survival of each tree, 1 for alive, 0 for dead
  real thresh;                      // dbh threshold separating first from second curve
  int cstart[n];                    // index for census first year
  int cend[n];                      // index for census last year
  int plot[n];                      //  plot identifier
  int Tind[T];                      //  year identifier
}

parameters {
  // Main model parameters
  real K_mu;

  // Non-centered priors (see transformed parameters block):
  vector[T] K_T_z;
  vector[n_plot] K_P_z;
  vector[n] eps_z;

  real<lower=0, upper=2> r1;
  real<lower=-2, upper=0> r2;
  real p2;    //<lower=max(dbh) * 1.2, upper= 2 * (max(dbh))>
  real p1;    //<lower=-1*thresh, upper=thresh*0.75>

  // varying effect parameters
  real<lower=0.25> sigma_K_mu;
  real<lower=0> sigma_K_T;
  real<lower=0> sigma_K_P;
}

transformed parameters {
  vector[T] K_T;
  vector[n_plot] K_P;
  vector[n] eps;

  // Non-centered priors (smuggle the hyperparameters out of the adaptive priors):
  K_T = K_T_z * sigma_K_T;
  K_P = K_P_z * sigma_K_P;
  eps = eps_z * sigma_K_mu;
}

model {
  real theta;
  real p;
  real K;
  real K_tmp;

  // Non-centered priors, to avoid having hyperparameters in the adaptive priors (more efficient MCMC sampling):
  K_T_z ~ normal(0, 1);
  K_P_z ~ normal(0, 1);
  eps_z ~ normal(0, 1);
  K_mu ~ normal(2, 1);

  // Priors on hyperparameters for random effects are here
  sigma_K_T ~ gamma(3, 4);
  sigma_K_P ~ gamma(3, 4);
  sigma_K_mu ~ gamma(3, 4);
  p2 ~ normal(max(dbh), 1);     // 1, here, represents one SD of the species' dbh (DBH is standardised)
  p1 ~ normal(min(dbh), 1);

  for(i in 1:n)
  {
    if(dbh[i] < thresh)
    {
      for(t in cstart[i]:cend[i])
      {
        K_tmp = K_mu + K_P[plot[i]] + K_T[Tind[t]] + eps[i];
        K = inv_logit(K_tmp);
        p = (K) / (1 + exp(-r1 * (dbh[i] - p1)));
        theta = p;
      }
    }
    else
    {
      for(t in cstart[i]:cend[i])
      {
        K_tmp = K_mu + K_P[plot[i]] + K_T[Tind[t]] + eps[i];
        K = inv_logit(K_tmp);
        p = (K) / (1 + exp(-r2 * (dbh[i] - p2)));
        theta = p;
      }
    }
    surv[i] ~ bernoulli(theta);
  }

}",
file = "CODE/MORTALITY_Year_Plot_DB_reparam.stan")

# Function to run the Stan model on one species, save the raw model output and 
# summarise the necessary survival parameter posteriors for further use:
# *****************************************************************************
# Beginning of run.stan
run.stan <- function(sp.df, surv.model = surv.model, n_iter = 2000, n_chains = 3, 
                     n_cores = 3, seed = 2020) {
  
  # 'Small' and 'big' individuals are considered based on the mean DBH of 
  # the species in the data (threshold):
  thresh <- mean(sp.df$dbh, na.rm = T) # 0, if dbh was centered
  tot.years <- max(sp.df$yrend, na.rm = T) - min(sp.df$yrstart, na.rm = T) + 1
  
  # Transform the data into a list format, as needed for the Stan model:
  # --------------------------------------------------------------------
  # compose_data creates a list of vectors and real values, where numeric variables are left unchanged
  # and factors and character strings are transformed into index variables, as needed in Stan. 
  # For the latter, an object called 'n_variableName' is created, with the number of indices of 'variableName'.
  surv.data <- compose_data(sp.df)
  
  surv.data <- purrr::list_modify(.x = surv.data,
                                  T = tot.years,
                                  thresh = thresh,
                                  cstart = sp.df$yrstart - min(sp.df$yrstart) + 1,
                                  cend = sp.df$yrend - min(sp.df$yrstart) + 1,
                                  Tind = seq(tot.years))
  
  # Run the survival model:
  surv.fit <- stan(fit = surv.model, data = surv.data, 
                   iter = n_iter, chains = n_chains, cores = n_cores, 
                   verbose = TRUE, init = 0, 
                   control = list(adapt_delta = 0.98, max_treedepth = 15),
                   seed = seed)
  
  # Save the full raw model output:
  sp_name <- sp.df$Sp[1]
  save(surv.fit, file = sprintf("OUTPUT/surv_fit_plot_year_%s.RData", sp_name))
  
  # Extract a sample of draws from the parameter posteriors (excluding warmups):
  # ----------------------------------------------------------------------------
  # tidybayes::spread_draws() generates for each level of K_P for example (i.e. each plot), 
  # a number of rows equal to the number of iterations (after removing the warmups), 
  # i.e. n_iter/2 by default, multiplied by the n_chains. 
  # tidybayes::recover_types() is a very handy function that back-transforms the index variables
  # used in Stan (i.e. plots and years) into their original values.
  
  extract_KP <- tidybayes::spread_draws(model = recover_types(surv.fit, sp.df),
                                        K_P[plot],  seed = seed) %>%
    dplyr::summarise(K_P = median(K_P))
  
  extract_KT <- tidybayes::spread_draws(model = surv.fit,
                                        K_T[T], seed = seed) %>%
    dplyr::summarise(K_T = median(K_T))
  
  # Convert 'T' back in corresponding years:
  extract_KT$T <- extract_KT$T + min(sp.df$year_0) - 1 
  
  extract_Kmu <- tidybayes::spread_draws(model = surv.fit,
                                         K_mu,  seed = seed) %>% 
    dplyr::summarise(K_mu = median(K_mu))
  
  extract_p1 <- tidybayes::spread_draws(model = surv.fit, 
                                        p1,  seed = seed) %>% 
    dplyr::summarise(p1 = median(p1))
  
  extract_r1 <- tidybayes::spread_draws(model = surv.fit,
                                        r1,  seed = seed) %>% 
    dplyr::summarise(r1 = median(r1))
  
  extract_p2 <- tidybayes::spread_draws(model = surv.fit,
                                        p2,  seed = seed) %>% 
    dplyr::summarise(p2 = median(p2))
  
  extract_r2 <- tidybayes::spread_draws(model = surv.fit,
                                        r2, seed = seed) %>% 
    dplyr::summarise(r2 = median(r2))
  
  extract_eps <- tidybayes::spread_draws(model = surv.fit,
                                         eps[n], seed = seed) %>% 
    dplyr::summarise(eps = median(eps))  
  
  sp.df <- data.table(sp.df)
  cyear <- sp.df[, list(plot = plot,
                        surv=surv,
                        dbh=dbh,
                        year= seq(year_0, year_1)),
                 by=1:nrow(sp.df) ]
  
  names(extract_KT) <- c('year','K_T')
  names(extract_eps) <- c('nrow','eps')
  
  cyear <- cyear %>% left_join(extract_KT, by='year')
  cyear <- cyear %>% left_join(extract_KP, by='plot')
  cyear <- cyear %>% left_join(extract_eps, by='nrow')
  
  Kmu <-  extract_Kmu$K_mu
  
  df_Ktmp <- data.table(species = rep(sp_name, sum(surv.data$cend - surv.data$cstart + 1)), 
                        obs = cyear$nrow, 
                        plot = cyear$plot, 
                        year = cyear$year, 
                        dbh = cyear$dbh, 
                        surv = cyear$surv,
                        K_tmp = NA,
                        K_mu = Kmu,
                        K_P = cyear$K_P,
                        K_T = cyear$K_T,
                        eps = cyear$eps, 
                        p1 = extract_p1$p1,
                        r1 = extract_r1$r1,
                        p2 = extract_p2$p2,
                        r2 = extract_r2$r2, 
                        theta = NA,
                        n = NA,
                        seed = NA)
  
  df_Ktmp$K_tmp <- df_Ktmp$K_mu + df_Ktmp$K_P + df_Ktmp$K_T + df_Ktmp$eps
  
  save(df_Ktmp, file = sprintf("OUTPUT/df_Ktmp_plot_year_%s.RData", sp_name))
  
  return(df_Ktmp)
}
# End of run.stan()
# **************************************************************************************** 

# Select species with enough observations:
# ----------------------------------------
data_summary <- data2 %>% 
  group_by(Sp) %>% 
  summarise(nb_obs = n(),
            nb_deaths = sum(surv == 0),
            nb_stems = n_distinct(stem),
            nb_plots = n_distinct(plot)) %>% 
  ungroup

big.sp <- unique(filter(data_summary, nb_obs >= 400)$Sp)
data.big <- subset(data2, Sp %in% big.sp)

# Create a list to save the df_Ktmp datasets generated for each species:
df_Ktmp_ls <- vector("list", length(unique(data.big$Sp)))
names(df_Ktmp_ls) <- unique(data.big$Sp)

# Size-dependent survival model parameters:
# -----------------------------------------
sp.df_list <- lapply(unique(data.big$Sp), function (x) subset(data.big, Sp == x))
names(sp.df_list) <- unique(data.big$Sp)

# Run the following chunk once only to build and save the stan model.
# Build the survival model on one species, instead of building it every time in run.stan():
# -----------------------------------------------------------------------------------------
tmp <- sp.df_list[[1]]
tmp <- dplyr::select(.data = tmp,
                     Sp, plot, year_0, year_1, yrstart, yrend, dbh, surv)
tmp$dbh <- as.numeric(scale(tmp$dbh))
thresh <- mean(tmp$dbh) 
tot.years <- max(tmp$yrend) - min(tmp$yrstart) + 1
surv.data <- compose_data(tmp)

surv.data <- list_modify(.x = surv.data,
                         T = tot.years,
                         thresh = thresh,
                         cstart = tmp$yrstart - min(tmp$yrstart) + 1,
                         cend = tmp$yrend - min(tmp$yrstart) + 1,
                         Tind = seq(tot.years))

surv.model <- stan(file = "MORTALITY_Year_Plot_DB_reparam.stan",
                   data = surv.data, chains = 0, save_dso = TRUE, verbose = TRUE)
save(surv.model, file = "surv_model_Year_Plot_DB_reparam.RData")

# load("CODE/surv_model_Year_Plot_DB_reparam.RData")

# Loop through the species in parallel:
# -------------------------------------
n_iter <- 3000
n_chains <- 3    # defines both the number of chains and cores used
seed <- 2020
# n_samples <- 100  # nb of samples from the posterior draws for each observation (to generate K_tmp)

doParallel::registerDoParallel(cores = 9) # 8 on clusters, 4 on laptop

df_Ktmp_ls <- foreach (a = sp.df_list, 
                       .packages = c("rstan", "tidyverse", 
                                     "tidybayes", "data.table")) %dopar%  {
                                       a <- dplyr::select(.data = a,
                                                          Sp, plot, year_0, year_1, yrstart, yrend, dbh, surv)
                                       a$dbh <- as.numeric(scale(a$dbh))
                                       
                                       df_Ktmp <- run.stan(sp.df = a, surv.model = surv.model,
                                                           n_iter = n_iter, n_chains = n_chains, n_cores = n_chains, 
                                                           seed = seed)
                                       return(df_Ktmp)
                                     }
foreach::registerDoSEQ()

names(df_Ktmp_ls) <- big.sp

save(df_Ktmp_ls, file = "OUTPUT/df_Ktmp_list_plot_year.RData")
# load('OUTPUT/df_Ktmp_list_plot_year.RData')

# CREATE THE DATASET FOR THE MODELS OF K_lat (here called K_tmp):
# ***********************************************
# ***********************************************

df_all_sp <- bind_rows(df_Ktmp_ls)

df_all_sp <- df_all_sp %>% 
  mutate(K_P_abs = K_mu + K_P,
         K_T_abs = K_mu + K_T) %>% 
  select(sp = species, obs, plot, year, dbh, surv, K_tmp, 
         K_mu:r2, K_P_abs, K_T_abs, theta)

View(df_all_sp)

# load("DATA/df_all_sp.RData")

# Add a column to differentiate species with and without trait data available:
# ----------------------------------------------------------------------------
load('DATA/datafinal_trait_u.RData')    # Species with trait data
sp_with_traits <- unique(datafinal_trait_u$code)

df_all_sp$trait_data <- sapply(df_all_sp$sp, 
                               function(x) ifelse(x %in% sp_with_traits, 1, 0))

# Add columns of nb_obs, nb_stems, nb_plots, nb_deaths to ease the filtering of species:
# --------------------------------------------------------------------------------------
df_all_sp <- df_all_sp %>% 
  left_join(data_summary %>% 
              rename(sp = Sp), by = "sp")

# Species with >=400obs (regardless of having trait data measured):
# -----------------------------------------------------------------
save(df_all_sp, file = "DATA/df_all_sp.RData")

# load("DATA/df_all_sp.RData")

# Add climate variables to 'df_all_sp':
# *************************************

# Climate covariates were built in a separate R script #

load("DATA/clim_final2.RData")

df_all_sp_covar <- df_all_sp %>% 
  left_join(clim_final, by = "plot")

# Add the trait data, when available (species mean traits):
# *********************************************************
df_all_sp_covar <- df_all_sp %>% 
  left_join(datafinal_trait_u %>% 
              select(sp = code, dbh_max:LMA, -c(Vcmax, Jmax)) %>% 
              rename(Vcmax = Vcmax_400op, Jmax = Jmax_1200op) %>% 
              group_by(sp) %>% 
              summarise_all(~mean(., na.rm = T)) %>% 
              ungroup,
            by = "sp")

save(df_all_sp_covar, file = "DATA/df_all_sp_covar.RData")
# load("DATA/df_all_sp_covar.RData")

## Visualisation of functional size-dependent survival curve per species:
# #######################################################################

# Create functions:
# -----------------
predict.surv <- function(x, params, time = 1) {
  K <- params[1]
  p <- params[2]
  r <- params[3]
  # Linear predictor parameters
  pred.y <-  K / (1 + exp(-(r * ((x - p) )))) ^ (time)
  return(pred.y)
}

get.curve.mat <- function(samps, subsample = TRUE, s.no = 500, thresh = 0,
                          mean.dbh, sd.dbh,
                          min.size = -2, max.size = 2, suppl.dbh = 1) {
  if (subsample == TRUE) {
    samp.id <- sample(seq(length(samps$p1)), s.no)
    K  <- boot::inv.logit(samps$K_mu[samp.id] + 
                            samps$K_P[samp.id] + 
                            samps$K_T[samp.id] + 
                            samps$eps[samp.id])
    p1 <- samps$p1[samp.id]
    r1 <- samps$r1[samp.id]
    p2 <- samps$p2[samp.id]
    r2 <- samps$r2[samp.id]
  } else {
    K  <- boot::inv.logit(samps$K_mu + samps$K_P + samps$K_T + samps$eps)
    p1 <- samps$p1
    r1 <- samps$r1
    p2 <- samps$p2
    r2 <- samps$r2
  }
  # combinations of 3 param.
  params.lo <- split(cbind(K, p1, r1), seq(length(K))) 
  # combinations of 3 param.
  params.hi <- split(cbind(K, p2, r2), seq(length(K))) 
  # dbh <- seq(from = min.size, to = max.size + suppl.dbh, length = 1000)
  # for same x-axis range in different sp.:
  dbh <- seq(from = min.size, to = (1600 - mean.dbh) / sd.dbh, length = 1000) 
  res.mat <- matrix(0, ncol = length(dbh), nrow = length(K))
  for(i in 1:nrow(res.mat)) {
    lo.prob <- predict.surv(x = dbh[dbh < thresh], params = params.lo[[i]])
    hi.prob <- predict.surv(x = dbh[dbh >= thresh], params = params.hi[[i]])
    res.mat[i, ] <- c(lo.prob, hi.prob)
  }
  return(res.mat)
}

make.surv.fig1 <- function(res.mat, min.dbh = -2, max.dbh = 2, suppl.dbh = 1,
                           mean.dbh, sd.dbh, species = "", 
                           probs= c(0.05, 0.25, 0.5, 0.75, 0.95)) {
  quants <- apply(res.mat, 2, quantile, probs = probs)
  rownames(quants) <- paste("q", str_replace(rownames(quants), "%", ""), sep = "")
  # dbh <- seq(from = min.dbh, to = max.dbh + suppl.dbh, length = 1000)
  # for same x-axis range in different sp.:
  dbh <- seq(from = min.dbh, to = (1600 - mean.dbh) / sd.dbh, length = 1000) 
  
  quants_t <- as.data.frame(cbind(t(quants), dbh))
  # mean.dbh and sd.dbh are used to back-transform the dbh to original scale:
  quants_t$dbh <- (quants_t$dbh * sd.dbh) + mean.dbh
  max <- (max.dbh * sd.dbh) + mean.dbh
  
  ggplot(data = quants_t, aes(x = dbh, y = q50)) +
    geom_line(colour = "red3") +
    geom_ribbon(aes(ymin = q5, ymax = q95), alpha = 0.3) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.5) +
    geom_vline(xintercept = max, linetype = 2) +
    # geom_text(data = tibble(x = max, y = 0), 
    #           aes(x, y), label = "max") +
    lims(y = c(0, 1)) +
    labs(x = "DBH (mm)",
         y = "Probability of survival",
         subtitle = sp_name,
         caption = "Median, 50%- and 90%-credibility intervals") +
    theme_classic()
}

# Load output from the stan model:
# --------------------------------
load("OUTPUT/surv_fit_plot_year_AcaCel.RData") 
sp_name <- "AcaCel"
load("OUTPUT/surv_fit_plot_year_AcrLae.RData") 
sp_name <- "AcrLae"
load("OUTPUT/surv_fit_plot_year_CarSub.RData") 
sp_name <- "CarSub"
load("OUTPUT/surv_fit_plot_year_FliBou.RData") 
sp_name <- "FliBou"

# Extract posterior draws:
# ------------------------
samps <- rstan::extract(surv.fit)

# Draw 's.no' combinations of K, p, r values from the model posterior,
# and generate predicted theta for dbh ranging from the min to the max dbh of sp.:

min_dbh <- min(data.big$dbh[which(data.big$Sp == sp_name)])
max_dbh <- max(data.big$dbh[which(data.big$Sp == sp_name)])

# To back-transform dbh to original scale, for the figure:
mean_dbh <- mean(data.big$dbh[which(data.big$Sp == sp_name)])
sd_dbh <- sd(data.big$dbh[which(data.big$Sp == sp_name)])

res.mat <- get.curve.mat(samps, s.no = 500, 
                         thresh = 0, 
                         min.size = (min_dbh - mean_dbh) / sd_dbh,
                         max.size = ((max_dbh - mean_dbh) / sd_dbh),
                         mean.dbh = mean_dbh,
                         sd.dbh = sd_dbh,
                         suppl.dbh = 1)

# Figure of theta across all plots and years (median, 50%-CI and 90%-CI):
make.surv.fig1(res.mat,
               min.dbh = (min_dbh - mean_dbh) / sd_dbh,
               max.dbh = ((max_dbh - mean_dbh) / sd_dbh),
               suppl.dbh = 1,
               mean.dbh = mean_dbh,
               sd.dbh = sd_dbh,
               species = sp_name)

# General features of the dataset:
# ********************************
tmp <- data.big %>% 
  select(plot, stem, sp = Sp, family, genus, taxon, year_0, year_1, 
         surv, dbh, elevation) %>% 
  mutate(start_end = paste(year_0, year_1, sep = "_"))

table_1 <- tmp %>% 
  group_by(taxon, family, sp) %>% 
  summarise(nb_stems = n_distinct(stem),
            nb_obs = n(),
            nb_plots = n_distinct(plot),
            year_0 = min(year_0),
            year_last = max(year_1),
            nb_deaths = length(which(surv == 0)),
            dbh_u = mean(dbh, na.rm = T),
            dbh_sd = sd(dbh, na.rm = T),
            dbh_min = min(dbh, na.rm = T),
            dbh_max = max(dbh, na.rm = T)) %>% 
  ungroup %>% 
  mutate(nb_years = year_last - year_0)

table_2 <- tmp %>% 
  group_by(sp, stem) %>% 
  summarise(nb_censuses = n_distinct(start_end) + 1,
            time_0 = min(year_0),
            time_last = max(year_1),
            time = time_last - time_0) %>% 
  group_by(sp) %>% 
  summarise(nb_cen_u = mean(nb_censuses),
            nb_cen_sd = sd(nb_censuses),
            nb_cen_min = min(nb_censuses),
            nb_cen_max = max(nb_censuses),
            time_u = mean(time),
            time_sd = sd(time),
            time_min = min(time),
            time_max = max(time)) %>% 
  ungroup

load("DATA/datafinal_trait_u.RData")
trait_sp <- unique(datafinal_trait_u$code)

table_final <- table_1 %>% 
  left_join(table_2) %>% 
  select(taxon:sp, nb_stems, nb_obs, nb_obs, nb_deaths, nb_plots, year_0, year_last, nb_years,
         everything()) %>% 
  mutate(trait_data = ifelse(sp %in% trait_sp, 1, 0))

write.table(table_final, file = "OUTPUT/Table_species_features.txt", sep = "\t")

# Illustration figure of the functional curve (Supplementary Figure) #
# ****************************************************************** #

predict.surv <- function(x, params, time = 1) {
  K <- params[1]
  p <- params[2]
  r <- params[3]
  # Linear predictor parameters
  pred.y <-  K / (1 + exp(-(r * ((x - p) )))) ^ (time)
  return(pred.y)
}

make.pred.sim <- function(K = 0.99, p1 = 20, r1 = 0.05, p2 = 1800, r2 = -0.02, 
                          thresh = 100, dbh.min = 10, dbh.max = 2000) {
  params.lo <- c(K, p1, r1)
  params.hi <- c(K, p2, r2)
  dbh <- seq(from = dbh.min, to = dbh.max, length = 1000)
  lo.prob <- predict.surv(x = dbh[dbh < thresh], params = params.lo)
  hi.prob <- predict.surv(x = dbh[dbh >= thresh], params = params.hi)
  pred <- c(lo.prob, hi.prob)
  df_fig <- as.data.frame(cbind(as.matrix(pred), dbh))
  colnames(df_fig) <- c("pred", "dbh")
  return(df_fig)
}

# Figure with one plot and one year, and one K_lat for this plot and year:
# ------------------------------------------------------------------------
# ------------------------------------------------------------------------

# Plot 1: 
K <- 0.98
p1 <- 20
r1 <- 0.05
p2 <- 1800
r2 <- -0.02
thresh <- 270

df_fig1 <- make.pred.sim(K, p1, r1, p2, r2, thresh)
df_fig1 <- cbind(df_fig1, K = rep(K, 1000))

# Year 1: 
K <- 0.86
df_fig2 <- make.pred.sim(K, p1, r1, p2, r2, thresh)
df_fig2 <- cbind(df_fig2, K = rep(K, 1000))

# Average (inv.logit(K_mu)): 
# Doesn't need be the mean(plots and years above) (supposes unrepresented groups)
K <- 0.95  
df_fig_u <- make.pred.sim(K, p1, r1, p2, r2, thresh)

# K_lat: # mean of the K chosen above
(K <- mean(c(0.95, 0.98, 0.86)))
boot::inv.logit(boot::logit(0.95) + boot::logit(0.03) - boot::logit(0.09))
df_fig_Klat <- make.pred.sim(K, p1, r1, p2, r2, thresh)
df_fig_Klat <- cbind(df_fig_Klat, K = rep(K, 1000))

df_fig_groups <- cbind(rbind(df_fig1, 
                             df_fig2,
                             df_fig_Klat), 
                       group = rep(c("plot1", "year1", "Klat"), each = 1000),
                       type = c(rep("plot", 1000), rep("year", 1000), rep("Klat", 1000)))

# Main figure:
(fig <- ggplot(data = df_fig_groups, aes(x = dbh, y = pred)) +
    geom_line(aes(lty = type, group = group, colour = K), size = 0.7) +
    scale_linetype_manual(values = c(6, 2, 3)) +
    scale_colour_gradient(expression(paste(italic(K))), low = "darkred", high = "blue",
                          limits = c(0.86, 0.98), breaks = seq(0.86, 0.98, 0.04)) +
    guides(linetype = F) +
    geom_line(data = df_fig_u, aes(dbh, pred), size = 1) +
    geom_vline(xintercept = thresh, colour = "grey60") +
    lims(y = c(0, 1),
         x = c(0, 2000)) +
    labs(x = "DBH (mm)",
         y = expression(paste(Survival~probability~per~year~(italic(theta))))) +
    theme_classic() +
    my_theme_big)

## Table of survival parameters per species:
# ******************************************

load("DATA/df_all_sp.RData")

surv_params_summary <- df_all_sp %>% 
  group_by(sp) %>% 
  summarise_at(vars(K_mu, r1, p1, r2, p2), ~mean(.)) %>% 
  ungroup %>% 
  left_join(datafinal_all_sp %>% 
              dplyr::select(sp = code, taxon, family) %>% 
              group_by(sp, taxon, family) %>% 
              summarise()) %>% 
  mutate(K_mu_inv = boot::inv.logit(K_mu),
         Mortality_risk = 1 - K_mu_inv) %>% 
  dplyr::select(taxon, family, sp, K_mu, K_mu_inv, Mortality_risk, everything())

write_tsv(surv_params_summary, "OUTPUT/surv_params_summary.txt")

#################################################################################
#################################################################################
### PART II: Analyses of temporal trend in K_lat and the underlying processes ###
#################################################################################
#################################################################################


## II.1. Model M1: Unconditional model ##
#########################################
#########################################

# Packages:
# ---------
library(tidyverse)
library(brms)

load("DATA/df_all_sp.RData")

# Nb real obs. and stems:
tmp <- data.big 
nrow(tmp)
length(unique(tmp$stem))

## Model:
# *******
# *******

# Model parameters:
iterations <- 6000
chains <- 4
seed <- 43

fit_brms <- brm(formula = K_tmp ~ 1 + (1 | plot) + (1 | sp) + (1 | year),
                data = df_all_sp,
                family = gaussian(),
                prior = c(prior(normal(4.4, 1), class = "Intercept"),
                          prior(exponential(1), class = "sigma"),
                          prior(exponential(1), class = "sd")),
                iter = iterations,
                warmup = 1000,
                chains = chains,
                cores = chains,
                seed = seed,
                control = list(adapt_delta = 0.99, max_treedepth = 15))

# save(fit_brms, file = "OUTPUT/fit_brms_80sp_unconditional_withBEK01_KtmpNotSTD.RData")

## Analysis of model output:
# **************************
# **************************

# Packages:
# ---------
library(tidyverse)
library(tidybayes)
library(viridis)

# Model output:
# -------------
load("OUTPUT/fit_brms_80sp_unconditional_withBEK01_KtmpNotSTD.RData")

# Summary of the model:
print(summary(fit_brms, prob = 0.9), digits = 4)

# Year-level mortality across species and plots:
# -------------------------------------------------
year_coeff <- fit_brms %>% 
  spread_draws(b_Intercept, r_year[year, ]) %>% 
  mutate(year_Intercept = -(b_Intercept + r_year),
         b_Intercept = -b_Intercept) %>% 
  dplyr::select(year, b_Intercept, year_Intercept) %>% 
  # Inv-logit mortality rate per year:
  mutate(year_Intercept_inv = boot::inv.logit(year_Intercept),
         b_Intercept_inv = boot::inv.logit(b_Intercept))

summ_year <- fit_brms %>% 
  spread_draws(b_Intercept, r_year[year, ]) %>% 
  mutate(year_Intercept = -(b_Intercept + r_year),
         b_Intercept = -b_Intercept) %>% 
  dplyr::select(year, b_Intercept, year_Intercept) %>% 
  mutate(year_Intercept_inv = boot::inv.logit(year_Intercept),
         b_Intercept_inv = boot::inv.logit(b_Intercept)) %>% 
  point_interval(.point = median, .interval = qi, .width = c(.9)) %>% 
  ungroup

# Figure:
summ_year %>% 
  ggplot(aes(year, year_Intercept_inv)) +
  geom_pointinterval(aes(ymin = year_Intercept_inv.lower, ymax = year_Intercept_inv.upper),
                     size = 0.7, alpha = 0.7) +
  labs(subtitle = "Tree mortality risk per year across sites and species",
       x = "Years",
       y = expression(paste(Mortality~risk~(yr**-1)))) +
  theme_bw() 

## II.2. Change-point analysis:
###############################
###############################

library(chngpt)

df <- df_all_sp %>% 
  mutate(year = year - 1971)

fit <- chngptm(formula.1 = K_tmp ~ 1, formula.2 = ~ year, 
               data = df,
               type = "segmented", family = "gaussian",
               est.method = "fastgrid", var.type = "bootstrap", save.boot = TRUE)
summary(fit)


save(fit, file = "OUTPUT/Changept_fit.RData")
load("OUTPUT/Changept_fit.RData")


## II.4.1. Model M3 (without climate covariates) - Not in manuscript
################################################
################################################

# Packages:
# ---------
library(tidyverse)
library(brms)

load("DATA/df_all_sp_covar.RData")

df_all_sp_covar2 <- df_all_sp_covar %>%
  filter(! plot %in% c("BEK01", "CBAY"), year >= 1984) %>%
  mutate_at(vars(K_tmp, year), ~scale(.))

# # If only for plots without cyclone intervals:
# load("DATA/df_all_sp_NoCyclCensus.RData")
# df_all_sp_2 <- df_all_sp %>% 
#   filter(year >= 1984) %>% 
#   mutate_at(vars(K_tmp, year), ~scale(.))
# 
# # If only for plots with no cyclone damage in 50 years:
# load("DATA/df_all_sp_NoCyclPlot.RData")
# df_all_sp_2 <- df_all_sp %>% 
#   filter(year >= 1984) %>% 
#   mutate_at(vars(K_tmp, year), ~scale(.))
# 
# # Nb real obs. and stems:
# tmp <- data.big %>%
#   filter(! plot %in% c("BEK01", "CBAY"), year_0 >= 1984)
# nrow(tmp)
# length(unique(tmp$stem))

## Model:
# *******
# *******

# Model parameters:
iterations <- 6000
chains <- 4
seed <- 42

fit_brms <- brm(formula = K_tmp ~ 1 + year + (1 + year | plot) + (1 | sp),
                data = df_all_sp_covar2,
                family = gaussian(),
                prior = c(prior(normal(0, 1), class = "Intercept"),
                          prior(normal(0, 1), class = "b"),
                          prior(exponential(1), class = "sigma"),
                          prior(exponential(1), class = "sd"),
                          prior(lkj(2), class = "cor")),
                iter = iterations,
                warmup = 1000,
                chains = chains,
                cores = chains,
                seed = seed,
                control = list(adapt_delta = 0.99, max_treedepth = 15))

save(fit_brms, file = "OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot.RData")
# save(fit_brms, file = "OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot_NoCyclCensus.RData")
# save(fit_brms, file = "OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot_NoCyclPlot.RData")

## Analysis of model output:
# **************************
# **************************

# Packages:
# ---------
library(tidyverse)
library(tidybayes)
library(viridis)
library(ggridges)

# Model output:
# -------------
load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot.RData")
# load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot_NoCyclCensus.RData")
# load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot_NoCyclPlot.RData")

# Summary of the model:
print(summary(fit_brms, prob = 0.95), digits = 4)


# Plot-level 'year' slope:
# ------------------------
plot <- fit_brms %>%
  spread_draws(b_year, r_plot[plot, var]) %>%
  filter(var == "year") %>%
  mutate(plot_coeff = -(b_year + r_plot),
         b_year = -b_year) %>%
  dplyr::select(-r_plot)

# Figure of credibility intervals (for composite figure):
plot %>%
  filter(! plot == "BEK01") %>% 
  point_interval(.point = median, .interval = hdi, .width = 0.9) %>%
  mutate(signif = ifelse(plot_coeff.lower < 0 & plot_coeff.upper > 0, 0, 1),
         colour = ifelse(plot_coeff < 0, "neg", "pos")) %>%
  ggplot(aes(plot_coeff, reorder(plot, plot_coeff))) +
  geom_pointinterval(aes(xmin = plot_coeff.lower, xmax = plot_coeff.upper,
                         alpha = signif, fill = colour), shape = 21, point_size = 1.5,
                     size = 1) +
  scale_fill_manual(values = c("red3", "blue")) +
  scale_alpha(range = c(0.3, 1)) +
  geom_vline(xintercept = 0, lty = 1, colour = "darkred") +
  geom_vline(aes(xintercept = b_year), lty = 1) +
  geom_vline(aes(xintercept = b_year.lower), lty = 2) +
  geom_vline(aes(xintercept = b_year.upper), lty = 2) +
  labs(y = "Plots",
       x = expression(paste(Slope~of~mortality~risk~change~per~year~(italic(beta[k])))),
       subtitle = "Plot-level temporal change of mortality (post-1984)") +
  my_theme_big +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  guides(fill = F, alpha = F)

# Figure of predictions from the model:
# *************************************

library(modelr)

# All plots with >2 censuses:
df_unstd <- df_all_sp_covar %>%
  filter(! plot %in% c("BEK01", "CBAY"), year >= 1984) %>%
  left_join(clim_final)

# Analysis based on plots without cyclone damage:
df_unstd <- df_all_sp %>% 
  filter(year >= 1984)

# Population-level prediction:
# ----------------------------
grid <- df_all_sp_covar2 %>%
  data_grid(year = seq_range(year, n = 50))

# NoCyclCensus or NoCyclPlots:
grid <- df_all_sp_2 %>%
  data_grid(year = seq_range(year, n = 50))

y_pred <- fit_brms %>%
  add_fitted_draws(newdata = grid,
                   allow_new_levels = FALSE,
                   re_formula = NA) %>%
  ungroup %>%
  # Back-transform:
  mutate(year = (year * sd(df_unstd$year)) + mean(df_unstd$year),
         .value = (.value * sd(df_unstd$K_tmp)) + mean(df_unstd$K_tmp)) %>%
  # Inverse-logit:
  mutate(.value_inv = 1 - boot::inv.logit(.value))

y_pred_summ <- y_pred %>%
  # Summarise posterior of predictions per covariate value:
  group_by(year) %>%
  point_interval(.point = median, .interval = qi, .width = 0.95) %>%
  ungroup

# Figure Mortality (grand 'year' slope):
y_pred_summ %>%
  ggplot(aes(year, .value_inv)) +
  geom_ribbon(aes(ymin = .value_inv.lower, ymax = .value_inv.upper),
              fill = "grey70", alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = c(1984, seq(1990, 2010, 10), 2019),
                     limits = c(1984, 2019)) +
  # scale_y_continuous(breaks = seq(0.01, 0.022, 0.004)) +
  # geom_line(aes(year, y = 0.00949), linetype = "dashed",
  #           colour = "red3", alpha = 0.4) +
  # geom_line(aes(year, y = 0.019573661), linetype = "dashed",
  #           colour = "red3", alpha = 0.4) +
  labs(y = expression(paste(Mortality~risk~(yr**-1)))) +
  guides(fill = FALSE, colour = FALSE) +
  # Specification for inset figure:
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(size = 0.5),
        axis.text = element_text(size = 23.10876),
        axis.title = element_text(size = 27.73051),
        axis.title.x = element_blank())

## II.4.2. Model M3: Plot-level analysis of the pattern, with climate
#####################################################################
#####################################################################

# Packages:
# ---------
library(tidyverse)
library(brms)

load("DATA/df_all_sp_covar.RData")
load("DATA/clim_final2.RData")

df_all_sp_covar2 <- df_all_sp_covar %>% 
  filter(! plot %in% c("BEK01", "CBAY"), year >= 1984) %>% 
  left_join(clim_final) %>%
  mutate_at(vars(K_tmp, year, Tmax_Q_up_u, vpd_Q_up_u, mcwd_up_u), ~scale(.))

df_sel_sp <- df_all_sp_covar2 

# Model:
# ******

# Model parameters:
iterations <- 2000
chains <- 4
seed <- 42

# VPD and MCWD:
fit_brms <- brm(formula = K_tmp ~ 1 + year + vpd_Q_up_u + mcwd_up_u + 
                  year:vpd_Q_up_u + year:mcwd_up_u + (1 + year | plot) + (1 | sp),
                data = df_sel_sp,
                family = gaussian(),
                prior = c(prior(student_t(2, 0, 1), class = "Intercept"),
                          prior(student_t(2, 0, 1), class = "b"),
                          prior(exponential(1), class = "sigma"),
                          prior(exponential(1), class = "sd"),
                          prior(lkj(2), class = "cor")),
                iter = iterations,
                warmup = 500,
                chains = chains,
                cores = chains,
                seed = seed,
                control = list(adapt_delta = 0.99, max_treedepth = 15))

# or Tmax and MCWD (to avoid collinearity issues between VPD and Tmax):

fit_brms <- brm(formula = K_tmp ~ 1 + year + Tmax_Q_up_u + mcwd_up_u + 
                  year:Tmax_Q_up_u + year:mcwd_up_u + (1 + year | plot) + (1 | sp),
                data = df_sel_sp,
                family = gaussian(),
                prior = c(prior(student_t(2, 0, 1), class = "Intercept"),
                          prior(student_t(2, 0, 1), class = "b"),
                          prior(exponential(1), class = "sigma"),
                          prior(exponential(1), class = "sd"),
                          prior(lkj(2), class = "cor")),
                iter = iterations,
                warmup = 500,
                chains = chains,
                cores = chains,
                seed = seed,
                control = list(adapt_delta = 0.99, max_treedepth = 15))

save(fit_brms, file = "output/fit_brms_effect_year_on_post84Ktmp_varsl_plot_level2_vpd_mcwd.RData")
save(fit_brms, file = "output/fit_brms_effect_year_on_post84Ktmp_varsl_plot_level2_Tmax_mcwd.RData")


## Analysis of model output:
# **************************
# **************************

# Packages:
# ---------
library(tidyverse)
library(tidybayes)
library(viridis)
library(ggridges)

# Model output:
# -------------
# load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot_level2_vpd_precip.RData")
load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot_level2_Tmax_mcwd.RData")
load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_plot_level2_vpd_mcwd.RData")

# Summary of the model:
print(summary(fit_brms, prob = 0.95), digits = 4)

# Plot-level 'year' coefficient:
# -----------------------
plot <- fit_brms %>% 
  spread_draws(b_year, r_plot[plot, var]) %>% 
  filter(var == "year") %>% 
  mutate(plot_coeff = b_year + r_plot) %>% 
  select(-r_plot) 

# Figure of credibility intervals: (Figure for Fig. 3)
plot %>% 
  point_interval(.point = median, .interval = hdi, .width = 0.95) %>% 
  mutate(signif = ifelse(plot_coeff.lower < 0 & plot_coeff.upper > 0, 0, 1),
         colour = ifelse(plot_coeff < 0, "neg", "pos")) %>% 
  ggplot(aes(plot_coeff, reorder(plot, plot_coeff))) +
  geom_pointinterval(aes(xmin = plot_coeff.lower, xmax = plot_coeff.upper,
                         alpha = signif, fill = colour), shape = 21, point_size = 1.5,
                     size = 1) +
  scale_fill_manual(values = c("red3", "blue")) +
  scale_alpha(range = c(0.3, 1)) +
  geom_vline(xintercept = 0, lty = 1, colour = "darkred") +
  geom_vline(aes(xintercept = b_year), lty = 1) +
  geom_vline(aes(xintercept = b_year.lower), lty = 2) +
  geom_vline(aes(xintercept = b_year.upper), lty = 2) +
  labs(y = "Plots",
       x = expression(paste(Slope~of~mortality~risk~change~per~year~(italic(beta[k])))),
       subtitle = "Plot-level temporal change of survivorship (post-1984)") +
  my_theme_big +
  guides(fill = F, alpha = F) +
  scale_x_continuous(limits = c(-0.2535, 0.6935), 
                     breaks = seq(-0.2, 0.6, 0.2))

# Effect of level 2 covariates (climate covariates):
# *****************************
# coeffs <- fit_brms %>%
#   gather_draws(b_year, b_Tmax_Q_up_u, b_mcwd_up_u,
#                `b_year:Tmax_Q_up_u`, `b_year:mcwd_up_u`) %>%
#   mutate(.value = .value * -1) %>% # express coeff. in terms of mortality risk
#   median_hdi(.width = 0.95)
# 
# coeffs <- coeffs %>%
#   mutate(.variable = c("MCWD", "Tmax", "year", "MCWD on year", "Tmax on year"))
# 
# ordered_covar <- rev(c("year", "Tmax", "MCWD", "Tmax on year", "MCWD on year"))

coeffs <- fit_brms %>%
  gather_draws(b_year, b_vpd_Q_up_u, b_mcwd_up_u,
               `b_year:vpd_Q_up_u`, `b_year:mcwd_up_u`) %>%
  mutate(.value = .value * -1) %>% # express coeff. in terms of mortality risk
  median_hdi(.width = 0.95)

coeffs <- coeffs %>%
  mutate(.variable = c("MCWD", "VPD", "year", "MCWD on year", "VPD on year"))

ordered_covar <- rev(c("year", "VPD", "MCWD", "VPD on year", "MCWD on year"))

coeffs <- coeffs %>%
  mutate(signif = ifelse(.lower < 0 & .upper > 0, 0, 1),
         sign = ifelse(.value > 0, "pos", "neg")) %>% 
  mutate(.variable = as_factor(fct_relevel(.variable, ordered_covar))) 

# Figure:
coeffs %>% 
  ggplot(aes(.value, .variable, xmin = .lower, xmax = .upper)) +
  geom_pointrange(aes(alpha = signif, fill = sign), shape = 21) +
  scale_alpha(range = c(0.3, 1)) +
  scale_fill_manual(values = c("red3", "blue")) +
  guides(alpha = F, fill = F) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = "Coefficient z-scores",
       y = "Covariates",
       subtitle = expression(paste(Climate~effects~on~average~mortality~risk~and~the~rate~of~risk~change~over~time))) +
  my_theme_big +
  theme(panel.grid = element_blank())


## II.3.1. Model M4: Post-84 temporal survivorship decrease per and across species:
#########################################################################
#########################################################################

# Packages:
# ---------
library(tidyverse)
library(brms)

load("DATA/df_all_sp_covar.RData")

# If only for plots without cyclone intervals:
load("DATA/df_all_sp_NoCyclCensus.RData")
df_all_sp_2 <- df_all_sp %>% 
  filter(year >= 1984) %>% 
  mutate_at(vars(K_tmp, year), ~scale(.))

# If only for plots with no cyclone damage in 50 years:
load("DATA/df_all_sp_NoCyclPlot.RData")
df_all_sp_2 <- df_all_sp %>% 
  filter(year >= 1984) %>% 
  mutate_at(vars(K_tmp, year), ~scale(.))

## Model:
# *******
# *******

df_all_sp_covar2 <- df_all_sp_covar %>% 
  filter(year >= 1984) %>% 
  mutate_at(vars(K_tmp, year), ~scale(.))

# Nb real obs. and stems:
tmp <- data.big %>% 
  filter(year_0 >= 1984)
nrow(tmp)
length(unique(tmp$stem))
nrow(df_all_sp_covar2)

# Model parameters:
iterations <- 6000
chains <- 4
seed <- 42

fit_brms <- brm(formula = K_tmp ~ 1 + year + (1 + year | sp) + (1 | plot),
                data = df_all_sp_covar2,
                family = gaussian(),
                prior = c(prior(normal(0, 0.5), class = "Intercept"),
                          prior(normal(0, 0.5), class = "b"),
                          prior(normal(0, 1), class = "sigma"),
                          prior(normal(0, 1), class = "sd"),
                          prior(lkj(2), class = "cor")),
                iter = iterations,
                warmup = 500,
                chains = chains,
                cores = chains,
                seed = seed,
                control = list(adapt_delta = 0.99, max_treedepth = 15))

save(fit_brms, file = "OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_sp_withBEK.RData")
# save(fit_brms, file = "OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_sp_noCyclDam.RData")

## Analysis of model output:
# **************************
# **************************

# Packages:
# ---------
library(tidyverse)
library(tidybayes)
library(viridis)
library(ggridges)

# Model output:
# -------------
load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_sp.RData")
# load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_sp_NoCyclCensus.RData")
# load("OUTPUT/fit_brms_effect_year_on_post84Ktmp_varsl_sp_NoCyclPlot.RData")

# Summary of the model:
print(summary(fit_brms, prob = 0.95), digits = 4)

# Species-level 'year' slope:
# ---------------------------
sp <- fit_brms %>% 
  spread_draws(b_year, r_sp[sp, var]) %>% 
  filter(var == "year") %>% 
  mutate(sp_coeff = -(b_year + r_sp),
         b_year = -b_year) %>% 
  dplyr::select(-r_sp) 

# Proportion of species with positive and negative slopes:
sptest <- sp %>% 
  ungroup %>% 
  dplyr::select(sp, sp_coeff) %>% 
  group_by(sp) %>% 
  median_qi(.width = 0.95) %>% 
  mutate(sign = ifelse(.lower < 0 & .upper > 0, 0, 1)) %>% 
  filter(sign == 1) %>% 
  mutate(posneg = ifelse(sp_coeff > 0, "positive", "negative"))

length(which(sptest$posneg == "positive")) / 81 ; length(which(sptest$posneg == "positive"))
length(which(sptest$posneg == "negative")) / 81 ; length(which(sptest$posneg == "negative")) 

# length(which(sptest$posneg == "positive")) / 63 ; length(which(sptest$posneg == "positive"))
# length(which(sptest$posneg == "negative")) / 63 ; length(which(sptest$posneg == "negative")) 

# Figure of credibility intervals (for composite figure):
sp %>% 
  point_interval(.point = median, .interval = hdi, .width = 0.95) %>% 
  mutate(signif = ifelse(sp_coeff.lower < 0 & sp_coeff.upper > 0, 0, 1),
         colour = ifelse(sp_coeff < 0, "neg", "pos")) %>% 
  ggplot(aes(sp_coeff, reorder(sp, sp_coeff))) +
  geom_pointinterval(aes(xmin = sp_coeff.lower, xmax = sp_coeff.upper,
                         alpha = signif, fill = colour), shape = 21, point_size = 1.5,
                     size = 1) +
  scale_fill_manual(values = c("red3", "blue")) +
  scale_alpha(range = c(0.3, 1)) +
  geom_vline(xintercept = 0, lty = 1, colour = "darkred") +
  geom_vline(aes(xintercept = b_year), lty = 1) +
  geom_vline(aes(xintercept = b_year.lower), lty = 2) +
  geom_vline(aes(xintercept = b_year.upper), lty = 2) +
  labs(y = "Species",
       x = expression(paste(Slope~of~mortality~risk~change~per~year~(italic(beta[j])))),
       subtitle = "Species-level temporal change of mortality (post-1984)") +
  my_theme_big +
  theme(axis.text.y = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(-0.2535, 0.6935), 
                     breaks = seq(-0.2, 0.6, 0.2)) +
  guides(fill = F, alpha = F) 

# Table of coefficients (Supplementary Material):
table_sp_level_year_effect <- sp %>% 
  point_interval(.point = median, .interval = hdi, .width = 0.95) %>% 
  ungroup %>% 
  select(sp, sp_coeff, sp_coeff.lower, sp_coeff.upper) %>% 
  left_join(datafinal_all_sp %>% 
              select(sp = code, species = taxon, genus, family) %>% 
              group_by(sp, species, genus, family) %>% 
              summarise()) %>% 
  select(species:family, sp, sp_coeff, sp_coeff.lower, sp_coeff.upper)

write_tsv(table_sp_level_year_effect, "OUTPUT/Table_sp_level_year_effect.txt")
write_tsv(table_sp_level_year_effect, "OUTPUT/Table_sp_level_year_effect_NoCyclCensus.txt")
write_tsv(table_sp_level_year_effect, "OUTPUT/Table_sp_level_year_effect_NoCyclPlot.txt")

# Figure of predictions from the model for a few species:
# *******************************************************

library(modelr)

df_unstd <- df_all_sp_covar %>% 
  filter(year >= 1984)

# NoCyclCensus or NoCyclPlot:
df_unstd <- df_all_sp %>%
  filter(year >= 1984) 

# Population-level prediction:
# ----------------------------
grid <- df_all_sp_covar2 %>% 
  data_grid(year = seq_range(year, n = 50))

# NoCyclCensus or NoCyclPlot:
grid <- df_all_sp_2 %>% 
  data_grid(year = seq_range(year, n = 50))

y_pred <- fit_brms %>% 
  add_fitted_draws(newdata = grid,
                   allow_new_levels = FALSE,
                   re_formula = NA) %>% 
  ungroup %>% 
  # Back-transform:
  mutate(year = (year * sd(df_unstd$year)) + mean(df_unstd$year),
         .value = (.value * sd(df_unstd$K_tmp)) + mean(df_unstd$K_tmp)) %>% 
  # Inverse-logit:
  mutate(.value_inv = 1 - boot::inv.logit(.value))

y_pred_summ <- y_pred %>% 
  # Summarise posterior of predictions per covariate value:
  group_by(year) %>% 
  point_interval(.point = median, .interval = qi, .width = 0.95) %>% 
  ungroup 

# Figure Mortality (population-level):
y_pred_summ %>%
  ggplot(aes(year, .value_inv)) +
  geom_ribbon(aes(ymin = .value_inv.lower, ymax = .value_inv.upper),
              fill = "grey70", alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = c(1984, seq(1990, 2010, 10), 2019),
                     limits = c(1984, 2019)) +
  # scale_y_continuous(breaks = seq(0.01, 0.022, 0.004)) +
  # geom_line(aes(year, y = 0.01084075), linetype = "dashed", 
  #           colour = "red3", alpha = 0.4) +
  # geom_line(aes(year, y = 0.01820285), linetype = "dashed", 
  #           colour = "red3", alpha = 0.4) +
  labs(y = expression(paste(Mortality~risk~(yr**-1)))) +
  guides(fill = FALSE, colour = FALSE) +
  # Specification for inset figure:
  my_theme_inset

## II.3.4. Model M5: VPD, Tmax and MCWD niche
##############################################
##############################################

## II.3.4.1. Post-1984 - Covariates = 'year' and VPD, Tmax and MCWD niche:
# ************************************************************************
# ************************************************************************

# Packages:
# ---------
library(tidyverse)
library(brms)

load("DATA/df_all_sp.RData")
load("DATA/df_niche_related_V5.RData")

niche_aus <- niche_aus_1km_thinned %>% 
  group_by(plot, sp, taxon) %>% 
  summarise_all(~mean(., na.rm = T)) %>% 
  ungroup

df_all_sp_covar2 <- df_all_sp %>% 
  filter(sp %in% niche_aus$sp) %>% 
  left_join(niche_aus, by = c("sp", "plot")) %>% 
  filter(nb_plots >= 4, year >= 1984)

df_all_sp_covar2_std <- df_all_sp_covar2 %>% 
  mutate_at(vars(K_tmp, year), ~scale(.)) %>%
  group_by(sp) %>%
  mutate_at(vars(Tmax_max_quant, vpd_max_quant, precip_an_quant, mcwd_quant), ~scale(.)) %>%
  ungroup  

# Nb real obs. and stems:
tmp <- data.big %>% 
  filter(year_0 >= 1984, Sp %in% unique(df_all_sp_covar2$sp))
nrow(tmp)
length(unique(tmp$stem))
nrow(df_all_sp_covar2)

# Model:
# ******

# Model parameters:
iterations <- 2000
chains <- 4
seed <- 42

# Varying slope on species:
fit_brms <- brm(formula = K_tmp ~ 1 + year + Tmax_max_quant + vpd_max_quant + mcwd_quant + 
                  (1 + year + Tmax_max_quant + vpd_max_quant + mcwd_quant | sp) + 
                  (1 | plot),
                data = df_all_sp_covar2_std,
                family = gaussian(),
                prior = c(prior(student_t(2, 0, 1), class = "Intercept"),
                          prior(student_t(2, 0, 1), class = "b"),
                          prior(exponential(1), class = "sigma"),
                          prior(exponential(1), class = "sd"),
                          prior(lkj(2), class = "cor")),
                iter = iterations,
                warmup = 500,
                chains = chains,
                cores = chains,
                seed = seed,
                control = list(adapt_delta = 0.99, max_treedepth = 15))

save(fit_brms, file = "output/fit_brms_effect_year_nich_Tmax_vpd_mcwd_STDperSP_on_post84_varsl_sp_ab4plots_1km_fast.RData")

## Analysis of model output:
# **************************
# **************************

# Packages:
# ---------
library(tidyverse)
library(tidybayes)
library(viridis)
# library(ggridges)

# Model output:
# -------------
load("OUTPUT/fit_brms_effect_year_nich_Tmax_vpd_mcwd_STDperSP_on_post84_varsl_sp_ab4plots_1km_fast.RData")
# load("OUTPUT/fit_brms_effect_year_nich_Tmax_vpd_mcwd_STDperSP_on_post84Ktmp_varsl_sp_fast.RData")

# Summary of the model:
print(summary(fit_brms, prob = 0.95), digits = 4)

# Sp-level coefficient:
# ---------------------
# Year:
sp <- fit_brms %>% 
  spread_draws(b_year, r_sp[sp, var]) %>% 
  filter(var == "year") %>% 
  mutate(sp_coeff = -(b_year + r_sp),
         b_year = -b_year) %>% 
  dplyr::select(-r_sp) 

# Proportion of species with positive and negative slopes:
sptest <- sp %>% 
  ungroup %>% 
  dplyr::select(sp, sp_coeff) %>% 
  group_by(sp) %>% 
  median_qi(.width = 0.9) %>% 
  mutate(sign = ifelse(.lower < 0 & .upper > 0, 0, 1)) %>% 
  filter(sign == 1) %>% 
  mutate(posneg = ifelse(sp_coeff > 0, "positive", "negative"))

length(which(sptest$posneg == "positive")) / 56 ; length(which(sptest$posneg == "positive"))
length(which(sptest$posneg == "negative")) / 56 ; length(which(sptest$posneg == "negative")) 

# VPD:
sp <- fit_brms %>% 
  spread_draws(b_vpd_max_quant, r_sp[sp, var]) %>% 
  filter(var == "vpd_max_quant") %>% 
  mutate(sp_coeff = -(b_vpd_max_quant + r_sp),
         b_vpd_max_quant = -b_vpd_max_quant) %>% 
  dplyr::select(-r_sp) 

# Proportion of species with positive and negative slopes:
sptest <- sp %>% 
  ungroup %>% 
  dplyr::select(sp, sp_coeff) %>% 
  group_by(sp) %>% 
  median_qi(.width = 0.95) %>% 
  mutate(sign = ifelse(.lower < 0 & .upper > 0, 0, 1)) %>% 
  filter(sign == 1) %>% 
  mutate(posneg = ifelse(sp_coeff > 0, "positive", "negative"))

length(which(sptest$posneg == "positive")) / 56 ; length(which(sptest$posneg == "positive"))
length(which(sptest$posneg == "negative")) / 56 ; length(which(sptest$posneg == "negative")) 

# Figure (not used in paper):
sp %>% 
  point_interval(.point = median, .interval = hdi, .width = 0.9) %>% 
  mutate(signif = ifelse(sp_coeff.lower < 0 & sp_coeff.upper > 0, 0, 1),
         colour = ifelse(sp_coeff < 0, "neg", "pos")) %>% 
  ggplot(aes(sp_coeff, reorder(sp, sp_coeff))) +
  geom_pointinterval(aes(xmin = sp_coeff.lower, xmax = sp_coeff.upper,
                         alpha = signif, fill = colour), shape = 21, point_size = 1.5,
                     size = 1) +
  scale_fill_manual(values = c("red3", "blue")) +
  scale_alpha(range = c(0.3, 1)) +
  geom_vline(xintercept = 0, lty = 1, colour = "darkred") +
  geom_vline(aes(xintercept = b_vpd_max_quant), lty = 1) +
  geom_vline(aes(xintercept = b_vpd_max_quant.lower), lty = 2) +
  geom_vline(aes(xintercept = b_vpd_max_quant.upper), lty = 2) +
  labs(y = "Species",
       x = expression(paste(Slope~of~mortality~risk~change~per~year~(italic(beta[j])))),
       subtitle = "Species-level temporal change of mortality (post-1984)") +
  my_theme_big +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  guides(fill = F, alpha = F) 

# Tmax:
sp <- fit_brms %>% 
  spread_draws(b_Tmax_max_quant, r_sp[sp, var]) %>% 
  filter(var == "Tmax_max_quant") %>% 
  mutate(sp_coeff = -(b_Tmax_max_quant + r_sp),
         b_Tmax_max_quant = -b_Tmax_max_quant) %>% 
  dplyr::select(-r_sp) 

# Proportion of species with positive and negative slopes:
sptest <- sp %>% 
  ungroup %>% 
  dplyr::select(sp, sp_coeff) %>% 
  group_by(sp) %>% 
  median_qi(.width = 0.9) %>% 
  mutate(sign = ifelse(.lower < 0 & .upper > 0, 0, 1)) %>% 
  filter(sign == 1) %>% 
  mutate(posneg = ifelse(sp_coeff > 0, "positive", "negative"))

length(which(sptest$posneg == "positive")) / 56 ; length(which(sptest$posneg == "positive"))
length(which(sptest$posneg == "negative")) / 56 ; length(which(sptest$posneg == "negative")) 

# Figure (not used in paper):
sp %>% 
  point_interval(.point = median, .interval = hdi, .width = 0.9) %>% 
  mutate(signif = ifelse(sp_coeff.lower < 0 & sp_coeff.upper > 0, 0, 1),
         colour = ifelse(sp_coeff < 0, "neg", "pos")) %>% 
  ggplot(aes(sp_coeff, reorder(sp, sp_coeff))) +
  geom_pointinterval(aes(xmin = sp_coeff.lower, xmax = sp_coeff.upper,
                         alpha = signif, fill = colour), shape = 21, point_size = 1.5,
                     size = 1) +
  scale_fill_manual(values = c("red3", "blue")) +
  scale_alpha(range = c(0.3, 1)) +
  geom_vline(xintercept = 0, lty = 1, colour = "darkred") +
  geom_vline(aes(xintercept = b_Tmax_max_quant), lty = 1) +
  geom_vline(aes(xintercept = b_Tmax_max_quant.lower), lty = 2) +
  geom_vline(aes(xintercept = b_Tmax_max_quant.upper), lty = 2) +
  labs(y = "Species",
       x = expression(paste(Slope~of~mortality~risk~change~per~year~(italic(beta[j])))),
       subtitle = "Species-level temporal change of mortality (post-1984)") +
  my_theme_big +
  theme(axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank()) +
  guides(fill = F, alpha = F) 

# See separate script for construction of Fig. 4

## Predicted mortality risk as a function of spatial changes in VPD niche:
# ------------------------------------------------------------------------

df_unstd <- df_all_sp_covar %>% 
  filter(year >= 1984)

# NoCyclCensus or NoCyclPlot:
df_unstd <- df_all_sp %>%
  filter(year >= 1984) 

# Population-level prediction:
# ----------------------------
grid <- df_all_sp_covar2 %>% 
  data_grid(year = seq_range(year, n = 50))

# NoCyclCensus or NoCyclPlot:
grid <- df_all_sp_2 %>% 
  data_grid(year = seq_range(year, n = 50))

y_pred <- fit_brms %>% 
  add_fitted_draws(newdata = grid,
                   allow_new_levels = FALSE,
                   re_formula = NA) %>% 
  ungroup %>% 
  # Back-transform:
  mutate(year = (year * sd(df_unstd$year)) + mean(df_unstd$year),
         .value = (.value * sd(df_unstd$K_tmp)) + mean(df_unstd$K_tmp)) %>% 
  # Inverse-logit:
  mutate(.value_inv = 1 - boot::inv.logit(.value))

y_pred_summ <- y_pred %>% 
  # Summarise posterior of predictions per covariate value:
  group_by(year) %>% 
  point_interval(.point = median, .interval = qi, .width = 0.95) %>% 
  ungroup 

# Figure Mortality (population-level):
y_pred_summ %>%
  ggplot(aes(year, .value_inv)) +
  geom_ribbon(aes(ymin = .value_inv.lower, ymax = .value_inv.upper),
              fill = "grey70", alpha = 0.4) +
  geom_line() +
  scale_x_continuous(breaks = c(1984, seq(1990, 2010, 10), 2019),
                     limits = c(1984, 2019)) +
  # scale_y_continuous(breaks = seq(0.01, 0.022, 0.004)) +
  # geom_line(aes(year, y = 0.01084075), linetype = "dashed", 
  #           colour = "red3", alpha = 0.4) +
  # geom_line(aes(year, y = 0.01820285), linetype = "dashed", 
  #           colour = "red3", alpha = 0.4) +
  labs(y = expression(paste(Mortality~risk~(yr**-1)))) +
  guides(fill = FALSE, colour = FALSE) +
  # Specification for inset figure:
  my_theme_inset


## II.3.2. Models M6: Species traits to explain intersp. differences in 'year' slope:
#########################################################################
#########################################################################

# Packages:
# ---------
library(tidyverse)
library(brms)

load("DATA/df_all_sp_covar.RData")

## Model:
# *******
# *******

df_all_sp_covar2 <- df_all_sp_covar %>% 
  filter(year >= 1984, trait_data == 1) %>% 
  dplyr::select(-leaf_d15N) %>% # remove unused traits
  mutate(leaf_d13C = leaf_d13C + abs(min(leaf_d13C) + min(leaf_d13C)/10)) %>%  
  mutate_at(vars(dbh_max:LMA), ~log(.)) %>%  # log-transformation of all traits across all plots
  mutate_at(vars(dbh_max:LMA), ~scale(.)) %>%
  mutate_at(vars(K_tmp, year), ~scale(.))

# Nb real obs. and stems:
tmp <- data.big %>% 
  filter(year_0 >= 1984, Sp %in% unique(df_all_sp_covar2$sp))
nrow(tmp)
length(unique(tmp$stem))
nrow(df_all_sp_covar2)

# Model parameters:
iterations <- 6000
chains <- 4
seed <- 42

# Example of model structure for leaf mass per area (LMA):
fit_brms <- brm(formula = K_tmp ~ 1 + year + LMA + year:LMA + (1 + year | sp) + (1 | plot),
                data = df_all_sp_covar2,
                family = gaussian(),
                prior = c(prior(normal(0, 1), class = "Intercept"),
                          prior(normal(0, 1), class = "b"),
                          prior(exponential(1), class = "sigma"),
                          prior(exponential(1), class = "sd"),
                          prior(lkj(2), class = "cor")),
                iter = iterations,
                warmup = 500,
                chains = chains,
                cores = chains,
                seed = seed,
                control = list(adapt_delta = 0.99, max_treedepth = 15))

save(fit_brms, file = "output/fit_brms_effect_year_on_post84Ktmp_varsl_sp_withBEK_trait_LMA.RData")

## Analysis of model output:
# **************************
# **************************

# Packages:
# ---------
library(tidyverse)
library(tidybayes)
library(viridis)
library(ggridges)

# Model output:
# -------------
# We extract the trait effect on species-level intercept and 'year' slope (cross-level interaction):

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_Asat.RData")
trait_effects <- fit_brms %>% 
  gather_draws(b_Asat, `b_year:Asat`) %>% 
  median_hdi(.width = 0.95) %>% 
  ungroup %>% 
  mutate(.variable = str_replace(.variable, "b_Asat", "Asat on K_lat")) %>% 
  mutate(.variable = str_replace(.variable, "b_year:Asat", "Asat on slope"))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_Amax.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_Amax, `b_year:Amax`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_Amax", "Amax on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:Amax", "Amax on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_gsat.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_gsat, `b_year:gsat`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_gsat", "gsat on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:gsat", "gsat on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_gmax.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_gmax, `b_year:gmax`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_gmax", "gmax on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:gmax", "gmax on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_Vcmax.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_Vcmax, `b_year:Vcmax`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_Vcmax", "Vcmax on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:Vcmax", "Vcmax on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_Jmax.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_Jmax, `b_year:Jmax`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_Jmax", "Jmax on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:Jmax", "Jmax on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_Rd.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_Rd, `b_year:Rd`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_Rd", "Rd on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:Rd", "Rd on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_d13C.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_leaf_d13C, `b_year:leaf_d13C`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_leaf_d13C", "d13C on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:leaf_d13C", "d13C on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_leaf_P.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_leaf_P, `b_year:leaf_P`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_leaf_P", "leaf_P on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:leaf_P", "leaf_P on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_leaf_N.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_leaf_N, `b_year:leaf_N`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_leaf_N", "leaf_N on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:leaf_N", "leaf_N on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_LMA.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_LMA, `b_year:LMA`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_LMA", "LMA on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:LMA", "LMA on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_LA.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_LA, `b_year:LA`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_LA", "LA on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:LA", "LA on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_WD.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>% 
              gather_draws(b_WD, `b_year:WD`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_WD", "WD on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:WD", "WD on slope")))

load("OUTPUT/tmp_withROB06_2019/fit_brms_effect_year_on_post84Ktmp_varsl_sp_trait_dbh_max.RData")
trait_effects <- trait_effects %>% 
  bind_rows(fit_brms %>%
              gather_draws(b_dbh_max, `b_year:dbh_max`) %>% 
              median_hdi(.width = 0.95) %>% 
              ungroup %>% 
              mutate(.variable = str_replace(.variable, "b_dbh_max", "DBHmax on K_lat")) %>% 
              mutate(.variable = str_replace(.variable, "b_year:dbh_max", "DBHmax on slope")))

# Summary of the model:
print(summary(fit_brms, prob = 0.95), digits = 4)

# Generate figure summarising the standardised effect sizes of trait-related effects:
ordered_covar <- rev(c("Asat on K_lat", "Amax on K_lat", "gsat on K_lat", "gmax on K_lat",
                       "Vcmax on K_lat", "Jmax on K_lat", "Rd on K_lat", "d13C on K_lat", 
                       "leaf_P on K_lat", "leaf_N on K_lat", "LMA on K_lat", "LA on K_lat", 
                       "WD on K_lat", "DBHmax on K_lat",
                       "Asat on slope", "Amax on slope", "gsat on slope", "gmax on slope",
                       "Vcmax on slope", "Jmax on slope", "Rd on slope", "d13C on slope", 
                       "leaf_P on slope", "leaf_N on slope", "LMA on slope", "LA on slope", 
                       "WD on slope", "DBHmax on slope"))

trait_effects2 <- trait_effects %>% 
  mutate(.value = .value * -1,
         .lower2 = .upper * -1,
         .upper2 = .lower * -1) %>% 
  dplyr::select(.variable, .value, .lower = .lower2, .upper = .upper2) %>% 
  mutate(signif = ifelse(.lower < 0 & .upper > 0, 0, 1),
         sign = ifelse(.value > 0, "pos", "neg")) %>% 
  mutate(.variable = as_factor(fct_relevel(.variable, ordered_covar))) %>% 
  mutate(type = ifelse(str_detect(.variable, "slope"), "slope", "intercept"))

# # Figure (all in one figure):
# trait_effects2 %>% 
#   ggplot(aes(.value, .variable, xmin = .lower, xmax = .upper)) +
#   geom_pointrange(aes(alpha = signif, fill = sign), shape = 21) +
#   scale_alpha(range = c(0.3, 1)) +
#   scale_fill_manual(values = c("red3", "blue")) +
#   guides(alpha = F, fill = F) +
#   geom_vline(xintercept = 0, lty = 2) +
#   labs(x = "Coefficient z-scores",
#        y = "Covariates",
#        subtitle = expression(paste(Trait~effects~on~average~italic(K_lat)~and~italic(K_lat)~temporal~decrease))) +
#   theme_classic()

# Figure (trait mediation of species-level intercept):
trait_effects2 %>% 
  filter(type == "intercept") %>% 
  mutate(.variable = str_remove(.variable, " on K_lat")) %>% 
  mutate(.variable = as_factor(fct_relevel(.variable, 
                                           rev(c("Asat", "Amax", "gsat", "gmax",
                                                 "Vcmax", "Jmax", "Rd", "d13C",
                                                 "leaf_P", "leaf_N", "LMA", "LA",
                                                 "WD", "DBHmax"))))) %>% 
  ggplot(aes(.value, .variable, xmin = .lower, xmax = .upper)) +
  geom_pointrange(aes(alpha = signif, fill = sign), shape = 21) +
  scale_alpha(range = c(0.3, 1)) +
  scale_fill_manual(values = c("red3", "blue")) +
  guides(alpha = F, fill = F) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = expression(paste(Coefficient~italic(z)*"-scores"~(alpha[1]))),
       y = "Covariates",
       subtitle = expression(paste(Trait~effects~on~average~mortality~risk))) +
  scale_x_continuous(limits = c(-0.4, 0.3)) +
  my_theme_big +
  theme(axis.title.y = element_blank())

# Figure (trait mediation of species-level 'year' slope):
trait_effects2 %>% 
  filter(type == "slope") %>% 
  mutate(.variable = str_remove(.variable, " on slope")) %>% 
  mutate(.variable = as_factor(fct_relevel(.variable, 
                                           rev(c("Asat", "Amax", "gsat", "gmax",
                                                 "Vcmax", "Jmax", "Rd", "d13C",
                                                 "leaf_P", "leaf_N", "LMA", "LA",
                                                 "WD", "DBHmax"))))) %>% 
  ggplot(aes(.value, .variable, xmin = .lower, xmax = .upper)) +
  geom_pointrange(aes(alpha = signif, fill = sign), shape = 21) +
  scale_alpha(range = c(0.1, 0.3)) +
  scale_fill_manual(values = c("red3", "blue")) +
  guides(alpha = F, fill = F) +
  geom_vline(xintercept = 0, lty = 2) +
  labs(x = expression(paste(Coefficient~italic(z)*"-scores"~(beta[1]))),
       y = "Covariates",
       subtitle = expression(paste(Trait~effects~on~year~slope))) +
  scale_x_continuous(limits = c(-0.4, 0.3)) +
  my_theme_big +
  theme(axis.title.y = element_blank())


## II.5. Models of climate over time:
#####################################

library(tidyverse)
library(brms)
library(modelr)

# Load and prepare data:
# ----------------------
load("DATA/clim.RData")

c <- clim %>% 
  select(plot, month, year, Tmax, vpd) %>% 
  mutate(ym = round(year + month/12, 2)) %>% 
  as_tibble

clim_for_mod <- c %>%
  group_by(plot) %>%
  mutate_at(vars(Tmax, vpd), ~scale(.)) %>%
  ungroup %>%
  mutate(ym = scale(ym))

fit_brms <- brm(vpd ~ 1 + s(ym, by = plot, bs = "bs") + (1 | plot), 
                data = clim_for_mod,
                chains = 4,
                cores = 4,
                iter = 1000)

# save(fit_brms, file = "output/fit_brms_vpd_over_time_GAM_varsl.RData")

# Analyse output:
# ---------------
load("OUTPUT/tmp_withROB06_2019/fit_brms_vpd_over_time_GAM.RData")

plot(conditional_smooths(fit_brms))

# fitted(fit_brms) %>% 
#   data.frame() %>% 
#   bind_cols(select(clim_for_mod, plot, vpd, Tmax, year)) %>% 
#   ggplot(aes(year, vpd, ymin = Q2.5, ymax = Q97.5)) +
#   # geom_hline(yintercept = fixef(gam_vpd2)[1, 1], color = "white", linetype = 2) +
#   geom_point(color = "#ffb7c5", alpha = 0.5) +
#   geom_ribbon(fill = "white", alpha = 0.6) +
#   geom_line(aes(year, Estimate), colour = "white") +
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         panel.background = element_rect(fill = "#4f455c")) +
#   facet_wrap(~reorder(plot, desc(vpd_u)), ncol = 6)

# Grid for prediction at the plot level:
# --------------------------------------

grid <- clim_for_mod %>% 
  data_grid(ym = seq_range(ym, n = 100))

# Predictions:
pred_clim <- fit_brms %>% 
  add_fitted_draws(newdata = grid,
                   allow_new_levels = TRUE,
                   re_formula = NA) %>% 
  mutate(ym = (ym * sd(c$ym)) + mean(c$ym)) %>% 
  # back-transform outcome:
  mutate(.value = (.value * sd(c$vpd)) + mean(c$vpd)) %>% 
  median_qi(.width = 0.9) %>% 
  ungroup 

# Figure:
# GAM with b-splines, and k = 3: (fit_brms_k3)
pred_clim %>% 
  ggplot(aes(ym, .value)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = "grey70", alpha = 0.6) +
  geom_line() +
  # add data:
  geom_point(data = c %>% 
               filter(plot %in% c("EP44")) %>% # EP9, 44
               group_by(plot, year) %>% 
               summarise_at(vars(Tmax, vpd), ~max(.)) %>% 
               ungroup %>% 
               rename(ym = year),
             aes(ym, vpd, fill = plot),
             shape = 21, alpha = 0.5) +
  labs(y = "Annual max VPD (hPa)") +
  scale_x_continuous(breaks = c(1971, seq(1980, 2010, 10), 2019)) +
  my_theme_big +
  theme(panel.grid = element_blank(),
        axis.title.x = element_blank())


### III. Model of absolute diameter growth rate as a function of tree DBH and 'year':
#####################################################################################

# Loading the data:
# *****************
load('DATA/datafinal_all_sp.RData')

# Species for which survival was modelled:
# ----------------------------------------
l_sp <- list.files(path = "OUTPUT/STAN_survival_model_outputs", 
                   pattern = "df_Ktmp_plot_year_")
# Isolate species code:
start <- stringr::str_locate(l_sp, "df_Ktmp_plot_year_")[, 2][1]
sp <- stringr::str_sub(l_sp, start + 1, start + 6) # sp name only

# Has AGR decreased / increased with time?

data <- datafinal_all_sp %>% 
  filter(dead == 0, outlier == 0, new_ID == 0) %>% 
  # Remove palm and fern:
  filter(! taxon %in% c("Normanbya normanbyi", "Cyathea cooperi")) %>% 
  filter(! plot %in% c("BEK01", "CBAY", "ROB06")) %>% # remove plots < 3 censuses 
  select(-c('stem2', 'new_ID', 'any_abnormal', 'abnormal_last', 'family', 'genus', 'date_0', 'date_1', 'nbdays', 
            'agr_ba', 'rgr_ba', 'lnba_0')) %>% 
  # Only keep species for which survival was modelled:
  filter(code %in% sp)

saveRDS(data, file = "DATA/data_81sp_for_model_AGR_over_time.rds")
data <- readRDS('DATA/data_81sp_for_model_AGR_over_time.rds')

min_agr <- min(data$agr_dbh, na.rm = T)

# AGR and year transformation:
data2 <- data %>% 
  mutate(agr_dbh = log(agr_dbh + abs(min_agr + min_agr/10))) %>% 
  mutate_at(vars(year_0, lndbh_0), ~scale(.))

# Hierarchical Bayesian model of AGR through time:
# ************************************************
# Priors:
prior_brms <- c(prior(normal(0, 1), class = "Intercept"), 
                prior(normal(0, 0.3), class = "b"),
                prior(normal(0, 1), class = "sd"),
                prior(lkj(2), class = "cor")) 

# Varying slope per species:
# --------------------------
fit_brms_AGR <- brm(formula = agr_dbh ~ 1 + lndbh_0 + year_0 + 
                      (1 + lndbh_0 + year_0 | code) + (1 | plot) +
                      (1 | stem),
                    data = data2,
                    family = gaussian(),
                    prior = prior_brms,
                    iter = 3000,
                    chains = 4,
                    cores = 4,
                    seed = 42,
                    control = list(adapt_delta = 0.98, max_treedepth = 12))

save(fit_brms_AGR, file = "OUTPUT/fit_brms_effect_year_AGR_varsl_sp.RData")

## Model outputs ##
###################

library(ggridges)

# Varying slope per species (81 species studied for survival):
# **************************
load("OUTPUT/fit_brms_effect_year_AGR_varsl_sp_81sp.RData")

# Summary:
print(summary(fit_brms_AGR, prob = 0.95), digits = 3)

# Supplementary Figure for the manuscript:
# ----------------------------------------

# 'year' hyperparameter (inset):
fig_hyperpar <- fit_brms_AGR %>% 
  gather_draws(b_year_0, b_lndbh_0) %>% 
  filter(.variable == "b_year_0") %>% 
  mutate(.variable = "Year") %>% 
  ggplot(aes(.value, .variable,
             fill = 0.5 - abs(0.5 - stat(ecdf)))) +
  stat_density_ridges(geom = "density_ridges_gradient", calc_ecdf = TRUE) +
  scale_fill_viridis_c(name = "Tail probability", direction = -1) +
  geom_vline(xintercept = 0, lty = 2) +
  theme_ridges() +
  labs(subtitle = expression(paste("Grand 'year' slope (" * italic(beta)["2, 0"]*")")),
       x = "Coefficient") +
  theme(axis.line = element_line(size = 0.5, colour = "black"),
        axis.text = element_text(size = 20, colour = "black"),
        axis.title = element_text(size = 24, colour = "black"),
        plot.subtitle = element_text(size = 24, colour = "black"),
        legend.text = element_text(size = 16, colour = "black"),
        axis.title.y = element_blank(),
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key.width = unit(1.2, "cm"),
        legend.title = element_text(size = 18, colour = "black",
                                    margin = margin(-20, 7, 0, 0, "pt")))

fig_hyperpar

# 'year' slope per species (main figure):
fit_brms_AGR %>% 
  spread_draws(b_year_0, r_code[code, variable]) %>% 
  filter(variable == "year_0") %>% 
  mutate(Year = b_year_0 + r_code) %>% 
  point_interval(.point = median, .interval = hdi, .width = 0.95) %>% 
  mutate(alpha = ifelse(Year.lower < 0 & Year.upper > 0, 0.3, 1)) %>% 
  ggplot(aes(Year, reorder(code, Year), alpha = alpha)) +
  geom_pointinterval(aes(xmin = Year.lower, xmax = Year.upper), 
                     point_size = 1, size = 0.5) +
  geom_vline(xintercept = 0, lty = 2, colour = "red", size = 1) +
  geom_vline(aes(xintercept = b_year_0), colour = "blue", size = 0.7, alpha = 0.7) +
  geom_vline(aes(xintercept = b_year_0.lower), colour = "blue", 
             size = 0.7, lty = 2, alpha = 0.7) +
  geom_vline(aes(xintercept = b_year_0.upper), colour = "blue", 
             size = 0.7, lty = 2, alpha = 0.7) +
  theme_classic() +
  labs(subtitle = expression(paste("Species-specific 'year' slope (" * italic(beta)["2"*italic("j")] * ")")),
       y = "Species",
       x = "Coefficient") +
  guides(alpha = F) +
  theme(axis.text.y = element_blank(),
        axis.line = element_line(size = 0.5, colour = "black"),
        axis.text = element_text(size = 12, colour = "black"),
        axis.title = element_text(size = 14, colour = "black"),
        plot.subtitle = element_text(size = 18, colour = "black"),
        legend.text = element_text(size = 16, colour = "black"))


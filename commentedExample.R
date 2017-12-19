# This file contains a number of examples that were used in my thesis.
# Hopefully the combination of code and comments is self evident enough
# so as to enable the use of the model in other applications.

# libs
library(rstan)
library(fda)
library(MonoPoly)
library(lattice)
library(latticeExtra)
library(gridExtra)
library(MASS)
library(parallel)

# other supporting functions
# source("supFuncs.R") # probably don't want to run through these tests
source("newSupportingFunctions.R")

# data set of interest is as ever, growth from the fda package
# which we will rescale. We will not exclude any observations here.
# Just
y <- as.vector(rescale(fda::growth$hgtm, scale = list(min = 0, max = 1)))
x <- rescale(fda::growth$age, scale = list(min = 0, max = 1))

dat_df <- data.frame(x = x, y = y)

# make a factor for each person

n_indiv <- dim(fda::growth$hgtm)[2]
n_obs_per_per <- dim(fda::growth$hgtm)[1]
id_vec <- c()

for (zz in 1:n_indiv) {
  id_vec <- c(id_vec, rep(zz, n_obs_per_per))
}

dat_df$person <- id_vec

# due to the lack of ragged data support in stan, we need to create a vector
# that contains the index of the first observation for each person, and also
# has last element equal to one more than the total number of obs in it.
# the following function should generate such a vector given a sorted
pos_vec <- position_generator(dat_df$person)


# model and mcmc parameters, choose carefully! high degree polynomials take
# longer and each mcmc sample is expensive
q <- 6
mcmc_iter <- 1000
n_random_effects <- 2
n_chains <- 4

# compute the degrees of the subpolynomials
if (q %% 2 == 0) {
  k_1 <- (q / 2) - 1
  k_2 <- (q / 2) - 1
} else {
    k_1 <- (q - 1) / 2
    k_2 <- k_1 - 1
}


## posterior predictive x values, these are the thing you are probably most
## interested in as other quantites are not always unique.
n_x_new <- 69
x_new <- seq(from = min(dat_df$x), to = max(dat_df$x), length.out = n_x_new)

# choose monotonicity direction, 1 = up, -1 = down
alpha <- 1

# stack it all in a named list to pass to stan
model_dat <- list(N_x = length(dat_df$x),
                  N_i = n_indiv,
                  x_vector = dat_df$x,
                  y_vector = dat_df$y,
                  position_vector = pos_vec,
                  q = q,
                  k_1 = k_1,
                  k_2 = k_2,
                  alpha = alpha,
                  rand_params = n_random_effects,
                  N_x_new = n_x_new,
                  x_new_vector = x_new)

# fit it
model_fit <- stan(file = "BMonPolREF_V_0_4.stan",
                  data = model_dat,
                  chains = n_chains,
                  cores = parallel::detectCores(),
                  iter = mcmc_iter,
                  refresh = mcmc_iter/100,
                  control = list(max_treedepth = 13, adapt_delta = 0.825))

# slim summary

pars_of_interest <- c("sd_y", "p1_rand_ef_mean",
                      "p1_rand_ef_sd", "p_1_fixed", "p_2",
                      "beta_zero_mean", "beta_zero_sd",
                      "mu_new_pop", # these last two lines are most of the parameters
                      "y_indiv_new")# but they should probably be the thing you care about

model_convergence_summary <- print(summary(model_fit,
                                           pars = pars_of_interest,
                                           probs = c(0.1, 0.9)))$summary

## now for all plots
back_plot <- xyplot(y ~ x|person, data = dat_df,
                    type = "p", pch = "+", cex = 1.2, strip = FALSE)

# use the mean from the ppd,
y_new_samples <- extract(model_fit, "y_indiv_new")[[1]]
y_new_fitted <- apply(y_new_samples, 2:3, median)

ppd_df <- data.frame(ppd_y_new = as.vector(t(y_new_fitted)), x = x_new)

ppd_per_vec <- c()
for (i in 1:n_indiv) {
  ppd_per_vec <- c(ppd_per_vec, rep(i, n_x_new))
}
ppd_df$person_vec <- ppd_per_vec

fit_plot_2 <- xyplot(ppd_y_new ~ x|person_vec, data = ppd_df, type = "l", col = "red", lty = 2)

back_plot + fit_plot_2

# now find the prediction intervals, mostly concerned with them, could look for
# credible intervals, check other code for those though.

ppd_lower_int <- as.vector(t(apply(y_new_samples, 2:3, function(X){quantile(X, 0.025)})))
ppd_df$lower <- ppd_lower_int

ppd_upper_int <- as.vector(t(apply(y_new_samples, 2:3, function(X){quantile(X, 0.975)})))
ppd_df$upper <- ppd_upper_int




fit_plot_3 <- xyplot(lower ~ x|person_vec, data = ppd_df, type = "l", col = "red", lty = 1)
fit_plot_4 <- xyplot(upper ~ x|person_vec, data = ppd_df, type = "l", col = "red", lty = 1)

back_plot + fit_plot_2 + fit_plot_3 + fit_plot_4


## can be hard to see individual behaviour on this scale. lets look at a subset
## of individuals.

indivs_of_interest <- c(1, 7, 20, 30, 34, 38)

dat_df_small <- subset(dat_df, subset = dat_df$person %in% indivs_of_interest)
dat_df_small$person <- as.factor(as.numeric(as.factor(dat_df_small$person)))

back_plot_small <- xyplot(y ~ x|person, data = dat_df_small,
                          type = "p", pch = "+", cex = 1.2, strip = FALSE)
ppd_df_small <- subset(ppd_df, subset = ppd_df$person_vec %in% indivs_of_interest)
ppd_df_small$person_vec <- as.factor(as.numeric(as.factor(ppd_df_small$person)))

fit_plot_2_small <- xyplot(ppd_y_new ~ x|person_vec, data = ppd_df_small, type = "l", col = "red", lty = 2)
fit_plot_3_small <- xyplot(lower ~ x|person_vec, data = ppd_df_small, type = "l", col = "red", lty = 1)
fit_plot_4_small <- xyplot(upper ~ x|person_vec, data = ppd_df_small, type = "l", col = "red", lty = 1)

back_plot_small + fit_plot_2_small + fit_plot_3_small + fit_plot_4_small


# see ScriptPriorTest.R for less well commented code, where I perform the interpolation
# testing.

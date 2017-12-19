# AAM, 22/11/16


library(rstan)
library(shinystan)
library(parallel)
library(MonoPoly) # for evalPoly
library(fda)
library(lattice)
library(latticeExtra)
library(MASS)
library(bayesplot)
library(gridExtra)

source("./supFuncs.R")
source("./newSupportingFunctions.R")

setwd("~/Dropbox/RefBMonPol/REFsims/REFsims/priorTest/")
####
# Interpolation question solution (will require adjusting the manner in which we plot the points.)


#### running parameters
# polynomial order / chain length / number of "Random effects"
q <- 8
# turn down for testing.
mcmc.iterations <- 5000
rand.param.count <- 3
n.chains <- 2

# individual of interest


# predictive params
n.x.new <- 69 # this is a very judicious choice, so that some of the
# new x values line up with the data points
# NB: for individual credible intervals  this doesn't
# actually matter as we use the variation in the Betas
# anyway.

# now this is the parameter of interest.

####
## get some temporary values to compile the model with
if (q %% 2 == 0) {
  k.1 <- (q / 2) - 1
  k.2 <- (q / 2) - 1

} else {
  k.1 <- (q - 1) / 2
  k.2 <- k.1 - 1

}

x.new <- newXvalGenerator(n.x.new)
x.new <- (x.new - min(x.new)) / diff(range(x.new))

temp.data <- dataRearranger(1,0)
temp.pos.vec <- posVecGenerator(0)

temp.data.in <- list(N_x = length(temp.data$x), N_i =  length(temp.pos.vec) - 1, x_vector = temp.data$x,
                     y_vector = temp.data$y, position_vector = temp.pos.vec, q = q,
                     k_1 = k.1, k_2 = k.2, alpha = 1, rand_params = 1,
                     N_x_new = n.x.new, x_new_vector = x.new)

# prefit to be used for all models.
# in the script version this doesn't actually matter as i'm
# not sure saving 40 seconds each iteration on an iteration
# that takes 5 hours is worht the time and effort for me to get right.
ref.mod.prefit  <- stan(file = "./BMonPolREF_V_0_4.stan", data = temp.data.in, chains = 0)



num.indiv <- dim(growth$hgtm)[2]


p_yes_no_vec <- c(0.6, 0.8, 1)
p_num_obs_vec <- c(0.75, 0.95)
max_to_remove <- 25

for (tt in p_yes_no_vec) {

  for (zz in p_num_obs_vec) {

    # if (tt == 0.6 & zz = 0.75) {
    #   next
    # }

    p_yes_no <- tt
    p_num_obs <- zz

    print(tt)
    print(zz)

    all_data_list <- dataProducer(p_yes_no, p_num_obs , max_to_remove )


    dropped_data <- all_data_list[[2]][[1]]
    dropped_pos_vec <- all_data_list[[2]][[2]]


    model_dat <- list(N_x = length(dropped_data$x), N_i = num.indiv, x_vector = dropped_data$x,
                      y_vector = dropped_data$y, position_vector = dropped_pos_vec,
                      q = q, k_1 = k.1, k_2 = k.2, alpha = 1, rand_params = rand.param.count, N_x_new = n.x.new,
                      x_new_vector = x.new)

    model.fit <- stan(fit = ref.mod.prefit, data = model_dat, chains = 2, cores = 2, iter = 250, warmup = 100,
                      refresh = 10, control = list(max_treedepth = 15, adapt_delta = 0.8))


    pars.of.interest <- c("sd_y", "p1_rand_ef_mean",
                          "p1_rand_ef_sd", "p_1_fixed", "p_2",
                          "beta_zero_mean", "beta_zero_sd")
    model.convergence.summary <- print(summary(model.fit, pars = pars.of.interest, probs = c(0.1, 0.9)))$summaryz

    write.csv(x = model.convergence.summary,
              file = paste("convSummary", as.character(tt),
                           "obsDropped", as.character(zz),
                           ".csv", sep = ""))

    # plot the originial data using lattice, using the drop column to colour the points.
    # plot the fitted curves to this data (use the estimated betas FIND THIS CODE SOMEWHERE)
    # calculate the prediction intervals for each person using the new y_new_indiv vector or whatever it
    # is called, this will be of a different length so will need a different df / plot the lines seperately
    # you have code to do this
    #
    # also plot the prediction intervals for the mean fitted curve (the old way just to compare?)
    full_plot_data <- all_data_list[[1]][[1]]

    back_plot <- xyplot(y ~ x|person_vec, data = full_plot_data, type = "p", groups = full_plot_data$was_dropped,
           pch = 21, cex = 0.82, xlab = "Scaled Age", ylab = "Scaled Height")#, strip = FALSE)
    back_plot


    ## lets do an equivelent 2x2
    indiv_of_interest <- c(2,4,6,12)

    small_plot_data <- subset(full_plot_data, full_plot_data$person_vec %in% indiv_of_interest)
    small_plot_data$person_vec <- as.factor(as.numeric(as.factor(small_plot_data$person_vec)))

    back_plot_small <- xyplot(y ~ x|person_vec, data = small_plot_data, type = "p",
                              groups = small_plot_data$was_dropped,
                              pch = 21, cex = 0.52, xlab = "Scaled Age", ylim = c(0,1.1),
                              ylab = "Scaled Height", strip = FALSE)
    back_plot_small
    ## fitted values

    beta_final_samples <- extract(model.fit, "beta_final")[[1]]
    beta_final_mat <- apply(beta_final_samples, 2:3, mean)

    fitted_vec <- c()
    fitted_vec_small <- c()

    for (ii in 1:num.indiv) {
      x_temp <- subset(dropped_data, person_vec == ii)$x


      fitted_values <- evalPol(x = x_temp, beta = beta_final_mat[ii, ])
      fitted_vec <- c(fitted_vec, fitted_values)

      if (ii %in% indiv_of_interest) {
        fitted_vec_small <- c(fitted_vec_small, fitted_values)
      }

    }
    dropped_data$fitted_value <- fitted_vec
    ## not sure i actually care about this, more care about ppd for 2x2

    fit_plot <- xyplot(fitted_value ~ x|person_vec, data = dropped_data, type = "l", col = "black", lwd = 1.2)
    back_plot + fit_plot

    ## credible intervals for the mean.

    n_samples <- dim(beta_final_samples)[1]

    big_list <- list()

    for (ii in 1:num.indiv) {
      x_temp <- subset(dropped_data, person_vec == ii)$x

      temp_fitted_mat <- matrix(0, nrow = n_samples, ncol = length(x_temp))

      for (qq in 1:n_samples) {
        y_fitted <- evalPol(x = x_temp, beta = beta_final_samples[qq, ii,])
        temp_fitted_mat[qq,] <- y_fitted

      }

      big_list[[ii]] <- temp_fitted_mat

    }

    dropped_data$lower <-  unlist(lapply(big_list, function(X){apply(X, 2, function(Z){quantile(Z, 0.025)})}))
    dropped_data$upper <-  unlist(lapply(big_list, function(X){apply(X, 2, function(Z){quantile(Z, 0.975)})}))

    fit_plot_5 <- xyplot(lower ~ x|person_vec, data = dropped_data, type = "l", col = "black", lwd = 1.2)
    fit_plot_6 <- xyplot(upper ~ x|person_vec, data = dropped_data, type = "l", col = "black", lwd = 1.2)

    back_plot + fit_plot + fit_plot_5 + fit_plot_6

    # ppd
    # This was stupid because it used x values way outside of the data range.
    y_new_samples <- extract(model.fit, "y_indiv_new")[[1]]
    y_new_fitted <- apply(y_new_samples, 2:3, median)

    ppd_df <- data.frame(ppd_y_new = as.vector(t(y_new_fitted)), x = x.new)

    ppd_per_vec <- c()
    for (i in 1:num.indiv) {
      ppd_per_vec <- c(ppd_per_vec, rep(i, n.x.new))
    }
    ppd_df$person_vec <- ppd_per_vec

    fit_plot_2 <- xyplot(ppd_y_new ~ x|person_vec, data = ppd_df, type = "l", col = "red", lty = 2)

    back_plot + fit_plot + fit_plot_2


    ## credible intervals

    ppd_lower_int <- as.vector(t(apply(y_new_samples, 2:3, function(X){quantile(X, 0.025)})))
    ppd_df$lower <- ppd_lower_int

    ppd_upper_int <- as.vector(t(apply(y_new_samples, 2:3, function(X){quantile(X, 0.975)})))
    ppd_df$upper <- ppd_upper_int




    fit_plot_3 <- xyplot(lower ~ x|person_vec, data = ppd_df, type = "l", col = "red", lty = 1)
    fit_plot_4 <- xyplot(upper ~ x|person_vec, data = ppd_df, type = "l", col = "red", lty = 1)

    back_plot + fit_plot_2 + fit_plot_3 + fit_plot_4

    ppd_df_small <- subset(ppd_df, ppd_df$person_vec %in% indiv_of_interest)
    ppd_df_small$person_vec <- as.factor(as.numeric(as.factor(ppd_df_small$person_vec)))

    fit_plot_2_small <- xyplot(ppd_y_new ~ x|person_vec, data = ppd_df_small,
                         type = "l", col = "red", lty = 2)
    fit_plot_3_small <- xyplot(lower ~ x|person_vec, data = ppd_df_small,
                               type = "l", col = "red", lty = 1)
    fit_plot_4_small <- xyplot(upper ~ x|person_vec, data = ppd_df_small,
                               type = "l", col = "red", lty = 1)

    back_plot_small + fit_plot_2_small + fit_plot_3_small + fit_plot_4_small

    pdf(file = "persons_2-4-6-12_all_dropped.pdf", width = 5.5, height = 4.5)
    print(
      back_plot_small + fit_plot_2_small + fit_plot_3_small + fit_plot_4_small
    )
    dev.off()
    ############################

    back_plot + fit_plot_2 + fit_plot_3 + fit_plot_4 +  fit_plot_5 + fit_plot_6


    ######## covereage by the predictive distribution.

    dropped_points_df <- subset(full_plot_data, subset = (was_dropped == TRUE))
    dropped_people <- unique(dropped_points_df$person_vec)

    num_dropped_obs <- length(dropped_points_df$x)
    was_predicted_vec <- c()

    for (qq in 1:num_dropped_obs) {
      temp_df <- dropped_points_df[qq, ]

      ppd_temp_df <- subset(ppd_df, subset = (x == temp_df$x & person_vec == temp_df$person_vec))

      res <- (temp_df$y < ppd_temp_df$upper) & (temp_df$y > ppd_temp_df$lower)

      was_predicted_vec <- c(was_predicted_vec, res)
    }

    dropped_points_df$was_predicted <- was_predicted_vec

    mean(dropped_points_df$was_predicted)

    # write out % data missing vs coverage

    write.csv(data.frame(missing_data_proportion = num_dropped_obs / length(full_plot_data$x),
                         coverage_proportion = mean(dropped_points_df$was_predicted)),
              file = paste0("./coverages/pYesNo", as.character(p_yes_no), "pNumObs", as.character(p_num_obs),
                            "maxRemoved", as.character(max_to_remove), ".csv"))
  }
}

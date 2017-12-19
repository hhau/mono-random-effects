# supporting functions for the interpolation question
library(fda)

dataProducer <- function(p_yes_no, p_num_obs, max_to_remove) {
  # this function produces two data frames (returns in a list)
  # first is the dataframe with the observations removed according to the method specified in the notes
  # second is the full data set with an extra column indicating if that particular observation was
  # removed or not.

  y <- as.vector(growth$hgtm)
  y <- (y - min(y))/diff(range(y))

  x <- growth$age
  x <- (x - min(x)) / diff(range(x))

  full_data <- data.frame(y = y, x = x)


  n_indiv <- dim(growth$hgtm)[2]
  n_obs_per_indiv <- length(growth$age)

  yes_no <- rbinom(n = n_indiv, size = 1, prob = p_yes_no)
  num_remove <- rbinom(n = n_indiv, size = max_to_remove, prob = p_num_obs)
  num_obs_to_remove <- yes_no * num_remove

  removal_to_start <- sample(x = 1:(n_obs_per_indiv - max_to_remove - 2), size = n_indiv, replace = TRUE) + 2
  pos_vec <- c(1, 1:(n_indiv) * n_obs_per_indiv + 1)

  drop_vec <- 1:length(y)

  person_vec <- c()
  for (i in 1:n_indiv) {
    person_vec <- c(person_vec, rep(i, 31))
  }

  full_data$person_vec <- person_vec

  for (ii in 1:n_indiv) {

    starting_index <- pos_vec[ii] + removal_to_start[ii] - 1
    number_to_drop <- num_obs_to_remove[ii]

    if (number_to_drop == 0) {
      next
    }

    stopping_index <- starting_index + num_obs_to_remove[ii] - 1

    drop_vec[starting_index:stopping_index] <- -drop_vec[starting_index:stopping_index]

  }

  full_data$drop_vector <- drop_vec
  full_data$was_dropped <- drop_vec < 0

  row_dropper <- as.numeric(full_data$was_dropped) * drop_vec
  dropped_data <- full_data[row_dropper, ]

  rownames(dropped_data) <- 1:nrow(dropped_data)

  new_pos_vec <- pos_vec - c(0, cumsum(num_obs_to_remove))

  res_list <- list(list(full_df = full_data), list(dropped_df = dropped_data, dropped_person_vec = new_pos_vec ))

  return(res_list)


}

# this should handle both vectors and matricies
rescale <- function(input,
                    scale = list(min = 0,
                                 max = 1),
                    ...) {

  res <- ((input - min(input)) * (scale$max - scale$min) / diff(range(input))) + scale$min
  res
}

# kinda assumes the id vector is sorted appropriately
# and starts at 1 and goes up to n_indivs
position_generator <- function(id_vec, ...) {
  n_obs <- length(id_vec)
  pos_vec <- c(1)
  for (ii in 1:(n_obs - 1)) {
    if (id_vec[ii] == (id_vec[ii + 1] - 1)) {
      pos_vec <- c(pos_vec, ii + 1)
    }
  }
  pos_vec <- c(pos_vec, n_obs + 1)
  return(pos_vec)
}

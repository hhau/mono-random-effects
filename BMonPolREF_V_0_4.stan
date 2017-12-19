// this is a stan file that should fit monotonic polynomials with random effects.
// for the moment this is only intended to fit polynomials on the [0, infty) space
// and only those that are monotonically increasing?.
// AAM 18/7/16

functions {
  // we will still need the convolve command here, and the
  // horner evaluation of the polynomials, however they might need adjusting
  // although we can probably do it without adjusting ? depends how smart
  // we want to be
  vector convolve(vector a, vector b){
    int na;
    int nb;
    int nc;
    vector[rows(b)] rev_b;
    vector[rows(a)+rows(b)-1] c;

    na = rows(a);
    nb = rows(b);
    nc = na+nb-1;
    for(i in 1:nb){
      rev_b[i] = b[nb-i+1];
    }
    for(j in 1:nc){
      int istart;
      int istop;
      istart = max(0, j-nb) + 1;
      istop = min(na, j);
      c[j] = sum(a[istart:istop] .* rev_b[(nb-j+istart):(nb-j+istop)]);
    }
    return c;
  }
   // we can make this horner function considerably smarter later
   // and probably avoid a loop or two.


  vector horner(vector x, vector beta){
    int nb;
    vector[rows(x)] res;

    nb = rows(beta);
    res = rep_vector(beta[nb], rows(x));
    for(i in 1:(nb-1)){
      res = res .* x + beta[nb-i];
    }
    return res;
  }
}

data {
  // we're going to take data in a vector, to allow for non-square data arrays.
  // total number of observations
  int <lower=0> N_x;

  // number of individuals
  int <lower=0> N_i;

  // vectors of x values and y values
  vector [N_x] x_vector;
  vector [N_x] y_vector;

  // vector of positions that has the starting points of each of the
  // individuals data runs. So if the first data run is observations
  // 1 to 31, the position_vector will have elements [1, 32, ...]
  // it will also have one more than the final value (i.e N_x + 1) in it
  int <lower=1>  position_vector[N_i + 1];

  // now we need to take the polynomial parameters
  // overall polynomial order
  int <lower=0> q;

  // sub polynomial order.
  int <lower=0> k_1;
  int <lower=0> k_2;

  // this is always going to be one for our current purpose.
  int alpha;

  // how many "bits" \ coefficients of one of the underlying vectors are
  // going to be random. if it is k + 1, that means all of the bits are random.
  // in which case all of the code is going to break.
  // maybe we have to have at least one fixed part?
  int <lower=0, upper=(k_1 + 1)> rand_params;

  // stuff for overall mean prediction  invervals
  int <lower=0> N_x_new;
  vector [N_x_new] x_new_vector;

}

transformed data{
  // we might need some temporary vectors to store intermediate conversations.
  // calculate the lengths for the temporary vectors?

  // // for the "final" gammas
  // // int <lower=0> gamma_1_length;
  // // int <lower=0> gamma_2_length;

  // // for the gamma after it is convolved with itself but before it is convolved with (0,1) ?
  // // actually as we can nest convolve commands i don't need this
  // int <lower=0> gamma_2_length;



  // construct vector equivalent of convolution with (x - a)
  vector[2] lower_bound_vector;
  lower_bound_vector[1] = 0;
  lower_bound_vector[2] = 1;

  // if(q % 2 == 0) {
  //   gamma_1_length = (2 * k_1) - 1;
  //   gamma_2_length = (2 * k_2) - 1;
  // } else {
  //     gamma_1_length = (2 * (k_1 + 1)) - 1;
  //     gamma_2_length = (2 * k_2) - 1;
  // }

}

parameters {
  // as is, this is still an overall quantity, might have to change it
  // to be on a per individual level in the future.
  real <lower=0> sd_y;

  // we are going to need some hierachical variances as well
  // also some hierachical means

  // underlying polynomials (p_1, p_2)

  // these could do with priors
  vector [N_i] beta_zero;
  real beta_zero_mean;
  real<lower=0> beta_zero_sd;


  matrix [N_i, rand_params] p_1_rand;
  vector[rand_params] p1_rand_ef_mean;
  vector <lower=0> [rand_params] p1_rand_ef_sd;

  // these are fixed parameters
  vector [(k_1 + 1) - rand_params] p_1_fixed;
  vector [k_2 + 1] p_2;
}

transformed parameters {
  // these two big quantities should now be matricies
  matrix[N_i, q + 1] beta_final;

  // These are always going to be different lengths.
  // check the other code for how to add the zero correctly
  matrix[N_i, (2*k_1) + 1] gamma_1;
  matrix[N_i, (2*k_2) + 2] gamma_2; //

  // when i add the above two variables,  i going to need to zero pad
  matrix[N_i, q] gamma;

  // going to need a mu vector to hold the means appropriately
  vector[N_x] mu;

  mu = rep_vector(0, N_x);
  // need container for the combo of the fixed and random parameters

  // there has to be some loop in here, at least until i figure out what
  // the vector wise operation is

  for(ii in 1:(N_i)) {
    // create appropriate p_1 for each indiviudal
    // might have to do this outside the loop and overwrite it each time
    vector[k_1 + 1] temp_p1;

    //vector[position_vector[ii + 1] - position_vector[ii]] temp_mu;
    int lower_loop;
    int upper_loop;

    lower_loop = position_vector[ii];
    upper_loop = position_vector[ii + 1] - 1;


    // set entires in matrix to zero (they aren't initalised to be)
    gamma_1[ii] = rep_row_vector(0, (2*k_1) + 1);
    gamma_2[ii] = rep_row_vector(0, (2*k_2) + 2);

    temp_p1 = append_row(p_1_rand[ii]', p_1_fixed);

    // put in correct rows of 'temporary' matrix gamma_i
    gamma_2[ii] = convolve(convolve(p_2, p_2), lower_bound_vector)';
    gamma_1[ii] = convolve(temp_p1, temp_p1)';

    // add the zero correctly (we append it to the end like before)
    if (q % 2 == 0) {
      gamma[ii] = (append_row(gamma_1[ii]',rep_row_vector(0,1)') + gamma_2[ii]')';
    } else {
        gamma[ii] = (append_row(gamma_2[ii]', rep_row_vector(0,1)') + gamma_1[ii]')';
    }

    beta_final[ii, 1] = beta_zero[ii];

    for(qq in 1:q) {
      beta_final[ii, qq + 1] = alpha *  (gamma[ii, qq] / qq);
    }

    mu[lower_loop:upper_loop] = horner(x_vector[lower_loop:upper_loop], beta_final[ii]');

  }

}

model {
  // loop over each individuals polynomial
  for (ii in 1:(N_i)) {
    // need to delcare y, mu will be declared elsewhere
    int loop_lower;
    int loop_upper;
    loop_lower = position_vector[ii];
    loop_upper = position_vector[ii + 1] - 1;
    y_vector[loop_lower:loop_upper] ~  normal(mu[loop_lower:loop_upper], sd_y);

    p_1_rand[ii] ~ normal(p1_rand_ef_mean, p1_rand_ef_sd);


  }

  beta_zero ~ normal(beta_zero_mean, beta_zero_sd);
  p1_rand_ef_sd ~ normal(0, 0.1);

  // 4/12/17 note
  // There are some quantites here i would really like to put priors on in retrospect.
  // Firstly a weakly informative prior on p1_rand_ef_mean, and beta_zero_mean

  p1_rand_ef_mean ~ normal(0, 100);
  beta_zero_mean ~ normal(0, 100);

  // Next, something with a bit less mass in the tails than the flat prior for
  // sd_y and beta_zero_sd, but only marginally less so as I'm really not sure
  // how sensitive the model is to these parameters.

  sd_y ~ normal(0,1) T[0,];
  beta_zero_sd ~ normal(0,1) T[0,];

  // Finally, some diffuse priors for computational reasons on the remaining
  // parameters

  p_1_fixed ~ normal(0, 100);
  p_2 ~ normal(0, 100);

  // also a note that the specificaition is non-unique, and hence we can have
  // chains that disagree considerably, with poor r^hat values, but all fit the
  // data reasonably well. Hence the convergence metrics of interest are really
  // the fitted values and the posterior predctive means. Maybe also the
  // variance paramter(s).

}

generated quantities {

  // containers needed for population level inference
  vector [q + 1] beta_final_pop;
  vector [(2*k_1) + 1] gamma_1_pop;
  vector [(2*k_2) + 2] gamma_2_pop; //
  vector[k_1 + 1] temp_p1_pop;
  vector[q] gamma_pop;
  vector[N_x_new] mu_new_pop;

  //vector[N_x_new] y_new_pop;

  matrix [N_i, N_x_new] mu_new;
  matrix [N_i, N_x_new] y_indiv_new;

  for (ii in 1:(N_i)) {
    mu_new[ii] = horner(x_new_vector, beta_final[ii]')';
    y_indiv_new[ii] = multi_normal_cholesky_rng(mu_new[ii]', diag_matrix(rep_vector(sd_y, N_x_new)))';
  }



  temp_p1_pop = append_row(p1_rand_ef_mean, p_1_fixed);
  gamma_2_pop = convolve(convolve(p_2, p_2), lower_bound_vector);
  gamma_1_pop = convolve(temp_p1_pop, temp_p1_pop);

  if (q % 2 == 0) {
    gamma_pop = (append_row(gamma_1_pop, rep_vector(0,1)) + gamma_2_pop);
  } else {
      gamma_pop = (append_row(gamma_2_pop, rep_vector(0,1)) + gamma_1_pop);
  }

  beta_final_pop[1] = beta_zero_mean;

  for(ii in 2:(q + 1)) {
    beta_final_pop[ii] = alpha * (gamma_pop[ii - 1] / (ii - 1));

  }

  mu_new_pop = horner(x_new_vector, beta_final_pop);

  // i'm fine with everything above, but if i want to use this formulation to
  // get a handle on the prediction intervals, then i will need the correct
  // variance at the population level, which may or may not be avaliable. Unsure


}



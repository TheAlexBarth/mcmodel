###
# Example testing script
##
# devtools::document()
# devtools::check()
# devtools::load_all()


y = mtcars$mpg
X = cbind(1, scale(mtcars$wt), scale(mtcars$hp)) # scaling is typically helpful


###
### Assuming a normal linear regression model
### With unknown variance and mean
### y_i ~ N(mu_i, sigma^2)
### mu_i = X_i %*% beta

# set up the data
data_list = list(
  y = y,
  X = X
)

# this is for indexing
const_list = list(
  K_x = ncol(X)
)

# needs initialization for model parameters
# initialization can be informed by data to ensure 
init_list = list(
  beta = c(20, rep(0, ncol(X)-1)),
  sigma = 6
)
init_list$mu = data_list$X %*% init_list$beta

# need to specify priors
# In this example:
# [beta] = N(mu_beta, sd_beta)
# [sigma^2] = IG(a_sig2, b_sig2)
prior_pars = list(
  mu_beta = c(10, 0, 0),
  sd_beta = c(3, 2, 2),
  a_sig2 = 2,
  b_sig2 = 2
)

# MH update for beta and Gibbs of sigma
# needs a prop_sd for beta
prop_sd = list()
prop_sd$beta = rep(1, const_list$K_x)

# this is an individual-wise update function
update_beta = function(prop_sd) {
  # updates for MH steps require a "acc" object to feed back into the algo.
  acc = rep(0, const_list$K_x)
  # initilalze acc for return on HM
  for(k in 1:K_x) {
    # in this example, minimal changes to beta are introduced
    # this is useful for large datasets to not recalc full matrix
    betak_star = rnorm(1, mean = beta[k], sd = prop_sd$beta[k])
    mu_star = mu + X[,k] * (betak_star - beta[k])

    mh_1 = sum(dnorm(y, mean = mu_star, sd = sigma, log = TRUE)) +
      dnorm(betak_star, mu_beta[k], sd_beta[k], log = TRUE)

    mh_2 = sum(dnorm(y, mean = mu, sd = sigma, log = TRUE)) +
      dnorm(beta[k], mu_beta[k], sd_beta[k], log = TRUE)

    # on update, etc
    if(log(runif(1)) < (mh_1-mh_2)) {
      beta[k] = betak_star
      acc[k] = 1
      mu = mu_star
    }
  }
  return(
    # the return structure of MH blocks MUST BE
    # a list with a `val` list and `acc` list
    list(
      val = list(
        # update to the global function environment
        beta = beta,
        mu = mu
      ),
      acc = list(
        beta = acc
      )
    )
  )
}

# conjugate prior is inverse gamma
# full conditional then  sigma^2 | a_sigma, b_sigma \propto IG(a+(n/2), b+0.5*SUM(y_i-mu_i))
update_sigma = function() {
  new_sig2 = 1/rgamma(1, a_sig2 + length(y)/2, b_sig2 + 0.5 * sum((y-mu)^2))
  return(
    list(
      val = list('sigma' = sqrt(new_sig2))
    )
  )
}

# mcmc_test = mcmc_run(
#   data_list = data_list,
#   const_list = const_list,
#   init_list = init_list,
#   prior_pars = prior_pars,
#   prop_sd = prop_sd,
#   update_functions = list(
#     'beta' = update_beta,
#     'sig' = update_sigma
#   ),
#   burn_adapt = TRUE,
#   save_names = c('beta','sigma'),
#   n_bur = 5000,
#   n_iter = 1e5,
#   thin = TRUE
# )


# system.time({mcmc_test1 = mcmc_run_internal(
#   data_list = data_list,
#   const_list = const_list,
#   init_list = init_list,
#   prior_pars = prior_pars,
#   prop_sd = prop_sd,
#   update_functions = list(
#     'beta' = update_beta,
#     'sig' = update_sigma
#   ),
#   burn_adapt = FALSE,
#   save_names = c('beta','sigma'),
#   n_burn = 5000,
#   n_iter = 1e5,
#   thin = 2
# )})



##########
# Check Diagnostics
########

# plot_post_chain(mcmc_test$samples, 'beta')

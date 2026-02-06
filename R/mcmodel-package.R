#' McModel
#' 
#' Tools for constructing Gibbs-MH samplers in base R with adaptive
#' MCMC algorithm support and basic diagnostics.
#' 
#' @docType package
#' @keywords internal
"_PACKAGE"


#' MCMC Run Object Structures
#' 
#' Extra details regarding the core arguments in 
#' \code{\link{mcmc_run}} construction. Examples here-in show a simple case of a GIBBS + MH sampler for
#' a normal regression model with two covariates and unknown mean. Additionally details are shown for
#' simulating derived quantities with a bayesian p-value based on MSE.
#' 
#' # Details 
#'
#' ##  `data_list`
#' A list containing all data used in the u/d functions that are non-mutable
#' but necessary for function runs (e.g., response and covariate data). List terms
#' MUST match the variable names used in fuctions
#' 
#' 
#' ##  `const_list`
#' A list object with named terms used in indexing.
#' 
#' Note that functionally terms provided in either const_list or data_list can be interchangable
#' the distinction is inspired by separating components of the model, as can be done in NIMBLE
#' or similarly through STAN.
#' 
#' 
#' ##  `init_list`
#' A list with ALL model parameters or latent states (both stochastic and deterministic nodes) which
#' must be accessible for update functions or saving. Note that if a term is definied within an u/d fn
#' it does not need to be accessible (however it then is not available for saving).
#' 
#' If running multiple chains, you must provide an init_generator function which mutates 
#' the centralized init_list. see below for details...
#' 
#' ## `prop_sd`
#' A list of values to be used in proposal functions. Typically assumes a normal proposal distribution
#' However theoritically possible to be used in an alternative approach.
#' 
#' @name mcmc_structures
#' @docType data
#' @keywords internal
#' 
NULL
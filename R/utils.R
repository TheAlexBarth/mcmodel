#############
#MARK: ADAPT TUNER
#############
#' Default proposal adaptation
#' 
#' By default, prop_sd will be updated for MH blocks based on a double exponential
#' targetting a 0.17 to 0.5 acceptance rate. Note this current set up allows for 
#' bulk updates of multiple parameters if any fall outside the goal range
#' 
#' @param term a character string for name of parameter to tune within acc_counter
#' @param acc_counter the list object with acceptance counts
#' @param prop_sd list of proposal sds
#' @param interval_length length of adaptation interval (default should be 100)
#' 
#' @return Updated proposal SD for this term
#' @keywords internal
adapt_tuner_double_exp <- function(term, acc_counter, prop_sd, interval_length) {
  acc_rate <- acc_counter[[term]] / interval_length
  if(any(acc_rate < 0.17 | acc_rate > 0.50)){
    acc_update <- which(acc_rate < 0.17 | acc_rate > 0.50)
    prop_sd[[term]][acc_update] <- prop_sd[[term]][acc_update] * exp(2*acc_rate[acc_update]-0.33)
    cat("\n", paste0('Updating acceptance rate for ',length(acc_update)," ",term, 's mean rate: ', mean(acc_rate)))
  }
  return(prop_sd[[term]])
}

#' No adaptation in proposal tuning
#' 
#' This is just to keep code smooth but allow for no_adaptation 
#' 
#' @inheritParams adapt_tuner_double_exp
#' @return original proposal sd for this term
#' @keywords internal
no_adapt = function(term, acc_counter, prop_sd, interval_length) {
    return(prop_sd[[term]])
}


#################
#MARK: Trim Chain
#################
#' Trim Chains
#' 
#' Takes a save object (or derv save) and trims to a specified index point
#' 
#' @param save_obj a list of paramters that are saved accross iterations
#' @param trim_point single index integer at which a trim should extend.
#' @param front Logical to keep front DEFAULT TRUE.
#' 
#' @return a list with the trimmed save object
#' @keywords internal
trim_chain = function(save_obj, trim_point, front = TRUE) {
    
    for(param in names(save_obj)) {
        max_len = dim(save_obj[[param]])[1]
        if(is.null(max_len)) max_len = length(save_obj[[param]])
        if(front) {
            from = 1
            to = trim_point
        } else {
            from = max_len-trim_point+1
            to = max_len
        }
        if(length(dim(save_obj[[param]])) != 2) {
           save_obj[[param]] = gafa(save_obj[[param]], from:to)
        } else {
            save_obj[[param]] = save_obj[[param]][from:to,]
        }
    }
    return(save_obj)
}

##################
#MARK: Assess Burn
##################
#' Assess Burnin
#' 
#' Loop through the save object to assess if the burnin meets 
#' some convergence criteria. Default is to use geweke_check but other
#' functions can be created.
#' 
#' @param save_obj a list of save objects that IS TRIMMED!
#' @param chain_check a function for checking chains - default should be geweke_check
#' 
#' @return boolean if it passess the checker
#' @keywords internal
assess_burnin = function(save_obj, chain_check) {
    passes = logical()
    for(param in names(save_obj)) {
        arr = save_obj[[param]]
        chain_mat = fafa(arr)
        for(j in 1:ncol(chain_mat)) {
            check = chain_check(chain_mat[,j])
            passes = c(passes, check)
        }
    }
    if(all(passes)){
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#' Geweke convergence diagnostic
#' 
#' A comparison of means between the first 10th of a chain and the latter half
#' The proper Geweke diagnostic uses a z-test of means but I used a t-test
#' 
#' @param chain a single vector of an MCMC chain
#' 
#' @examples
#' geweke_check(rnorm(1000))
#' 
#' @return boolean
#' 
#' @references
#' Geweke, J. 1992. Evaluting the accuracy of sampling-based approaches to the calculations 
#' of posterior moments. Bayesian statistics. 4: 641-649.
#' 
#' @importFrom stats t.test
#' 
#' @export
geweke_check = function(chain) {
    first_quant = chain[1:round(length(chain)*0.1)]
    latter_half = chain[round(length(chain)*0.5):length(chain)]
    if(length(first_quant) < 10 | length(latter_half) < 10) {
        return(FALSE) # insufficient chain lengths
    }
    if(var(first_quant) == 0 | var(latter_half) == 0) {
        return(FALSE)
    }
    p = t.test(first_quant, latter_half)$p.value
    if(p > 0.05) {
        return(TRUE)
    } else {
        return(FALSE)
    }
}

#########################
#MARK: Assign Along First Axis
#########################
#' Along First Axis functions
#' 
#' These series of helper functions provide an extremely 
#' valuable approach to interact with the first axis of any shaped array.
#' For MCMC applications where the first axis is iterations of a saved value.
#' A- assign (replace values); G- get (index values); F- flatten (return set of chains)
#' 
#' @name afas
NULL

#' Asign Along First Axis
#' 
#' 
#' @param arr array or vector
#' @param idx index to get on first axis
#' @param val value to assign at index
#' 
#' @return array with changes at idx
#' 
#' @rdname afas
aafa = function(arr, idx, val) {
    d = dim(arr)
    if(is.null(d)) {
        arr[idx] = val
    } else {
        ldx = c(list(idx), rep(list(TRUE), length(dim(arr))-1))
        arr = do.call('[<-', c(list(arr), ldx, list(val)))
    }
    return(arr)
}

#' Get Along First Axis
#' 
#' @param arr array or vector
#' @param idx index to get on first axis
#' 
#' @rdname afas
gafa = function(arr, idx) {
    d = dim(arr)
    if(is.null(d)) {
        arr = arr[idx]
    } else {
        ldx = c(list(idx), rep(list(TRUE), length(d)-1)) # list index
        arr = do.call(
            '[', c(list(arr), ldx, list(drop = FALSE))
        )
    }
    return(arr)
}

#' Flatten Along First Axis
#' 
#' @param arr array or vector
#' @rdname afas
fafa = function(arr) {
    d = dim(arr)
    if(!is.null(d)) {
        arr = matrix(arr, nrow = d[1])
    }
    return(as.matrix(arr))
}


############
# MARK: PSRF 
############
#' Potential Scale Reduction Factor (Rhat)
#' 
#' Functions to calculate rhat. The default used in \link{mcmc_run} is the
#' rank normalized, split chain R-hat described by Vehtari et al (2021). However,
#' there is also the ability to calculate the classic Gelman-Rubin Diagnostic.
#' 
#' @references
#' Vehtari A. Gelman A. Simpson D. Carpenter B. B ̈urkner P-C. 2021.
#'   Rank-Normalization, Folding, and Localization:  An Improved ̂R for Assessing Convergence of  MCMC.
#'   Bayesian Analysis 16(2): 667-718. 10.1214/20-BA1221.
#' Gelman, A. and Rubin, D. B. (1992). 
#'   “Inference from iterative simulation using multiple sequences (with discussion).” 
#'   Statistical Science, 7(4): 457–511. 667, 668, 671, 675
#' 
#' @name rhat
NULL

#' Gelman-Rubin Rhat
#' 
#' Classic approach which REQUIRES a chain list with equal length chains.
#' This is not as built out as the veharti approach and thus cannot support multidimensional output
#' 
#' @param chain_list a list of MCMC Chains
#' 
#' @examples
#' n = 1000
#' gelman_rubin_psrf(list(rnorm(n), rnorm(n))) # should be close to 1
#' gelman_rubin_psrf(list(rnorm(n), rnorm(n, 100))) #should be large
#' 
#' @returns numeric score of rhat
#' 
#' @importFrom stats var
#' 
#' @rdname rhat
#' @export
gelman_rubin_psrf = function(chain_list) {
    ns = sapply(chain_list, length) |>
        unique()

    if(length(ns) > 1) stop('Gelman-Rubin Stat requires all chains equal length')
    chain_means = sapply(chain_list, mean)
    mean_of_mean = mean(chain_means)
    B = var(chain_means)
    W = mean(sapply(chain_list, var))
    
    Rhat = ((ns-1)/ns * W + (1/ns)*B) / W

    return(Rhat)
}

#' Split Chain, Rank Normalized R-hat
#' 
#' Based on Vehtari et al 2021; 10.1214/20-BA1221. This version supports both vector chains or multi-dimension arrays (assuming iterations along first axis)
#' This script is a wrapper which calls \link{vsp_internal}
#' 
#' @param chains either a list of chains or vector of one chain (which will be split)
#' @param rank_norm boolean for rank normalization
#' 
#' @examples
#' 
#' n = 1000
#' vehtari_split_psrf(rnorm(n))
#' vehtari_split_psrf(list(rnorm(n), rnorm(n)), rank_norm = FALSE)
#' 
#' @returns vector of rhat score(s)
#' 
#' @export
#' @rdname rhat
vehtari_split_psrf = function(chains, rank_norm = TRUE){
    #need to check dimension sizing
    if(is.list(chains)) {
        dims = unique(sapply(chains, function(x) length(dim(x))))
        chains = lapply(chains, fafa)
        mat_len = unique(sapply(chains, function(x) dim(x)[2]))
        chains = lapply(1:mat_len, function(x) lapply(chains, '[', ,x))
    } else {
        chains = fafa(chains)
        dims = dim(chains)
        mat_len = dim(chains)[2]
        chains = lapply(1:mat_len, function(x) chains[,x])
    }

    split_r = sapply(1:mat_len, function(x) vsp_internal(chains[[x]], rank_norm = rank_norm))
    return(split_r)
}

#' Interal features for Vehtari function
#' 
#' @param chains a single chains (e.g., like the gelman rubin)
#' @param rank_norm boolean to ranknoramlize for not
#' 
#' @importFrom stats qnorm
#' 
#' @returns numeric of rhat values
#' @keywords internal
vsp_internal = function(chains, rank_norm) {
    if(is.list(chains)) {
        split_chains = chains |> lapply(
            function(x) split(x, cut(seq_along(x), 2, labels = FALSE))
        ) |>
        do.call(what = c,)
    } else {
        split_chains = split(chains, cut(seq_along(chains), 2, labels = FALSE))
    }

    M = length(split_chains)
    N = unique(sapply(split_chains, length))
    if(rank_norm) {
        r = rank(unlist(split_chains, use.names = FALSE))
        z = qnorm((r-(3/8))/(N*M + 1/4))
        split_chains = split(z, cut(seq_along(z), M, labels = FALSE))
    }
    # naming references chain-indexing
    dot_m = sapply(split_chains, mean)
    dot_dot = mean(dot_m)
    B = (N/(M-1)) * sum((dot_m - dot_dot)^2) #between var
    
    s_sq_m = sapply(split_chains, function(x) (1/(N-1))*sum((x-mean(x))^2))
    W = (1/M) * sum(s_sq_m) #within var

    varhat = ((N-1)/N)*W + B/N

    split_r = sqrt(varhat / W)
    return(split_r)
}

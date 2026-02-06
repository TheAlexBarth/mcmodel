require(coda)

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
#' 
#' @return a list with the trimmed save object
#' @keywords internal
trim_chain = function(save_obj, trim_point) {
    for(param in names(save_obj)) {
       save_obj[[param]] = gafa(save_obj[[param]], 1:trim_point)
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
#' @return boolean
#' @keywords internal
geweke_check = function(chain) {
    first_quant = chain[1:round(length(chain)*0.1)]
    latter_half = chain[round(length(chain)*0.5):length(chain)]
    if(length(first_quant) < 10 | length(latter_half) < 10) {
        return(FALSE) # insufficient chain lengths
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
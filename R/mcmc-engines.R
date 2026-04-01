#' Internal dispatch for a single chain run
#' 
#' @inheritParams mcmc_run
#' @param prop_adapt integer interval for proposal adaptation. Default is 100
#' @param adapt_tuner function for proposal adaptation. Default is \link{adapt_tuner_double_exp}. no_adapt is also a built-in option
#' @param chain_check function to check chain convergence within adaptive burnin feaure.
#' 
#' @return list object with mcmc features
#' 
#' @importFrom utils txtProgressBar setTxtProgressBar
#' 
#' @keywords internal
mcmc_run_internal = function(
    data_list,
    const_list,
    init_list,
    prior_pars,
    prop_sd,
    update_functions,
    save_names,
    derv_quants,
    derived_functions,
    n_iter,
    thin,
    n_burn,
    burn_adapt,
    min_burn,
    prop_adapt = 100,
    adapt_tuner = adapt_tuner_double_exp,
    chain_check = geweke_check
) {
   #region \- load params -----------------------
    list2env(data_list, envir = environment())
    list2env(const_list, envir = environment())
    list2env(init_list, envir = environment())
    list2env(prior_pars, envir = environment())
    for(i in 1:length(update_functions)) {
        environment(update_functions[[i]]) = environment()
    }
    environment(format_output) = environment()
    environment(check_function_format) = environment()
    #endregion ----------------------------------

    #region \- Save construct ------------------------------
    if(is.null(save_names)) save_names = names(init_list)
    max_save_size = (n_iter)/thin
    # for speed on indexing/assignment, every thing is stored as a matrix
    # then reshpaed if needed later
    plen = list()
    pd = list()
    save = list()
    if(burn_adapt) run_save = list()
    for(param in save_names) {
        pd[[param]] = dim(init_list[[param]])
        plen[[param]] = length(init_list[[param]])

        save[[param]] = array(NA, dim = c(max_save_size, plen[[param]]))
        if(burn_adapt) run_save[[param]] = array(NA, dim = c(min_burn+1, plen[[param]]))
    }
    #endregion---------------------------------------------

    #region \- d.quants ----------------------------------------
    # currently only supports one-dimensional quantities
    if(!is.null(derv_quants)) {
        if(is.null(derived_functions)) {
            stop('Need to specify a derived_functions list')
        } else {
            for(i in 1:length(derived_functions)) {
                environment(derived_functions[[i]]) = environment()
            }
        }
        derv_save = list()
        for(dname in derv_quants) {
            assign(dname, 1L)
            derv_save[[dname]] = array(NA, dim = max_save_size)
        }
    }
    #endregion -------------------------------------------

    #region \- set acc counter if needed ----------------
    if(!is.null(prop_sd)) {
        acc_counter = reset_acc_counter(prop_sd)
    }
    #endregion ---------------------------------------

    #region \- Check update/dev functions --------------
    for(fn in update_functions) {
        check_function_format(fn)
    }
    if(!is.null(derived_functions)) {
        for(fn in derived_functions) {
            check_function_format(fn)
        }
    }
    #endregion -------------------------------------------

    # run trackers
    pb = txtProgressBar(0, n_iter, style = 3)
    pbiter = 0
    save_count = 0
    cat('\n', 'Launching Run')
    in_burnin = TRUE
    for(iter in 1:n_burn) {
        pbiter = pbiter + 1
        setTxtProgressBar(pb, pbiter)

        #region \- Run update block -------------------
        for(i in 1:length(update_functions)) {
            update_fn = update_functions[[i]]
                if('prop_sd' %in% names(formals(update_fn))) {
                    res = update_fn(prop_sd)
                    # only need to update here if in MH
                    for(cpar in names(res$acc)) {
                        acc_counter[[cpar]] = acc_counter[[cpar]] + res$acc[[cpar]]
                    }
                } else {
                    # for gibbs functions with no prop_sd
                    res = update_fn()
                }
                # this loop hits regardless of if it were an MH or Gibbs
                for(par in names(res$val)) {
                    assign(par, res$val[[par]])
            }
        }
        #endregion----------------------------------


        #region \- check burn ------------------
        # burn-in behavior is different
        # if adatpvive burnin
        
        #region \- prop tune --------------
        if(iter %% prop_adapt == 0 ) {
            if(!is.null(prop_sd)) {
                mean_acc_rate = sapply(
                    names(acc_counter),
                    function(param) acc_counter[[param]] / prop_adapt
                ) |>
                    unlist() |>
                    mean()
                acc_good = (mean_acc_rate > 0.17 & mean_acc_rate < 0.50)
                #update prop_sd
                for(param in names(acc_counter)) {
                    prop_sd[[param]] = adapt_tuner(param, acc_counter, prop_sd, prop_adapt)
                }
                acc_counter = reset_acc_counter(prop_sd)
            }
        }
        #endregion -------------------------------

        if(burn_adapt) {
            if(iter <= min_burn) {
                for(param in names(run_save)) {
                    run_save[[param]][iter, ] = as.vector(get(param))
                }
            } else {
                # streaming save of length min_burn
                for(param in names(run_save)) {
                    pdim = dim(run_save[[param]])
                    pdx = rep(list(TRUE), length(pdim)-1L)
                    run_save[[param]][1:min_burn, ] = run_save[[param]][2:(min_burn+1), ]
                    run_save[[param]][min_burn, ] = as.vector(get(param))
                }
                if(iter %% prop_adapt == 0) {                 
                    check_obj = trim_chain(run_save, min_burn)
                    if(assess_burnin(check_obj, chain_check) & acc_good) {
                        in_burnin = FALSE
                        acc_idx = 0
                        cat('\n', paste0("Exiting Burnin After ", iter, ' iterations..'))
                        cat('\n', 'Launching Save Intvs')
                        pbiter = n_burn # update progress bar to reflect iters
                        iter = n_burn #jump to end of burnin
                        cat('\n', paste0('jumping to ', iter, ' iterations...'))
                    }
                }
                if(iter >= n_burn) {
                    in_burnin = FALSE
                    acc_idx = 0
                    cat('\n', 'Burn-in exited at set maximum n_burn')
                }
            }
        } else {
            if(iter > n_burn) {
                in_burnin = FALSE
                acc_idx = 0
                cat('\n', 'Burnin ended at set length \n starting saving')
            }
        }
    }
    acc_idx = 0
    for (iter in (n_burn + 1):n_iter) {
        setTxtProgressBar(pb, pbiter)
        #region \- Run update block -------------------
        for(i in 1:length(update_functions)) {
            update_fn = update_functions[[i]]
                if('prop_sd' %in% names(formals(update_fn))) {
                    res = update_fn(prop_sd)
                    # only need to update here if in MH
                    for(cpar in names(res$acc)) {
                        acc_counter[[cpar]] = acc_counter[[cpar]] + res$acc[[cpar]]
                    }
                } else {
                    # for gibbs functions with no prop_sd
                    res = update_fn()
                }
                # this loop hits regardless of if it were an MH or Gibbs
                for(par in names(res$val)) {
                    assign(par, res$val[[par]])
            }
        }
        #endregion----------------------------------
    
        
        #region \- run derv fn ----------------
        if(!is.null(derv_quants)) {
            for(i in 1:length(derived_functions)) {
                derv_fn = derived_functions[[i]]
                res = derv_fn()
                for(par in names(res$val)) {
                    assign(par, res$val[[par]])
                }
            }
        }
        #endregion
        
        #region \- track acceptance rate ---------------------------
        # this keeps track of all iterations regardless of thinning
        # but only records every prop_adapt for memory sake
        if(iter %% prop_adapt == 0) {
            if(!is.null(prop_sd)) {
                acc_idx = acc_idx + prop_adapt
                mean_acc_rate = sapply(
                    names(acc_counter),
                    function(param) acc_counter[[param]] / acc_idx
                ) |>
                    unlist() |>
                    mean()
            }
        }
        #endregion ------------------------------

        #region \- save --------------
        if((iter%%thin) == 0) {
            save_count = save_count + 1
            for(param in save_names) {
                save[[param]][save_count, ] = as.vector(get(param))
            }

            if(!is.null(derv_quants)) {
                for(dname in derv_quants) {
                    derv_save[[dname]][save_count] = get(dname)
                }
            } 
        }
        #endregion --------------------------
    }
    cat("\n Done with MCMC, Saving....")

    # reformatting to original structure
    for(param in save_names) {
        ns = save_count
        d = pd[[param]]
        kp = plen[[param]]

        if(is.null(d)) {
            if(kp==1) {
                save[[param]] = save[[param]][1:ns, 1]
            } else {
                save[[param]] = save[[param]][1:ns,,drop = FALSE]
            }
        } else {
            save[[param]] = array(save[[param]], dim = c(ns, d))
        }
    }

    save = trim_chain(save, save_count)
    if(is.null(derv_quants)) derv_save = NULL else derv_save = trim_chain(derv_save, save_count)
    out = format_output(save, derv_save)
    return(out)
}


################
#MARK: Output
################
#' format outputs for mcmc algos
#' 
#' @param save_item save object as a list with named parameters
#' @param derv_save derv save object as list with named parameters OR NULL
#' 
#' @return formatted output for mcmc engines
#' @keywords internal
#' notes to self 
#' - remove drop front in non-adaptive stopping.
#' - add option for info calculation to not run if done in multi-chain cases.
format_output = function(save_item, derv_save) {
    #chop of the pre-thinned items
    if(!is.null(derv_save)){
        out = c(save_item, derv_save)
    } else {
        out = save_item
    }

    
    chain_info = list(
        'size' = save_count
    )
    if(exists('mean_acc_rate')) {
      chain_info$mean_acc = get("mean_acc_rate")
    }

    return(
        list(
            samples = out,
            info = chain_info
        )
    )
}

###############
#MARK: Check Function Format
################
#' Check update and derived function formats
#' 
#' This makes sure that the update/derived functions are properly formatted
#' and will return a list with val and acc (if MH).
#' Also it will check formatting between val names and prop_sd names (if needed).
#' 
#' @param fn function to check
#' 
#' @keywords internal
check_function_format = function(fn) { 
    needs_acc = 'prop_sd' %in% names(formals(fn))
    if(needs_acc) {
        temp_res = fn(prop_sd)
    } else {
        temp_res = fn()
    }
    if(!is.list(temp_res)) {
        stop('Update/Derived function must return a list')
    } else {
        # check naming format
        names_res = names(temp_res)
        if(!('val' %in% names_res)) {
            stop('Update/Derived function must return a list with a val list')
        } else if(needs_acc & !('acc' %in% names_res)) {
            stop('MH Update function must return a list with an acc list')
        }

        # check actual structure:
        for(par in names(temp_res$val)) {
            if(!exists(par)) {
                stop(paste0('Function returns ', par, ' but does not exist in scope'))
            } else {
                real_par = get(par)
                if(!(length(real_par) == length(temp_res$val[[par]]))) {
                    stop(paste0('Returned value for ', par, ' does not match length of parameter in scope'))
                }
            }
        }
        if(needs_acc) {
            for(cpar in names(temp_res$acc)) {
                if(!cpar %in% names(acc_counter)) {
                    stop(paste0('Function returns acc for ', cpar, ' but no prop_sd provided for this parameter'))
                } else {
                    if(!(length(acc_counter[[cpar]]) == length(temp_res$acc[[cpar]]))) {
                        stop(paste0('Returned acc for ', cpar, ' does not match length of prop_sd for this parameter'))
                    }
                }
            } 
        }
    }
}


####################
#MARK: Reset Acc Counter
####################

reset_acc_counter = function(prop_sd) {
    acc_list = list()
    for(param in names(prop_sd)) {
        # wrap to accomidate non-array entries
        if(is.null(dim(prop_sd[[param]]))) {
            pdim = length(prop_sd[[param]])
        } else {
            pdim = dim(prop_sd[[param]])
        }
        acc_list[[param]] = array(0, dim = pdim)
    }
    return(acc_list)

}

#' Internal dispatch for a single chain run
#' 
#' @inheritParams mcmc_run
#' @param n_burn length of burnin period if burn_adapt = FALSE, else minimum burnin
#' @param ... additional args passed to SEE LINK
#' 
mcmc_run_internal = function(
    data_list,
    const_list,
    init_list,
    prior_pars,
    update_functions,
    n_burn = 5000,
    n_iter = 1e5,
    thin = 1,
    prop_adapt = 100,
    burn_adapt = FALSE,
    min_burn = 1000,
    prop_sd = NULL,
    save_names = NULL,
    derv_quants = NULL,
    derived_functions = NULL,
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
    for(iter in 1:n_iter) {
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
        if(in_burnin) {
            #region \- prop tune --------------
            if(iter %% prop_adapt == 0 ) {
                if(exists("acc_counter")) {
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
                        if(assess_burnin(check_obj, chain_check)) {
                            in_burnin = FALSE
                            cat('\n', paste0("Exiting Burnin After ", iter, ' iterations..'))
                            cat('\n', 'Launching Save Intvs')
                            pbiter = n_burn
                        }
                    }
                    if(iter >= round(n_iter*0.9)) {
                        stop(paste0('Burn in failed to converge after ', iter, ' iterations'))
                    }
                }
            } else {
                if(iter > n_burn) {
                    in_burnin = FALSE
                    cat('\n', 'Burnin ended at set length \n starting saving')
                }
            }
        #endregion ----------------------------
        } else {
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
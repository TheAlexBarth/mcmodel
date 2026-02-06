###################
#MARK: ADAPTIVE
##################
#' MCMC Run
#' 
#' A generic mcmc engine for flexible model building
#' 
#' This function is the main driver of mcmc algos in a base R framework,
#' it is designed for maximum flexibility in mind, not necessarily speed.
#' However performance largely depends
#' 
#' There are several useful features including an adapative-burnin, thinning, and tuning
#' for metropolis hastings algorithms. Users should carefully note that the burnin and thinning
#' adaptation algorithms rely on parameters which are saved ONLY. Thus, if you are running a complex
#' model with multiple variables but only saving a few, there may be hidden issues.
#' Also be mindful of default max/min args, these can create a big memory constraint if
#' saving many paramters. Finally, the adaptive checking process for thinning, exit, and burnin
#' are extremely costly with many parameters. They really only offer memory benefits, but rarely improve
#' total computte time unless memory constraint is a large consideration.
#' 
#' @param min_burn integer of minimum number of burnin interations before saving
#' @param max_iter the maximum number of iterations
#' @param prop_adapt adaptation interval for both proposal tunning and thinnin tuning
#' @param max_thin maximum thinning interval
#' @param goal_ess target ess for a chain
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom stats acf quantile t.test
#' @importFrom coda effectiveSize
#' 
#' @return a list of sample objects and info.
#' 
#' @keywords internal
adaptive_mcmc_run = function(
    data_list,
    const_list,
    init_list,
    prior_pars,
    update_functions,
    prop_sd = NULL,
    save_names = NULL,
    derv_quants = NULL,
    derived_functions = NULL,
    adapt_tuner = adapt_tuner_double_exp,
    chain_check = geweke_check,
    adaptive_thin = FALSE,
    min_burn = 5000,
    max_iter = 1e6,
    prop_adapt = 100,
    adapt_intv = 1000,
    max_thin = 10,
    goal_ess = 100
) {
    #region \- load inputs
    list2env(data_list, envir = environment())
    list2env(const_list, envir = environment())
    list2env(init_list, envir = environment())
    list2env(prior_pars, envir = environment())
    for(i in 1:length(update_functions)) {
        environment(update_functions[[i]]) = environment()
    }
    environment(format_output) = environment()
    environment(check_function_format) = environment()
    
    #region \- structure save array
    if(is.null(save_names)) save_names = names(init_list)
    max_save = max_iter - min_burn
    # auto set max burnin to half of max iterations
    save = list()
    run_save = list() # for burn-in assessment
    for(param in save_names) {
        param_dims = dim(init_list[[param]])
        if(is.null(param_dims)) {
            param_dims = length(init_list[[param]])
        }
        save[[param]] = array(NA, dim = c(max_save, param_dims))
        run_save[[param]] = array(NA, dim = c(min_burn+1, param_dims))
    }

    # set derived quanities
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
            derv_save[[dname]] = array(NA, dim = max_save)
        }
    }

    #region \- set acc counter if needed ---------
    if(!is.null(prop_sd)) {
        acc_counter = reset_acc_counter(prop_sd)
    }
    #endregion -----------------------------------

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


    #region \- set up running trackers
    pb1 = txtProgressBar(0, max_iter/2, style = 3)
    in_burnin = TRUE
    save_count = 0
    pb2 = txtProgressBar(0, max_iter, style = 3)
    thin_intv = 1
    #MARK: LOOP
    cat('\n', "Launching Run - Burning in")
    for(iter in 1:max_iter) {
        #region \- parameter updates -------------------
        # update functions should be written to have scoping 
        # in mind for accessing items in the input lists (data, const, init, prior)
        # update functions must return a list with a val object and if MH an accepted list
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

        #region \- Adapt Tuning ----------------------------
        if(iter %% prop_adapt == 0 & in_burnin) {
            if(exists("acc_counter")) {
                for(param in names(acc_counter)) {
                    prop_sd[[param]] = adapt_tuner(param, acc_counter, prop_sd, prop_adapt)
                }
                acc_counter = reset_acc_counter(prop_sd)
            }
        }

        #region \- Assess Burnin ----------------------------
        if(in_burnin) {
            setTxtProgressBar(pb1, iter)
            if(iter <= min_burn) {
                for(param in names(run_save)) {
                    run_save[[param]] = aafa(run_save[[param]], iter, get(param))
                }
            } else {
                if(iter %% prop_adapt == 0) {
                    for(param in names(run_save)) {
                        pdim = dim(run_save[[param]])
                        pdx = rep(list(TRUE), length(pdim)-1L)

                        run_save[[param]] = aafa(
                            run_save[[param]],
                            1:min_burn,
                            gafa(run_save[[param]], 2:(min_burn+1))
                        )
                    }
                    check_obj = trim_chain(run_save, min_burn)
                    if(assess_burnin(check_obj, chain_check)){
                        in_burnin = FALSE
                        setTxtProgressBar(pb1, max_iter/2)
                        flush.console()
                        cat('\n', paste0("Exiting Burnin After ", iter, ' iterations..'))
                        cat('\n','Launching Save Intvs')
                    }
                    if(iter >= max_iter/2) {
                        in_burnin = FALSE
                        flush.console()
                        warning('\n', 'Maximum Burn In Hit Before Stabilizing.....')
                    }
                }
            }
        } else {
            setTxtProgressBar(pb2, iter)
            #region \- run derv functions ------------
            if(!is.null(derv_quants)) {
                for(i in 1:length(derived_functions)) {
                    derv_fn = derived_functions[[i]]
                    res = derv_fn()
                    for(par in names(res$val)) {
                        assign(par, res$val[[par]])
                    }
                }
            }

            #region \- Saving  --------------------
            if((iter %% thin_intv) == 0) {
                save_count = save_count + 1
                for(param in save_names) {
                    save[[param]] = aafa(save[[param]], save_count, get(param))
                }

                if(!is.null(derv_quants)) {
                    for(dname in derv_quants) {
                        derv_save[[dname]][save_count] = get(dname)
                    }
                }
            }

            #region \- Assess Exit  --------------------
            if((save_count %% adapt_intv) == 0 & save_count > goal_ess) {
                all_good = NULL # tracker to see if they are all good length
                cur_save = trim_chain(save, save_count)
                if(adaptive_thin) sugg_thin = NULL #only for adaptive thinning
                for(param in save_names) {
                    chain_mat = fafa(cur_save[[param]])
                    all_good = c(all_good, mean(effectiveSize(chain_mat)) > goal_ess)
                    if(adaptive_thin) {
                        # thining updates
                        for(j in 1:ncol(chain_mat)) {
                            ac = acf(chain_mat[,j], plot = FALSE)$acf
                            min_thin = which(ac[-1] <0.5)[1]
                            if(!is.na(min_thin)) sugg_thin = c(sugg_thin, min_thin)
                        }
                    }
                }
                if(!is.null(all_good) & all(all_good)) {
                    cat('\n', paste0('Exiting at ', iter, ' iterations'))
                    if(is.null(derv_quants)) derv_save = NULL else derv_save = trim_chain(derv_save, save_count)
                    out = format_output(cur_save, derv_save)
                    return(out)
                }
                if(adaptive_thin) {
                    if(length(sugg_thin) > 0) {
                        new_thin = round(mean(sugg_thin))
                        new_thin = max(1, min(new_thin, max_thin))
                        if(new_thin != thin_intv) {
                            thin_intv = new_thin
                            cat('\n', paste0('Setting new thin to ', thin_intv))
                        }
                    }
                }
            }
        }
    }

    if(is.null(derv_quants)) derv_save = NULL
    out = format_output(save, derv_save)
    return(out)
}


###################
#MARK: Simple
###################
#' MCMC Run Simple
#' 
#' Non-adaptive (except proposal tunning) mcmc run
#' 
#' @param burnin length of burnin interval
#' @param n_iter total numer of iterations
#' @param thin thinning interval
#' @param prop_adapt adaptation interval for both proposal tunning and thinnin tuning
#' 
#' @return list of samples and info from MCMC
#' @keywords internal
simple_mcmc_run = function(
    data_list,
    const_list,
    init_list,
    prior_pars,
    update_functions,
    burnin = 5000,
    n_iter = 1e5,
    thin = 1,
    prop_adapt = 100,
    prop_sd = NULL,
    save_names = NULL,
    derv_quants = NULL,
    derived_functions = NULL,
    adapt_tuner = adapt_tuner_double_exp
) {
    # load params
    list2env(data_list, envir = environment())
    list2env(const_list, envir = environment())
    list2env(init_list, envir = environment())
    list2env(prior_pars, envir = environment())
    for(i in 1:length(update_functions)) {
        environment(update_functions[[i]]) = environment()
    }
    environment(format_output) = environment()
    environment(check_function_format) = environment()

    # make save array
    if(is.null(save_names)) save_names = names(init_list)
    save_size = (n_iter-burnin)/thin
    save = list()
    for(param in save_names) {
        param_dims = dim(init_list[[param]])
        if(is.null(param_dims)) {
            param_dims = length(init_list[[param]])
        }
        save[[param]] = array(NA, dim = c(save_size, param_dims))
    }

    # set derv qs
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
            derv_save[[dname]] = array(NA, dim = save_size)
        }
    }

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
    save_count = 0
    cat('\n', 'Launching Run')
    for(iter in 1:n_iter) {
        setTxtProgressBar(pb, iter)
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

        #region \- check for burn -----------------------
        if(iter <= burnin) {
            # tuning
            if(exists('acc_counter') & (iter %% prop_adapt) == 0) {
                for(param in names(acc_counter)) {
                    prop_sd[[param]] = adapt_tuner(param, acc_counter, prop_sd, prop_adapt)
                }
                acc_counter = reset_acc_counter(prop_sd)
            }
        } else {
            #region \- saving
            #region \- run derv functions ------------
            if(!is.null(derv_quants)) {
                for(i in 1:length(derived_functions)) {
                    derv_fn = derived_functions[[i]]
                    res = derv_fn()
                    for(par in names(res$val)) {
                        assign(par, res$val[[par]])
                    }
                }
            }

            #region \- Saving  --------------------
            if((iter %% thin) == 0) {
                save_count = save_count + 1
                for(param in save_names) {
                    save[[param]] = aafa(save[[param]], save_count, get(param))
                }

                if(!is.null(derv_quants)) {
                    for(dname in derv_quants) {
                        derv_save[[dname]][save_count] = get(dname)
                    }
                }
            }
        }
    }

    if(is.null(derv_quants)) derv_save = NULL
    out = format_output(save, derv_save, drop_front = FALSE)
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
format_output = function(save_item, derv_save, drop_front = TRUE) {
    #chop of the pre-thinned items
    if(!is.null(derv_save)){
        out = c(save_item, derv_save)
    } else {
        out = save_item
    }

    if(drop_front) {
        for(item in names(out)) {
            out[[item]] = gafa(out[[item]], -c(1:goal_ess))
        }
        chain_info = list(
            'ESS' = sapply(out, function(x) effectiveSize(fafa(x))),
            'size' = save_count-goal_ess
        )
    } else {
        chain_info = list(
            'ESS' = sapply(out, function(x) effectiveSize(fafa(x))),
            'size' = save_count
        )
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
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
    derived_functions,
    save_names,
    derv_quants,
    n_iter,
    thin,
    burn_adapt,
    min_burn,
    n_burn,
    prop_adapt,
    chkpt_path,
    chkpt_freq,
    resume = TRUE,
    reconstruct = TRUE,
    chkpt_clean = TRUE,
    adapt_tuner = adapt_tuner_double_exp,
    chain_check = geweke_check
) {
    
    # resume logic:
    if(resume) {
        # add resume logic here
        # find resumable files
        # load state
        state = load_recent(chkpt_path)
        # reset state parameters to match new inputs
        if(is.null(state)) {
            cat('\n', 'No checkpoints or blocks found... starting new', '\n')
            resume = FALSE
        } else {
            cat('\n', "Loading from saved point", '\n')
            state$n_iter = n_iter
            if(state$burn$min < min_burn) {
                if(state$phase == 'burnin') {
                    state$burn$min = min_burn
                } else {
                    warning('Loaded state past burnin with shorter burn period')
                }
            }
        }
    } 
    if(!resume)  {
        # define state:
        state = make_state(
            data_list = data_list,
            const_list = const_list,
            init_list = init_list,
            prior_pars = prior_pars,
            prop_sd = prop_sd,
            update_functions = update_functions,
            derived_functions = derived_functions,
            save_names = save_names,
            derv_quants = derv_quants,
            n_iter = n_iter,
            thin = thin,
            burn_adapt = burn_adapt,
            min_burn = min_burn,
            n_burn = n_burn,
            prop_adapt = prop_adapt,
            adapt_tuner = adapt_tuner,
            chain_check = chain_check,
            chkpt_path = chkpt_path,
            chkpt_freq = chkpt_freq
        )
    }
    
    # run burnin
    state = run_burn(state)

    remaining = ((state$n_iter - state$iter) %/% state$block_size)
    #loop through 
    if(remaining > 0) {
        for(block in 1:remaining) {
            state = run_samples(state, state$block_size)
        }
    } else {
        cat('\n', 'Loading complete blocks!', '\n')
    }


    cat("\n Done with MCMC, Saving....")

    if(reconstruct) {
        out = load_all_blks(state)
    } else {
        cat('\n', 'Runs completed and saved as blocks')
        return(":)")
        #note it returns here just to avoid possible clean up below..
    }
    if(chkpt_clean) {
        unlink(state$chk_path, TRUE)
    }
    return(out)
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
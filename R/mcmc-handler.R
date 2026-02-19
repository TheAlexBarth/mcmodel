
#' Run MCMC
#'
#' High-level wrapper for running simple or adaptive MCMC chains with
#' options for single or multichain runs (in parallel). See \code{\link{mcmc_structures}} for full detail
#' for details and examples on properly constructing runs. Algorithmic development was inspired by Hooten & Hobbs 2019, however the code
#' allows for substantial flexibility. The structure of data inputs were largely inspired by NIMBLE, but intentially without a bulky object-oriented
#' structure.
#' 
#' @param data_list list object with non-changing data
#' @param const_list list object with relevant indexing used in updates/derived fns
#' @param init_list list object with ALL parameters used in update/derived, see \link{mcmc_structures}
#' @param init_generator (optional) IF n_chains >1, a list of functions to change init_list points. refer to \link{mcmc_structures}
#' @param prior_pars list of prior values used in update functions
#' @param prop_sd IF any Metropolis-Hastings algos are used with, prop_sd much be provided.
#' @param update_functions a list of update functions used in blocks for MCMC application
#' @param save_names (optional) list of parameters in init_list to save out. Defaults to init_list
#' @param derv_quants Character vector specifying derived quantities to track.
#' @param derived_functions List of functions for computing derived quantities.
#' @param n_iter integer for number of interations
#' @param thin integer for thinning interval. Default is 1.
#' @param n_burn integer for maximum burn-in length
#' @param burn_adapt logical to decide to use adaptive burn-in procedure.
#' @param min_burn integer used if burn_adapt=TRUE, how long minimum burn in needs to be
#' @param n_chains number of chains to run (if > 1 will be done in parallel)
#' @param n_cores optionally specify number of cores for parallel run. If don't want parallel execution, set n_cores=1
#' @param seeds optional vector of seeds (equal to n_chains)
#' @param log_files logical, should chain outputs be recorded if run in parallel? Default TRUE
#' @param log_dir if log_files, where to write a log directory (note one log will be made for each chain)
#' @param delete_logs logical, delete log files once completed. default is TRUE
#' @param ... arguments passed to the internal MCMC engines, see \link{mcmc_run_internal}
#'
#' @details
#' There are several features in this default MCMC algorithm which are worth understanding to correct set parameters. First, users
#' should be familar with how the burn-in procedure is selected. If adaptive-burn in is used, the burn-in is set to be 
#' AT LEAST min_burn and AT MOST n_burn. Thus, if using burn_adapt, it is recommended to increase n_burn to ensure that sufficient
#' burn-in is provided as the default settings are likely too low in most real cases. If using adaptive burn-in the default exit
#' is based on \code{\link{geweke_check}}. If MH-steps are used, early exit in adaptive burn-in also requires all acceptance rates to fall
#' between 0.17 and 0.50. Finally, if working in parallel, it is possible that some chains will end sooner than others in adaptive burn-in
#' To allow for convergence checks, I automatically trucate longer chains to the shortest one - effectivelly imposing a similar burn-in across samples.
#' 
#' In all cases, there is an adaptive proposal adjustment for MH-updates. The function used to update the prop_sd list is provided with adapt_tuner.
#' By default this is \code{\link{adapt_tuner_double_exp}}, although user definied functions are possible. To avoid proposal adaptation, set adapt_tuner=no_adapt
#' While acceptance rate per parameter is not recorded, the overall mean acceptance rate is reported if MH steps are present.
#' 
#' A full working example is available with package source code OR can be produced throgugh following the set-up with \link{mcmc_structures}
#' 
#' 
#' @return A list of MCMC output with two items, samples list and info list
#'
#' @examples
#' \dontrun{
#' mcmc_out <- mcmc_run(
#'  data_list = data_list,
#'  const_list = const_list,
#'  init_list = init_list,
#'  init_generator = init_generator,
#'  prior_pars = my_prior,
#'  prop_sd = prop_sd,
#'  update_functions = update_fns,
#'  save_names = c('beta', 'sigma'),
#'  derv_quants = c('p_mse'),
#'  derived_functions = list('val' = calc_pmse),
#'  n_iter = 1e4,
#'  thin = 1,
#'  n_burn = 5000,
#'  burn_adapt = TRUE,
#'  min_burn = 1000,
#'  n_chains = 10
#' )
#' }
#' 
#' @importFrom abind abind
#' @importFrom coda effectiveSize
#' 
#' @references
#' Hooten, M.B. and T.J. Hefley. (2019). Bringing Bayesian Models to Life. Chapman and Hall/CRC.
#' NIMBLE Development Team. 2026. NIMBLE: MCMC, Particle Filtering, and Programmable Hierarchical Modeling.
#'    doi: 10.5281/zenodo.1211190. R package version 1.4.1, https://cran.r-project.org/package=nimble.
#'
#' @export
mcmc_run = function(
    data_list,
    const_list,
    init_list,
    init_generator = NULL,
    prior_pars,
    prop_sd = NULL,
    update_functions,
    save_names = NULL,
    derv_quants = NULL,
    derived_functions = NULL,
    n_iter = 1e5,
    thin = 1,
    n_burn = 5e3,
    burn_adapt = FALSE,
    min_burn = 1e3,
    n_chains = 1,
    n_cores = NULL,
    seeds = NULL,
    log_files = TRUE,
    log_dir = '.',
    delete_logs = TRUE,
    ...
) {
    dots = list(...)
    # MARK: Checks
    if(!is.null(seeds) & length(seeds) != n_chains) {
        stop('if providing seeds - it must be same length as chains')
    }

    if(n_chains > 1) {
        #MARK: Multi-chain run
        init_lists = list()
        if(is.null(init_generator)) {
            warning('No init_generator provided, all inits will be same point')
            for(i in 1:n_chains) {
                init_lists[[i]] = init_list
            }
        } else {
            #region \- initgen
            for(i in 1:n_chains) {
                chain_init = init_list
                for(par in names(init_generator)) {
                    if(!(par %in% names(init_list))) {
                        stop(paste0(par, " not found in init_list but provided in init_gen"))
                    } else {
                        chain_init[[par]] = init_generator[[par]](chain_init[[par]])
                    }
                }
                init_lists[[i]] = chain_init
            }            
        }
        #region \- parallel run ---------------------------
        # set up cores
        if(is.null(n_cores)) {
            # auto select number of cores
            n_cores = min(round(future::availableCores()*0.75), n_chains)
        } else {
            if(n_cores > future::availableCores()-1) warning("Machine does not support requested cores")
            n_cores = max(1, min(n_cores, future::availableCores()-1))
        }
        # set up plan
        # per copilot - make it posssible for user to overwrite their own
        # future plan (e.g., they could set a multicore or cluster plan)
        current_plan = future::plan()
        if(inherits(current_plan, 'sequential') & n_cores > 1) {
            future::plan(future::multisession, workers = n_cores)
        }
        message("Running ", n_chains, " chains using future with ", n_cores, " workers....")
        #endregion -----------------------------------------------------------

        #region \- worker function --------------------------
        run_chain = function(
            chain_id,
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
            extra_args = dots
        ) {
            if(log_files) {
                log_dir = file.path(log_dir, 'logs')
                if(!dir.exists(log_dir)) {
                    dir.create(log_dir, recursive = TRUE)
                }
                log_file = file.path(
                    log_dir, paste0(chain_id,".log")
                )
                zz = file(log_file, open = 'wt')
                sink(zz)
                sink(zz, type = "message")
            }
            base_args = list(
                    data_list = data_list,
                    const_list = const_list,
                    init_list = init_list, # note this is confusing because here it is a worker-list, but it is informed from init_lists[[chain_id]]
                    prior_pars = prior_pars,
                    prop_sd = prop_sd,
                    update_functions = update_functions,
                    save_names = save_names,
                    derv_quants = derv_quants,
                    derived_functions = derived_functions,
                    n_iter = n_iter,
                    thin = thin,
                    n_burn = n_burn,
                    burn_adapt = burn_adapt,
                    min_burn = min_burn
            )
            on.exit(
                {
                    if(!is.null(log_file)) {
                        sink(type = 'message')
                        sink()
                        close(zz)

                        if(delete_logs) {
                            file.remove(log_file)
                        }
                    }
                },
                add = TRUE
            )
            all_args = c(base_args, extra_args)
            worker_out = do.call(mcmc_run_internal, all_args)
            return(worker_out)
      
        }

        chain_res = future.apply::future_lapply(
            X = seq_len(n_chains),
            FUN = function(chain_id) {
                run_chain(
                    chain_id = chain_id,
                    data_list = data_list,
                    const_list = const_list,
                    init_list = init_lists[[chain_id]],
                    prior_pars = prior_pars,
                    prop_sd = prop_sd,
                    update_functions = update_functions,
                    save_names = save_names,
                    derv_quants = derv_quants,
                    derived_functions = derived_functions,
                    n_iter = n_iter,
                    thin = thin,
                    n_burn = n_burn,
                    burn_adapt = burn_adapt,
                    min_burn = min_burn,
                    extra_args = dots
                )
            },
            future.seed = if (is.null(seeds)) TRUE else seeds,
            future.globals = list(
                run_chain = run_chain,
                data_list = data_list,
                const_list = const_list,
                init_lists = init_lists,
                prior_pars = prior_pars,
                prop_sd = prop_sd,
                update_functions = update_functions,
                save_names = save_names,
                derv_quants = derv_quants,
                derived_functions = derived_functions,
                n_iter = n_iter,
                thin = thin,
                n_burn = n_burn,
                burn_adapt = burn_adapt,
                min_burn = min_burn,
                log_files = log_files,
                log_dir = log_dir,
                delete_logs = delete_logs,
                dots = dots
            )
        )
        #endregion -----------------------------------------------------------------

        # could build optional to not run in parallel

        #region \- post-loop formatting -------------
        # get samples
        
        min_length = chain_res |>
            lapply('[[', 'info') |>
            sapply('[[', 'size') |>
            min()

        chain_samples = truncate_samples(chain_res, min_length)

        # need to flip samples to vector form
        if(!is.null(derv_quants)) {
            keep_names = c(save_names, derv_quants)
        } else {
            keep_names = save_names
        }
        chain_samples = lapply(
            keep_names,
            function(par) lapply(chain_samples, '[[', par)
        )
        names(chain_samples) = keep_names
        split_r = lapply(chain_samples, vehtari_split_psrf)

        combn_samples = lapply(
            chain_samples,
            abind,
            along = 1
        )
        out = list(
            samples = combn_samples,
            info = list(
                "n_chains" = n_chains,
                'rhat' = split_r,
                "size" = min_length*n_chains,
                "ESS" = sapply(combn_samples, function(x) effectiveSize(fafa(x))),
                "chain_idx" = rep(1:n_chains, each = min_length)
            )
        )
      
        if(!is.null(chain_res[[1]]$info$mean_acc)) {
            out$info$mean_acc = sapply(chain_res, function(x) x[['info']][['mean_acc']]) |> 
                mean()
        }

        #endregion ------------------------------------
        

    } else {
        #MARK: Single Run
        if(!is.null(seeds)) set.seed(seeds[1])
        out = do.call(
            mcmc_run_internal,
            c(
                list(        
                    data_list = data_list,
                    const_list = const_list,
                    init_list = init_list,
                    prior_pars = prior_pars,
                    prop_sd = prop_sd,
                    update_functions = update_functions,
                    save_names = save_names,
                    derv_quants = derv_quants,
                    derived_functions = derived_functions,
                    n_iter = n_iter,
                    thin = thin,
                    n_burn = n_burn,
                    burn_adapt = burn_adapt,
                    min_burn = min_burn
                ),
                dots
            )
        )
        
        out$info$ESS = sapply(out$samples, function(x) effectiveSize(fafa(x)))
        out$info$rhat = lapply(out$samples, vehtari_split_psrf)
        out$info$n_chains = 1
    }
    return(out)
}

########
#MARK: Formatting
##########
#' Truncate to the shortest chain length
#' 
#' for multi-chain runs, chains must be similar length
#' in order to calculate Rhat. This discards some small number of 
#' points. Future extensions of this package may consider tracking/warning
#' if chains are substantially different lengths.
#' 
#' @param chain_res result from future loop
#' @param min_length the minimum lenght of all chains
#' 
#' @return a list of truncated samples
#' @keywords internal
truncate_samples = function(chain_res, min_length) {
    #could use trim chain here
    trunc_samples = chain_res |>
        lapply('[[', 'samples') |>
        lapply(trim_chain, min_length, front = FALSE)

    return(trunc_samples)
}

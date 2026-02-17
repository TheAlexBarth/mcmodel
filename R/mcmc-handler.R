
#' Run MCMC
#'
#' High-level wrapper for running simple or adaptive MCMC chains with
#' options for single or multichain runs (in parallel). See \code{\link{mcmc_structure}} for full detail
#' for details and examples on properly constructing runs. 
#' 
#' @param data_list list object with non-changing data
#' @param const_list list object with relevant indexing used in updates/derived fns
#' @param init_list list object with ALL parameters used in update/derived, see \link{mcmc_structure}
#' @param init_generator (optional) IF n_chains >1, a list of functions to change init_list points. refer to \link{mcmc_structure}
#' @param prior_pars list of prior values used in update functions
#' @param update_functions a list of update functions used in blocks for MCMC application
#' @param adaptive logical - if true adaptive burnin and thinning will be used. Adaptive tunning occurs in both unless adapt_tuner=no_adapt
#' @param save_names (optional) list of parameters in init_list to save out. Defaults to init_list
#' @param prop_sd IF any Metropolis-Hastings algos are used with, prop_sd much be provided.
#' @param derv_quants Character vector specifying derived quantities to track.
#' @param derived_functions List of functions for computing derived quantities.
#' @param n_chains number of chains to run (if > 1 will be done in parallel)
#' @param n_cores optionally specify number of cores for parallel run. If don't want parallel execution, set n_cores=1
#' @param seeds optional vector of seeds (equal to n_chains)
#' @param adapt_tuner function for tuning prop_sd in MH applications
#' @param chain_check function for assess_burnin if run in adaptive burnin.
#' @param ... arguments passed to the internal MCMC engines, see
#'
#' @return A list of MCMC output with two items;
#' \itemize{
#'     \item{samples}{list of posterior samplews}
#'     \item{info}{chain diagnostics}
#' }
#'
#' @examples
#' \dontrun{
#' mcmc_out <- mcmc_run(
#'  data_list = my_data,
#'  const_list = my_consts,
#'  init_list = my_init,
#'  prior_pars = my_prior,
#'  update_functions = update_fns,
#'  chains = 1
#' )
#' }
#'
#' @export
mcmc_run = function(
    data_list,
    const_list,
    init_list,
    init_generator = NULL,
    prior_pars,
    update_functions,
    burn_adapt = FALSE,
    save_names = NULL,
    prop_sd = NULL,
    derv_quants = NULL,
    derived_functions = NULL,
    n_chains = 1,
    n_cores = NULL,
    seeds = NULL,
    adapt_tuner = adapt_tuner_double_exp,
    chain_check = geweke_check,
    log_dir = '.',
    log_files = TRUE,
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
        stop('Parallel Chains not yet implemented')

        #region \- worker function --------------------------
        run_chain = function(
            chain_id,
            data_list,
            const_list,
            init_list,
            prior_pars,
            update_functions,
            burn_adapt,
            save_names,
            prop_sd,
            derv_quants,
            derived_functions,
            adapt_tuner,
            chain_check,
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
                init_list = init_list, # note this is a worker-list and will be passed from lists
                prior_pars = prior_pars,
                update_functions = update_functions,
                burn_adapt = burn_adapt,
                save_names = save_names,
                prop_sd = prop_sd,
                derv_quants = derv_quants,
                derived_functions = derived_functions,
                adapt_tuner = adapt_tuner,
                chain_check = chain_check,
            )
            all_args = c(base_args, extra_args)
            worker_out = do.call(mcmodel:::mcmc_run_internal, all_args)
            return(worker_out)
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
                    update_functions = update_functions,
                    burn_adapt = burn_adapt,
                    save_names = save_names,
                    prop_sd = prop_sd,
                    derv_quants = derv_quants,
                    derived_functions = derived_functions,
                    adapt_tuner = adapt_tuner,
                    chain_check = chain_check,
                    extra_args = dots
                )
            },
            future.seed = if (is.null(seeds)) TRUE else seeds
        )

        run_chain = function(
            chain_id,
            data_list,
            const_list,
            init_list,
            prior_pars,
            update_functions,
            burn_adapt,
            save_names,
            prop_sd,
            derv_quants,
            derived_functions,
            adapt_tuner,
            chain_check
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
                init_list = init_list, # note this is a worker-list and will be passed from lists
                prior_pars = prior_pars,
                update_functions = update_functions,
                burn_adapt = burn_adapt,
                save_names = save_names,
                prop_sd = prop_sd,
                derv_quants = derv_quants,
                derived_functions = derived_functions,
                adapt_tuner = adapt_tuner,
                chain_check = chain_check
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
            all_args = c(base_args)
            worker_out = do.call(mcmodel:::mcmc_run_internal, all_args)
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
                    update_functions = update_functions,
                    burn_adapt = burn_adapt,
                    save_names = save_names,
                    prop_sd = prop_sd,
                    derv_quants = derv_quants,
                    derived_functions = derived_functions,
                    adapt_tuner = adapt_tuner,
                    chain_check = chain_check
                )
            },
            future.seed = if (is.null(seeds)) TRUE else seeds,
            future.globals = list(
                run_chain = run_chain,
                data_list = data_list,
                const_list = const_list,
                init_lists = init_lists,
                prior_pars = prior_pars,
                update_functions = update_functions,
                burn_adapt = burn_adapt,
                save_names = save_names,
                prop_sd = prop_sd,
                derv_quants = derv_quants,
                derived_functions = derived_functions,
                adapt_tuner = adapt_tuner,
                chain_check = chain_check,
                log_files = log_files,
                log_dir = log_dir,
                delete_logs = delete_logs
            )
        )


        #region \- post-loop formatting -------------
        # get samples
        chain_samples = chain_res |> lapply('[[', 'samples')
        out_samples = list()
        for(param in save_names) {
            out_samples[[param]] = lapply(chain_samples, '[[', param)
        }

        

    } else {
        if(!is.null(seeds)) set.seed(seeds[1])
        out = do.call(
            mcmc_run_internal,
            c(
                list(        
                    data_list = data_list,
                    const_list = const_list,
                    init_list = init_list,
                    prior_pars = prior_pars,
                    update_functions = update_functions,
                    burn_adapt = burn_adapt,
                    save_names = save_names,
                    prop_sd = prop_sd,
                    derv_quants = derv_quants,
                    derived_functions = derived_functions,
                    adapt_tuner = adapt_tuner,
                    chain_check = chain_check
                ),
                dots
            )
        )
        return(out)
    }
}

########
# Formatting
##########

psrf_from_arr = function(chain_list, psrf_fx = gelman_rubin_psrf) {
    # need to format this into a bigger fx
}

#' Truncate to the shortest chain length
#' 
#' for multi-chain runs, chains must be similar length
#' in order to calculate Rhat.
#' 
#' 
auto_trucate = function(chain_res) {
    min_length = chain_res |>
        lapply('[[', 'info') |>
        sapply('[[', 'size') |>
        min()

    #could use trim chain here
    
}

##############
#MARK: STATE
##############

#' Construct a new state for mcmc
#' 
#' McModel States are the core interal object for mcmc runs
#' They allow for checkpoint (chkpt) saving which can allow for in-memory action
#' 
#' @inheritParams mcmc_run
#' 
#' @return state object
#' 
#' @keywords internal
make_state = function(
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
    max_burn,
    prop_adapt,
    adapt_tuner,
    chain_check,
    chkpt_path,
    chkpt_freq
) {
  env = new.env(parent = parent.frame())
  list2env(c(data_list, const_list, prior_pars, init_list), env)
  
  for(i in 1:length(update_functions)) {
    environment(update_functions[[i]]) = env
  }
  if (!is.null(derived_functions)) {
    for(i in 1:length(derived_functions)) {
      environment(derived_functions[[i]]) = env
    }
  }

  block_size = chkpt_freq %/% thin # how big to write chunks

  state = list(
    phase = 'burnin',
    iter = 0, #for tracking progress in larger scale
    env = env,
    pars = init_list,
    prop_sd = prop_sd,
    acc_counter = if (!is.null(prop_sd)) reset_acc_counter(prop_sd) else NULL,
    acc_good = TRUE,
    updates = update_functions,
    dervs = derived_functions,
    burn = list(
      adapt = burn_adapt,
      min = min_burn,
      max = max_burn,
      prop_adapt = prop_adapt,
      tuner = adapt_tuner,
      check = chain_check,
      run_save = if(burn_adapt) make_save_block(init_list, save_names, min_burn) else NULL
    ),
    save = make_save_block(init_list, save_names, block_size),
    derv_save = if (!is.null(derv_quants)) make_save_block(NULL, derv_quants, block_size) else NULL,
    save_names = save_names,
    derv_quants = derv_quants,
    # these are specific to checkpointing
    use_chkpts = if (chkpt_freq > 0) TRUE else FALSE,
    chk_path = chkpt_path,
    block_size = block_size,
    blk_id = 0,
    blk_count = 0, #internal iterations
    chkpt_freq = chkpt_freq,
    pb = txtProgressBar(0, n_iter, style = 3)
  )
}


#' Constructor for checkpoint paths
#' 
#' @param init_list a list for initialization
#' @param save_names pars to keep track of 
#' @param alloc how big each block needs to be
#' 
#' @returns block list strucutre
#' 
#' @keyworks internal
make_save_block = function(init_list, save_names, alloc) {
  blk = list()
  if(!is.null(init_list)) {
    for(p in save_names) {
      blk[[p]] = array(NA, c(alloc, length(init_list[[p]])))
    }
  } else {
    for(p in save_names) {
      blk[[p]] = array(NA, alloc)
    }
  }
  return(blk)
}



##############
#MARK: IO
##############

#' Write out block point
#' 
#' @param state a mcmodel state
#' 
#' @returns a flushed state AND writes a block out
#' @keywords internal
write_block = function(state) {
  blk_id = state$blk_id + 1

  file = file.path(
    state$chk_path,
    paste0('block_',state$blk_id,'.rds')
  )

  #save out chunk
  saveRDS(state, file)

  #flush save carrier
  state$save = make_save_block(state$pars, state$save_names, state$block_size)
  if(!is.null(state$derv_save)) {
    state$derv_save = make_save_block(NULL, state$derv_quants, state$block_size)
  }
  state$blk_count = 0
  return(state)
}

#' Burn Checkponit
#' 
#' @param state mcmodel state obj
#' @returns nothing - just writes current state
#' @keywords internal
write_check = function(state) {
  file = file.path(state$chk_path, paste0("chkpt_", state$iter,".rds"))
  saveRDS(state, file)
}


###########
#MARK: STEPPERS
###########

#' Update Step
#' 
#' @param state a mcmodel state list
#' 
#' @returns a state object
#' @keyword internal
step_update = function(state) {
  for(fn in state$updates) {
    if('prop_sd' %in% names(formals(fn))) {
      res = fn(state$prop_sd)
      for(cpar in names(res$acc)) {
        state$acc_counter[[cpar]] = state$acc_counter[[cpar]] + res$acc[[cpar]]
      }
    } else {
      res = fn()
    }
    for (p in names(res$val)) {
      state$pars[[p]] = res$val[[p]]
      assign(p, res$val[[p]], state$env) #assign's to the states env
    }
  }
  return(state)
}

#' Tune proposal distribution
#' 
#' @param state mcmodel state object
#' 
#' @returns state object
#' @keywords internal
step_adapt = function(state) {
  if(is.null(state$prop_sd)) return(state)
  if(state$iter %% state$burn$prop_adapt != 0) return(state)
  
  mean_acc_rate = lapply(
    state$acc_counter, function(x) {
      x/state$burn$prop_adapt
    }
  ) |>
    unlist() |>
    mean()

  state$acc_good <- mean_acc_rate > 0.17 && mean_acc_rate < 0.50
  for(p in names(state$ac_counter)) {
    state$prop_sd[[p]] = state$burn$tuner(
      p, state$acc_counter, state$prop_sd, state$burn$prop_adapt
    )
  }
  state$acc_counter = reset_acc_counter(state$prop_sd)
  return(state)
}

#' Streaming burnin save
#' 
#' @param state a mcmodel state list
#' 
#' @returns a state object
#' @keywords internal
step_burn_stream = function(state) {
  rs = state$burn$run_save
  i = state$iter
  min = state$burn$min

  if(i <= min) {
    for(p in names(rs)) {
      rs[[p]][i, ] = state$pars[[p]]
    }
  } else {
    for(p in names(rs)) {
      rs[[p]][1:min, ] = rs[[p]][2:(min+1),]
      rs[[p]][min, ] = state$pars[[p]]
    }
  }
  state$burn$run_save = rs
  return(state)
}

#' Check Burn exit 
#' 
#' @param state a mcmodel state
#' 
#' @returns boolean exit
#' @keywords internal
step_check_burn = function(state) {
  if(!state$burn$adapt) return(FALSE)
  if(state$iter < state$burn$min) return(FALSE)
  if(!state$acc_good) return(FALSE)
  
  status = trim_chain(state$burn$run_save, state$burn$min) |>
    assess_burnin(state$burn$check)

  return(status)
}


#' Save Step
#' 
#' @param state a mcmodel state
#' @returns a state
#' @keywords internal
step_save = function(state) {
  if(state$phase != 'sample') return(state)
  if(state$iter %% state$thin != 0) return(state)
  
  state$blk_count = state$blk_count + 1
  i = state$blk_count
  for(p in state$save_names) {
    state$save[[p]][i, ] = state$pars[[p]]
  }

  if(!is.null(state$derv_quants)) {
    for(d in state$derv_quants) {
      state$save[[d]][i] = get(d, state$env)
    }
  }
  if(state$use_chkpts == TRUE & i >= state$block_size) {
    state = write_block(state)
  }
  return(state)
}

##############
#MARK: RUNNERS
##############

#' Intenral burn in loop
#' 
#' @param state mcmodel state object
#' @returns state after burnin
#' @keywords internal
run_burn = function(state) {
  for(i in 1:state$burn$max) {
    state$iter = state$iter + 1
    state = step_update(state)
    state = step_adapt(state)
    state = step_burn_stream(state)
    if(state$iter %% state$chkpt_freq == 0) {
      write_check(state)
    }
    if(step_check_burn(state)) {
      break
    }
  }
  state$phase = 'sample'
  state$burn$run_save = NULL
  return(state)
}

#' Intenral burn in loop
#' 
#' @param state mcmodel state object
#' @returns state at process
#' @keywords internal
run_samples = function(state, n_steps) {
  for(i in 1:n_steps) {
    state$iter = state$iter + 1
    state = step_update(state)
    state = step_save(state)
  }
  return(state)
}


###############
#MARK: SECRET HELP
###############


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

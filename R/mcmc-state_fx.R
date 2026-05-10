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
    n_burn,
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
  if(block_size == 0) {
    block_size = (n_iter - n_burn) / thin
  }

  state = structure(
    list(
      phase = 'burnin',
      n_iter = n_iter,
      thin = thin,
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
        max = n_burn,
        prop_adapt = prop_adapt,
        tuner = adapt_tuner,
        check = chain_check,
        run_save = if(burn_adapt) make_save_block(init_list, save_names, min_burn+1) else NULL
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
    ),
    class = 'mcmod.state'
  )

  # ensure that environment for needed things are all there
  environment(state$burn$tuner) = state$env
  environment(state$burn$check) = state$env

  if(state$use_chkpts) {
    dir.create(state$chk_path, showWarnings = FALSE, recursive = TRUE)
  }
  return(state)
}


#' Constructor for checkpoint paths
#' 
#' @param init_list a list for initialization
#' @param save_names pars to keep track of 
#' @param alloc how big each block needs to be
#' 
#' @returns block list strucutre
#' 
#' @keywords internal
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

####################
#MARK: State Methods 
####################

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
  state$blk_id = state$blk_id + 1
  state$blk_count = 0

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


#' Find Most Recent Save Point
#' 
#' @param chk_path path to look in
#' @returns mcmodel state objs
#' @keywords internal
load_recent = function(chk_path) {
  files = dir(chk_path)
  save_state = NULL
  # look for files
  while(length(files) > 0) {
    read_file = find_file(chk_path, files)
    if(is.null(read_file)){
      return(NULL)
    } else {
      save_state = tryCatch({
        readRDS(file.path(chk_path, read_file))
      }, error = function(e) {
        message("Corrupt File: ", read_file)
        NULL
      })
      if(!is.null(save_state)) {
        break
      } else {
        files = files[-which(files == read_file)]
      }
    }
  }
  if(is.null(save_state)) {
    return(NULL)
  }
  # empty out save or run blocks
  if(save_state$phase =='sample') {
    save_state$save = make_save_block(save_state$pars, save_state$save_names, save_state$block_size)
    if (!is.null(save_state$derv_quants)) {
      save_state$derv_save = make_save_block(NULL, save_state$derv_quants, save_state$block_size)  
    }
  }
  # set iters to 
  return(save_state)
}

#' Load Recent helper:
#' 
#' @param chk_path the path to look for file
#' @param files the file list to try and readfrom
#' @returns file path
#' @keywords internal
find_file = function(chk_path, files) {
  if(any(grepl('block_', files))) {
    nums = gsub('.rds', "", gsub('block_', "", files))
    max_num = max(suppressWarnings(as.numeric(nums)),na.rm = TRUE)
    read_file = paste0('block_', max_num, ".rds")
    # restart save block
  } else if (any(grepl('chkpt_', files))) {
    nums = gsub('.rds', "", gsub('chkpt_', "", files))
    max_num = max(suppressWarnings(as.numeric(nums)),na.rm = TRUE)
    read_file = paste0('chkpt_', max_num, ".rds")
  } else {
    return(NULL)
  }
  return(read_file)
}

#' Reconstruct states
#' 
#' @param state path to check
#' @returns completed set of chains
#' @keywords internal
load_all_blks = function(state) {
  if(state$use_chkpts) {
    files = dir(state$chk_path)
    blks = grepl('block_', files)
    if(sum(blks) == 0) {
      stop('Error - no blocks produced')
    }
    #load up save blocks
    all_blks = lapply(files[blks], function(x) readRDS(file.path(state$chk_path, x)))  
  } else {
    all_blks = list('bad_code' = state)
  }
  save_blks = all_blks |>
    lapply('[[', 'save')
  save = lapply(state$save_names, function(p) {
    lapply(save_blks, '[[', p) |>
      do.call(what = rbind, )
  })
  names(save) = state$save_names
  
  # reformat to original shape
  ns = state$block_size * length(all_blks) #num samples
  for(p in state$save_names) {
    d = dim(state$pars[[p]]) #dimensions
    kp = length(state$pars[[p]]) # K-par
    
    if(is.null(d)) {
      if(kp==1) {
          save[[p]] = save[[p]][1:ns, 1]
      } else {
          save[[p]] = save[[p]][1:ns,,drop = FALSE]
      }
    } else {
        save[[p]] = array(save[[p]], dim = c(ns, d))
    }
  }
  
  #optionally load derv blocks
  if(!is.null(state$derv_save)) {
    drv_blks = all_blks |>
      lapply('[[', 'derv_save')
    drv = lapply(state$derv_quants, function(p) {
      lapply(drv_blks, '[[', p) |>
        unlist()
    })
  } else {
    drv = NULL
  }
  names(drv) = state$derv_quants
  return(
    list(
      samples = c(save, drv),
      info = list('size' = ns)
    )
  )
}

###########
#MARK: STEPPERS
###########

#' Update Step
#' 
#' @param state a mcmodel state list
#' 
#' @returns a state object
#' @keywords internal
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
  # if in sample, update derived functions
  if(state$phase == 'sample') {
    for(fn in state$dervs) {
      res = fn()
      for(d in names(res$val)) {
        assign(d, res$val[[d]],state$env)
      }
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
  for(p in names(state$acc_counter)) {
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
      rs[[p]][i, ] = get(p, state$env)
    }
  } else {
    for(p in names(rs)) {
      rs[[p]][1:min, ] = rs[[p]][2:(min+1),]
      rs[[p]][min, ] = get(p, state$env)
    }
  }
  state$burn$run_save = rs
  return(state)
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
    state$save[[p]][i, ] = get(p, state$env)
  }

  if(!is.null(state$derv_quants)) {
    for(d in state$derv_quants) {
      state$derv_save[[d]][i] = get(d, state$env)
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
  if(state$phase == 'burnin') {
    for(i in state$iter:state$burn$max) {
      state$iter = state$iter + 1
      state = step_update(state)
      state = step_adapt(state)
      state = step_burn_stream(state)
      if(state$use_chkpts & state$iter %% state$chkpt_freq == 0) {
        write_check(state)
      }
      if(state$iter %% state$burn$prop_adapt == 0 & step_check_burn(state)) {
        break
      }
      setTxtProgressBar(state$pb, state$iter)
      flush.console()
    }
    state$phase = 'sample'
    cat('\n', 'Burn-in Complete: ', state$iter, 'iterations.',"\n")
    state$iter = state$burn$max
    setTxtProgressBar(state$pb, state$iter)
    flush.console()
    state$burn$run_save = NULL
  }
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
    setTxtProgressBar(state$pb, state$iter)
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

#' Interal print method
#' 
#' @param x the state object
#' @param ... addition args not used but for method matching
#' 
#' @returns Prints overview information from state
#' @export
#' @method print mcmod.state
print.mcmod.state = function(x, ...) {
  cat('\n', 'McModel Internal State Object', '\n')
  cat('\n', 'Phase: ', x$phase)
  cat('\n', 'Iteration: ', x$iter, ' out of ', x$n_iter)
  cat('\n', 'Block ', x$blk_id, ' out of ', (x$n_iter-x$burn$max) %/% x$block_size)
}
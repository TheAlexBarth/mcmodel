# mcmodel
A lightweight base-R framework for constructing flexible, Gibbs-MH MCMC sampling algorithms.
Rather than relying on external sampling algorithms and probabilistic coding languages, some researchers might prefer to write "homemade"
MCMC algorithms for fitting bayesian models. This is particularly the case when dealing with high dimensional, large datasets, where it is beneficial to consider some steps around memory management and sampling design. However, these homemade algorithms typically come at a cost of losing several attractive benefits from established workflows (e.g., STAN, JAGS, NIMBLE). McModel aims to bridge that gap by offering a suite of tools that make parallelization and back-up points much

# General Workflow
A brief example for using McModel is avaialble in /R/example.R, listed throughout the documentation.

## Set-up:
The design of the code is largely inspired by NIMBLE framework, where users present a series of lists that each correspond to aspects of the model run. Thus, the core set up requires all designs to specify a "data_list", "const_list", "init_list", "prior_pars". These all just get pushed into state environment (see below), so just like NIMBLE, *technically* anything provided to these lists are accessible for interior functions. While the distinction between "data" and "constants" are somewhat arbitrary, users should keep to the expected strucutre given some of the default behavior, particularly with the init_list which are used to set up structures for saving containers. Similar to other bayesian algorithm-generating software, users also must specify n_iter, n_burn, and thin. 

The big difference, and advantage, of McModel is that users must write the sampling functions for their algorithms. Inspired by STAN, sampling functions are specified as either "update_functions", which provide updates and sampling logic to the MCMC, or "derived_functions" which calculate output values in-line with MCMC processing after the burnin period. Correspondingly users must specific vectors of "save_names" and "derv_quants" respectively which are tracked throught this system. The large advantage here is that not all important variables to the model must be recorded - which is particularly useful if models are parameter intensive. When writing update or derived functions, users must treat all objects as accessible in the environment (because they will be). This may be a little uncomfortable for some programmers as the code will specify variables not accessible in the writer's environment. Additionally, all functions must return a list object with at least a "val" list which returns the name of parameter used in the model AND for MH steps, an "acc" list which matches the prop_sd list...

When writing MH algorithms, users must provide a prop_sd list, which has the standard deviations used for normal proposals. As described in the previous section, users writing MH functions, must include prop_sd as argument in the function and return acceptance outcomes! Note if users are using non-normal proposal distributions they must avoid specifying these in prop_sd, but it could work with some thoughtful tweaks. 

Finally, users can specify additional details associated with algorithm flow (e.g., checkpoint frequency, sampling adaptation, etc.)

## State Blocks
Under the hood, all the set-up objects are put into a "state" object. This is fundamentally a list that has a specified environment, in which model parameters are all accessible. A core feature of McModel design is that users can specify a "work_dir" where the algorithm will periodically write out save-points ('chkpt' during burnin and 'block' during sampling). These are advantageous because they offer the ability to start-up an MCMC run if the sampler crashes due to external reasons (e.g., server time-out). Also, for very memory intensive models, running blocks reduces the size of a posterior chain that get bulky in memory (this is substantial but can be influential in some cases). By default, McModel will write out save points at a frequency of 10\% of the total n_iter. These are then collected and deleted once the sampler finishes. However, users can keep these objects and avoid reconstruction if desired. This behavior can be turned off by setting chkpt_freq = 0.

## Adaptive sampling structure
Two default and unique advantages of McModel are an adaptive burnin phase and sampler tuning for MH.

## Thoughts
This software was written with Gibbs-MH MCMC algorithms for fitting bayesian models in mind. However, similar to NIMBLE and STAN, this software could be coereced to accomplish other models - so long as the update functions match necessary format. I'd be interested to see how others might use this code so please reach out if you are using it!


# Change log:
## Version 0.2.0
Date: 2026-03-10
### Core Updates
 - Shifted to state-block format to function with save points

### Minor updates
 - Documentation updates and changes to example code to match new format

### Update Author
Dr. Alex Barth

## Version 0.1.0
Date: 2026-03-02
### Minor Updates:
 - fixed R-hat calculation issues.
### Update Author
Dr. Alex Barth

## Version 0.1.0
Date: 2026-02-19
### Core Updates
 - Initial format launch
 - Deleted the original format which included adaptive thinning at a massive computational cost

### Update Author
Dr. Alex Barth
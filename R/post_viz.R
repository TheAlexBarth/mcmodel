
#############
#MARK: Post Chain
##############

#' Plot posterior chain trace plots
#' 
#' @param run_obj a list of named samples
#' @param par_name the name to access in samples
#' @param h_lines if known, a constant to plot across traces (typically only for simulation cases)
#' 
#' @examples
#' \dontrun{
#' mcmc_out = mcmc_run(...)
#' 
#' plot_post_chain(mcmc_out$samples, 'parname')
#' }
#' 
#' @importFrom graphics abline lines par plot.new plot.window title axis
#' 
#' @return plot in active window
#' @export 
plot_trace <- function(run_obj, par_name, h_lines = NULL) {
    dat = run_obj$samples[[par_name]] |> fafa()
    par(mfrow = c(ncol(dat), 1))
    for(c in 1:ncol(dat)) {
        plot.new()
        plot.window(
            xlim = c(0, run_obj$info$size/run_obj$info$n_chains), 
            ylim = c(min(dat[,c]), max(dat[,c]))
        )
        axis(1); axis(2)
        title(ylab = paste0(par_name, c))
        if(run_obj$info$n_chains == 1) {
            lines(dat[,c])
        } else {
            for(i in 1:run_obj$info$n_chains) {
                lines(dat[run_obj$info$chain_idx == i, c], col = i)
            }
        }
        abline(h = h_lines[c])
    }
}
 
#' Posterior Summary
#' 
#' @param run_obj output from \link{mcmc_run}
#' 
#' @examples
#' \dontrun{
#'  mcmc_out = run_mcmc(...)
#'  post_summary(run_obj)
#' }
#' 
#' @return a matrix of summary values
#' @export
post_summary = function(run_obj) {
    
    parnames = names(run_obj$samples)
    run_obj$samples[['beta']]
    core = parnames |> 
        lapply(function(x) {
            mat = fafa(run_obj$samples[[x]]) |> 
                apply(2, summarize_chain) |> 
                rbind(
                    'ESS' = run_obj$info$ESS[[x]],
                    'rhat' = run_obj$info$rhat[[x]]
                ) |> 
                t()
            rownames(mat) = paste0(x, 1:length(run_obj$info$ESS[[x]]))
            return(mat)
        }) 
    
    out = do.call(rbind, core)
    
    return(out)
}

#MARK: Summary
#' Create Chain Summary
#' 
#' @param chain a chain of single values
#' 
#' @examples
#' summarize_chain(rnorm(10000))
#' 
#' @return a vector of summary values, slightly more detailed than summarize.
#' 
#' @importFrom stats median quantile
#' 
#' @export
summarize_chain = function(chain) {
    quants = quantile(chain, probs = c(0.25, 0.75, 0.125, 0.875, 0.025,0.975))
    names(quants) = c('low.50', "high.50", "low.75", "high.75", "low.95", "high.95")
    out = c(
        'mean' = mean(chain),
        'median' = median(chain),
        quants
    )
    return(out)
}

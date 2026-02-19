
#############
#MARK: Post Chain
##############

#' Plot posterior chain trace plots
#' 
#' @param samples a list of named samples
#' @param par_name the name to access in samples
#' @param real_par if known, a constant to plot (typically only for simulation cases)
#' 
#' @examples
#' \dontrun{
#' mcmc_out = mcmc_run(...)
#' 
#' plot_post_chain(mcmc_out$samples, 'parname')
#' }
#' 
#' @importFrom graphics abline
#' 
#' @return produces a plot in active window
#' @export 
plot_trace <- function(run_obj, par_name, h_lines = NULL) {
    dat = run_obj$samples[[par_name]]
    par(mfrow = c(ncol(dat), 1))
    for(c in 1:ncol(dat)) {
        plot.new()
        plot.window(
            xlim = c(0, run_obj$info$size/run_obj$info$n_chains), 
            ylim = c(min(dat[,c]), max(dat[,c]))
        )
        axis(1); axis(2)
        if(run_obj$info$n_chains == 1) {
            lines(dat[,c])
        } else {
            for(i in 1:run_obj$info$n_chains) {
                lines(dat[run_obj$info$chain_idx == i, c], col = i)
            }
        }
    }

}
plot_post_chain <- function(samples, par_name, real_par = NULL) {
    par(mfrow = c(length(par_name), 1))
    for(par in par_name) {
        if(is.matrix(samples[[par]])) {
          par(mfrow = c(ncol(samples[[par]]), 1))
          for(c in 1:ncol(samples[[par]])) {
            ylim = c(min(samples[[par]][,c]), max(samples[[par]][,c]))
            credI = quantile(samples[[par]][,c], probs = c(0.025,0.975)) |> sapply(round, 2)
            if(!is.null(real_par[c])) {
                if(real_par[c] < ylim[1]) {
                    ylim[1] = real_par[c]
                } else if (real_par[c] > ylim[2]) {
                    ylim[2] = real_par[c]
                }
                plot(
                    samples[[par]][,c], type = 'l', 
                    ylim = ylim, ylab = par, xlab = "",
                    sub = paste0('(', credI[1], ", ", credI[2],"); ", real_par[c])
                )
                abline(h = real_par[c])    
            } else {
                plot(
                    samples[[par]][,c], type = 'l', 
                    ylim = ylim, ylab = par, xlab = "",
                    sub = paste0('(', credI[1], ", ", credI[2],"); ", real_par)
                )
            }
          }
        } else {
            ylim = c(min(samples[[par]]), max(samples[[par]]))
            credI = quantile(samples[[par]], probs = c(0.025,0.975)) |> sapply(round, 2)
            if(!is.null(real_par)) {
                if(real_par < ylim[1]) {
                    ylim[1] = real_par
                } else if (real_par > ylim[2]) {
                    ylim[2] = real_par
                }
                plot(
                    samples[[par]], type = 'l', 
                    ylim = ylim, ylab = par, xlab = "",
                    sub = paste0('(', credI[1], ", ", credI[2],"); ", real_par)
                )
                abline(h = real_par)    
            } else {
                plot(
                    samples[[par]], type = 'l', 
                    ylim = ylim, ylab = par, xlab = "",
                    sub = paste0('(', credI[1], ", ", credI[2],"); ", real_par)
                )
            }
        }
    }
    par(mfrow = c(1,1))
}
 
#' Posterior Summary
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

#' Create Chain Summary
#' 
#' 
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


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

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

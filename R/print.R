##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title Summary for BMGARCH object
##' @param object 
##' @param CrI 
##' @return Summary object
##' @author philippe
summary.bmgarch = function(object, CrI = c(0.05, 0.95) ){
    ## depending on parameteriztion, different parameters will be returned:
    if(object$param == 'CCC') {
        model_summary = rstan::summary(object$model_fit,
                          pars = c('c_h', 'a_h', 'b_h', 'R', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'),
                          probs = CrI)$summary} else {
    if(object$param == 'DCC') {
       model_summary = rstan::summary(object$model_fit,
                          pars = c('a_q', 'b_q', 'c_h', 'a_h', 'b_h', 'S', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'),
                          probs = CrI)$summary } else {
   if(object$param == 'BEKK') {
        model_bekk = rstan::summary(object$model_fit,
                          pars = c('Cnst', 'A', 'B', 'corC', 'b0', 'b1', 'b2', 'b3', 'b4', 'lp__'),
                          probs = CrI)$summary }
   }}
    cat( "Model:", paste0(object$param, "-MGARCH\n"))
    cat("---\n")
    cat("Iterations: ", object$iter, "\n")
    cat("Chains: ", object$chains, "\n")
    cat("Date: ", object$date, "\n")
    cat("Elapsed time (min):",
    round((max(object$elapsed_time[,1]) + max(object$elapsed_time[,2]))/60, 2) , "\n\n")
    print(round( model_summary, 2) )
}

#' Optimize over a function using the TAO POUNDERs optimization library.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param fn A function to be minimized (or maximized), with first argument 
#'        the vector of parameters over which minimization is to take place. 
#'        It should return a scalar result.
#' @param gr A function to return the gradient, if using a gradient-based
#'        optimization method.
#' @param hs A function to return the hessian, if using an algorithm which
#'        uses the hessian
#' @param method The method to be used. See ‘Details’.
#' @param control A list of control parameters. See ‘Details’.
#'
#' @examples
#' objfun = function(x) c((x[1] - 3), (x[2] + 1))
#' ret = tao.optim(c(1, 2), 
#'                 objfun,
#'                 method = "pounders",
#'                 control = list(tao_pounders_delta = "0.1"))
tao.optim = function(par, fn, gr = NULL, hs = NULL,
                     method = c("nm", "pounders", "lmvm", "blmvm"),
                     control = list()) {
    
    funclist = list(objFun = fn);
    if (!is.null(gr)) {
        funclist = rbind(funclist, list(jacFun = gr))
    }
    if (!is.null(hs)) {
        funclist = rbind(funclist, list(hesFun = hs))
    }
    
    ret = tao(functions = funclist,
              startValues = par,
              method = method,
              options = control)

}
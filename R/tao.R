#' Optimize over a function using the TAO optimization library.
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param fn A function to be minimized (or maximized), with first argument 
#'        the vector of parameters over which minimization is to take place. 
#'        It should return a scalar result.
#' @param gr A function to return the gradient, if using a gradient-based
#'        optimization method.
#' @param hs A function to return the hessian, if using an algorithm which
#'        uses the hessian
#' @param method The method to be used. See 'Details'.
#' @param control A list of control parameters. See 'Details'.
#' @param n The number of elements of objfun.
#' @return A list with final parameter values, the objective function, and
#'        information on why the optimizer stopped.
#'
#' @examples
#' # Gradient-free method
#' objfun = function(x) c((x[1] - 3), (x[2] + 1))
#' ret = tao.optim(c(1, 2), 
#'                 objfun,
#'                 method = "pounders",
#'                 control = list(),
#'                 n = 2)
#' 
#' # Gradient-based method
#' objfun = function(x) sum(c((x[1] - 3)^2, (x[2] + 1))^2)
#' grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
#'     
#' ret = tao.optim(c(1, 2), 
#'                 objfun,
#'                 gr = grafun,
#'                 method = "lmvm")
tao.optim = function(par, fn, gr = NULL, hs = NULL,
                     method = c("nm", "pounders", "lmvm", "blmvm"),
                     control = list(),
                     n = 1) {
    
    funclist = list(objfun = fn);
    if (!is.null(gr)) {
        funclist = c(funclist, grafun = gr)
    }
    if (!is.null(hs)) {
        funclist = c(funclist, hesfun = hs)
    }
    
    ret = tao(functions = funclist,
              startValues = par,
              method = method,
              options = control,
              n)

}
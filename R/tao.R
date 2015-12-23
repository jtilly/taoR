#  taoR -- Toolkit for Advanced Optimization (TAO) Bindings for R
#
#  Copyright (C) 2015  Jan Tilly <jtilly@econ.upenn.edu>
#                      Nick Janetos <njanetos@econ.upenn.edu>
#
#  This file is part of taoR.
#
#  taoR is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 2 of the License, or
#  (at your option) any later version.
#
#  taoR is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.

#' R bindings for the TAO optimization library.
#' 
#' Various optimization routines from the TAO optimization library. See
#' the TAO documentation for a complete listing. 
#'
#' @param par Initial values for the parameters to be optimized over.
#' @param fn A function to be minimized (or maximized), with first argument 
#'        the vector of parameters over which minimization is to take place. 
#'        It should return a scalar result.
#' @param gr A function to return the gradient, if using a gradient-based
#'        optimization method.
#' @param hs A function to return the hessian, if using an algorithm which
#'        uses the hessian.
#' @param method The method to be used. See 'Details'.
#' @param control A list of control parameters. See 'Details'.
#' @param lb A vector with lower variable bounds (optional) 
#' @param ub A vector with upper variable bounds (optional)
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
#'                 control = list(tao_pounders_delta="0.1"),
#'                 n = 2)
#' ret$x
#' 
#' # Gradient-based method: Limited memory variable metric method
#' objfun = function(x) (x[1] - 3)^2 + (x[2] + 1)^2
#' grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
#'     
#' ret = tao.optim(c(1, 2), 
#'                 objfun,
#'                 gr = grafun,
#'                 method = "lmvm")
#' ret$x
#' 
#' # Gradient-based method: Limited memory variable metric method with bounds
#' objfun = function(x) (x[1] - 3)^2 + (x[2] + 1)^2
#' grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
#' inequal = function(x) c(x[1] - 2, x[2] - 2)
#'     
#' ret = tao.optim(c(1, 2), 
#'                 objfun,
#'                 gr = grafun,
#'                 inequal = inequal,
#'                 method = "blmvm")
#' ret$x
#' 
#' # Hessian (Newton Trust Region)
#' objfun = function(x) (x[1] - 3)^2 + (x[2] + 1)^2
#' grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
#' hesfun = function(x) matrix(c(2, 0, 0, 2), nrow = 2, ncol = 2)
#'     
#' ret = tao.optim(c(1, 2), 
#'                 objfun,
#'                 gr = grafun,
#'                 hs = hesfun,
#'                 method = "ntr")
#' ret$x

tao.optim = function(par, fn, gr = NULL, hs = NULL,
                     method = c("lmvm", "nls", "ntr", "ntl", 
                                "cg", "tron", "blmvm", "gpcg",
                                "nm", "pounders"),
                     control = list(),
                     n = 1, 
                     lb = NULL, 
                     ub = NULL) {
    
    funclist = list(objfun = fn);
    
    if (!is.null(gr)) {
        funclist = c(funclist, grafun = gr)
    }
    
    if (!is.null(hs)) {
        funclist = c(funclist, hesfun = hs)
    }
    
    if (!is.null(inequal)) {
      funclist = c(funclist, inequal = inequal)
    }
    

    # if method requires gradient and none was provided, make sure
    # that tao_fd_gradient is set, i.e. that finite differences are
    # computed
    if (is.null(gr) && method %in% c("lmvm", "nls", "ntr", "ntl", 
                                     "cg", "tron", "blmvm", "gpcg")) {
        if(!("tao_fd_gradient" %in% names(control))) {
            control = c(control, list("tao_fd_gradient"="true"))
        }
    }
    
    # if method doesn't use gradient, but one was provided, throw warning
    if (!is.null(gr) && !(method %in% c("lmvm", "nls", "ntr", "ntl", 
                                        "cg", "tron", "blmvm", "gpcg"))) {
        warning("method ", method, " does not make use of user-defined gradient.")
    }
    
    # if method requires hessian and none was provided, use finite differences
    if (is.null(hs) && method %in% c("nls", "ntr", "ntl", "tron", "gpcg")) {
        stop("method ", method, " requires user-defined hessian, but non was provided.")
    }
    
    # if method doesn't use hessian, but one was provided, throw warning
    if (!is.null(hs) && !(method %in% c("nls", "ntr", "ntl", "tron", "gpcg"))) {
        warning("method ", method, " does not make use user-defined hessian.")
    }
    
    if(!is.null(lb) & length(lb) != length(par)) {
        stop("If set, the vector with lower bounds lb must have the same length as the vector par.")
    }
    
    if(!is.null(ub) & length(ub) != length(par)) {
        stop("If set, the vector with upper bounds ub must have the same length as the vector par.")
    }
    
    if(is.null(ub)) {
        ub = rep(1e16, length(par))
    }
    
    if(is.null(lb)) {
        lb = rep(-1e16, length(par))
    }
    
    ret = tao(functions = funclist,
              start_values = par,
              method = method,
              options = control,
              n, lb, ub)

}
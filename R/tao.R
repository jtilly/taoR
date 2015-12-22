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
#'                 control = list(tao_pounders_delta="0.1"),
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
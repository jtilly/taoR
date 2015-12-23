//  taoR -- Toolkit for Advanced Optimization (TAO) Bindings for R
//
//  Copyright (C) 2015  Jan Tilly <jtilly@econ.upenn.edu>
//                      Nick Janetos <njanetos@econ.upenn.edu>
//
//  This file is part of taoR.
//
//  taoR is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 2 of the License, or
//  (at your option) any later version.
//
//  taoR is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.

//
// This program is based on 
// http://www.mcs.anl.gov/petsc/petsc-current/src/tao/leastsquares/examples/tutorials/chwirut1.c.html
//
// Copyright (c) 1991-2014, UChicago Argonne, LLC and the PETSc Development Team
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
// 
// * Redistributions of source code must retain the above copyright notice, this
//   list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice, this
//   list of conditions and the following disclaimer in the documentation and/or
//   other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
// ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
// ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include <taoR.h>
#include "utils.h"
#include "evaluate.h"

//' Use TAO to minimize an objective function
//'
//' @param functions is a list of Rcpp functions. The first is always the objective 
//'        function. The second and third are optionally the Jacobian and the Hessian 
//'        functions.
//' @param start_values is a vector containing the starting values of the parameters.
//' @param method is a string that determines the type of optimizer to be used.
//' @param options is a list containing option values for the optimizer
//' @param n is the number of elements in the objective function.
//' @param lower_bounds is a vector with lower bounds
//' @param upper_bounds is a vector with upper bounds
//' @return a list with the objective function and the final parameter values
//' @examples
//' # use pounders
//' objfun = function(x) c((x[1] - 3), (x[2] + 1))
//' ret = tao(functions = list(objfun = objfun), 
//'               start_values = c(1, 2), 
//'               method = "pounders", 
//'               options = list(), 
//'               n = 2,
//'               lower_bounds = c(-2, -2),
//'               upper_bounds = c(5, 5))
//' ret$x
//'     
//' # use Nelder-Mead
//' objfun = function(x) sum(c((x[1] - 3)^2, (x[2] + 1))^2)
//' ret = tao(functions = list(objfun = objfun), 
//'                   start_values = c(1, 2), 
//'                   method = "nm", 
//'                   options = list(),
//'                   n = 1,
//'                   lower_bounds = c(-2, -2),
//'                   upper_bounds = c(5, 5))
//' ret$x
// [[Rcpp::export]]
List tao(List functions,
         NumericVector start_values, 
         String method, 
         List options, 
         int n, 
         NumericVector lower_bounds,
         NumericVector upper_bounds) {

    // Redirect output to the R console
    PetscVFPrintf = print_to_rcout;
    
    // Initialize PETSc
    initialize(options);
    
    // Problem-defined work context 
    Problem problem; 
    
    // Read in problem dimensions
    problem.n = n;
    problem.k = start_values.size();
    
    if(method != "pounders") {
        if(n > 1)  {
            stop("n must be equal 1 unless you are using Pounders.");
        }
    }
    
    // Read in the objective function
    // Add it to problem context
    Function objfun = functions["objfun"];
    problem.objfun = &objfun;
    
    // Check whether we need to read in the jacobian
    // to the problem context.
    Function grafun = functions["objfun"];
    if (functions.containsElementNamed("grafun")) {
        grafun = functions["grafun"];
        problem.grafun = &grafun;
    }
    
    // Check whether we need to read in the hessian
    // to the problem context.
    Function hesfun = functions["objfun"];
    if (functions.containsElementNamed("hesfun")) {
        hesfun = functions["hesfun"];
        problem.hesfun = &hesfun;
    }
    
    Vec x, f; // solution, function
    Vec ub, lb, ci; // upper and lower bounds
    Tao tao_context; // Tao solver context 
    PetscReal fc, gnorm, cnorm, xdiff;
    PetscInt its;

    // Allocate vectors
    catch_error(VecCreateSeq(MPI_COMM_SELF, start_values.size(), &x));
    catch_error(VecCreateSeq(MPI_COMM_SELF, lower_bounds.size(), &lb));
    catch_error(VecCreateSeq(MPI_COMM_SELF, upper_bounds.size(), &ub));
    catch_error(VecCreateSeq(MPI_COMM_SELF, n, &f));
    
    // Create TAO solver
    catch_error(TaoCreate(PETSC_COMM_SELF, &tao_context));
    catch_error(TaoSetType(tao_context, method.get_cstring()));
    
    // Form starting values and define functions
    catch_error(createVec(x, start_values));
    catch_error(TaoSetInitialVector(tao_context, x));
    
    // Form lower, upper bounds vectors
    catch_error(createVec(lb, lower_bounds));
    catch_error(createVec(ub, upper_bounds));
    catch_error(createVec(ci, problem.k));
    
    // Check whether lower / upper bounds inequalities have been set
    Function inequal = functions["inequal"];
    if (functions.containsElementNamed("inequal")) {
        problem.inequal = &inequal;
        catch_error(TaoSetInequalityConstraintsRoutine(tao_context, ci, evaluate_inequalities, &problem));
    }
    
    // Check whether lower / upper bounds equalities have been set
    Function equal = functions["equal"];
    if (functions.containsElementNamed("equal")) {
        problem.equal = &equal;
        catch_error(TaoSetConstraintsRoutine(tao_context, ci, evaluate_equalities, &problem));
    }
    
    // Create a matrix to hold hessians
    Mat H;
    MatCreate(PETSC_COMM_WORLD, &H);
    MatSetSizes(H, PETSC_DECIDE, PETSC_DECIDE, problem.k, problem.k);
    MatSetUp(H);
    
    // Define objective functions and gradients
    if(method == "pounders") {
        catch_error(TaoSetSeparableObjectiveRoutine(tao_context, f, evaluate_objective_separable, (void*)&problem));
    } else {
        catch_error(TaoSetObjectiveRoutine(tao_context, evaluate_objective, (void*)&problem));
    }
    
    if (functions.containsElementNamed("grafun")) {
        catch_error(TaoSetGradientRoutine(tao_context, evaluate_gradient, (void*)&problem));
    }
    
    if (functions.containsElementNamed("hesfun")) {
        catch_error(TaoSetHessianRoutine(tao_context, H, H, evaluate_hessian, (void*)&problem));
    }
    
    // Set variable bounds
    catch_error(TaoSetVariableBounds(tao_context, lb, ub));
    
    // Define monitor
    catch_error(TaoSetMonitor(tao_context, my_monitor, &problem, NULL));
    
    // Check for any TAO command line arguments 
    catch_error(TaoSetFromOptions(tao_context));
    
    // Perform the Solve
    catch_error(TaoSolve(tao_context));
    catch_error(TaoView(tao_context, PETSC_VIEWER_STDOUT_SELF));
    catch_error(TaoGetSolutionStatus(tao_context, &its, &fc, &gnorm, &cnorm, &xdiff, 0));
    
    // Free TAO data structures
    catch_error(TaoDestroy(&tao_context));
    NumericVector xVec(start_values.size());
    xVec = get_vec(x, start_values.size());
    catch_error(VecDestroy(&x));
    
    NumericVector fVec(n);
    if(method == "pounders") {
        fVec = get_vec(f, n);
        catch_error(VecDestroy(&f));
    } else {
        fVec[0] = fc;
    }
    
    return List::create( 
        Named("x")  = xVec,
        Named("f")  = fVec,
        Named("iterations")  = its,
        Named("gnorm")  = gnorm,
        Named("cnorm")  = cnorm,
        Named("xdiff")  = xdiff
    );
    
}
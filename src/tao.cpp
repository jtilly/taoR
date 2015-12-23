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

#include "utils.h"
#include "evaluate.h"

// Forward declarations
PetscErrorCode my_monitor(Tao, void*);
PetscErrorCode print_to_rcout(FILE*, const char*, va_list);
NumericVector get_vec(Vec, int);

//' Use TAO to minimize an objective function
//'
//' @param functions is a list of Rcpp functions. The first is always the objective 
//'        function. The second and third are optionally the Jacobian and the Hessian 
//'        functions.
//' @param start_values is a vector containing the starting values of the parameters.
//' @param method is a string that determines the type of optimizer to be used.
//' @param options is a list containing option values for the optimizer
//' @param n is the number of elements in the objective function.
//' @param lb is a vector with lower bounds
//' @param ub is a vector with upper bounds
//' @return a list with the objective function and the final parameter values
//' @examples
//' # use pounders
//' objfun = function(x) c((x[1] - 3), (x[2] + 1))
//' ret = tao(functions = list(objfun = objfun), 
//'               start_values = c(1, 2), 
//'               method = "pounders", 
//'               options = list(), 
//'               n = 2)
//' ret$x
//'     
//' # use Nelder-Mead
//' objfun = function(x) sum(c((x[1] - 3)^2, (x[2] + 1))^2)
//' ret = tao(functions = list(objfun = objfun), 
//'                   start_values = c(1, 2), 
//'                   method = "nm", 
//'                   options = list())
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
    
    PetscErrorCode error_code; // used to check for functions returning nonzeros 
    Vec x, f; // solution, function
    Vec ub, lb; // upper and lower bounds
    Tao tao_context; // Tao solver context 
    PetscReal fc, gnorm, cnorm, xdiff;
    PetscInt its;

    // Allocate vectors
    error_code = VecCreateSeq(MPI_COMM_SELF, start_values.size(), &x); CHKERRQ(error_code);
    error_code = VecCreateSeq(MPI_COMM_SELF, lower_bounds.size(), &lb); CHKERRQ(error_code);
    error_code = VecCreateSeq(MPI_COMM_SELF, upper_bounds.size(), &ub); CHKERRQ(error_code);
    error_code = VecCreateSeq(MPI_COMM_SELF, n, &f); CHKERRQ(error_code);
    
    // Create TAO solver
    error_code = TaoCreate(PETSC_COMM_SELF, &tao_context); CHKERRQ(error_code);
    error_code = TaoSetType(tao_context, method.get_cstring()); CHKERRQ(error_code);
    
    // Define starting values and define functions
    error_code = createVec(x, start_values); CHKERRQ(error_code);
    error_code = TaoSetInitialVector(tao_context, x); CHKERRQ(error_code);
    
    // Define lower bounds
    error_code = createVec(lb, lower_bounds); CHKERRQ(error_code);
    error_code = createVec(ub, upper_bounds); CHKERRQ(error_code);
    
    // Create a matrix to hold hessians
    Mat H;
    MatCreate(PETSC_COMM_WORLD, &H);
    MatSetSizes(H, PETSC_DECIDE, PETSC_DECIDE, problem.k, problem.k);
    MatSetUp(H);
    
    // Define objective functions and gradients
    if(method == "pounders") {
        error_code = TaoSetSeparableObjectiveRoutine(tao_context, f, evaluate_objective_separable, (void*)&problem); CHKERRQ(error_code);
    } else {
        error_code = TaoSetObjectiveRoutine(tao_context, evaluate_objective, (void*)&problem); CHKERRQ(error_code);
    }
    
    if (functions.containsElementNamed("grafun")) {
        error_code = TaoSetGradientRoutine(tao_context, evaluate_gradient, (void*)&problem); CHKERRQ(error_code);
    }
    
    if (functions.containsElementNamed("hesfun")) {
        error_code = TaoSetHessianRoutine(tao_context, H, H, evaluate_hessian, (void*)&problem); CHKERRQ(error_code);
    }
    
    // set variable bounds
    error_code = TaoSetVariableBounds(tao_context, lb, ub); CHKERRQ(error_code);
    
    // Define monitor
    error_code = TaoSetMonitor(tao_context, my_monitor, &problem, NULL); CHKERRQ(error_code);
    
    // Check for any TAO command line arguments 
    error_code = TaoSetFromOptions(tao_context); CHKERRQ(error_code);
    
    // Perform the Solve
    error_code = TaoSolve(tao_context); CHKERRQ(error_code);
    error_code = TaoView(tao_context, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(error_code);
    error_code = TaoGetSolutionStatus(tao_context, &its, &fc, &gnorm, &cnorm, &xdiff, 0);
    
    // Free TAO data structures
    error_code = TaoDestroy(&tao_context); CHKERRQ(error_code);
    NumericVector xVec(start_values.size());
    xVec = get_vec(x, start_values.size());
    error_code = VecDestroy(&x); CHKERRQ(error_code);
    
    NumericVector fVec(n);
    if(method == "pounders") {
        fVec = get_vec(f, n);
        error_code = VecDestroy(&f); CHKERRQ(error_code);
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
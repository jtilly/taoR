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

// problem structure
typedef struct {
    Rcpp::Function *objfun;
    Rcpp::Function *grafun;
    Rcpp::Function *hesfun;
    int k;
    int n;
} Problem;

// Forward declarations
PetscErrorCode FormStartingPoint(Vec, Rcpp::NumericVector);
PetscErrorCode EvaluateSeparableFunction(Tao, Vec, Vec, void *);
PetscErrorCode EvaluateFunction(Tao, Vec, PetscReal*, void *);
PetscErrorCode EvaluateGradient(Tao, Vec, Vec, void *);
PetscErrorCode MyMonitor(Tao, void*);
Rcpp::NumericVector getVec(Vec, int);
Rcpp::Environment base("package:base");

//' Use Pounders to minimize a non-linear sum of squares problem 
//'
//' @param functions is a list of Rcpp functions. The first is always the objective 
//'        function. The second and third are optionally the Jacobian and the Hessian 
//'        functions.
//' @param startValues is a vector containing the starting values of the parameters.
//' @param method is a string that determines the type of optimizer to be used.
//' @param options is a list containing option values for the optimizer
//' @param n is the number of elements in the objective function.
//' @return a list with the objective function and the final parameter values
//' @examples
//' # use pounders
//' objfun = function(x) c((x[1] - 3), (x[2] + 1))
//'     ret = tao(functions = list(objfun = objfun), 
//'               startValues = c(1, 2), 
//'               method = "pounders", 
//'               options = list(tao_pounders_npmax = "1", tao_pounders_delta = "0.2"), 
//'               n = 2)
//'     ret$x
//'     
//' # use Nelder-Mead
//'     objfun = function(x) sum(c((x[1] - 3)^2, (x[2] + 1))^2)
//'         ret = tao(functions = list(objfun = objfun), 
//'                   startValues = c(1, 2), 
//'                   method = "nm", 
//'                   options = list())
//'         ret$x
// [[Rcpp::export]]
Rcpp::List tao(Rcpp::List functions,
               Rcpp::NumericVector startValues, 
               std::string method, 
               Rcpp::List options, 
               int n) {

    // Initialize PETSc
    petscInitialize(options);
    
    // Problem-defined work context 
    Problem problem; 
    
    // Read in problem dimensions
    problem.n = n;
    problem.k = startValues.size();
    
    if(method != "pounders") {
        if(n > 1)  {
            Rcpp::stop("n must be equal 1 unless you are using Pounders.");
        }
    }
    
    // Read in the objective function
    // Add it to problem context
    Rcpp::Function objfun = functions["objfun"];
    problem.objfun = &objfun;
    
    // Check whether we need to read in the jacobian
    // to the problem context.
    Rcpp::Function grafun = base["identity"]; 
    if (functions.containsElementNamed("grafun")) {
        grafun = functions["grafun"];
        problem.grafun = &grafun;
    }
    
    // Check whether we need to read in the hessian
    // to the problem context.
    Rcpp::Function hesfun = base["identity"];
    if (functions.containsElementNamed("hesfun")) {
        hesfun = functions["hesfun"];
        problem.hesfun = &hesfun;
    }
    
    PetscErrorCode ierr; // used to check for functions returning nonzeros 
    Vec x, f; // solution, function 
    Tao tao; // Tao solver context 
    PetscReal fc, gnorm, cnorm, xdiff;
    PetscInt its;

    // allocate vectors
    ierr = VecCreateSeq(MPI_COMM_SELF, startValues.size(), &x); CHKERRQ(ierr);
    ierr = VecCreateSeq(MPI_COMM_SELF, n, &f); CHKERRQ(ierr);
    
    // Create TAO solver
    ierr = TaoCreate(PETSC_COMM_SELF, &tao); CHKERRQ(ierr);
    ierr = TaoSetType(tao, method.c_str()); CHKERRQ(ierr);
    
    // Define starting values and define functions
    ierr = FormStartingPoint(x, startValues); CHKERRQ(ierr);
    ierr = TaoSetInitialVector(tao, x); CHKERRQ(ierr);
    
    // Define objective functions and gradients
    ierr = TaoSetSeparableObjectiveRoutine(tao, f, EvaluateSeparableFunction, (void*)&problem); CHKERRQ(ierr);
    ierr = TaoSetObjectiveRoutine(tao, EvaluateFunction, (void*)&problem); CHKERRQ(ierr);
    ierr = TaoSetGradientRoutine(tao, EvaluateGradient, (void*)&problem); CHKERRQ(ierr);

    // Define monitor
    ierr = TaoSetMonitor(tao, MyMonitor, &problem, NULL); CHKERRQ(ierr);
    
    // Check for any TAO command line arguments 
    ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);
    
    // Perform the Solve
    ierr = TaoSolve(tao); CHKERRQ(ierr);
    ierr = TaoView(tao, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
    ierr = TaoGetSolutionStatus(tao, &its, &fc, &gnorm, &cnorm, &xdiff, 0);
    
    // Free TAO data structures
    ierr = TaoDestroy(&tao); CHKERRQ(ierr);
    
    Rcpp::NumericVector xVec(startValues.size());
    xVec = getVec(x, startValues.size());
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    
    Rcpp::NumericVector fVec(n);
    if(method == "pounders") {
        fVec = getVec(f, n);
        ierr = VecDestroy(&f); CHKERRQ(ierr);
    } else {
        fVec[0] = fc;
    }
    
    // PetscFinalize();
    
    return Rcpp::List::create( 
        Rcpp::Named("x")  = xVec,
        Rcpp::Named("f")  = fVec,
        Rcpp::Named("iterations")  = its,
        Rcpp::Named("gnorm")  = gnorm,
        Rcpp::Named("cnorm")  = cnorm,
        Rcpp::Named("xdiff")  = xdiff
    );
    
}

// this function transforms a vector of type Vec to a vector of type
// Rcpp::NumericVector
Rcpp::NumericVector getVec(Vec X, int k) {
    PetscInt i;
    PetscReal *x;
    VecGetArray(X, &x);
    Rcpp::NumericVector xVec(k);
    for (i=0; i < k; i++) {
        xVec[i] = x[i];
    }
    return xVec;
}

// this function evaluates the separable objective function
PetscErrorCode EvaluateSeparableFunction(Tao tao, Vec X, Vec F, void *ptr) {
    
    Problem *problem = (Problem *)ptr;
    PetscInt i;
    PetscReal *x,*f;
    PetscErrorCode ierr;
    Rcpp::Function objfun = *(problem->objfun);
    int n = problem->n;
    int k = problem->k;
    
    PetscFunctionBegin;
    ierr = VecGetArray(X, &x); CHKERRQ(ierr);
    ierr = VecGetArray(F, &f); CHKERRQ(ierr);
    
    Rcpp::NumericVector xVec(k);
    Rcpp::NumericVector fVec(n);
    
    for (i=0; i < k; i++) {
        xVec[i] = x[i];
    }
    
    fVec = objfun(xVec);
    
    for (i=0; i < n; i++) {
        f[i] = fVec[i];
    }
    
    ierr = VecRestoreArray(X, &x); CHKERRQ(ierr);
    ierr = VecRestoreArray(F, &f); CHKERRQ(ierr);
    
    PetscFunctionReturn(0);
}

// this function evaluates the objective function
PetscErrorCode EvaluateFunction(Tao tao, Vec X, PetscReal *f, void *ptr) {

    Problem *problem = (Problem *)ptr;
    PetscInt i;
    PetscReal *x;
    PetscErrorCode ierr;
    Rcpp::Function objfun = *(problem->objfun);
    int k = problem->k;
    
    PetscFunctionBegin;
    ierr = VecGetArray(X, &x); CHKERRQ(ierr);
    
    Rcpp::NumericVector xVec(k);
    Rcpp::NumericVector fVec(1);
    
    for (i=0; i < k; i++) {
        xVec[i] = x[i];
    }
    
    fVec = objfun(xVec);
    *f = fVec[0];
    
    ierr = VecRestoreArray(X, &x); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}


// this function evaluates the gradient
PetscErrorCode EvaluateGradient(Tao tao, Vec X, Vec G, void *ptr) {
    
    Problem *problem = (Problem *)ptr;
    PetscInt i;
    PetscReal *x;
    PetscReal *g;
    PetscErrorCode ierr;
    Rcpp::Function grafun = *(problem->grafun);
    int k = problem->k;
    
    PetscFunctionBegin;
    ierr = VecGetArray(X, &x); CHKERRQ(ierr);
    ierr = VecGetArray(G, &g); CHKERRQ(ierr);
    
    Rcpp::NumericVector xVec(k);
    Rcpp::NumericVector gVec(k);
    
    for (i=0; i < k; i++) {
        xVec[i] = x[i];
    }
    
    gVec = grafun(xVec);

    for (i=0; i < k; i++) {
        g[i] = gVec[i];
    }
    
    ierr = VecRestoreArray(X, &x); CHKERRQ(ierr);
    ierr = VecRestoreArray(G, &g); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

// this function set the starting value
PetscErrorCode FormStartingPoint(Vec X, Rcpp::NumericVector startValues) {
    
    PetscReal *x;
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    ierr = VecGetArray(X,&x); CHKERRQ(ierr);
    for(int iX=0; iX<startValues.size(); iX++) {
        x[iX] = startValues[iX];
    }
    VecRestoreArray(X,&x); CHKERRQ(ierr);
    PetscFunctionReturn(0);
}

// this function reports the progress of the optimizer
PetscErrorCode MyMonitor(Tao tao, void *ptr) {
    
    PetscReal fc, gnorm;
    PetscInt its;
    PetscViewer viewer = PETSC_VIEWER_STDOUT_SELF;
    PetscErrorCode ierr;
    
    PetscFunctionBegin;
    ierr = TaoGetSolutionStatus(tao, &its, &fc, &gnorm, 0, 0, 0);
    ierr = PetscViewerASCIIPrintf(viewer, "iter = %3D,", its); CHKERRQ(ierr);
    ierr = PetscViewerASCIIPrintf(viewer, " Function value %g,", (double) fc); CHKERRQ(ierr);
    if (gnorm > 1.e-6) {
        ierr = PetscViewerASCIIPrintf(viewer, " Residual: %g \n", (double) gnorm); CHKERRQ(ierr);
    } else if (gnorm > 1.e-11) {
        ierr = PetscViewerASCIIPrintf(viewer, " Residual: < 1.0e-6 \n"); CHKERRQ(ierr);
    } else {
        ierr = PetscViewerASCIIPrintf(viewer, " Residual: < 1.0e-11 \n"); CHKERRQ(ierr);
    }
    PetscFunctionReturn(0);
}

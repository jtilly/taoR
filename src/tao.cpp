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
    Function *objfun;
    Function *grafun;
    Function *hesfun;
    int k;
    int n;
} Problem;

// Forward declarations
PetscErrorCode form_starting_point(Vec, NumericVector);
PetscErrorCode evaluate_objective_separable(Tao, Vec, Vec, void *);
PetscErrorCode evaluate_objective(Tao, Vec, PetscReal*, void *);
PetscErrorCode evaluate_gradient(Tao, Vec, Vec, void *);
PetscErrorCode evaluate_hessian(Tao, Vec, Mat, Mat, void*);
PetscErrorCode my_monitor(Tao, void*);
PetscErrorCode print_to_rcout(FILE*, const char*, va_list);
NumericVector get_vec(Vec, int);
Environment base("package:base");

//' Use TAO to minimize an objective function
//'
//' @param functions is a list of Rcpp functions. The first is always the objective 
//'        function. The second and third are optionally the Jacobian and the Hessian 
//'        functions.
//' @param start_values is a vector containing the starting values of the parameters.
//' @param method is a string that determines the type of optimizer to be used.
//' @param options is a list containing option values for the optimizer
//' @param n is the number of elements in the objective function.
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
               int n = 1) {

    // Redirect output to the R console
    PetscVFPrintf = print_to_rcout;
    
    // Initialize PETSc
    petsc_initialize(options);
    
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
    Function grafun = base["identity"]; 
    if (functions.containsElementNamed("grafun")) {
        grafun = functions["grafun"];
        problem.grafun = &grafun;
    }
    
    // Check whether we need to read in the hessian
    // to the problem context.
    Function hesfun = base["identity"];
    if (functions.containsElementNamed("hesfun")) {
        hesfun = functions["hesfun"];
        problem.hesfun = &hesfun;
    }
    
    PetscErrorCode error_code; // used to check for functions returning nonzeros 
    Vec x, f; // solution, function
    Tao tao_context; // Tao solver context 
    PetscReal fc, gnorm, cnorm, xdiff;
    PetscInt its;

    // Allocate vectors
    error_code = VecCreateSeq(MPI_COMM_SELF, start_values.size(), &x); CHKERRQ(error_code);
    error_code = VecCreateSeq(MPI_COMM_SELF, n, &f); CHKERRQ(error_code);
    
    // Create TAO solver
    error_code = TaoCreate(PETSC_COMM_SELF, &tao_context); CHKERRQ(error_code);
    error_code = TaoSetType(tao_context, method.get_cstring()); CHKERRQ(error_code);
    
    // Define starting values and define functions
    error_code = form_starting_point(x, start_values); CHKERRQ(error_code);
    error_code = TaoSetInitialVector(tao_context, x); CHKERRQ(error_code);
    
    // Create a matrix to hold hessians
    Mat H;
    MatCreate(PETSC_COMM_WORLD, &H);
    MatSetSizes(H, PETSC_DECIDE, PETSC_DECIDE, problem.k, problem.k);
    MatSetUp(H);
    
    // Define objective functions and gradients
    error_code = TaoSetSeparableObjectiveRoutine(tao_context, f, evaluate_objective_separable, (void*)&problem); CHKERRQ(error_code);
    error_code = TaoSetObjectiveRoutine(tao_context, evaluate_objective, (void*)&problem); CHKERRQ(error_code);
    error_code = TaoSetGradientRoutine(tao_context, evaluate_gradient, (void*)&problem); CHKERRQ(error_code);
    error_code = TaoSetHessianRoutine(tao_context, H, H, evaluate_hessian, (void*)&problem); CHKERRQ(error_code);

    // Define monitor
    error_code = TaoSetMonitor(tao_context, my_monitor, &problem, NULL); CHKERRQ(error_code);
    
    // Check for any TAO command line arguments 
    error_code = TaoSetFromOptions(tao_context); CHKERRQ(error_code);
    
    // Perform the Solve
    // The Solve is here
    // Yes yes here is the Solve rejoice!
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
    
    // PetscFinalize();
    
    return List::create( 
        Named("x")  = xVec,
        Named("f")  = fVec,
        Named("iterations")  = its,
        Named("gnorm")  = gnorm,
        Named("cnorm")  = cnorm,
        Named("xdiff")  = xdiff
    );
    
}

// this function transforms a vector of type Vec to a vector of type
// NumericVector
NumericVector get_vec(Vec X, int k) {
    PetscErrorCode error_code;
    PetscReal *x;
    VecGetArray(X, &x);
    NumericVector xVec(k);
    for (int i = 0; i < k; i++) {
        xVec[i] = x[i];
    }
    error_code = VecRestoreArray(X, &x); CHKERRQ(error_code);
    return xVec;
}

// this function evaluates the separable objective function
PetscErrorCode evaluate_objective_separable(Tao tao_context, Vec X, Vec F, void *ptr) {
    Problem *problem = (Problem *)ptr;
    PetscReal *x,*f;
    PetscErrorCode error_code;
    Function objfun = *(problem->objfun);
    int n = problem->n;
    int k = problem->k;
    
    PetscFunctionBegin;
    error_code = VecGetArray(X, &x); CHKERRQ(error_code);
    error_code = VecGetArray(F, &f); CHKERRQ(error_code);
    
    NumericVector xVec(k);
    NumericVector fVec(n);
    
    for (int i = 0; i < k; i++) {
        xVec[i] = x[i];
    }
    
    fVec = objfun(xVec);
    
    for (int i = 0; i < n; i++) {
        f[i] = fVec[i];
    }
    
    error_code = VecRestoreArray(X, &x); CHKERRQ(error_code);
    error_code = VecRestoreArray(F, &f); CHKERRQ(error_code);
    
    PetscFunctionReturn(0);
}

// this function evaluates the objective function
PetscErrorCode evaluate_objective(Tao tao_context, Vec X, PetscReal *f, void *ptr) {

    Problem *problem = (Problem *)ptr;
    PetscReal *x;
    PetscErrorCode error_code;
    Function objfun = *(problem->objfun);
    int k = problem->k;
    
    PetscFunctionBegin;
    error_code = VecGetArray(X, &x); CHKERRQ(error_code);
    
    NumericVector xVec = get_vec(X, k);
    NumericVector fVec(1);
    
    for (int i = 0; i < k; i++) {
        xVec[i] = x[i];
    }
    
    fVec = objfun(xVec);
    *f = fVec[0];
    
    error_code = VecRestoreArray(X, &x); CHKERRQ(error_code);
    PetscFunctionReturn(0);
}


// this function evaluates the gradient
PetscErrorCode evaluate_gradient(Tao tao_context, Vec X, Vec G, void *ptr) {
    
    Problem *problem = (Problem *)ptr;
    PetscReal *x;
    PetscReal *g;
    PetscErrorCode error_code;
    Function grafun = *(problem->grafun);
    int k = problem->k;
    
    PetscFunctionBegin;
    error_code = VecGetArray(X, &x); CHKERRQ(error_code);
    error_code = VecGetArray(G, &g); CHKERRQ(error_code);
    
    NumericVector xVec(k);
    NumericVector gVec(k);
    
    for (int i = 0; i < k; i++) {
        xVec[i] = x[i];
    }
    
    gVec = grafun(xVec);

    for (int i = 0; i < k; i++) {
        g[i] = gVec[i];
    }
    
    error_code = VecRestoreArray(X, &x); CHKERRQ(error_code);
    error_code = VecRestoreArray(G, &g); CHKERRQ(error_code);
    PetscFunctionReturn(0);
}

// this function evaluates the hessian
PetscErrorCode evaluate_hessian(Tao tao_context, Vec X, Mat H, Mat Hpre, void *ptr) {
    
    Problem *problem = (Problem *)ptr;
    PetscReal *x;
    PetscErrorCode error_code;
    Function hesfun = *(problem->hesfun);
    int k = problem->k;
    
    PetscFunctionBegin;
    
    error_code = VecGetArray(X, &x); CHKERRQ(error_code);

    NumericVector xVec(k);
    NumericMatrix hMat(k, k);
    
    for (int i = 0; i < k; i++) {
        xVec[i] = x[i];
    }
    
    hMat = hesfun(xVec);
    
    // Assemble the matrix
    for (int row = 0; row < k; ++row) {
        for (int col = 0; col < k; ++col) {
            MatSetValues(H, 1, &row, 1, &col, &(hMat(row, col)), INSERT_VALUES);
        }
    }
    
    error_code =  MatAssemblyBegin(H, MAT_FINAL_ASSEMBLY); CHKERRQ(error_code);
    error_code =  MatAssemblyEnd(H, MAT_FINAL_ASSEMBLY); CHKERRQ(error_code);
    error_code = VecRestoreArray(X, &x); CHKERRQ(error_code);

    PetscFunctionReturn(0);
}


// this function set the starting value
PetscErrorCode form_starting_point(Vec X, NumericVector start_values) {
    
    PetscReal *x;
    PetscErrorCode error_code;
    
    PetscFunctionBegin;
    error_code = VecGetArray(X,&x); CHKERRQ(error_code);
    for(int iX=0; iX<start_values.size(); iX++) {
        x[iX] = start_values[iX];
    }
    VecRestoreArray(X,&x); CHKERRQ(error_code);
    PetscFunctionReturn(0);
}

// this function reports the progress of the optimizer
PetscErrorCode my_monitor(Tao tao_context, void *ptr) {
    
    PetscReal fc, gnorm;
    PetscInt its;
    PetscViewer viewer = PETSC_VIEWER_STDOUT_SELF;
    PetscErrorCode error_code;
    
    PetscFunctionBegin;
    error_code = TaoGetSolutionStatus(tao_context, &its, &fc, &gnorm, 0, 0, 0);
    error_code = PetscViewerASCIIPrintf(viewer, "iter = %3D,", its); CHKERRQ(error_code);
    error_code = PetscViewerASCIIPrintf(viewer, " Function value %g,", (double) fc); CHKERRQ(error_code);
    if (gnorm > 1.e-6) {
        error_code = PetscViewerASCIIPrintf(viewer, " Residual: %g \n", (double) gnorm); CHKERRQ(error_code);
    } else if (gnorm > 1.e-11) {
        error_code = PetscViewerASCIIPrintf(viewer, " Residual: < 1.0e-6 \n"); CHKERRQ(error_code);
    } else {
        error_code = PetscViewerASCIIPrintf(viewer, " Residual: < 1.0e-11 \n"); CHKERRQ(error_code);
    }
    PetscFunctionReturn(0);
}

// Checks if output is going to stdout or stderr, if so, redirects to Rcout or Rcerr.
// Overrides PetscVFPrintf.
PetscErrorCode print_to_rcout(FILE *file, const char format[], va_list argp) {
    PetscErrorCode error_code;
    
    PetscFunctionBegin;
    if (file != stdout && file != stderr) {
        error_code = PetscVFPrintfDefault(file, format, argp); CHKERRQ(error_code);
    } else if (file == stdout) {
        char buff[1024];
        size_t length;
        error_code = PetscVSNPrintf(buff, 1024, format, &length, argp); CHKERRQ(error_code);
        Rcout << buff;
    } else if (file == stderr) {
        char buff[1024];
        size_t length;
        error_code = PetscVSNPrintf(buff, 1024, format, &length, argp); CHKERRQ(error_code);
        Rcerr << buff;
    }
    PetscFunctionReturn(0);
}
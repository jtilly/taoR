#include <Rcpp.h>
#include <petsctao.h>

// problem structure
typedef struct {
    Rcpp::Function *objFun;
    int k;
    int n;
} Problem;

// Forward declarations
PetscErrorCode FormStartingPoint(Vec, Rcpp::NumericVector);
PetscErrorCode EvaluateFunction(Tao, Vec, Vec, void *);
PetscErrorCode MyMonitor(Tao, void*);
Rcpp::NumericVector getVec(Vec, int);

//' Use Pounders to minimize a non-linear sum of squares problem 
//'
//' @param objFun is an R objective function that maps k parameters into n equations.
//' @param startValues is a vector with k elements
//' @param k is the number of parameters
//' @param n is the number of elements in the objective function
//' @return a list with the objective function and the final parameter values
//' @examples
//' objfun = function(x) c((x[1] - 3), (x[2] + 1))
//' ret = pounders(objfun, c(1,2), 2, 2)
//' ret$x
// [[Rcpp::export]]
Rcpp::List pounders(Rcpp::Function objFun, Rcpp::NumericVector startValues, int k, int n) {
    
    // create command line arguments
    char* dummy_args[] = {NULL};
    int argc = sizeof(dummy_args) / sizeof(dummy_args[0]) - 1;
    char** argv = dummy_args;
    
    PetscErrorCode ierr; // used to check for functions returning nonzeros 
    Vec x, f; // solution, function 
    Tao tao; // Tao solver context 
    PetscInt i; // iteration information 
    Problem problem; // problem-defined work context 
    
    PetscInitialize(&argc, &argv, (char *)0, (char *)0);
    
    // allocate vectors
    ierr = VecCreateSeq(MPI_COMM_SELF, k, &x); CHKERRQ(ierr);
    ierr = VecCreateSeq(MPI_COMM_SELF, n, &f); CHKERRQ(ierr);
    
    // add objective function to problem
    problem.objFun = &objFun;
    problem.n = n;
    problem.k = k;
    
    // Create TAO solver
    ierr = TaoCreate(PETSC_COMM_SELF, &tao); CHKERRQ(ierr);
    ierr = TaoSetType(tao, TAOPOUNDERS); CHKERRQ(ierr);
    
    // Define starting values and define functions
    ierr = FormStartingPoint(x, startValues); CHKERRQ(ierr);
    ierr = TaoSetInitialVector(tao, x); CHKERRQ(ierr);
    ierr = TaoSetSeparableObjectiveRoutine(tao, f, EvaluateFunction, (void*)&problem); CHKERRQ(ierr);
    
    // define monitor
    ierr = TaoSetMonitor(tao, MyMonitor, &problem, NULL); CHKERRQ(ierr);
    
    // Check for any TAO command line arguments 
    ierr = TaoSetFromOptions(tao); CHKERRQ(ierr);
    
    // Perform the Solve
    ierr = TaoSolve(tao); CHKERRQ(ierr);
    ierr = TaoView(tao, PETSC_VIEWER_STDOUT_SELF); CHKERRQ(ierr);
    
    // Free TAO data structures
    ierr = TaoDestroy(&tao); CHKERRQ(ierr);
    
    Rcpp::NumericVector xVec(k);
    xVec = getVec(x, k);
    Rcpp::NumericVector fVec(n);
    fVec = getVec(f, n);
    
    // Free PETSc data structures
    ierr = VecDestroy(&x); CHKERRQ(ierr);
    ierr = VecDestroy(&f); CHKERRQ(ierr);
    
    //PetscFinalize();
    return Rcpp::List::create( 
        Rcpp::Named("x")  = xVec,
        Rcpp::Named("f")  = fVec
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

// this function evaluates the objective function
PetscErrorCode EvaluateFunction(Tao tao, Vec X, Vec F, void *ptr) {
    
    Problem *problem = (Problem *)ptr;
    PetscInt i;
    PetscReal *x,*f;
    PetscErrorCode ierr;
    Rcpp::Function objFun = *(problem->objFun);
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
    
    fVec = objFun(xVec);
    
    for (i=0; i < n; i++) {
        f[i] = fVec[i];
    }
    
    ierr = VecRestoreArray(X, &x); CHKERRQ(ierr);
    ierr = VecRestoreArray(F, &f); CHKERRQ(ierr);
    
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
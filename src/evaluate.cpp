#include "evaluate.h"

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
PetscErrorCode createVec(Vec X, NumericVector start_values) {
    
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

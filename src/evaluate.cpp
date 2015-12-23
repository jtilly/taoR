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
PetscErrorCode createVec(Vec X, NumericVector y) {
    
    PetscReal *x;
    PetscErrorCode error_code;
    
    PetscFunctionBegin;
    error_code = VecGetArray(X,&x); CHKERRQ(error_code);
    for(int iX=0; iX<y.size(); iX++) {
        x[iX] = y[iX];
    }
    VecRestoreArray(X,&x); CHKERRQ(error_code);
    PetscFunctionReturn(0);
}

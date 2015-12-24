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

// this function evaluates the separable objective function
PetscErrorCode evaluate_objective_separable(Tao tao_context, Vec X, Vec F, void *ptr) {
    Problem *problem = (Problem *)ptr;
    Function objfun = *(problem->objfun);
    int n = problem->n;
    int k = problem->k;
    
    return evaluate_function(X, F, &objfun, k, n);
}

// this function evaluates the objective function
PetscErrorCode evaluate_objective(Tao tao_context, Vec X, PetscReal *f, void *ptr) {
    
    Problem *problem = (Problem *)ptr;
    Function objfun = *(problem->objfun);
    int k = problem->k;
    
    return evaluate_function(X, f, &objfun, k);
}

// this function evaluates the gradient
PetscErrorCode evaluate_gradient(Tao tao_context, Vec X, Vec G, void *ptr) {
    Problem *problem = (Problem *)ptr;
    Function grafun = *(problem->grafun);
    int k = problem->k;

    return evaluate_function(X, G, &grafun, k);
}

// this function evaluates the hessian
PetscErrorCode evaluate_hessian(Tao tao_context, Vec X, Mat H, Mat Hpre, void *ptr) {
    Problem *problem = (Problem *)ptr;
    Function hesfun = *(problem->hesfun);
    int k = problem->k;
    
    return evaluate_function(X, H, &hesfun, k);
}

// this function evaluates the vector of inequalities
PetscErrorCode evaluate_inequalities(Tao tao_context, Vec X, Vec Ci, void *ptr) {
    Problem *problem = (Problem *)ptr;
    Function inequal = *(problem->inequal);
    int k = problem->k;
    
    return evaluate_function(X, Ci, &inequal, k);
}


// this function evaluates the vector of equalities
PetscErrorCode evaluate_equalities(Tao tao_context, Vec X, Vec Ci, void *ptr) {
    Problem *problem = (Problem *)ptr;
    Function equal = *(problem->equal);
    int k = problem->k;
    
    return evaluate_function(X, Ci, &equal, k);
}

// this function set the starting value
PetscErrorCode create_vec(Vec X, NumericVector y) {
    
    PetscReal *x;

    PetscFunctionBegin;
    catch_error(VecGetArray(X,&x));
    for(int iX=0; iX<y.size(); iX++) {
        x[iX] = y[iX];
    }
    catch_error(VecRestoreArray(X,&x));
    PetscFunctionReturn(0);
}

// this function set the starting value to zero
PetscErrorCode create_vec(Vec X, int k) {
  
  PetscReal *x;
  
  PetscFunctionBegin;
  catch_error(VecGetArray(X, &x));
  for(int iX = 0; iX < k; iX++) {
    x[iX] = 0;
  }
  catch_error(VecRestoreArray(X, &x));
  PetscFunctionReturn(0);
}

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

/*
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
*/
// this function set the starting value
PetscErrorCode create_vec(Vec X, NumericVector y) {
    
    PetscReal *x;

    PetscFunctionBegin;
    catch_error(VecGetArray(X, &x));
    for(int iX = 0; iX < y.size(); iX++) {
        x[iX] = y[iX];
    }
    catch_error(VecRestoreArray(X, &x));
    PetscFunctionReturn(0);
}

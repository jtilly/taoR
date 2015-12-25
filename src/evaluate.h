#ifndef evaluate_h
#define evaluate_h

#include "taoR.h"

// Evaluates a separable function of the form f1(X), f2(X)... Used for
// Pounders, which minimizes f1(X)^2 + f2(X)^2 + ...
//
// @param tao_context The tao context.
// @param X The parameter values to evaluate the function at.
// @param F The vector to write f1, f2,... into.
// @param ptr User-defined problem context.
// @return Error code.
PetscErrorCode evaluate_objective_separable(Tao tao_context, Vec X, Vec F, void *ptr);

// Evaluates an objective function of the form f(X).
//
// @param tao_context The tao context.
// @param X The parameter values to evaluate the function at.
// @param f The location to write the value of the objective function.
// @param ptr User-defined problem context.
// @return Error code.
PetscErrorCode evaluate_objective(Tao tao_context, Vec X, PetscReal *f, void *ptr);

// Evaluates the gradient function.
//
// @param tao_context The tao context.
// @param X The parameter values to evaluate the function at.
// @param G The vector to write the gradient to.
// @param ptr User-defined problem context.
// @return Error code.
PetscErrorCode evaluate_gradient(Tao tao_context, Vec X, Vec G, void *ptr);

// Evaluates the Hessian matrix.
//
// @param tao_context The tao context.
// @param X The parameter values to evaluate the function at.
// @param H The matrix to which to write the Hessian.
// @param ptr User-defined problem context.
// @return Error code.
PetscErrorCode evaluate_hessian(Tao tao_context, Vec X, Mat H, Mat Hpre, void *ptr);

// If the problem has constraints of the form g1(X) >=0, g2(X) >= 0, ...
// this evaluates the functions g1, g2, ...
//
// @param tao_context The tao context.
// @param X The parameter values to evaluate the function at.
// @param Ci The vector to which to write the values of g1, g2, ...
// @param ptr User-defined problem context.
// @return Error code.
PetscErrorCode evaluate_inequalities(Tao tao_context, Vec X, Vec Ci, void *ptr);

// If the problem has constraints of the form g1(X) = 0, g2(X) = 0, ...
// this evaluates the functions g1, g2, ...
//
// @param tao_context The tao context.
// @param X The parameter values to evaluate the function at.
// @param Ci The vector to which to write the values of g1, g2, ...
// @param ptr User-defined problem context.
// @return Error code.
PetscErrorCode evaluate_equalities(Tao tao_context, Vec X, Vec Ci, void *ptr);

// Initializes a new PETSc vector from an Rcpp NumericVector.
//
// @param X The vector to initialize..
// @param y The vector to copy values from.
// @return Error code.
PetscErrorCode create_vec(Vec X, NumericVector y);

#endif

#ifndef utils_h
#define utils_h

#include "taoR.h"

// Initializes petsc. Wraps messy command line parsing. Call me but once.
//
// @param options is the list to read. The column names are the flags, the values are
//        passed in as strings.
void initialize(List options);

// Returns the values of a Petsc vector written into an Rcpp vector.
//
// @param X A Petsc vector to get the values of.
// @param k The length of the vector.
// @returns An Rcpp vector with the entries of X.
NumericVector get_vec(Vec X, int k);

// Prints the progress of TAO optimization algorithms to the console.
//
// @param tao_context The TAO context.
// @param ptr The user problem context.
// @returns Error code checked with catch_error.
PetscErrorCode my_monitor(Tao tao_context, void *ptr);

// Re-directs all Petsc output from stdout to Rcpp:Rcout.
PetscErrorCode print_to_rcout(FILE *file, const char format[], va_list argp);

// Evaluates an Rcpp function of the form f(X).
//
// @param X Vector to evalute function on.
// @param y Stores the result of the function.
// @param f The function to evaluate.
// @returns Error code checked with catch_error.
PetscErrorCode evaluate_function(Vec X, PetscReal *y, Function *f, int k);

// Evaluates an Rcpp function which maps R^k to R^k.
//
// @param X k-vector to evalute function on.
// @param Y k-vector to store result.
// @param f The function to evaluate.
// @returns Error code checked with catch_error.
PetscErrorCode evaluate_function(Vec X, Vec Y, Function *f, int k);

// Evaluates an Rcpp function which maps R^k to R^n.
//
// @param X k-vector to evalute function on.
// @param Y n-vector to store result.
// @param f The function to evaluate.
// @returns Error code checked with catch_error.
PetscErrorCode evaluate_function(Vec X, Vec Y, Function *f, int k, int n);

// Evaluates an Rcpp function which maps R^k to R^(k^2).
//
// @param X k-vector to evalute function on.
// @param Y kxk matrix to store result.
// @param f The function to evaluate. Should return a kxk matrix.
// @returns Error code checked with catch_error.
PetscErrorCode evaluate_function(Vec X, Mat Y, Function *f, int k);

// Evaluates an Rcpp function which maps R^k to R^(n^2).
//
// @param X k-vector to evalute function on.
// @param Y nxn matrix to store result.
// @param f The function to evaluate. Should return an nxn matrix.
// @returns Error code checked with catch_error.
PetscErrorCode evaluate_function(Vec X, Mat Y, Function *f, int k, int n);

#endif

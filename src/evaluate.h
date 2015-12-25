#ifndef evaluate_h
#define evaluate_h

// Evaluates a separable objective function of the form f1(X), f2(X)...
//
// @param tao_context The tao context.
// @param X The parameter values to evaluate the function at.
// @param F The vector to write f1, f2,... into.
// @param ptr User-defined problem context.
// @return Error code.
PetscErrorCode evaluate_objective_separable(Tao tao_context, Vec X, Vec F, void *ptr);
PetscErrorCode evaluate_objective(Tao tao_context, Vec X, PetscReal *f, void *ptr) ;
PetscErrorCode evaluate_gradient(Tao tao_context, Vec X, Vec G, void *ptr);
PetscErrorCode evaluate_hessian(Tao tao_context, Vec X, Mat H, Mat Hpre, void *ptr);
PetscErrorCode evaluate_inequalities(Tao tao_context, Vec X, Vec Ci, void *ptr);
PetscErrorCode evaluate_equalities(Tao tao_context, Vec X, Vec Ci, void *ptr);

// Initializes a new PETSc vector from an Rcpp NumericVector.
//
// @param X The vector to initialize..
// @param y The vector to copy values from.
// @return Error code.
PetscErrorCode create_vec(Vec X, NumericVector y);

#endif

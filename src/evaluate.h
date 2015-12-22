#ifndef evaluate_h
#define evaluate_h

#include "utils.h"

PetscErrorCode form_starting_point(Vec, NumericVector);
PetscErrorCode evaluate_objective_separable(Tao, Vec, Vec, void *);
PetscErrorCode evaluate_objective(Tao, Vec, PetscReal*, void *);
PetscErrorCode evaluate_gradient(Tao, Vec, Vec, void *);
PetscErrorCode evaluate_hessian(Tao, Vec, Mat, Mat, void*);

#endif
#ifndef evaluate_h
#define evaluate_h

#include "common.h"
#include "utils.h"

PetscErrorCode createVec(Vec, NumericVector);
PetscErrorCode evaluate_objective_separable(Tao, Vec, Vec, void *);
PetscErrorCode evaluate_objective(Tao, Vec, PetscReal*, void *);
PetscErrorCode evaluate_gradient(Tao, Vec, Vec, void *);
PetscErrorCode evaluate_hessian(Tao, Vec, Mat, Mat, void*);

#endif
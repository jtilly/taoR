#ifndef evaluate_h
#define evaluate_h

PetscErrorCode create_vec(Vec, NumericVector);
PetscErrorCode evaluate_objective_separable(Tao, Vec, Vec, void *);
PetscErrorCode evaluate_objective(Tao, Vec, PetscReal*, void *);
PetscErrorCode evaluate_gradient(Tao, Vec, Vec, void *);
PetscErrorCode evaluate_hessian(Tao, Vec, Mat, Mat, void*);
PetscErrorCode evaluate_inequalities(Tao, Vec, Vec, void*);
PetscErrorCode evaluate_equalities(Tao, Vec, Vec, void*);

#endif
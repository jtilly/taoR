#ifndef taoR_h
#define taoR_h

#include <Rcpp.h>
#include <petsctao.h>

using namespace Rcpp;
using namespace std;

// problem structure
typedef struct {
  Function *objfun;
  Function *grafun;
  Function *hesfun;
  Function *inequal;
  Function *equal;
  int k;
  int n;
} Problem;

#define catch_error(operation) do { PetscErrorCode error_code = operation; CHKERRQ(error_code); } while (0)

#endif

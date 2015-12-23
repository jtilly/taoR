#ifndef common_h
#define common_h

#define catch_error(operation) do { PetscErrorCode error_code = operation; CHKERRQ(error_code); } while (0)

#include <Rcpp.h>
#include <petsctao.h>

using namespace Rcpp;
using namespace std;

// problem structure
typedef struct {
  Function *objfun;
  Function *grafun;
  Function *hesfun;
  int k;
  int n;
} Problem;

#endif
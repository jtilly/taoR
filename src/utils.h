#ifndef utils_h
#define utils_h

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

//' Initializes petsc. Wraps messy command line parsing. Call me but once.
//'
//' @param options is the list to read. The column names are the flags, the values are
//'        passed in as strings.
void initialize(List options);

NumericVector get_vec(Vec X, int k);

PetscErrorCode my_monitor(Tao tao_context, void *ptr);

PetscErrorCode print_to_rcout(FILE *file, const char format[], va_list argp);

#endif

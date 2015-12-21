#ifndef galeshapley_h
#define galeshapley_h

#include <Rcpp.h>
#include <petsctao.h>

//' Initializes petsc. Wraps messy command line parsing. Call me but once.
//'
//' @param options is the list to read. The column names are the flags, the values are
//'        passed in as strings.
void petscInitialize(Rcpp::List options);

#endif

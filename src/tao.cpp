#include <Rcpp.h>
using namespace Rcpp;

extern "C" {
    int chwirut1(int argc,char **argv);
}

// [[Rcpp::export]]
int test() {
    int argc;
    char **argv;
    chwirut1(argc, argv);
    return 0;
}

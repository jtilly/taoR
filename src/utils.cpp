#include "utils.h"

static std::vector<char*> parseCommandLineArguments(Rcpp::List options) {
    std::vector<char*> retVal(0);
    
    Rcpp::CharacterVector names = options.names();
    
    for (int i = 0; i < names.size(); ++i) {
        Rcpp::Rcout << names[i];
    }
}

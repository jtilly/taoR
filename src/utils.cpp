#include "utils.h"

void petscInitialize(Rcpp::List options) {
    
    // Values will be read into this vector
    std::vector<char*> args;
    
    // Get the list of the column names
    Rcpp::CharacterVector names = options.names();
    
    // Read List into vector of char arrays
    for (int i = 0; i < names.size(); ++i) {
        std::string flag = "-" + names[i];
        char *arg = new char[flag.size() + 1];
        std::copy(flag.begin(), flag.end(), arg);
        arg[flag.size()] = '\0';
        args.push_back(arg);
        
        std::string val = options[i];
        char *argVal = new char[val.size() + 1];
        std::copy(val.begin(), val.end(), argVal);
        argVal[val.size()] = '\0';
        args.push_back(argVal);
    }
    args.push_back(0);
    
    // Read vector into char array
    int argc = args.size() - 1;
    char** argv = &args[0];
    
    // Initialize petsc
    PetscInitialize(&argc, &argv, (char *)0, (char *)0);
    
    // Delete command line options
    for (size_t i = 0 ; i < args.size(); i++) {
        delete[] args[i];
    }
}
#include "utils.h"

void petscInitialize(Rcpp::List options) {
    
    // Values will be read into this vector
    std::vector<char*> args;
    int argc = 1;
    char** argv;
    
    if(options.size() > 0) {
        
        // Get the list of the column names
        Rcpp::CharacterVector names = options.names();
        
        // Read in the name
        char *name = new char[1];
        name[0] = '\0';
        args.push_back(name);
        
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
        
        // Read vector into char array
        argc = args.size();
        argv = &args[0u];
        
    } 
    
    // Check if already initialized
    PetscBool isInitialized;
    PetscInitialized(&isInitialized);
    
    if (isInitialized == PETSC_FALSE) {
        // Initialize PETSc
        PetscInitialize(&argc, &argv, (char *)0, (char *)0);
    } else {
        // Re-set the options
        PetscOptionsClear();
        PetscOptionsInsert(&argc, &argv, (char *)0);
    }
    
    // Delete command line options
    if(options.size() > 0) {
        for (size_t i = 0 ; i < args.size(); i++) {
            delete[] args[i];
        }
    }
}
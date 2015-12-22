#include "utils.h"

// [[Rcpp::export]]
int tao_init() {
    initialize(Rcpp::List::create());
    return 0;
}

// [[Rcpp::export]]
int tao_finalize() {
    PetscFinalize();
    return 0; 
}

void initialize(Rcpp::List options) {
    
    int argc = 1 + 2 * options.size();
    char** argv = new char*[argc];
    
    // Read in the name
    std::string name = "\0";
    argv[0] = new char[name.size() + 1];
    strcpy(argv[0], name.c_str());
    
    if(options.size() > 0) {
        
        // Get the list of the column names
        Rcpp::CharacterVector names = options.names();
        
        // internal counter
        int counter = 1;
        
        // Read List into vector of char arrays
        for (int i = 0; i < names.size(); ++i) {
            std::string flag = "-" + names[i];
            argv[counter] = new char[flag.size() + 1];
            strcpy(argv[counter++], flag.c_str());
            
            std::string val = options[i];
            argv[counter] = new char[val.size() + 1];
            strcpy(argv[counter++], val.c_str());
        }
        
    }
    
    // Check if already initialized
    PetscBool isInitialized;
    PetscInitialized(&isInitialized);
    
    if (isInitialized == PETSC_FALSE) {
        // Initialize PETSc
        PetscInitialize(&argc, &argv, (char *)0, (char *) 0);
    } else { 
        // Re-set the options
        PetscOptionsClear();
        PetscOptionsInsert(&argc, &argv, (char *) 0);
    }
    
    for(int i = 0; i < argc; ++i)
        delete[] argv[i];
    
    delete[] argv;

}
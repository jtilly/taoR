#include "utils.h"

//' Initialize TAO
//' 
//' This function is called automatically when the package is loaded.
// [[Rcpp::export]]
void tao_init() {
    initialize(Rcpp::List::create());
}

//' Finalize TAO
//' 
//' This function is called automatically when the package is unloaded.
// [[Rcpp::export]]
void tao_finalize() {
    PetscFinalize();
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

// this function transforms a vector of type Vec to a vector of type
// NumericVector
NumericVector get_vec(Vec X, int k) {
    PetscErrorCode error_code;
    PetscReal *x;
    VecGetArray(X, &x);
    NumericVector xVec(k);
    for (int i = 0; i < k; i++) {
        xVec[i] = x[i];
    }
    error_code = VecRestoreArray(X, &x); CHKERRQ(error_code);
    return xVec;
}

// this function reports the progress of the optimizer
PetscErrorCode my_monitor(Tao tao_context, void *ptr) {
    
    PetscReal fc, gnorm;
    PetscInt its;
    PetscViewer viewer = PETSC_VIEWER_STDOUT_SELF;
    PetscErrorCode error_code;
    
    PetscFunctionBegin;
    error_code = TaoGetSolutionStatus(tao_context, &its, &fc, &gnorm, 0, 0, 0);
    error_code = PetscViewerASCIIPrintf(viewer, "iter = %3D,", its); CHKERRQ(error_code);
    error_code = PetscViewerASCIIPrintf(viewer, " Function value %g,", (double) fc); CHKERRQ(error_code);
    if (gnorm > 1.e-6) {
        error_code = PetscViewerASCIIPrintf(viewer, " Residual: %g \n", (double) gnorm); CHKERRQ(error_code);
    } else if (gnorm > 1.e-11) {
        error_code = PetscViewerASCIIPrintf(viewer, " Residual: < 1.0e-6 \n"); CHKERRQ(error_code);
    } else {
        error_code = PetscViewerASCIIPrintf(viewer, " Residual: < 1.0e-11 \n"); CHKERRQ(error_code);
    }
    PetscFunctionReturn(0);
}

// Checks if output is going to stdout or stderr, if so, redirects to Rcout or Rcerr.
// Overrides PetscVFPrintf.
PetscErrorCode print_to_rcout(FILE *file, const char format[], va_list argp) {
    PetscErrorCode error_code;
    
    PetscFunctionBegin;
    if (file != stdout && file != stderr) {
        error_code = PetscVFPrintfDefault(file, format, argp); CHKERRQ(error_code);
    } else if (file == stdout) {
        char buff[1024];
        size_t length;
        error_code = PetscVSNPrintf(buff, 1024, format, &length, argp); CHKERRQ(error_code);
        Rcout << buff;
    } else if (file == stderr) {
        char buff[1024];
        size_t length;
        error_code = PetscVSNPrintf(buff, 1024, format, &length, argp); CHKERRQ(error_code);
        Rcerr << buff;
    }
    PetscFunctionReturn(0);
}

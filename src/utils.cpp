#include <taoR.h>
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
    string name = "\0";
    argv[0] = new char[name.size() + 1];
    strcpy(argv[0], name.c_str());
    
    if(options.size() > 0) {
        
        // Get the list of the column names
        CharacterVector names = options.names();
        
        // internal counter
        int counter = 1;
        
        // Read List into vector of char arrays
        for (int i = 0; i < names.size(); ++i) {
            string flag = "-" + names[i];
            argv[counter] = new char[flag.size() + 1];
            strcpy(argv[counter++], flag.c_str());
            
            string val = options[i];
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
    PetscReal *x;
    VecGetArray(X, &x);
    NumericVector xVec(k);
    for (int i = 0; i < k; i++) {
        xVec[i] = x[i];
    }
    catch_error(VecRestoreArray(X, &x));
    return xVec;
}

// this function reports the progress of the optimizer
PetscErrorCode my_monitor(Tao tao_context, void *ptr) {
    
    PetscReal fc, gnorm;
    PetscInt its;
    PetscViewer viewer = PETSC_VIEWER_STDOUT_SELF;
    
    PetscFunctionBegin;
    catch_error(TaoGetSolutionStatus(tao_context, &its, &fc, &gnorm, 0, 0, 0));
    catch_error(PetscViewerASCIIPrintf(viewer, "iter = %3D,", its));
    catch_error(PetscViewerASCIIPrintf(viewer, " Function value %g,", (double) fc));
    if (gnorm > 1.e-6) {
        catch_error(PetscViewerASCIIPrintf(viewer, " Residual: %g \n", (double) gnorm));
    } else if (gnorm > 1.e-11) {
        catch_error(PetscViewerASCIIPrintf(viewer, " Residual: < 1.0e-6 \n"));
    } else {
        catch_error(PetscViewerASCIIPrintf(viewer, " Residual: < 1.0e-11 \n"));
    }
    PetscFunctionReturn(0);
}

// Checks if output is going to stdout or stderr, if so, redirects to Rcout or Rcerr.
// Overrides PetscVFPrintf.
PetscErrorCode print_to_rcout(FILE *file, const char format[], va_list argp) {
        
    PetscFunctionBegin;
    if (file != stdout && file != stderr) {
        catch_error(PetscVFPrintfDefault(file, format, argp));
    } else if (file == stdout) {
        char buff[1024];
        size_t length;
        catch_error(PetscVSNPrintf(buff, 1024, format, &length, argp));
        Rcout << buff;
    } else if (file == stderr) {
        char buff[1024];
        size_t length;
        catch_error(PetscVSNPrintf(buff, 1024, format, &length, argp));
        Rcerr << buff;
    }
    PetscFunctionReturn(0);
}

PetscErrorCode evaluate_function(Vec X, PetscReal *y, Function *f, int k) {
    
    PetscReal *x;

    PetscFunctionBegin;
    
    // Read out arrays
    catch_error(VecGetArray(X, &x));

    // Write into Rcpp vector and evaluate
    NumericVector xVec = get_vec(X, k);
    NumericVector yVec = (*f)(xVec);
    
    y = &yVec[0];
    
    catch_error(VecRestoreArray(X, &x));
    PetscFunctionReturn(0);
    
}

PetscErrorCode evaluate_function(Vec X, Vec Y, Function *f, int k) {
    return evaluate_function(X, Y, f, k, k);
}

PetscErrorCode evaluate_function(Vec X, Vec Y, Function *f, int k, int n) {
  
    PetscReal *x;
    PetscReal *y;
  
    PetscFunctionBegin;
    
    // Read out arrays
    catch_error(VecGetArray(X, &x));
    catch_error(VecGetArray(Y, &y));
    
    // Write into Rcpp vector and evaluate
    NumericVector xVec = get_vec(X, k);
    NumericVector yVec = (*f)(xVec);
    
    // Write back into array
    
    for (int i = 0; i < n; ++i) {
      y[i] = yVec[i];
    }
    
    catch_error(VecRestoreArray(X, &x));
    catch_error(VecRestoreArray(Y, &y));
    PetscFunctionReturn(0);
  
}

PetscErrorCode evaluate_function(Vec X, Mat Y, Function *f, int k) {
    return evaluate_function(X, Y, f, k, k);
}

PetscErrorCode evaluate_function(Vec X, Mat Y, Function *f, int k, int n) {
  
    PetscReal *x;
    
    PetscFunctionBegin;
    
    // Read out arrays
    catch_error(VecGetArray(X, &x));

    // Write into Rcpp vector and evaluate
    NumericVector xVec = get_vec(X, k);
    NumericMatrix yMat = (*f)(xVec);
    
    // Assemble the matrix
    for (int row = 0; row < n; ++row) {
        for (int col = 0; col < n; ++col) {
            MatSetValues(Y, 1, &row, 1, &col, &(yMat(row, col)), INSERT_VALUES);
        }
    }
    
    catch_error(MatAssemblyBegin(Y, MAT_FINAL_ASSEMBLY));
    catch_error(MatAssemblyEnd(Y, MAT_FINAL_ASSEMBLY));
    catch_error(VecRestoreArray(X, &x));
    PetscFunctionReturn(0);
  
}
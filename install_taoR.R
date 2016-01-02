######################################################
# install_taoR.R
# 
# taoR requires a working installation of PETSc. There
# are three different ways to get this done:
# - download_binaries = TRUE will download PETSc binaries
#   that may or may not work on Linux and OSX depending
#   on the specification.
# - compile_binaries = TRUE will download the PETSc source
#   and then compile the binaries. This will take a while
#   but should always work
# - if you already have PETSc on your system, then set
#   PETSC_DIR and -- if needed -- PETSC_ARCH to point
#   to your installation of PETSc.
#
# Usage:
#   install_taoR(download_binaries = TRUE)
#   install_taoR(compile_binaries = TRUE)
#   install_taoR(PETSC_DIR = "~/petsc-3.6.3", PETSC_ARCH = "arch-linux2-c-opt")
#
install_taoR = function(download_binaries = TRUE, 
                        compile_binaries = FALSE,
                        PETSC_DIR = NULL,
                        PETSC_ARCH = NULL) {
    
    # get install_github function
    source("http://jtilly.io/install_github/install_github.R")
    
    if(! download_binaries && ! compile_binaries && is.null(PETSC_DIR)) {
        stop("Please choose a method to install taoR.")
    }
    
    if (!is.null(PETSC_DIR)) {
        if (!dir.exists(file.path(PETSC_DIR, PETSC_ARCH))) {
            stop("PETSC not found in ", file.path(PETSC_DIR, PETSC_ARCH))
        } 
        
        Sys.setenv("DOWNLOAD_PETSC_BINARIES"=0)
        Sys.setenv("PETSC_DIR" = path.expand(PETSC_DIR))
        if (!is.null(PETSC_ARCH)) {
            Sys.setenv("PETSC_ARCH" = PETSC_ARCH)
        }
    } 
    
    else if (compile_binaries) {
        Sys.setenv("PETSC_DIR"="")
        Sys.setenv("PETSC_ARCH"="")
        Sys.setenv("DOWNLOAD_PETSC_BINARIES"=0)
    } 
    
    else if (download_binaries) {
        Sys.setenv("DOWNLOAD_PETSC_BINARIES"=1)
        Sys.setenv("PETSC_DIR"="")
        Sys.setenv("PETSC_ARCH"="")
    }
    

    install_github("jtilly/taoR")

    # clean up
    Sys.unsetenv("DOWNLOAD_PETSC_BINARIES")
    Sys.unsetenv("PETSC_DIR")
    Sys.unsetenv("PETSC_ARCH")
}
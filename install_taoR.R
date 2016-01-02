#' install_github installs R packages directly from GitHub
#' 
#' install_github is a very lightweight function with no dependencies
#' that installs R packages from GitHub.
#' 
#' @param repo identifies the GitHub repository, i.e. username/repo
#' @param branch identifies the branch
#' @param dependencies is a vector that indicates which type of additional 
#'   packages are also to be installed. It can include any of the values in 
#'   \code{c("Depends", "Imports", "LinkingTo", "Suggests", "Enhances")}. A 
#'   value of \code{TRUE} corresponds to \code{c("Depends", "Imports", 
#'   "LinkingTo")}. If \code{dependencies} equals \code{FALSE} no further 
#'   packages will be installed.
#' @param method specifies which method \code{download.file} uses to download 
#'   the zip archive from GitHub. Can be set to any of the values in 
#'   \code{c("auto", "wget", "libcurl", "curl", "wininet")}. When set to 
#'   \code{auto}, this function will make an educated guess.
#' @examples
#' \dontrun{install_github("jtilly/matchingR")}
install_github = function(repo, branch = "master", dependencies = TRUE, method = "auto") {
    
    # unix system?
    unix = Sys.info()["sysname"]=="Darwin" || Sys.info()["sysname"]=="Linux"
    
    # select a method to download files over https
    if(method == "auto") {
        
        if(unix) {
            if(capabilities("libcurl")) {
                method = "libcurl"
            } else if (file.exists(Sys.which("wget"))) {
                method = "wget"
            } else if (file.exists(Sys.which("curl"))) {
                method = "curl"   
            } else {
                stop("Could not find method to download files over https.")
            }
        } else {
            # there's probably a way to fix this...
            if(getRversion() < "3.2") {
                stop("Need R 3.2 or higher to download files over https.")
            }
            method = "wininet"
        }
    }
    
    # assemble URL
    url = paste0("https://github.com/", repo, "/archive/", branch, ".zip")
    
    # create a temporary directory
    tmpdir = tempdir()
    pkg.zip = file.path(tmpdir, paste0(branch, ".zip"))
    
    # download the zip archive from GitHub
    message(paste("Downloading", url, "using", method))
    if (method == "curl") {
        # follow redirects with curl requires passing in the -L option
        download.file(url = url, destfile = pkg.zip, method = method, extra = "-L")
    } else {
        download.file(url = url, destfile = pkg.zip, method = method)
    }
    message(paste("Saved as", pkg.zip))
    if (!file.exists(pkg.zip)) {
        stop("Download failed.")
    }
    
    # extract the archive
    pkg.root = file.path(tmpdir, branch)
    unzip(pkg.zip, overwrite = TRUE, exdir = pkg.root)
    message(paste("Extracted package to", pkg.root))
    
    # assemble path to downloaded package
    pkg.dir = file.path(pkg.root, list.files(pkg.root)[1])
    
    # default dependencies
    if(dependencies == TRUE) {
        dependencies = c("Depends", "Imports", "LinkingTo")
    }
    
    # Get dependencies from the package's DESCRIPTION
    if(!is.null(dependencies)) {
        
        # http://stackoverflow.com/questions/30223957/elegantly-extract-r-package-dependencies-of-a-package-not-listed-on-cran
        dcf = read.dcf(file.path(pkg.dir, "DESCRIPTION"))
        deps = unique(gsub("\\s.*", "", trimws(unlist(strsplit(dcf[, intersect(dependencies, colnames(dcf))], ","), use.names = FALSE))))
        
        # install dependencies that aren't already installed
        deps = deps[!(deps %in% c("R", rownames(installed.packages())))]
        if(length(deps) > 0) {
            message(paste("Also installing:", paste(deps)))
            install.packages(deps)
        }
    }
    
    # check if the package has a configure / cleanup script
    if(unix) {
        if(file.exists(file.path(pkg.dir, "configure"))) {
            Sys.chmod(file.path(pkg.dir, "configure"), mode = "0755", use_umask = TRUE)
        }
        
        if(file.exists(file.path(pkg.dir, "cleanup"))) {
            Sys.chmod(file.path(pkg.dir, "cleanup"), mode = "0755", use_umask = TRUE)
        }
    }
    
    # install the package
    install.packages(pkg.dir, repos = NULL, type = "source")
    
    # clean up
    unlink(pkg.root, recursive = TRUE, force = TRUE)
    unlink(pkg.zip, force = TRUE)
}

#' install_taoR installs taoR from GitHub
#' 
#' taoR requires a working installation of PETSc. There
#' are three different ways to get this done:
#' @param download_binaries: Setting this to TRUE will download PETSc binaries 
#'   that may or may not work on Linux and OSX depending on the specification.
#' @param compile_binaries: Setting this to true will download the PETSc source 
#'   and then compile the binaries. This will take a while but should always
#'   work
#' @param PETSC_DIR: if you already have PETSc on your system, then set 
#'   PETSC_DIR and -- if needed -- PETSC_ARCH to point to your installation of
#'   PETSc.
#' @param PETSC_ARCH: if you already have PETSc on your system, then set 
#'   PETSC_DIR and -- if needed -- PETSC_ARCH to point to your installation of
#'   PETSc.
#' @examples
#'   install_taoR(download_binaries = TRUE)
#'   install_taoR(compile_binaries = TRUE)
#'   install_taoR(PETSC_DIR = "~/petsc-3.6.3", PETSC_ARCH = "arch-linux2-c-opt")
install_taoR = function(download_binaries = TRUE, 
                        compile_binaries = FALSE,
                        PETSC_DIR = NULL,
                        PETSC_ARCH = NULL) {
    
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

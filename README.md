# taoR [![Build Status](https://travis-ci.org/jtilly/taoR.svg?branch=master)](https://travis-ci.org/jtilly/taoR)

Toolkit for Advanced Optimization ([TAO](http://www.mcs.anl.gov/petsc/petsc-current/docs/tao_manual.pdf)) bindings for R. TAO includes a range of different optimizers. Among them is [Pounders](http://www.mcs.anl.gov/papers/P5120-0414.pdf), a local derivative-free optimizer for non-linear least squares problems. Pounders can be a useful tool for economists who estimate structural models using indirect inference. In situations, where one would commonly use Nelder-Mead, Pounders may be a better choice if the objective function takes on the form of a non-linear least squares problem. In contrast to Nelder-Mead, Pounders exploits the specific shape of the objective function and therefore requires far fewer iterations to arrive at a local optimum than Nelder-Mead. 


TAO is part of [PETSc](http://www.mcs.anl.gov/petsc/), which needs to be built and installed before installing this package.

This package is still at the "proof of concept" stage.

## Install

#### Install [PETSc](http://www.mcs.anl.gov/petsc/)
Our preferred way of installing PETSc is `python-pip`, which works on both Linux and Mac OS:
```{bash}
pip install petsc --allow-external petsc
```
Alternatively, you can install the PETSc libraries [by hand](http://www.mcs.anl.gov/petsc/documentation/installation.html) or use your system's package manager. On Debian-based systems, you can use `apt-get` and install from [sid] (https://packages.debian.org/sid/libpetsc3.6). On Mac OS, you can use [homebrew](http://brew.sh/): `brew install petsc`.

#### Install this package
You can install this package using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) from inside R:
```
# install.packages("devtools")
install_github("jtilly/taoR")
```
There are three environmental variables that ensure that R can find your installation of PETSc. You may have to set some of them by hand, *before* running `install_github()`.
* **PETSC_DIR**: This variable points to your PETSc installation. To change it, run `Sys.setenv("PETSC_DIR"="/where/is/petsc")`
* **PETSC_ARCH**: In case you compiled PETSc by hand, then this is the name of the directory where all the PETSc binaries are installed. To change it, run `Sys.setenv("PETSC_ARCH"="linux-debug-c")`
* **MPI_INCLUDE**: If you have an MPI library on your system, you may need to tell R where to look for the header file `mpi.h`. To change it, run `Sys.setenv("MPI_INCLUDE"="/where/is/mpi")`


## Example
We minimize the objective function `(x[1] - 3) ^ 2 + (x[2] + 1) ^ 2` with respect to `x`. 
```{r}
library("taoR")

# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2 
# with solution vector c(3, -1)

# use pounders
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao(functions = list(objfun = objfun), 
          start_values = c(1, 2), 
          method = "pounders", 
          options = list(tao_pounders_delta = "0.2"), 
          n = 2)
ret$x
# [1]  3 -1
ret$iterations
# [1] 5

    
# use Nelder-Mead
objfun = function(x) sum(c((x[1] - 3), (x[2] + 1))^2)
ret = tao(functions = list(objfun = objfun), 
          start_values = c(1, 2), 
          method = "nm", 
          options = list())
ret$x
# [1]  3.005468 -1.004479
ret$iterations
# [1] 20
```

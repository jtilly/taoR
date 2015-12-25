# Toolkit for Advanced Optimization (TAO) Bindings for R

[![Build Status](https://travis-ci.org/jtilly/taoR.svg?branch=master)](https://travis-ci.org/jtilly/taoR) [![Coverage Status](https://coveralls.io/repos/jtilly/taoR/badge.svg?branch=master&service=github)](https://coveralls.io/github/jtilly/taoR?branch=master) [![Project Status: Wip - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://img.shields.io/badge/status-WIP-yellow.svg)](https://img.shields.io/badge/status-WIP-yellow.svg)

taoR is an R package which lets you use the [TAO library](http://www.mcs.anl.gov/petsc/petsc-current/docs/tao_manual.pdf) from R. TAO is a library of optimization algorithms. Among them is [Pounders](http://www.mcs.anl.gov/papers/P5120-0414.pdf). Pounders can be a useful tool for economists who estimate structural models using indirect inference, because unlike commonly used algorithms such as Nelder-Mead, Pounders is tailored for minimizing a non-linear sum of squares objective function, and therefore may require fewer iterations to arrive at a local optimum than Nelder-Mead. For more details, see [here](http://arxiv.org/pdf/1406.5464.pdf) and [here](http://ftp.iza.org/dp8548.pdf).

TAO is part of [PETSc](http://www.mcs.anl.gov/petsc/). If you want to use taoR, you need to build and install PETSc first. The PETSc website contains detailed [installation instructions](http://www.mcs.anl.gov/petsc/documentation/installation.html), you can also see the section below for help getting started.

Please note that this package is currently a work in progress.

## Install

#### Install [PETSc](http://www.mcs.anl.gov/petsc/)
Our preferred way to install PETSc is to use the Python installer, `python-pip`, which works on both Linux and Mac OS:
```{bash}
pip install petsc --allow-external petsc
```
Alternatively, you can install the PETSc libraries [by hand](http://www.mcs.anl.gov/petsc/documentation/installation.html) or use your system's package manager. On Mac OS, you can use [homebrew](http://brew.sh/): `brew install petsc` (recommended). On Debian-based systems, you can use `apt-get` and install from [sid] (https://packages.debian.org/sid/libpetsc3.6) (not recommended). 

#### Install this package
You can install this package using [devtools](https://cran.r-project.org/web/packages/devtools/index.html) from inside R:
```{r}
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
ret = tao(par = c(1, 2),
                fn = objfun, 
                method = "pounders", 
                control = list(tao_pounders_delta = "0.2"), 
                n = 2)
ret$x
# [1]  3 -1
ret$iterations
# [1] 5

    
# use Nelder-Mead
objfun = function(x) sum(c((x[1] - 3), (x[2] + 1))^2)
ret = tao(par = c(1, 2),
                fn = objfun, 
                method = "nm", 
                control = list())
ret$x
# [1]  3.005468 -1.004479
ret$iterations
# [1] 20
```

## Available Optimization Algorithms
The parameter `method` can be set to one of the following optimizers.
* `nls`: Newton's method with line search for unconstrained minimization
* `ntr`: Newton's method with trust region for unconstrained minimization
* `ntl`: Newton's method with trust region, line search for unconstrained minimization
* `lmvm`: Limited memory variable metric method for unconstrained minimization
* `cg`: Nonlinear conjugate gradient method for unconstrained minimization
* `nm`: Nelder-Mead algorithm for derivate-free unconstrained minimization
* `tron`: Newton Trust Region method for bound constrained minimization
* `gpcg`: Newton Trust Region method for quadratic bound constrained minimization
* `blmvm`: Limited memory variable metric method for bound constrained minimization
* `pounders`: Derivate-free model-based algorithm for nonlinear least squares

# Toolkit for Advanced Optimization (TAO) Bindings for R

[![Build Status](https://travis-ci.org/jtilly/taoR.svg?branch=master)](https://travis-ci.org/jtilly/taoR) [![Coverage Status](https://coveralls.io/repos/jtilly/taoR/badge.svg?branch=master&service=github)](https://coveralls.io/github/jtilly/taoR?branch=master) [![Project Status: Wip - Initial development is in progress, but there has not yet been a stable, usable release suitable for the public.](https://img.shields.io/badge/status-WIP-yellow.svg)](https://img.shields.io/badge/status-WIP-yellow.svg)

taoR is an R package which lets you use the [TAO library](http://www.mcs.anl.gov/petsc/petsc-current/docs/tao_manual.pdf) from R. TAO is a library of optimization algorithms. Among them is [Pounders](http://www.mcs.anl.gov/papers/P5120-0414.pdf). Pounders can be a useful tool for economists who estimate structural models using indirect inference, because unlike commonly used algorithms such as Nelder-Mead, Pounders is tailored for minimizing a non-linear sum of squares objective function, and therefore may require fewer iterations to arrive at a local optimum than Nelder-Mead. For more details, see [here](http://arxiv.org/pdf/1406.5464.pdf) and [here](http://ftp.iza.org/dp8548.pdf).

Please note that this package is currently a work in progress.

## Install

TAO is part of [PETSc](http://www.mcs.anl.gov/petsc/). If you do not have PETSc installed on your system, then taoR will attempt to install PETSc for you. There are two different ways to install taoR.

#### 1. Build PETSc binaries as part of package installation
PETSc will be compiled using the same set of compilers that R uses. Building the PETSc binaries will take several minutes.
```{r}
source("http://jtilly.io/taoR/install_taoR.R")
install_taoR(compile_binaries = TRUE)
```

#### 2. Use existing PETSc installation
See [here](http://www.mcs.anl.gov/petsc/documentation/installation.html) for detailed instructions on how to install PETSc.
```{r}
source("http://jtilly.io/taoR/install_taoR.R")
install_taoR(PETSC_DIR = "/path/to/petsc", PETSC_ARCH = "...")
```
Note that taoR copies the PETSc binary into the R package directory. Installing the package this way may require manual adjustment of the flags in [`src/Makevars.in`](https://github.com/jtilly/taoR/blob/master/src/Makevars.in) before installation.

## Example
We minimize the objective function `(x[1] - 3) ^ 2 + (x[2] + 1) ^ 2` with respect to `x`. The syntax is similar to R's `optim` function.

```{r}
library("taoR")

# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2 
# with solution vector c(3, -1)

# use pounders
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao(par = c(1, 2),
                fn = objfun, 
                method = "pounders", 
                control = list("tao_pounders_gqt" = ""))
ret$x
# [1]  3 -1
ret$iterations
# [1] 7

    
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
* `pounders`: Derivative-free model-based algorithm for nonlinear least squares

# taoR [![Build Status](https://travis-ci.org/jtilly/taoR.svg?branch=master)](https://travis-ci.org/jtilly/taoR) ![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/taoR)
Toolkit for Advanced Optimization ([TAO](http://www.mcs.anl.gov/petsc/petsc-current/docs/tao_manual.pdf)) bindings for R. TAO includes a range of different optimizers. Among them is [Pounders](http://www.mcs.anl.gov/papers/P5120-0414.pdf), a local derivative-free optimizer for non-linear least squares problems. Pounders can be a useful tool for economists who estimate structural models using indirect inference. In situtations, where one would commonly use Nelder-Mead, Pounders may be a better choice if the objective function takes on the form of a non-linear least squares problem. In contrast to Nelder-Mead, Pounders exploits the specific shape of the objective function and therefore requires far fewer iterations to arrive at a local optimum than Nelder-Mead. 

TAO is part of [PETSc](http://www.mcs.anl.gov/petsc/), which needs to be built and installed before installing this package.

This package is still at the "proof of concept" stage.

## Install

#### Install [PETSc](http://www.mcs.anl.gov/petsc/)
```{bash}
# sudo apt-get install python-pip
pip install petsc --allow-external petsc
```
#### Install this package 
```{r}
Sys.setenv("PETSC_DIR" = "/where/is/petsc")
# install.packages("devtools")
devtools::install_github("jtilly/taoR")
```

## Example
We minimize the objective function `(x[1] - 3) ^ 2 + (x[2] + 1) ^ 2` with respect to `x`. 
```{r}
library("taoR")

# use Pounders
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao(objfun, startValues = c(1, 2), optimizer = "pounders", k = 2, n = 2)
ret$x
#> [1]  3 -1
ret$iterations
#> [1] 6

# use Nelder-Mead
objfun = function(x) sum(c((x[1] - 3), (x[2] + 1))^2)
ret = tao(objfun, startValues = c(1, 2), optimizer = "nm", k = 2)
ret$x
#> [1]  3.005468 -1.004479
ret$iterations
#> [1] 20
```

## Docker Container

[Docker](https://docs.docker.com/) is a lightweight virtual machine that can be run on almost any architecture. If you want to try `taoR` without installing the package and its dependencies on your machine, you can run R and `taoR` in a Docker container. The Docker container for this application is based on an Ubuntu image. For details, see the [Dockerfile](https://github.com/jtilly/taoR/blob/docker/Dockerfile) that creates the container. To run `taoR` in a Docker container, do the following:
```{bash}
docker run -ti jtilly/taor R
library("taoR")
demo("test")
```

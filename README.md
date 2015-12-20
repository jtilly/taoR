# taoR [![Build Status](https://travis-ci.org/jtilly/taoR.svg?branch=master)](https://travis-ci.org/jtilly/taoR)
Toolkit for Advanced Optimization (TAO) Bindings for R

This package is a proof of concept.

* Install [PETSc](http://www.mcs.anl.gov/petsc/): e.g. `pip install petsc --allow-external petsc`
* Set environmental variable: `PETSC_DIR=/where/is/petsc`
* Install this package: `devtools::install_github("jtilly/taoR")`

## Example

```{r}
# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2 
# with solution vector c(3, -1)

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

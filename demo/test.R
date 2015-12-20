# test.R
library("taoR")

# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2
# with solution vector [3, -1]

# use pounders
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao(objfun, c(1,2), "pounders", 2, 2)
ret$x

# use Nelder-Mead
objfun = function(x) sum(c((x[1] - 3), (x[2] + 1))^2)
ret = tao(objfun, c(1,2), "nm", 2)
ret$x
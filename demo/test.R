library("taoR")
# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2 
# with solution vector c(3, -1)

# use pounders
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao(par = c(1, 2), fn = objfun, method = "pounders", control = list(tao_pounders_delta="0.1"), n = 2)
ret$x
    
# use Nelder-Mead
objfun = function(x) sum(c((x[1] - 3), (x[2] + 1))^2)
ret = tao(par = c(1, 2), fn = objfun, method = "nm")
ret$x

# with gradient
objfun = function(x) sum(c((x[1] - 3)^2, (x[2] + 1))^2)
grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
     
ret = tao(par = c(1, 2), fn = objfun, gr = grafun, method = "lmvm")

ret$x

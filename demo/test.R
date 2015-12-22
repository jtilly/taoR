# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2 
# with solution vector c(3, -1)

# use pounders
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao(functions = list(objFun = objfun), 
          startValues = c(1, 2), 
          method = "pounders", 
          options = list(), 
          n = 2)
ret$x
    
# use Nelder-Mead
objfun = function(x) sum(c((x[1] - 3), (x[2] + 1))^2)
ret = tao(functions = list(objFun = objfun), 
          startValues = c(1, 2), 
          method = "nm", 
          options = list())
ret$x

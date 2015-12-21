# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2 
# with solution vector c(3, -1)

# use pounders
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao(functions = list(objFun = objfun), 
          startValues = c(1, 2), 
          method = "pounders", 
          options = list(tao_pounders_npmax = "1", tao_pounders_delta = "0.2"), 
          n = 2)
ret$x
    
# use Nelder-Mead
objfun = function(x) sum(c((x[1] - 3), (x[2] + 1))^2)
ret = tao(functions = list(objFun = objfun), 
          startValues = c(1, 2), 
          method = "nm", 
          options = list(),
          n = 1)
ret$x

# with gradient
objfun = function(x) sum(c((x[1] - 3)^2, (x[2] + 1))^2)
grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
     
ret = tao.optim(c(1, 2), 
                 objfun,
                 gr = grafun,
                 method = "nls",
                 n = 2)
 
ret
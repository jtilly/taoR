library("taoR")
# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2 
# with solution vector c(3, -1)

# use pounders
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao(functions = list(objfun = objfun), 
          start_values = c(1, 2), 
          method = "pounders", 
          options = list(tao_pounders_delta="0.1"), 
          n = 2)
ret$x
    
# use Nelder-Mead
objfun = function(x) sum(c((x[1] - 3), (x[2] + 1))^2)
ret = tao(functions = list(objfun = objfun), 
          start_values = c(1, 2), 
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
                 method = "lmvm",
                 n = 1)
ret$x
 
ret = tao(functions = list("objfun" = objfun, "grafun" = grafun),
          start_values = c(1, 2),
          method = "lmvm",
          options = list(),
          n = 1)
ret$x
# test.R

# the objective function is (x[1] - 3) ^ 2 + (x[2] + 1) ^2
# with solution vector [3, -1]
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = pounders(objfun, c(1,2), 2, 2)
ret$x
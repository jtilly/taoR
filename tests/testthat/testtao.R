# Gradient-free method
objfun = function(x) c((x[1] - 3), (x[2] + 1))
ret = tao.optim(c(1, 2), 
                objfun,
                method = "pounders",
                control = list(tao_pounders_delta="0.1"),
                n = 2)
expect_equal(sum(abs(c(3, -1) - ret$x)) < 0.01, TRUE)

# Gradient-based method
objfun = function(x) sum(c((x[1] - 3)^2, (x[2] + 1))^2)
grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
 
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                method = "lmvm")

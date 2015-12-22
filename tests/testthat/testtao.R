# TAOLMVM
objfun = function(x) (x[1] - 3)^2 + (x[2] + 1)^2
grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))

ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                method = "lmvm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAOCG
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                method = "cg")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAOOWLQN
# Returns [2.5, -0.5]
# ret = tao.optim(c(2, 0), 
#                 objfun,
#                 gr = grafun,
#                 method = "owlqn")
# expect_equal(sum(abs(c(3, -1) - ret$x)) < 0.01, TRUE)

# TAOCG
# Returns [2.2, -1.3]
#ret = tao.optim(c(1, 2), 
#                objfun,
#                gr = grafun,
#                method = "bmrm")
#expect_equal(sum(abs(c(3, -1) - ret$x)) < 0.01, TRUE)

# TAOBLMVM
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                method = "blmvm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAONM
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                method = "nm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# POUNDers
objfun = function(x) c(x[1] - 3, x[2] + 1)
ret = tao.optim(c(1, 2), 
                objfun,
                method = "pounders",
                control = list(),
                n = 2)
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAONM
# LCL Solver requires an initial state index set -- use TaoSetStateIS()
# ret = tao.optim(c(1, 2), 
#                objfun,
#                gr = grafun,
#                method = "lcl")
#expect_equal(sum(abs(c(3, -1) - ret$x)) < 0.01, TRUE)
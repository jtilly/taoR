objfun = function(x) (x[1] - 3)^2 + (x[2] + 1)^2
grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
hesfun = function(x) matrix(c(2, 0, 0, 2), nrow = 2, ncol = 2)

# TAOLMVM
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                method = "lmvm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAONLS
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "nls")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAONTR
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "ntr")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAONTL
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "ntl")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAOCG
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                method = "cg")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAOTRON
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "tron")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAOBLMVM
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                method = "blmvm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAOGPCG
ret = tao.optim(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "gpcg")
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
expect_equal(sum(objfun(ret$x)) < 0.01, TRUE)

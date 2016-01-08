library("taoR")
library("testthat")

objfun = function(x) (x[1] - 3)^2 + (x[2] + 1)^2
grafun = function(x) c(2*(x[1] - 3), 2*(x[2] + 1))
hesfun = function(x) matrix(c(2, 0, 0, 2), nrow = 2, ncol = 2)

# TAOLMVM
ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                method = "lmvm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                method = "lmvm")
expect_equal(objfun(ret$x) < 0.01, TRUE)


# TAONLS
ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "nls")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                hs = hesfun,
                method = "nls")
expect_equal(objfun(ret$x) < 0.01, TRUE)

expect_error(tao(c(1, 2), 
                objfun,
                method = "nls"))

# TAONTR
ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "ntr")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                hs = hesfun,
                method = "ntr")
expect_equal(objfun(ret$x) < 0.01, TRUE)

expect_error(tao(c(1, 2), 
                objfun,
                gr = grafun,
                method = "ntr"))

# TAONTL
ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "ntl")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                hs = hesfun,
                method = "ntl")
expect_equal(objfun(ret$x) < 0.01, TRUE)

expect_error(tao(c(1, 2), 
                objfun,
                gr = grafun,
                method = "ntl"))

# TAOCG
ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                method = "cg")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                method = "cg")
expect_equal(objfun(ret$x) < 0.01, TRUE)

# TAOTRON
ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "tron")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                hs = hesfun,
                method = "tron")
expect_equal(objfun(ret$x) < 0.01, TRUE)

expect_error(tao(c(1, 2), 
                objfun,
                gr = grafun,
                method = "tron"))

ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "tron", 
                lb = c(0, 0), ub = c(5, 5))
expect_equal(ret$x, c(3, 0))

# TAOBLMVM
ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                method = "blmvm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                method = "blmvm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                method = "blmvm",
                ub = c(5, 5), 
                lb = c(0, 0))
expect_equal(ret$x, c(3, 0))

# TAOGPCG
ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "gpcg")
expect_equal(objfun(ret$x) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                hs = hesfun,
                method = "gpcg")
expect_equal(objfun(ret$x) < 0.01, TRUE)

expect_error(tao(c(1, 2), 
                objfun,
                gr = grafun,
                method = "gpcg"))

ret = tao(c(1, 2), 
                objfun,
                gr = grafun,
                hs = hesfun,
                method = "gpcg", 
                lb = c(0, 0), ub = c(5, 5))
expect_equal(ret$x, c(3, 0))

# TAONM
ret = tao(c(1, 2), 
          objfun,
          method = "nm")
expect_equal(objfun(ret$x) < 0.01, TRUE)

expect_warning(tao(c(1, 2), 
          objfun,
          method = "nm",
          gr = grafun))

expect_warning(tao(c(1, 2), 
          objfun,
          method = "nm",
          gr = grafun,
          hs = hesfun))

expect_error(tao(c(1, 2), 
          objfun,
          method = "nm", n = 2))

# POUNDers
objfun = function(x) c(x[1] - 3, x[2] + 1)

ret = tao(c(1, 2), 
          objfun,
          method = "pounders",
          control = list())

expect_equal(sum(objfun(ret$x)) < 0.01, TRUE)

objfun = function(x) c(x[1] - 3, x[2] + 1)

ret = tao(c(1, 2), 
          objfun,
          method = "pounders",
          control = list(),
          n = 2)

expect_equal(sum(objfun(ret$x)) < 0.01, TRUE)

# POUNDers
objfun = function(x) c(x[1] - 3, x[2] + 1)

ret = tao(c(1, 2), 
                objfun,
                method = "pounders",
                control = list(),
                n = 2)

expect_equal(sum(objfun(ret$x)) < 0.01, TRUE)

ret = tao(c(1, 2), 
                objfun,
                method = "pounders",
                control = list(),
                n = 2,
                lb = c(0, 0), 
                ub = c(5, 5))

expect_equal(ret$x, c(3, 0))

expect_error(tao(c(1, 2), 
                 objfun,
                 method = "pounders",
                 control = list(),
                 n = 2,
                 lb = c(0, 0), 
                 ub = c(5)))

expect_error(tao(c(1, 2), 
                 objfun,
                 method = "pounders",
                 control = list(),
                 n = 2,
                 lb = c(0), 
                 ub = c(5, 5)))

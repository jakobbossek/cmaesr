context("CMA-ES run")

test_that("CMA-ES finds optimum of some BBOB functions", {
	max.iters = 50L
	lambda = 50L
	sigma = 1.5
  dims = c(2, 3, 4, 5)

  # accepted tolerance value for parameter and fitness values
  tol = 0.05
  fun.generators = c(makeSphereFunction, makeAckleyFunction, makeDoubleSumFunction)

  for (generator in fun.generators) {
    for (dim in dims) {
      fn = do.call(generator, list(dim))
      par.set = getParamSet(fn)
      opt = getGlobalOptimum(fn)
      lb = getLower(par.set)[1L]; ub = getUpper(par.set)[1L]

      res = runCMAES(
        fn,
        start.point = runif(dim, min = lb, max = ub),
        max.iter = max.iters,
        monitor = NULL,
        control = list(lambda = lambda * dim, sigma = sigma)
      )

      expect_true(abs(res$best.fitness - opt$value) < tol,
        info = sprintf("Desired fitness level not reached for dim = %i and function '%s'", dim, getName(fn)))
      expect_true(sum((res$best.param - opt$param)^2) < tol,
        info = sprintf("Desired parameter approximation not reached for dim = %i and function '%s'", dim, getName(fn)))
    } # dims
  } # fun.generators
})

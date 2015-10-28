context("CMA-ES run")

test_that("CMA-ES finds optimum of sphere function", {
	max.iters = 50L
	population.size = 50L
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
        population.size = population.size * dim,
        max.iter = max.iters,
        sigma = sigma,
        monitor = NULL
      )

      expect_true(abs(res$best.fitness - opt$value) < tol,
        info = sprintf("Desired fitness level not reached for dim = %i and function '%s'", dim, getName(fn)))
      expect_true(sum((res$best.param - opt$param)^2) < tol,
        info = sprintf("Desired parameter approximation not reached for dim = %i and function '%s'", dim, getName(fn)))
      #expect_equal(res$n.iters, max.iters)
      #expect_equal(res$convergence, 0L) # 0 for iter > max.iter
      #expect_equal(res$message, cma:::getTerminationMessage(0L))
    } # dims
  } # fun.generators
})


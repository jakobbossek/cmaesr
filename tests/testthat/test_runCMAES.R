context("CMA-ES run")

test_that("CMA-ES finds optimum of sphere function", {
	max.iters = 1000L
	population.size = 50L
	sigma = 0.5
  dims = c(2, 3, 4, 5)

  for (dim in dims) {
    objective.fun = makeSphereFunction(dim)
    par.set = getParamSet(objective.fun)
    lb = getLower(par.set)[1L]; ub = getUpper(par.set)[1L]
    res = runCMAES(
      objective.fun,
      start.point = runif(dim, min = lb, max = ub),
      population.size = population.size,
      max.iter = max.iters,
      sigma = sigma,
      monitor = NULL
    )
    print(res$best.fitness)

    expect_less_than(res$best.fitness, 0.01, info = sprintf("Desired fitness level not reached for dim = %i", dim))
    expect_true(all(res$best.param <= 0.01), info = sprintf("Desired parameter approximation not reached for dim = %i", dim))
    expect_equal(res$n.iters, max.iters)
    expect_equal(res$convergence, 0L) # 0 for iter > max.iter
    expect_equal(res$message, cma:::getTerminationMessage(0L))
  }
})


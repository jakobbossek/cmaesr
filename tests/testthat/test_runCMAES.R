context("CMA-ES run")

test_that("CMA-ES finds optimum of sphere function", {
	max.iters = 100L
	population.size = 50L
	sigma = 0.5
	objective.fun = makeSingleObjectiveFunction(
		name = "Sphere fun",
		fn = function(x) sum(x^2),
		par.set = makeParamSet(
			makeNumericParam("x1", lower = -10, upper = 10),
			makeNumericParam("x2", lower = -10, upper = 10)
		)
	)

	res = runCMAES(objective.fun,
		start.point = c(7, 7),
		population.size = population.size,
		max.iter = max.iters,
		sigma = sigma,
		monitor = NULL)

	expect_less_than(res$best.fitness, 0.01)
	expect_equal(res$n.iters, max.iters)
	expect_equal(res$convergence, 0L) # 0 for iter > max.iter
	expect_equal(res$message, cma:::getTerminationMessage(0L))
})


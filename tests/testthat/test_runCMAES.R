context("CMA-ES run")

test_that("CMA-ES finds optimum of sphere function", {
	objective.fun = makeSingleObjectiveFunction(
		name = "Sphere fun",
		fn = function(x) sum(x^2),
		par.set = makeParamSet(
			makeNumericParam("x1", lower = -10, upper = 10),
			makeNumericParam("x2", lower = -10, upper = 10)
		)
	)

	res = runCMAES(objective.fun, start.point = c(7, 7), population.size = 50L, max.iter = 100L, sigma = 0.5, monitor = NULL)
	expect_true(res$best.fitness < 0.01)
})
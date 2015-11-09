context("CMA-ES monitoring")

test_that("CMA-ES monitoring works well", {
	fn = makeSphereFunction(2L)
  res = cmaes(fn, monitor = makeSimpleMonitor())
  expect_true(is.numeric(res$best.fitness))
})

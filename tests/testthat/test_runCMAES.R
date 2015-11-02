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

test_that("CMA-ES works on Sphere with default parameters", {
  # accepted tolerance value for parameter and fitness values
  tol = 0.05

  for (dim in c(2, 3, 5)) {
    fn = makeSphereFunction(dim)
    res = runCMAES(
      fn,
      max.iter = 50L,
      monitor = NULL,
      control = list(sigma = 1, lambda = dim * 2 * 10)
    )
    expect_true(is.numeric(res$best.fitness))
    expect_true(all(is.numeric(res$best.param)))
    expect_true(res$best.fitness < tol, info = sprintf("For '%s' the desired fitness level was not reached.", getName(fn)))
  }
})

test_that("CMA-ES stops on invalid input", {
  control = list(sigma = 1)

  # multi-objective functions not supported
  fn = makeZDT1Function(2L)
  expect_error(runCMAES(fn, control = control))

  # noisy functions not allowed
  fn = makeSphereFunction(2L)
  attr(fn, "noisy") = TRUE
  expect_error(runCMAES(fn, control = control))

  # mixed functions
  fn = makeSingleObjectiveFunction(
    name = "Mixed",
    fn = function(x) {
      return(x$x^2 + as.numeric(x$y == "a"))
    },
    par.set = makeParamSet(
      makeNumericParam("x", lower = -10, upper = 10),
      makeDiscreteParam("y", values = c("a", "b"))
    )
  )
  expect_error(runCMAES(fn, control = control))
})

test_that("CMA-ES computes reasonanable results on noiseless 2D BBOB test set", {
  # check all functions
  fids = 1:24
  dims = 2
  lambda = 250L
  tol = 0.05

  for (fid in fids) {
    for (dim in dims) {
      # skip the hardest (very multimodal) functions
      if (fid %in% c(4, 5, 16, 23, 24)) {
        next
      }
      #lambda2 =  ifelse (fid %in% c(4, 5, 16, 23, 24), lambda * 4, lambda)
      fn = makeBBOBFunction(fid = fid, iid = 1L, dimension = dim)
      par.set = getParamSet(fn)
      opt = getGlobalOptimum(fn)
      lb = getLower(par.set)[1L]; ub = getUpper(par.set)[1L]
      control = list(sigma = (ub - lb) / 2, lambda = lambda)

      res = runCMAES(fn, control = control, monitor = NULL, max.iter = 150L)
      expect_true(is.numeric(res$best.fitness))
      expect_true(abs(res$best.fitness - opt$value) < tol,
        info = sprintf("Desired fitness level not reached for dim = %i and function '%s'", dim, getName(fn)))
      expect_true(sum((res$best.param - opt$param)^2) < tol,
        info = sprintf("Desired parameter approximation not reached for dim = %i and function '%s'", dim, getName(fn)))
    }
  }
})

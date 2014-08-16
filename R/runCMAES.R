#' Covariance-Matrix-Adaption
#'
#' Performs non-linear, non-convex optimization by means of the Covariance
#' Matrix Adaption Evolutionary Strategy by Hansen.
#'
#' @param objective.fun [\code{\link[otf]{otf_function}}]\cr
#'   Numerical objective function of type \code{\link[otf]{otf_function}}. The function
#'   must expect a list of numerical values and return a scaler numerical value.
#' @param start.point [\code{numeric}]\cr
#'   Initial solution vector.
#' @param population.size [\code{integer(1)}]\cr
#'   Population size.
#' @param sigma [\code{numeric(1)}]\cr
#'   Initial step-size, i. e., standard deviation in each coordinate direction.
#' @param max.iter [\code{integer(1)}]\cr
#'   Maximal number of sequential iterations.
#' @return [\code{CMAES_result}] Result object.
runCMAES = function(objective.fun, start.point, population.size = NULL, sigma, max.iter = 10L) {
	assertClass(objective.fun, "otf_function")

	# extract relevant data
	par.set = getParamSet(objective.fun)
	n = getNumberOfParameters(objective.fun)

	# sanity checks
	if (isNoisy(objective.fun)) {
		stopf("Noisy optimization is not supported at the moment.")
	}

	if (!isNumeric(par.set, include.int = FALSE)) {
		stopf("CMA-ES only works for objective functions with numeric parameters.")
	}

	if (isMultiobjective(objective.fun)) {
		stopf("CMA-ES can only handle single-objective functions.")
	}

	assertNumeric(start.point, len = n, any.missing = FALSE)
	assertNumber(sigma, lower = 0L, finite = TRUE)
	assertCount(max.iter, positive = TRUE)

	## SELECTION AN RECOMBINATION
	if (is.null(population.size)) {
		population.size = 4L + floor(3 * log(n))
	} else {		
		assertInt(population, lower = 4L)	
	}

	# offspring size
	lambda = population.size

	# parent size
	mu = floor(lambda / 2)

	# initialize recombination weights
	weights = log(mu + 0.5) - log(1:mu)

	# normalize weight vector
	weights = weights / sum(weights)

	# variance-effectiveness / variance effective selection mass of sum w_i x_i
	mu.eff = 1 / sum(weights^2)


	## STEP-SIZE CONTROL
	c.sigma = (me.eff + 2) / (n + mu.eff + 5)
	d.sigma = 1 + 2 * max(0, sqrt((mu.eff - 1) / (n + 1))) + c.sigma

	## COVARIANCE MATRIX ADAPTION
	c.c = (4 + mu.eff / n) / (n + 4 + 2 * mu.eff / n)
	c.1 = 2 / ((n + 1.3)^2 + mu.eff)
	alpha.mu = 2L
	c.mu = min(1 - c.1, alpha.mu * (mu.eff - 2 + 1/mu.eff) / ((n + 2)^2 + alpha.mu * mu.eff / 2))

	# path for covariance matrix C and stepsize sigma
	path.c = path.sigma = rep(0, n)
	D = diag(n)
	C = diag(n)
	C.inv = diag(n)

	eigen.eval = 0L
	count.eval = 0L

	# best individual
	best = Inf

	iter = 1L
	repeat {
		catf("Starting iteration %i.", iter)

		if (iter >= max.iter)
			break
		iter = iter + 1L
	}

}
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
	n.params = getNumberOfParameters(objective.fun)

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

	assertNumeric(start.point, len = n.params, any.missing = FALSE)
	assertNumber(sigma, lower = 0L, finite = TRUE)
	if (is.null(population.size)) {
		population.size = 4L + floor(3 * log(n.params))
	} else {		
		assertInt(population, lower = 4L)	
	}

	assertCount(max.iter, positive = TRUE)

	# offspring size
	lambda = population.size

	# parent size
	mu = floor(lambda / 2)

	# initialize recombination weights
	weights = log(mu + 0.5) - log(1:mu)

	# normalize weight vector
	weights = weights / sum(weights)

	# variance-effectiveness of sum w_i x_i
	mueff = sum(weights)^2 / sum(weights^2)

	iter = 1L
	repeat {
		catf("Starting iteration %i.", iter)

		if (iter >= max.iter)
			break
		iter = iter + 1L
	}

}
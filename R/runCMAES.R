#' Covariance-Matrix-Adaption
#'
#' Performs non-linear, non-convex optimization by means of the Covariance
#' Matrix Adaption Evolutionary Strategy by Hansen.
#'
#' @param objective.fun [\code{\link[otf]{otf_function}}]\cr
#'   Numerical objective function of type \code{\link[otf]{otf_function}}. The function
#'   must expect a list of numerical values and return a scaler numerical value.
#' @param start.point [\code{numeric}]\cr
#'   Initial solution vector. If \code{NULL}, one is generated randomly within the 
#'   box constraints offered by the paramter set of the objective function.
#' @param population.size [\code{integer(1)}]\cr
#'   Population size.
#' @param sigma [\code{numeric(1)}]\cr
#'   Initial step-size, i. e., standard deviation in each coordinate direction.
#' @param max.iter [\code{integer(1)}]\cr
#'   Maximal number of sequential iterations.
#' @param monitor [\code{cma_monitor}]\cr
#'   Monitoring object.
#' @return e[\code{CMAES_result}] Result object.
#' @export
runCMAES = function(objective.fun, start.point = NULL, population.size = NULL, sigma, max.iter = 10L, monitor = makeSimpleMonitor()) {
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

	if (!is.null(start.point)) {
		assertNumeric(start.point, len = n, any.missing = FALSE)	
	} else {
		if (!hasFiniteBoxConstraints(par.set)) {
			stopf("No start point provided. Cannot generate one, because parameter set cannot sample with Inf bounds!")
		}
		start.point = unlist(samplePoint(par.set))
	}
	assertNumber(sigma, lower = 0L, finite = TRUE)
	assertCount(max.iter, positive = TRUE)
	if (!is.null(monitor)) {		
		assertClass(monitor, "cma_monitor")	
	}

	## SELECTION AN RECOMBINATION
	if (is.null(population.size)) {
		population.size = 4L + floor(3 * log(n))
	} else {		
		assertInt(population.size, lower = 4L)	
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
	c.sigma = (mu.eff + 2) / (n + mu.eff + 5)
	d.sigma = 1 + 2 * max(0, sqrt((mu.eff - 1) / (n + 1))) + c.sigma

	## COVARIANCE MATRIX ADAPTION
	c.c = (4 + mu.eff / n) / (n + 4 + 2 * mu.eff / n)
	c.1 = 2 / ((n + 1.3)^2 + mu.eff)
	alpha.mu = 2L
	c.mu = min(1 - c.1, alpha.mu * (mu.eff - 2 + 1/mu.eff) / ((n + 2)^2 + alpha.mu * mu.eff / 2))

	# path for covariance matrix C and stepsize sigma
	path.c = path.sigma = rep(0, n)
	B = diag(n)
	D = diag(n)
	C = diag(n)
	C.invsqrt = diag(n)
	chi.n = sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))

	eigen.eval = 0L
	count.eval = 0L

	# best individual
	best.param = rep(NA, n)
	best.fitness = Inf

	## INITIALIZE BOOKKEEPING VARIABLES
	fitness = numeric(lambda)
	population = matrix(0, nrow = lambda, ncol = n)
	x.mean = start.point

	iter = 1L

	doMonitor = function(monitor, moment, ...) {
		if (!is.null(monitor))
			monitor[[moment]](...)
	}

	doMonitor(monitor, "before", iter, best.param, best.fitness)

	repeat {
		#FIXME: this is ugly as sin

		y = matrix(0, nrow = lambda, ncol = n)

		for (i in 1:lambda) {
			y[i, ] = B %*% D %*% t(rmvnorm(n = 1, sigma = diag(n)))
			population[i, ] = x.mean + sigma * y[i, ]
			fitness[i] = objective.fun(population[i, ])
			count.eval = count.eval + 1L
		}

		fitness.order = order(fitness)
		x.mean = colSums(population[fitness.order[1:mu], ] * weights)

		#FIXME: ugly
		y.w = colSums(y[fitness.order[1:mu], ] * weights)

		best.param = population[fitness.order[1], ]
		best.fitness = fitness[fitness.order[1]]


		# Update evolution path
		path.sigma = (1 - c.sigma) * path.sigma + sqrt(c.sigma * (2 - c.sigma) * mu.eff) * C.invsqrt %*% y.w
		h.sigma = sum(path.sigma) / sqrt(1 - (1 - c.sigma)^(2*(iter + 1))) < chi.n * (1.4 + 2 / (n + 1))
		path.sigma.norm = sqrt(sum(path.sigma^2))
		sigma = sigma * exp(c.sigma / d.sigma * ((path.sigma.norm / chi.n) - 1))

		# catf("h.sigma: %i", as.integer(h.sigma))
		# catf("path.sigma norm: %f", path.sigma.norm)
		# catf("SIGMA: %f", sigma)

		# Update covariance matrix
		path.c = (1 - c.c) * path.c + h.sigma * sqrt(c.c * (2 - c.c) * mu.eff) * y.w
		C = (1 - c.1 - c.mu) * C + c.1 * (path.c %*% t(path.c) + (1 - h.sigma) * c.sigma * (2 - c.c) * C)
		tmp = 0
		for (i in 1:mu) {
			yi = y[fitness.order[i], ]
			tmp = weights[i] * yi %*% t(yi)
		}
		C = C + c.mu * tmp

		doMonitor(monitor, "step", iter, best.param, best.fitness, population)
		
		#FIXME: write helpers getTerminationCode, getTerminationMessage
		if (iter >= max.iter)
			break
		iter = iter + 1L
	}
	doMonitor(monitor, "after", iter, best.param, best.fitness)
	makeS3Obj(
		best.param = best.param,
		best.fitness = best.fitness,
		convergence = 0L,
		classes = "cma_result"
	)
}
#' @title Covariance-Matrix-Adaption
#'
#' @description
#' Performs non-linear, non-convex optimization by means of the Covariance
#' Matrix Adaption Evolutionary Strategy by Hansen.
#'
#' @param objective.fun [\code{smoof_function}]\cr
#'   Numerical objective function of type \code{smoof_function}. The function
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
#'   Default is 10.
#' @param max.evals [\code{integer(1)}]\cr
#'   Maximal number of function evaluations.
#'   Default is \code{Inf}.
#' @param max.time [\code{integer(1)}]\cr
#'   Maximal time budget in seconds.
#'   Default is \code{Inf}.
#' @param monitor [\code{cma_monitor}]\cr
#'   Monitoring object.
#' @return [\code{CMAES_result}] Result object.
#' @export
#FIXME: add handling of noisy functions. See Hansen et al 2009, A Method for Handling Uncertainty in Evolutionary Optimization...
#FIXME: add restart options
runCMAES = function(objective.fun, start.point = NULL,
	population.size = NULL, sigma,
	max.iter = 10L, max.evals = Inf, max.time = Inf,
	monitor = makeSimpleMonitor()) {
	assertClass(objective.fun, "smoof_function")

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
		start.point = unlist(sampleValue(par.set))
	}
	assertNumber(sigma, lower = 0L, finite = TRUE)
	assertCount(max.iter, positive = TRUE)

	if (!is.infinite(max.evals)) {
		assertNumber(max.evals, lower = 0L, na.ok = FALSE)
	}

	if (!is.infinite(max.time)) {
		assertNumber(max.time, lower = 0L, na.ok = FALSE)
	}

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
  #FIXME: make sure, that weights are decreasing, i.e., w_1 >= w_2 >= ... >= w_n and sum w_i = 1
	weights = weights / sum(weights)

	# variance-effectiveness / variance effective selection mass of sum w_i x_i
	mu.eff = sum(weights)^2 / sum(weights^2) # chosen such that mu.eff ~ lambda/4


	## STEP-SIZE CONTROL
	c.sigma = (mu.eff + 2) / (n + mu.eff + 5)
	d.sigma = 1 + 2 * max(0, sqrt((mu.eff - 1) / (n + 1)) - 1) + c.sigma

	## COVARIANCE MATRIX ADAPTION
	c.c = (4 + mu.eff / n) / (n + 4 + 2 * mu.eff / n)
	c.1 = 2 / ((n + 1.3)^2 + mu.eff)
	alpha.mu = 2L
	c.mu = min(1 - c.1, alpha.mu * (mu.eff - 2 + 1/mu.eff) / ((n + 2)^2 + mu.eff))

	# path for covariance matrix C and stepsize sigma
	p.c = rep(0, n)
  p.sigma = rep(0, n)
	B = diag(n)
	D = diag(n)
	BD = B %*% D
  C = BD %*% t(BD) # C = B D^2 B^T = B B^T, since D equals I_n
  Cinvsqrt = B %*% diag(1 / sqrt(diag(D))) %*% t(B)

  # Precompute E||N(0,I)||
	chi.n = sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))

	# bookkeep best individual
	best.param = rep(NA, n)
	best.fitness = Inf

  # set initial distribution mean
	m = start.point

  # init some termination criteria stuff
	iter = 0L
  n.evals = 0L
	start.time = Sys.time()

	monitor$before()

	repeat {
    iter = iter + 1L

    # create new population of search points
		z = matrix(rnorm(n * lambda), ncol = lambda)
    y = BD %*% z # ~ N(0, C)
    x = m + sigma * y # ~ N(m, sigma^2 C)

    # compute fitness values (each idividual is a column of x)
    fitn = apply(x, 2L, function(x) objective.fun(x))

    # update evaluation
		n.evals = n.evals + lambda

    # order fitness values
    fitn.ordered.idx = order(fitn, decreasing = FALSE)
    fitn.ordered = fitn[fitn.ordered.idx]

    # lambda best individuals
    fitn.best = fitn.ordered[1:mu]

    # update best solution so far
    if (fitn.ordered[1L] < best.fitness) {
      best.fitness = fitn.ordered[1L]
      best.param = x[, fitn.ordered.idx[1L], drop = TRUE]
    }

    # update mean value / center of mass
    new.pop.idx = fitn.ordered.idx[1:mu]
    x.best = x[, new.pop.idx]
    m = drop(x.best %*% weights)

    #FIXME: do we really need y.w and z.w as variables?
    y.best = y[, new.pop.idx]
    y.w = drop(y.best %*% weights)
    z.best = z[, new.pop.idx]
    z.w = drop(z.best %*% weights)

		# Update evolution path with cumulative step-size adaption (CSA) / path length control
    # For an explanation of the last factor see appendix A in https://www.lri.fr/~hansen/cmatutorial.pdf
    p.sigma = (1 - c.sigma) * p.sigma + sqrt(c.sigma * (2 - c.sigma) * mu.eff) * (Cinvsqrt %*% y.w)
		h.sigma = as.integer(norm(p.sigma) / sqrt(1 - (1 - c.sigma)^(2 * (iter + 1))) < chi.n * (1.4 + 2 / (n + 1)))

		# Update covariance matrix
    p.c = (1 - c.c) * p.c + h.sigma * sqrt(c.c * (2 - c.c) * mu.eff) * y.w
    y = BD %*% z.best
    delta.h.sigma = as.numeric((1 - h.sigma) * c.c * (2 - c.c) <= 1)
		C = (1 - c.1 - c.mu) * C + c.1 * (p.c %*% t(p.c) + delta.h.sigma * C) + c.mu * y %*% diag(weights) %*% t(y)

    # Update step-size sigma
    sigma = sigma * exp(c.sigma / d.sigma * ((norm(p.sigma) / chi.n) - 1))

    # Finally do decomposition C = B D^2 B^T
    e = eigen(C, symmetric = TRUE)
    B = e$vectors
    D = diag(sqrt(e$values))
    BD = B %*% D
    C = BD %*% t(BD)
    Cinvsqrt = B %*% diag(1/diag(D)) %*% t(B) # update C^-1/2

		monitor$step()

    # now check if covariance matrix is positive definite
    if (any(e$values <= sqrt(.Machine$double.eps) * abs(max(e$values)))) {
      messagef("Covariance matrix is not numerically positive definite.")
      break
    }

		# check if we have to stop
		termination.code = getTerminationCode()
		if (termination.code > -1L) {
			break
		}
		iter = iter + 1L
	}
	monitor$after()
	makeS3Obj(
		par.set = par.set,
		best.param = best.param,
		best.fitness = best.fitness,
		convergence = termination.code,
		n.evals = n.evals,
		past.time = as.integer(difftime(Sys.time(), start.time, units = "secs")),
		n.iters = iter,
		message = getTerminationMessage(termination.code),
		classes = "cma_result"
	)
}

#' @export
print.cma_result = function(x, ...) {
	par.set = x$par.set
	par.names = getParamIds(par.set, repeated = TRUE, with.nr = TRUE)
	best.param = as.list(x$best.param)
	names(best.param) = par.names
	catf("Best parameter      : %s", paramValueToString(par.set, best.param))
	catf("Best fitness value  : %.6g", x$best.fitness)
	catf("Termination         : %s", x$message)
	catf("  #Iterations       : %i", x$n.iters)
	catf("  #Evaluations      : %i", x$n.evals)
	catf("  Time (in seconds) : %i", x$past.time)
}

getTerminationCode = function(envir = parent.frame()) {
	if (envir$iter > envir$max.iter)
		return(0L)
	if (envir$n.evals > envir$max.evals)
		return(1L)
	time.diff = difftime(Sys.time(), envir$start.time, units = "secs")
	if (time.diff > envir$max.time)
		return(2L)
	return(-1L)
}

getTerminationMessage = function(termination.code) {
	if (termination.code == 0L)
		return("Max iterations reached.")
	if (termination.code == 1L)
		return("Max functions evaluations exceeded.")
	if (termination.code == 2L)
		return("Time budget exceeded.")
	stopf("Unknown termination code '%i'", termination.code)
}

doMonitor = function(monitor, moment, ...) {
	if (!is.null(monitor))
		monitor[[moment]](...)
}

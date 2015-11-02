#' @title Covariance-Matrix-Adaption
#'
#' @description
#' Performs non-linear, non-convex optimization by means of the Covariance
#' Matrix Adaption Evolution Strategy (CMA-ES).
#'
#' @details
#' You may pass additional parameters to the CMA-ES via the \code{control} argument.
#' This must be a named list. The following elements will be considered by the
#' algorithm:
#' \describe{
#'   \item{lambda}{Number of offspring generaded in each generation.}
#'   \item{mu}{Number of individuals in each population. Defaults to \eqn{\lfloor \lambda / 2\rfloor}.}
#'   \item{weights}{Numeric vector of positive weights.}
#'   \item{sigma}{Initial step-size.}
#'   \item{do.restart}{Logical value indicating whether restarts should be triggered after certain
#'   stopping conditions fired. If \code{TRUE}, IPOP-CMA-ES is executed.}
#'   \item{restart.multiplier}{Factor which is used to increase the population size after restart.}
#'   \item{opt.value}{Scalar numeric known optimum.}
#'   \item{opt.param}{Numeric vector of target parameters. Algorithm stops if the euclidean distance to
#'     \code{opt.param} is lower then \code{tol.value}.}
#'   \code{tol.value}{Scalar numeric tolerance value. See \code{opt.value} and \code{opt.param} control
#'     parameters.}
#' }
#'
#' @references
#' [1] Auger and Hansen (2005). A Restart CMA Evolution Strategy With Increasing
#' Population Size. In IEEE Congress on Evolutionary Computation, CEC 2005, Proceedings,
#' pp. 1769-1776.
#' [2] N. Hansen (2006). The CMA Evolution Strategy: A Comparing Review. In J.A. Lozano,
#' P. Larranaga, I. Inza and E. Bengoetxea (Eds.). Towards a new evolutionary computation.
#' Advances in estimation of distribution algorithms. Springer, pp. 75-102.
#' [3] Hansen and Ostermeier (1996). Adapting arbitrary normal mutation distributions in evolution
#' strategies: The covariance matrix adaptation. In Proceedings of the 1996 IEEE
#' International Conference on Evolutionary Computation, pp. 312-317.
#'
#' @note
#' The restart variant is not yet implemented. Hence, setting \code{do.restart}
#' in \code{control} has no effect.
#'
#' @param objective.fun [\code{smoof_function}]\cr
#'   Numerical objective function of type \code{smoof_function}. The function
#'   must expect a list of numerical values and return a scaler numerical value.
#' @param start.point [\code{numeric}]\cr
#'   Initial solution vector. If \code{NULL}, one is generated randomly within the
#'   box constraints offered by the paramter set of the objective function.
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
#' @param control [\code{list}]\cr
#'   Futher paramters for the CMA-ES. See the details section for more in-depth
#'   information.
#' @return [\code{CMAES_result}] Result object.
#'
#' @examples
#' # generate objective function from smoof package
#' fn = makeRosenbrockFunction(dimensions = 2L)
#' res = runCMAES(fn, max.iter = 100L, monitor = NULL, control = list(sigma = 1.5, lambda = 40))
#' print(res)
#'
#' @export
runCMAES = function(objective.fun, start.point = NULL,
	max.iter = 10L, max.evals = Inf, max.time = Inf,
	monitor = makeSimpleMonitor(),
  control = list()) {
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

  # population and offspring size
  lambda = getCMAESParameter(control, "lambda", 4L + floor(3 * log(n)))
  assertInt(lambda, lower = 4)
  mu = getCMAESParameter(control, "mu", floor(lambda / 2))
  assertInt(mu)

  # initialize recombination weights
  weights = getCMAESParameter(control, "weights", log(mu + 0.5) - log(1:mu))
  if (any(weights < 0)) {
    stopf("All weights need to be positive, but there are %i negative ones.", sum(which(weights < 0)))
  }
  weights = weights / sum(weights)
  if (!(sum(weights) - 1.0) < .Machine$double.eps) {
    stopf("All 'weights' need to sum up to 1, but actually the sum is %f", sum(weights))
  }

  tol.value = getCMAESParameter(control, "tol.value", 0.05)
  opt.value = getCMAESParameter(control, "opt.value", NULL)
  opt.param = getCMAESParameter(control, "opt.param", NULL)

  #FIXME: default value should be derived from bounds
  sigma = getCMAESParameter(control, "sigma", 0.5)
  assertNumber(sigma, lower = 0L, finite = TRUE)

	# variance-effectiveness / variance effective selection mass of sum w_i x_i
	mu.eff = sum(weights)^2 / sum(weights^2) # chosen such that mu.eff ~ lambda/4

	# step-size control
	c.sigma = (mu.eff + 2) / (n + mu.eff + 5)
	d.sigma = 1 + 2 * max(0, sqrt((mu.eff - 1) / (n + 1)) - 1) + c.sigma

	# covariance matrix adaption parameters
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

  # tolerance values
  # (see Auger & Hansen, 2005 for details)
  tol.x = 10^(-12) * sigma # stop if all standard devitions are below this
  tol.fun = 10^(-12) # stop if range of best objective values in last generations is below this value
  tol.cond = 10^14 # stop if condition number of C exceeds this

	# bookkeep best individual
	best.param = rep(NA, n)
	best.fitness = Inf

  # set initial distribution mean
	m = start.point

  # init some termination criteria stuff
	iter = 0L
  n.evals = 0L
	start.time = Sys.time()

	callMonitor(monitor, "before")

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
    m.old = m
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

		callMonitor(monitor, "step")

    # CHECK STOPPING CONDITIONS
    # =========================
    stop.msg = NA

    # is covariance matrix not positive definite anymore?
    if (any(e$values <= sqrt(.Machine$double.eps) * abs(max(e$values)))) {
      stop.msg = "Covariance matrix is not numerically positive definite."
      break
    }

    # generations/iterations
    if (iter >= max.iter) {
      stop.msg = "Maximum number of iterations reached."
      break
    }

    # close enough to optimal parameters?
    if (!is.null(opt.param)) {
      if (sqrt(sum(best.param - opt.param)^2) < tol.value) {
        stop.msg = "Close enough to optimal parameter values."
        break
      }
    }

    # approximated optimal value?
    if (!is.null(opt.value)) {
      if (abs(best.fitness - opt.value) < tol.value) {
        stop.msg = "Found optimal value."
        break
      }
    }

    # running time
    if (difftime(Sys.time(), start.time, units = "secs") > max.time) {
      stop.msg = "Time limit reached."
      break
    }

    # function evalutions
    if (n.evals >= max.evals) {
      stop.msg = "Number of functions evaluations reached."
      break
    }

    # is the standard deviation below tolerance value in all coordinates?
    if (all(D < tol.x) && all((sigma * p.c) < tol.x)) {
      stop.msg = "All standard deviations below tolerance value."
      break
    }

    # Does addition of 0.1*sigma in a principal axis direction change m?
    # called 'noeffectaxis' (see Auger & Hansen, 2005)
    # [experimental]
    ii = (iter %% n) + 1L
    ui = e$vectors[, ii]
    lambdai = sqrt(e$values[ii])
    if (sum((m.old - (m.old + 0.1 * sigma * lambdai * ui))^2) < .Machine$double.eps) {
      stop.msg = "Adding fraction of standard deviation in principal axis direction does not change m."
      break
    }

    # Does addition of 0.2*sigma to each coordinate of m change m?
    if (sum((m.old - (m.old + 0.2 * sigma))^2) < .Machine$double.eps) {
      stop.msg = "Adding fraction of standard deviation to each coordinate of m does not change m."
      break
    }
	}

  callMonitor(monitor, "after")

	makeS3Obj(
		par.set = par.set,
		best.param = best.param,
		best.fitness = best.fitness,
		n.evals = n.evals,
		past.time = as.integer(difftime(Sys.time(), start.time, units = "secs")),
		n.iters = iter - 1L,
		message = stop.msg,
		classes = "cma_result"
	)
}

#' @export
print.cma_result = function(x, ...) {
	best.param = list(x$best.param)
	names(best.param) = getParamIds(x$par.set)
	catf("Best parameter      : %s", paramValueToString(x$par.set, best.param))
	catf("Best fitness value  : %.6g", x$best.fitness)
	catf("Termination         : %s", x$message)
	catf("  #Iterations       : %i", x$n.iters)
	catf("  #Evaluations      : %i", x$n.evals)
	catf("  Time              : %i (secs)", x$past.time)
}

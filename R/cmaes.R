#' @title Covariance-Matrix-Adaption
#'
#' @description
#' Performs non-linear, non-convex optimization by means of the Covariance
#' Matrix Adaption - Evolution Strategy (CMA-ES).
#'
#' @details
#' This a pure R implementation of the popular CMA-ES optimizer for numeric
#' black box optimization [2, 3]. It features a flexible system of stopping conditions
#' and enables restarts [1], which can be triggered by arbitrary stopping conditions.
#'
#' You may pass additional parameters to the CMA-ES via the \code{control} argument.
#' This argument must be a named list. The following control elements will be considered
#' by the CMA-ES implementation:
#' \describe{
#'   \item{lambda [\code{integer(1)}]}{Number of offspring generaded in each generation.}
#'   \item{mu [\code{integer(1)}]}{Number of individuals in each population. Defaults to \eqn{\lfloor \lambda / 2\rfloor}.}
#'   \item{weights [\code{numeric}]}{Numeric vector of positive weights.}
#'   \item{sigma [\code{numeric(1)}]}{Initial step-size.}
#'   \item{do.restart [\code{logical(1)}]}{Logical value indicating whether restarts should be triggered after certain
#'   stopping conditions are fired. If \code{TRUE} the CMA-ES is restarted with a
#'   bigger population size (IPOP-CMA-ES).}
#'   \item{restart.triggers [\code{character}]}{List of stopping condition codes / short names (see
#'   \code{\link{makeStoppingCondition}}). All stopping conditions which are placed in this vector do trigger a restart
#'   instead of leaving the main loop. Default is the empty character vector, i.e., restart is not triggered.}
#'   \item{max.restarts [\code{integer(1)}]}{Maximal number of restarts. Default is 0.}
#'   \item{restart.multiplier [\code{numeric(1)}]}{Factor which is used to increase the population size after restart.}
#'   \item{stop.ons [\code{list}]}{List of stopping conditions. The default is to stop after 10 iterations or after a
#'   kind of a stagnation (see \code{\link{getDefaultStoppingConditions}})}.
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
#' @keywords optimize
#'
#' @param objective.fun [\code{smoof_function}]\cr
#'   Numerical objective function of type \code{smoof_function}. The function
#'   must expect a list of numerical values and return a scaler numerical value.
#' @param start.point [\code{numeric}]\cr
#'   Initial solution vector. If \code{NULL}, one is generated randomly within the
#'   box constraints offered by the paramter set of the objective function.
#'   Default is \code{NULL}.
#' @param monitor [\code{cma_monitor}]\cr
#'   Monitoring object.
#'   Default is \code{\link{makeSimpleMonitor}}.
#' @param control [\code{list}]\cr
#'   Futher paramters for the CMA-ES. See the details section for more in-depth
#'   information. Stopping conditions are also defined here.
#'   By default only some stopping conditions are passed. See \code{\link{getDefaultStoppingConditions}}.
#' @return [\code{CMAES_result}] Result object.
#'
#' @examples
#' # generate objective function from smoof package
#' fn = makeRosenbrockFunction(dimensions = 2L)
#' res = cmaes(
#'   fn,
#'   monitor = NULL,
#'   control = list(
#'     sigma = 1.5, lambda = 40,
#'     stop.ons = c(list(stopOnMaxIters(100L)), getDefaultStoppingConditions())
#'   )
#' )
#' print(res)
#'
#' @export
cmaes = function(
  objective.fun,
  start.point = NULL,
	monitor = makeSimpleMonitor(),
  control = list(
    stop.ons = c(
      getDefaultStoppingConditions()
    )
  )) {
	assertClass(objective.fun, "smoof_function")

	# extract relevant data
	par.set = getParamSet(objective.fun)
  lb = getLower(par.set); ub = getUpper(par.set)
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

	if (!is.null(monitor)) {
		assertClass(monitor, "cma_monitor")
	}

  # get stopping conditions
  stop.ons = getCMAESParameter(control, "stop.ons", NULL)
  if (is.null(stop.ons)) {
    stopf("There must be at least one stopping condition!")
  }
  assertList(stop.ons, min.len = 1L, types = "cma_stopping_condition")

  # restart mechanism (IPOP-CMA-ES)
  restart.triggers = getCMAESParameter(control, "restart.triggers", character(0L))
  stop.ons.names = sapply(stop.ons, function(stop.on) stop.on$code)
  if (!isSubset(restart.triggers, stop.ons.names)) {
    stopf("Only codes / short names of active stopping conditions allowed as restart trigger, but '%s' are no stopping conditions.", collapse(setdiff(restart.triggers, stop.ons.names), sep = ", "))
  }
  restart.multiplier = getCMAESParameter(control, "restart.multiplier", 2)
  assertNumber(restart.multiplier, lower = 1, na.ok = FALSE, finite = TRUE)
  max.restarts = getCMAESParameter(control, "max.restarts", 0L)
  assertInt(max.restarts)

  #FIXME: default value should be derived from bounds
  sigma = getCMAESParameter(control, "sigma", 0.5)
  assertNumber(sigma, lower = 0L, finite = TRUE)

	# path for covariance matrix C and stepsize sigma
	p.c = rep(0, n)
  p.sigma = rep(0, n)

  # Precompute E||N(0,I)||
	chi.n = sqrt(n) * (1 - 1 / (4 * n) + 1 / (21 * n^2))

	# bookkeep best individual
	best.param = rep(NA, n)
	best.fitness = Inf

  # set initial distribution mean
	m = start.point

  # logs
  population.trace = list()

  # init some termination criteria stuff
	iter = 0L
  n.evals = 0L
	start.time = Sys.time()

	callMonitor(monitor, "before")

  # somehow dirty trick to "really quit" if stopping condition is met and
  # now more restart should be triggered.
  do.terminate = FALSE

  for (run in 0:max.restarts) {
    # population and offspring size
    if (run == 0) {
      lambda = getCMAESParameter(control, "lambda", 4L + floor(3 * log(n)))
      assertInt(lambda, lower = 4)
      mu = getCMAESParameter(control, "mu", floor(lambda / 2))
      assertInt(mu)
    } else {
      lambda = getCMAESParameter(control, "lambda", 4L + floor(3 * log(n)))
      # increase population size (IPOP-CMA-ES)
      lambda = ceiling(restart.multiplier^run * lambda)
      mu = floor(lambda / 2)
    }

    # initialize recombination weights
    weights = getCMAESParameter(control, "weights", log(mu + 0.5) - log(1:mu))
    if (any(weights < 0)) {
      stopf("All weights need to be positive, but there are %i negative ones.", sum(which(weights < 0)))
    }
    weights = weights / sum(weights)
    if (!(sum(weights) - 1.0) < .Machine$double.eps) {
      stopf("All 'weights' need to sum up to 1, but actually the sum is %f", sum(weights))
    }

    # variance-effectiveness / variance effective selection mass of sum w_i x_i
    mu.eff = sum(weights)^2 / sum(weights^2) # chosen such that mu.eff ~ lambda/4

    # step-size control
    c.sigma = (mu.eff + 2) / (n + mu.eff + 5)
    damps = 1 + 2 * max(0, sqrt((mu.eff - 1) / (n + 1)) - 1) + c.sigma

    # covariance matrix adaption parameters
    c.c = (4 + mu.eff / n) / (n + 4 + 2 * mu.eff / n)
    c.1 = 2 / ((n + 1.3)^2 + mu.eff)
    alpha.mu = 2L
    c.mu = min(1 - c.1, alpha.mu * (mu.eff - 2 + 1/mu.eff) / ((n + 2)^2 + mu.eff))

    # covariance matrix
    sigma = getCMAESParameter(control, "sigma", 0.5)
    B = diag(n)
    D = diag(n)
    BD = B %*% D
    C = BD %*% t(BD) # C = B D^2 B^T = B B^T, since D equals I_n
    Cinvsqrt = B %*% diag(1 / sqrt(diag(D))) %*% t(B)

    # no restart trigger fired until now
    restarting = FALSE

    # break inner loop if terminating stopping condition active or
    # restart triggered
  	while (!restarting) {
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

      # log population
      population.trace[[iter]] = z.best

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
      sigma = sigma * exp(c.sigma / damps * ((norm(p.sigma) / chi.n) - 1))

      # Finally do decomposition C = B D^2 B^T
      e = eigen(C, symmetric = TRUE)
      B = e$vectors
      D = diag(sqrt(e$values))
      BD = B %*% D
      C = BD %*% t(BD)
      Cinvsqrt = B %*% diag(1 / diag(D)) %*% t(B) # update C^-1/2

      callMonitor(monitor, "step")

      # escape flat fitness values
      if (fitn.ordered[1L] == fitn.ordered[ceiling(0.7 * lambda)]) {
        sigma = sigma * exp(0.2 + c.sigma / damps)
        warningf("Flat fitness values; increasing mutation step-size. Consider reformulating the objective!")
      }

      # CHECK STOPPING CONDITIONS
      # =========================
      stop.obj = checkStoppingConditions(stop.ons)

      n.stop.codes = length(stop.obj$codes)
      if (max.restarts > 0L && any(stop.obj$codes %in% restart.triggers)) {
        messagef("Restart tigger fired! Restarting!!!")
        n.stop.codes = sum(!(stop.obj$codes %in% restart.triggers))
        restarting = TRUE
      }

      # check if CMA-ES should really quit, i.e., is there a stopping condition,
      # that is active and does not trigger a restart?
      if (!restarting && (n.stop.codes > 0L)) {
        do.terminate = TRUE
        break
      }
  	}

    # really quit without more restarts
    if (do.terminate) {
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
    n.restarts = run,
    population.trace = population.trace,
		message = stop.obj$stop.msgs,
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

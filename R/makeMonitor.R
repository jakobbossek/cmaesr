#' Factory method for monitor objects.
#'
#' @param before [\code{function}]\cr
#'   Function called one time after initialization of the EA.
#' @param step [\code{function}]\cr
#'   Function applied after each iteration of the algorithm.
#' @param after [\code{function}]\cr
#'   Function applied after the EA terminated.
#' @param ... [\code{any}]\cr
#'   Not used.
#' @return [\code{cma_monitor}]
#'   Monitor object.
#' @export
makeMonitor = function(before = NULL, step = NULL, after = NULL, ...) {
  if (!is.null(before)) assertFunction(before)
  if (!is.null(step)) assertFunction(step)
  if (!is.null(after)) assertFunction(after)
  dummy = function(...) {}
  structure(
    list(
      before = coalesce(before, dummy),
      step = coalesce(step, dummy),
      after = coalesce(after, dummy)
    ),
    class = "cma_monitor")
}


#' Generator for simple monitor.
#'
#' The simple monitor prints the iteration, current best parameter values and best fitness
#' to the standard output.
#'
#' @return [\code{cma_monitor}]
#' @export
makeSimpleMonitor = function() {
	makeMonitor(
		before = function(envir = parent.frame()) {
      catf("Starting optimization.")
    },
		step = function(envir = parent.frame()) {
      best.param = as.list(envir$best.param)
      names(best.param) = getParamIds(envir$par.set)
      pars = paramValueToString(envir$par.set, x = best.param, num.format = "%.4f")
			catf("Iteration %i: %s, y = %.4f", envir$iter, pars, envir$best.fitness)
		},
		after = function(envir = parent.frame()) {
      catf("Optimization terminated.")
    }
	)
}

#' Generator for visualizing monitor.
#'
#' This generator visualizes the optimization process for two-dimensional functions
#' by means of \pkg{ggplot2}.
#'
#' @return [\code{cma_monitor}]
#' @export
makeVisualizingMonitor = function() {
  makeMonitor(
    before = function(envir = parent.frame()) {},
    step = function(envir = parent.frame()) {
      # get the population and mean/center
      x = envir$x
      m = envir$m

      df = as.data.frame(t(cbind(x, m)))
      df$Type = "Population"
      df[nrow(df), "Type"] = "Mean"
      colnames(df) = c("x1", "x2", "Type")
      df$Type = as.factor(df$Type)

      # use smoof's autoplot function to generate the contour plot
      pl = autoplot(envir$objective.fun)
      # ... and decorate with the points
      pl = pl + geom_point(data = df, aes(x = x1, y = x2, colour = Type))
      pl = pl + theme(legend.position = "bottom")
      print(pl)
      pause()
    },
    after = function(envir = parent.frame()) {}
  )
}


callMonitor = function(monitor, step, envir = parent.frame()) {
  if (!is.null(monitor)) {
    monitor[[step]](envir = envir)
  }
}

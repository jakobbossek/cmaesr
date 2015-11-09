#' @title Factory method for monitor objects.
#'
#' @description Monitors can be pluged in the main \code{\link{cmaes}} function.
#' They have full access to the environment of the optimization routine and can
#' be used to write/log/visualize relevant data in each iteration.
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
#' @seealso \code{\link{makeSimpleMonitor}}, \code{\link{makeVisualizingMonitor}}
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

#' @title Generator for simple monitor.
#'
#' @description The simple monitor prints the iteration, current best parameter values and best fitness
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
      best.param = list(envir$best.param)
      names(best.param) = getParamIds(envir$par.set)
      pars = paramValueToString(envir$par.set, x = best.param, num.format = "%.4f")
			catf("Iteration %i: %s, y = %.4f", envir$iter, pars, envir$best.fitness)
		},
		after = function(envir = parent.frame()) {
      catf("Optimization terminated.")
    }
	)
}

#' @title Generator for visualizing monitor.
#'
#' @description This generator visualizes the optimization process for two-dimensional functions
#' by means of \pkg{ggplot2}.
#'
#' @details The plot contains points representing the current population, the center
#' of mass or mean value of the population respectively. Optionally an ellipsis
#' represneting the normal distribution of the points can be depicted.
#'
#' @param show.last [\code{logical(1)}]\cr
#'   Should the last population be visualized as well?
#'   Default is \code{FALSE}.
#' @param show.distribution [\code{logical(1)}]\cr
#'   Should an ellipsis of the normal distribution be plotted?
#'   Default is \code{TRUE}.
#' @return [\code{cma_monitor}]
#' @export
makeVisualizingMonitor = function(show.last = FALSE, show.distribution = TRUE) {
  assertFlag(show.last, na.ok = FALSE)
  assertFlag(show.distribution, na.ok = FALSE)

  # store last population here
  last.x = NULL
  force(last.x)
  force(show.last)
  force(show.distribution)

  makeMonitor(
    before = function(envir = parent.frame()) {},
    step = function(envir = parent.frame()) {
      # get the population and mean/center
      x = envir$x
      m = envir$m.old

      # visualization only applicable for the 2D case
      if (length(m) != 2L) {
        invisible(NULL)
      }

      #FIXME: the following lines are ugly as sin, but refactor later.
      df = as.data.frame(t(cbind(x, m)))
      df$Type = "Current population"
      df[nrow(df), "Type"] = "Mean"
      colnames(df) = c("x1", "x2", "Type")

      # if last population is available, append
      if (!is.null(last.x) && show.last) {
        df2 = as.data.frame(t(last.x))
        df2$Type = "Last population"
        colnames(df2) = c("x1", "x2", "Type")
        df = rbind(df, df2)
      }

      # type needs to be factor in order to use ggplot
      df$Type = as.factor(df$Type)
      rownames(df) = NULL

      # use smoof's autoplot function to generate the contour plot
      pl = autoplot(envir$objective.fun)
      # ... and decorate with the points
      pl = pl + geom_point(data = df, aes_string(x = "x1", y = "x2", colour = "Type"))
      pl = pl + theme(legend.position = "bottom")

      # show ellipsis of normal distribution
      if (show.distribution) {
        pop.idx = which(grepl("population", as.character(df$Type)))
        pl = pl + stat_ellipse(data = df[pop.idx, , drop = FALSE], aes_string(colour = "Type"),
          linetype = "dashed", type = "norm")
      }

      # update last population
      last.x <<- x
      print(pl)
      pause()
    },
    after = function(envir = parent.frame()) {}
  )
}

#' @title Helper to call certain step function of a monitor.
#'
#' @description This funtions serves to call a specific monitor step.
#'
#' @param monitor [\code{CMAES_monitor}]\cr
#'   Monitor.
#' @param step [\code{character(1)}]\cr
#'   One of before, step, after.
#' @param envir [\code{environment}]\cr
#'   The environment to pass.
callMonitor = function(monitor, step, envir = parent.frame()) {
  if (!is.null(monitor)) {
    monitor[[step]](envir = envir)
  }
}

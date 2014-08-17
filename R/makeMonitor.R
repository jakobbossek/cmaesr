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
#FIXME: this is a copy of ecr::makeMonitor. Since it does not use any specific
#       interfaces it could be exported to a monitoring package?
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
		before = function(...) catf("Starting optimization."),
		step = function(iter, best.param, best.fitness) {
			#FIXME: best.param may have more than 2 components
			catf("Iteration %i: x1 = %f, x2 = %f, y = %f", iter, best.param[1], best.param[2], best.fitness)
		},
		after = function(...) catf("Optimization terminated.")
	)
}

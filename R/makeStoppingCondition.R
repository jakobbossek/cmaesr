#' @title Generate a stopping condition object.
#'
#' @param name [\code{character(1)}]\cr
#'   Name of the stopping condition.
#' @param message [\code{character(1)}]\cr
#'   Message returned if the stopping conditions is active.
#' @param stop.fun [\code{function}]\cr
#'   Function which expects an environment \code{envir} as its only argument and
#'   returns a single logical value.
#' @param code [\code{character(1)}]\cr
#'   Internal code, i.e., short name used to potentially trigger restarts.
#'   Default is \code{name}.
#' @param control [\code{list}]\cr
#'   Control params.
#' @return [\code{cma_stopping_condition}] Stopping condition object.
#' @export
makeStoppingCondition = function(name, message, stop.fun, code = name, control = list()) {
  assertString(name, na.ok = FALSE)
  assertString(message, na.ok = FALSE)
  assertFunction(stop.fun, args = "envir")
  assertString(code, na.ok = FALSE)
  assertList(control)
  makeS3Obj(
    name = name,
    message = message,
    stop.fun = stop.fun,
    code = code,
    control = control,
    classes = "cma_stopping_condition"
  )
}

shouldStop = function(x, envir) {
  UseMethod("shouldStop")
}

shouldStop.cma_stopping_condition = function(x, envir) {
  return(x$stop.fun(envir))
}

checkStoppingConditions = function(stop.ons, envir = parent.frame()) {
  assertList(stop.ons, min.len = 1L, types = "cma_stopping_condition")
  stop.msgs = character(0L)
  codes = character(0L)
  for (stop.on in stop.ons) {
    if (shouldStop(stop.on, envir = envir)) {
      stop.msgs = c(stop.msgs, stop.on$message)
      codes = c(codes, stop.on$code)
    }
  }
  return(list(stop.msgs = stop.msgs, codes = codes))
}

#' @title Return list of default stopping conditions.
#'
#' @return [\code{list}]
#' @export
getDefaultStoppingConditions = function() {
  return(
    list(
      stopOnIndefCovMat(),
      stopOnTolX(),
      stopOnNoEffectAxis(),
      stopOnNoEffectCoord(),
      stopOnCondCov()
    )
  )
}

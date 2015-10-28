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

# cmaesr: Covariance Matrix Adaption - Evolution Strategy in R

[![CRAN Status Badge](http://www.r-pkg.org/badges/version/cmaesr)](http://cran.r-project.org/web/packages/cmaesr)
[![Build Status](https://travis-ci.org/jakobbossek/cmaesr.svg?branch=master)](https://travis-ci.org/jakobbossek/cmaesr)
[![Build status](https://ci.appveyor.com/api/projects/status/eu0nns2dsgocwntw/branch/master?svg=true)](https://ci.appveyor.com/project/jakobbossek/cmaesr/branch/master)
[![Coverage Status](https://coveralls.io/repos/jakobbossek/cmaesr/badge.svg)](https://coveralls.io/r/jakobbossek/cmaesr)

## Description

The *cmaesr* package implements the popular [Covariance Matrix Adaption - Evolution Strategy](https://www.lri.fr/~hansen/cmatutorial.pdf) optimizer for numerical optimization problems in pure R. The main features of the package are:
* Extensible S3 based system for stopping conditions
* Possibility to enable [restarting](https://www.google.de/url?sa=t&rct=j&q=&esrc=s&source=web&cd=4&cad=rja&uact=8&ved=0CDgQFjADahUKEwiHyr2B3-fIAhVEOBoKHZFPBgs&url=https%3A%2F%2Fwww.lri.fr%2F~hansen%2Fcec2005ipopcmaes.pdf&usg=AFQjCNGwtYnwiRizaVZzbrfeXZjj-DYLtg&sig2=kMpEze_3Qe965UZ08wl-sw&bvm=bv.106130839,d.bGg) (flexible declaration of restart-triggering stopping conditions)
* Nice visualization of 2D optimization runs for, e.g., teaching purposes

## Installation Instructions

The package will be available in a first version at [CRAN](http://cran.r-project.org) soon. If you are interested in trying out and playing around with the current github developer version use the [devtools](https://github.com/hadley/devtools) package and type the following command in R:

```splus
devtools::install_github("jakobbossek/cmaesr")
```

## Example

Assume we want to minimize the 2D [Ackeley Function](http://www.sfu.ca/~ssurjano/ackley.html). To accomplish this task with *cmaesr* we need to define the objective function as a [smoof](https://github.com/jakobbossek/smoof) function. This function is then passed with some control arguments to the main function of the packge.

```splus
library(cmaesr)

# first generate the objective smoof function
fn = makeAckleyFunction(dimensions = 2L)
res = cmaes(
    fn, 
    monitor = makeSimpleMonitor(),
    control = list(
        sigma = 1.5, # initial step size
        lambda = 50, # number of offspring
        stop.ons = c(
            list(stopOnMaxIters(100)), # stop after 100 iteration ...
            getDefaultStoppingConditions() # ... or after some default stopping conditions
        )
    )
)
print(res)
```

For 2D functions a monitor for visualization is included.
```splus
library(cmaesr)

# generate the objective function
fn = makeSphereFunction(dimensions = 2L)
res = cmaes(
    fn,
    monitor = makeVisualizingMonitor(
        show.distribution = TRUE, show.last = TRUE
    ),
    control = list(
        sigma = 1, lambda = 100,
        stop.ons = list(stopOnMaxIters(15L))
    )
)
```

## Contact

Please address questions and missing features about the **cmaesr package** to the author Jakob Bossek <j.bossek@gmail.com>. Found some nasty bugs? Please use the [issue tracker](https://github.com/jakobbossek/cmaesr/issues) for this. Pay attention to explain the problem as good as possible. At its best you provide an example, so I can reproduce your problem.

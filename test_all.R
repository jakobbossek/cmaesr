library(methods)
library(devtools)
library(testthat)

if (interactive()) {
  load_all(".")
} else {
  library(cmaesr)
}

test_dir("tests/testthat")

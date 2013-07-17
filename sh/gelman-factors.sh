#!/bin/bash

if [  $# -eq 0 ]
then
  if [ -d out ]
    Rscript ../R/gelman-factors.r out
  else
    echo ERROR: must specify directory containing mcmc output.
  fi
elif [ -d $1 ]
then
  Rscript ../R/gelman-factors.r $1
else
  echo ERROR: directory, $1, does not exist.
fi
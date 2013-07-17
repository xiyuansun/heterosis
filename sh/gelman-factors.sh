#!/bin/bash

if [  $# -eq 0 ]
then
  if [ -d out ]
  then
    Rscript ../R/gelman-factors.r out
  else
    echo ERROR: must specify directory containing mcmc output.
    echo Default output directory, out/, not found.
  fi
elif [ -d $1 ]
then
  Rscript ../R/gelman-factors.r $1
else
  echo ERROR: directory, $1, does not exist.
fi
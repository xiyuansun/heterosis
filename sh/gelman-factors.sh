#!/bin/bash

if [  $# -eq 0 ]
then
  if [ -d out ]
  then
    Rscript ../R/gelman-factors.r out
  else
    echo ERROR: default output directory, out/, not found.
    echo Please input directory containing mcmc output.
  fi
elif [ -d $1 ]
then
  Rscript ../R/gelman-factors.r $1
else
  echo ERROR: directory, $1, does not exist.
fi
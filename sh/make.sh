#!/bin/bash

function cpu {
  echo Making CPU version.

  CC=gcc
  CFLAGS="-c -Wall -pedantic -I../include"
  LDFLAGS=-lm

  DEP=(printArrays)
  DEP+=(config getopts printConfig freeConfig)
  DEP+=(mySampleIntHost readGrp readData)
  DEP+=(allocChainHost newChainHost printChain freeChainHost)
  DEP+=(mu uniformHost normalHost gammaHost betaHost)
  DEP+=(cHost sigCHost epsHost etaHost dHost tauHost)
  DEP+=(phiHost alpHost delHost phiAlpDelJointHost phiAlpDel)
  DEP+=(thePhiHost theAlpHost theDelHost)
  DEP+=(sigPhiHost sigAlpHost sigDelHost)
  DEP+=(piAlpHost piDelHost)
  DEP+=(runChain oneChain summarizeChainHost)
  DEP+=(printProbs printRates printHyper printParms)
  DEP+=(main)

  OBJ=()

  if [ ! -d ../obj ]
  then
    mkdir ../obj
  fi

  if [ ! -d ../bin ]
  then
    mkdir ../bin
  fi

  for dep in ${DEP[@]}
  do
    OBJ+=(../obj/${dep}.o)
    ${CC} ../src/${dep}.c -o ../obj/${dep}.o ${CFLAGS} 
  done

  $CC ${OBJ[@]} -o ../bin/mcmc ${LDFLAGS}
}

function gpu {
  echo Coming soon...
}

function clean {

  if [ -d ../obj ]
  then
    rm -rf ../obj
  fi

  if [ -d ../bin ]
  then
    rm -rf ../bin
  fi
}


if [ $# -eq 0 ]
then
  cpu
elif [[ $1 =~ [cC][pP][uU] ]]
then
  cpu
elif [[ $1 =~ [gG][pP][uU] ]]
then
  gpu
elif [[ $1 =~ [cC][lL][eE][aA][nN] ]]
then
  clean
fi
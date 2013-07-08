#!/bin/bash

function mkdirs {
  
  if [ ! -d ../obj ]
  then
    mkdir ../obj
  fi

  if [ ! -d ../bin ]
  then
    mkdir ../bin
  fi

}


function cpu {
  echo Making CPU version...

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
  DEP+=(runChain oneChain summarizeChain)
  DEP+=(printProbsHost printRatesHost printHyperHost printParmsHost)
  DEP+=(main)

  OBJ=()

  for dep in ${DEP[@]}
  do
    OBJ+=(../obj/${dep}.o)
    ${CC} ../src/${dep}.c -o ../obj/${dep}.o ${CFLAGS} 
  done

  $CC ${OBJ[@]} -o ../bin/mcmc ${LDFLAGS}
}

function gpu {
  
  echo Making GPU version...

  CC=nvcc
  CFLAGS="-I../include -c -Wall -pedantic "
  LDFLAGS=-lm 

  DEP=(printArrays)
  DEP+=(config getopts printConfig freeConfig)
  DEP+=(readGrp readData)
  DEP+=(main)

  OBJ=()

  for dep in ${DEP[@]}
  do
    OBJ+=(../obj/${dep}.o)
    ${CC} ../src/${dep}.c -o ../obj/${dep}.o ${CFLAGS} 
  done

  $CC ${OBJ[@]} -o ../bin/gpu_mcmc ${LDFLAGS}

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

mkdirs

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
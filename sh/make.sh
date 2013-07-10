#!/bin/bash

function mkdirs {

  if [ ! -d ../bin ]
  then
    mkdir ../bin
  fi

  if [ ! -d ../obj ]
  then
    mkdir ../obj
  fi

  if [ ! -d ../obj/cpu ]
  then
    mkdir ../obj/cpu
  fi

  if [ ! -d ../obj/gpu ]
  then
    mkdir ../obj/gpu
  fi

}


function cpu {
  echo Making CPU version...

  CC=gcc
  CFLAGS="-c -Wall -pedantic -I../include/cpu"
  LDFLAGS=-lm 

  DEP=(printArrays)
  DEP+=(config getopts printConfig freeConfig)
  DEP+=(mySampleInt readGrp readData)
  DEP+=(allocChain newChain printChain freeChain)
  DEP+=(mu runiform rnormal rgamma rbeta)
  DEP+=(c sigC eps eta d tau)
  DEP+=(phi alp del phiAlpDelJoint phiAlpDel)
  DEP+=(thePhi theAlp theDel)
  DEP+=(sigPhi sigAlp sigDel)
  DEP+=(piAlp piDel)
  DEP+=(runChain oneChain summarizeChain)
  DEP+=(printProbs printRates printHyper printParms)
  DEP+=(main)

  OBJ=()

  for dep in ${DEP[@]}
  do
    OBJ+=(../obj/cpu/${dep}.o)
    ${CC} ../src/cpu/${dep}.c -o ../obj/cpu/${dep}.o ${CFLAGS} 
  done

  $CC ${OBJ[@]} -o ../bin/mcmc ${LDFLAGS}
}

function gpu {
  
  echo Making GPU version...

  CC=nvcc
  CFLAGS="-c -I../include/gpu -arch=sm_20"
  LDFLAGS=-lm 

  DEP=(printArrays)
  DEP+=(config getopts printConfig freeConfig)
  DEP+=(mySampleInt readGrp readData)
  DEP+=(allocChain chainDeviceToHost newChain printChain freeChain) 
  DEP+=(runiform rnormal rgamma rbeta)
  DEP+=(c tau piAlp piDel d) # sigC eps eta)
  DEP+=(thePhi) # theAlp theDel)
#  DEP+=(phi alp del phiAlpDelJoint phiAlpDel)

#  DEP+=(sigPhi sigAlp sigDel)
  DEP+=(runChain oneChain summarizeChain)
  DEP+=(printProbs printRates printHyper printParms)
  DEP+=(main)

  OBJ=()
 
  for dep in ${DEP[@]}
  do
    OBJ+=(../obj/gpu/${dep}.o)
    ${CC} ../src/gpu/${dep}.cu -o ../obj/gpu/${dep}.o ${CFLAGS} 
  done

  $CC ${OBJ[@]} -o ../bin/gpumcmc ${LDFLAGS}

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

  if [ -d ../out ]
  then
    rm -rf ../out
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
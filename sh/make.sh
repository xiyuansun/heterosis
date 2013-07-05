#!/bin/bash

function cpu {
  echo Making CPU version.

  CC=gcc
  CFLAGS="-c -Wall -pedantic -ansi -I../include"
  LDFLAGS=-lm

  DEP=(betaHost config freeConfig gammaHost getopts normalHost printConfig test uniformHost)
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

  $CC ${OBJ[@]} -o ../bin/test ${LDFLAGS}
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
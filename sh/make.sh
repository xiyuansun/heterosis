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
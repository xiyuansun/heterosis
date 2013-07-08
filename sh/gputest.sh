if [ ! -d ../out ]
then
  mkdir ../out
fi

make gpu
../bin/gpumcmc --data ../data/smallData.txt --group ../data/smallGroup.txt
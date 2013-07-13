git pull origin master
make -C .. gpu

../bin/gpumcmc --data ../data/mediumData.txt --group ../data/mediumGroup.txt -r -h -p -j
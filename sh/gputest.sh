git pull origin master
make -C .. gpu

../bin/gpumcmc --data ../data/data.txt --group ../data/group.txt -r -h -p
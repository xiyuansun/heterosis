make -C .. cpu
valgrind ../bin/mcmc --data ../data/mediumData.txt --group ../data/mediumGroup.txt -r -h -p --chains 1 -v
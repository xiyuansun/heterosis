make -C .. cpu
valgrind --track-origins=yes ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h -p --chains 1 -v
../bin/mcmc --data ../data/test/largeData.txt --group ../data/test/largeGroup.txt -r -h -p --chains 1 -v
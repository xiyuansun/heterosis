make -C .. cpu
valgrind --track-origins=yes ../bin/mcmc --data ../data/smallData.txt --group ../data/smallGroup.txt -r -h -p --chains 1 -v
make -C .. cpu

#valgrind --track-origins=yes ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h -p --chains 1 -v -t

../bin/mcmc --data ../data/test/mediumData.txt --group ../data/test/mediumGroup.txt -p --chains 1 -v -h -r -t 
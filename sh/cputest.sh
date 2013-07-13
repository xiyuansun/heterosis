make -C .. cpu

#valgrind --track-origins=yes ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h -p --chains 1 -v -t

../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -p --chains 1 -v -h -r -t -M 10000
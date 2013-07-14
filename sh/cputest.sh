make -C .. cpu

#valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt  -h -p --chains 1 -v  -j 

../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -p --chains 1 -v -h -r -t -M 10 
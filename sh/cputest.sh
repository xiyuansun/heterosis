make -C .. cpu

#valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -h -p --chains 3 -r -b 5 --diagnostics

#valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -h -p --chains 3 -r -b 5 

../bin/mcmc --data ../data/test/largeData.txt --group ../data/test/largeGroup.txt --chains 1 -M 100 -v -r -h -t --diagnostics
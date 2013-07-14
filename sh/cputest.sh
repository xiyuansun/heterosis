make -C .. cpu

#valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -h -p --chains 3 -r -b 5 --diagnostics

#valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -h -p --chains 3 -r -b 5 

 ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --chains 20 -M 10 --diagnostics 
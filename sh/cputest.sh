make -C .. cpu

#valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -h -p --chains 3 -r -b 5 --diagnostics

#valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -h -p --chains 3 -r -b 5 

valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --chains 3 -M 10 --diagnostics --tau 4

#valgrind --track-origins=yes ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt --chains 3 -M 10 --diagnostics 
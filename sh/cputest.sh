make -C .. cpu

valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 100 -c 100 -v --diagnostics

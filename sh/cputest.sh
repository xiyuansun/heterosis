make -C .. cpu

valgrind ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 20 -c 3 -h -p -v -P  -r -t --dic
make -C .. cpu

valgrind --track-origins=yes ../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 20 -c 3 -h -p -v --gelman -P --dic -r -t

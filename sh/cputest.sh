make -C .. cpu

../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 20 -c 1 -h --gelman

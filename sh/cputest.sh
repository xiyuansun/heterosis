make -C .. cpu

../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 20 -c 5 -v -p --gelman --dic

make -C .. gpu

../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 10 -c 10 -b 4 -h --dic
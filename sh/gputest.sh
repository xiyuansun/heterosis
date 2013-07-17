make -C .. gpu

../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroupDe.txt -M 20 -c 3 -h -p -P  -r -t --dic 

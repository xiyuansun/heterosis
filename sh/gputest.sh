make -C .. gpu

../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 3 -c 2 -h -p -P  -r -t --dic --phi-prior 4 --alpha-prior 5 --delta-prior 6

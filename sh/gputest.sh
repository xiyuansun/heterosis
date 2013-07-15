make -C .. gpu

../bin/gpu-mcmc --data ../data/test/largeData.txt --group ../data/test/largeGroup.txt -M 20 -c 1 -v -b 10 --diagnostics 

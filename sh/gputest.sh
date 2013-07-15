make -C .. gpu

../bin/gpu-mcmc --data ../data/test/largeData.txt --group ../data/test/largeGroup.txt -M 100 -c 10 -v -b 1000 -h -p
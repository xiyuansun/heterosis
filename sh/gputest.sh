make -C .. gpu

../bin/gpu-mcmc --data ../data/test/largeData.txt --group ../data/test/largeGroup.txt -h -p -P --dic -r -t -v -M 1000 -c 2
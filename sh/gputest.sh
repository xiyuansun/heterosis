make -C .. gpu

#cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h -p --chains 1 -v -t

../bin/gpu-mcmc --data ../data/test/largeData.txt --group ../data/test/largeGroup.txt -M 10 -c 3 -p -v -r -h -t --diagnostics

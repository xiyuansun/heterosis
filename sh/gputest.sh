make -C .. gpu

cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 20 -c 1 -v -b 10 --diagnostics

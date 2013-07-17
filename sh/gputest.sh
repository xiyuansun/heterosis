make -C .. gpu

cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 20 -c 3 -h -p -P  -r -t --dic

valgrind ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 20 -c 3 -h -p -P  -r -t --dic
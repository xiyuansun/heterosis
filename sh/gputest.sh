make -C .. clean
make -C .. gpu
cuda-memcheck ../bin/gpu-mcmc --data ../data/mediumData.txt --group ../data/mediumGroup.txt -r -h -p --chains 1 -v
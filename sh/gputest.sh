make -C .. gpu
cuda-memcheck ../bin/gpu-mcmc --data ../data/smallData.txt --group ../data/smallGroup.txt -r -h -p --chains 1 -v
../bin/gpu-mcmc --data ../data/data.txt --group ../data/group.txt -r -h -p --chains 1 -v
make -C .. gpu
# valgrind ../bin/gpu-mcmc --data ../data/smallData.txt --group ../data/smallGroup.txt -r -h -p --chains 1 -v
cuda-memcheck --show-backtrace=yes ../bin/gpu-mcmc --data ../data/smallData.txt --group ../data/smallGroup.txt -r -h -p --chains 1 -v
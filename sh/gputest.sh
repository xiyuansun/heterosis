make -C .. gpu

#cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h -p --chains 1 -v -t

../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 100 -c 1 -v -r -h -t 
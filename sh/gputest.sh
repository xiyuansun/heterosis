git pull origin master
make -C .. gpu

#cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h -p --chains 1 -v -t

#cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h -p --chains 1 -v -t -j

#cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt  -p --chains 1 -t

#cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h  --chains 1 -v -t

#../bin/gpu-mcmc --data ../data/test/mediumData.txt --group ../data/test/mediumGroup.txt -r -h -p --chains 1 -v -t -M 100

../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroupDe.txt -r -h -p --chains 1 -v -t -M 10
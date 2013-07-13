git pull origin master
make -C .. gpu

#valgrind ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -r -h -p --chains 1 -v -t

../bin/gpu-mcmc --data ../data/test/mediumData.txt --group ../data/test/mediumGroup.txt -r -h -p --chains 1 -v -t -M 100
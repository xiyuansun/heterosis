git pull origin master
make -C .. gpu

cuda-memcheck ../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 10 -c 3 -b 4 -p -P -r -t -s 100 -h -r --dic
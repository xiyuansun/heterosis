make -C .. gpu

../bin/gpu-mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -h -p -P --dic -r -t -j -o out2j

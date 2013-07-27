make -C .. cpu
rm -rf out

../bin/mcmc --data ../data/test/smallData.txt --group ../data/test/smallGroup.txt -M 15 -c 3 -h -p -P  -r -t --dic -b 11 --debug 3
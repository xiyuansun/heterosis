if [ ! -d ../out ]
then
  mkdir ../out
fi

make gpu

../bin/gpumcmc --data ../data/smallData.txt --group ../data/smallGroup.txt --probs ../out/probs.txt --hyper ../out/hyper.txt --rates ../out/rates.txt --all-parms ../out/allparms.txt --some-parms ../out/someparms.txt

diff ../out/allparms.txt ../goodout/allparms.txt
diff ../out/hyper.txt ../goodout/hyper.txt
diff ../out/probs.txt ../goodout/probs.txt
diff ../out/rates.txt ../goodout/rates.txt
diff ../out/someparms.txt ../goodout/someparms.txt
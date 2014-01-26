../bin/heterosis --data data.txt --group group.txt \
  --iter 2000000 --chains 3 --seed 0 \
  --hyper --verbose --burnin 1950000 \
  --sigma-c0 2 \
  --d0 100 \
  --a-tau 2 \
  --a-alpha 0.5 \
  --a-delta 0.5 \
  --b-tau 2 \
  --b-alpha 0.5 \
  --b-delta 0.5 \
  --gamma-phi 2\
  --gamma-alpha 2 \
  --gamma-delta 2 \
  --sigma-phi0 2 \
  --sigma-alpha0 2 \
  --sigma-delta0 2 \
./gelman.diag.sh
R CMD BATCH compareParms.r
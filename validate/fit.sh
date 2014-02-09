../bin/heterosis --data data.txt --group group.txt \
  --iter 5000 --chains 3 --seed 2 --verbose --hyper --parms \
  --sigma-c0 5 \
  --d0 10 \
  --a-tau 2 \
  --a-alpha 0.5 \
  --a-delta 0.5 \
  --b-tau 2 \
  --b-alpha 0.5 \
  --b-delta 0.5 \
  --gamma-phi 2 \
  --gamma-alpha 2 \
  --gamma-delta 2 \
  --sigma-alpha0 5 \
  --sigma-delta0 5 \
  --sigma-c 1.548 \
  --d 27.764 \
  --tau 0.9297 \
  --theta-alpha 0.6627 \
  --theta-delta 3.3358 \
  --sigma-alpha 1.148 \
  --sigma-delta 1.008 \
  --pi-alpha 0.175 \
  --pi-delta 0.3937 \
  --sigma-phi 1.025 \
  --theta-phi 0.6846 \
  
./gelman-diag.sh
R CMD BATCH compareParms.r
R CMD BATCH phi-alp-del.r
../bin/heterosis --data data.txt --group group.txt \
  --iter 500000 --chains 3 --seed 2 --verbose --hyper \
  --sigma-c 0.000000009272379 \
#  --d0 10 \
#  --a-tau 2 \
#  --a-alpha 0.5 \
#  --a-delta 0.5 \
#  --b-tau 2 \
#  --b-alpha 0.5 \
#  --b-delta 0.5 \
#  --gamma-phi 2 \
#  --gamma-alpha 2 \
#  --gamma-delta 2 \
#  --sigma-alpha0 5 \
#  --sigma-delta0 5 \
#  --d 5.157904 \
#  --tau 0.6852717 \
#  --theta-phi 5.018923 \
#  --theta-alpha 0.0001 \
#  --theta-delta 0.0001 \
#  --sigma-phi 1.666049 \
#  --sigma-alpha 1.829988 \
#  --sigma-delta 0.8901228 \
#  --pi-alpha 0.1101401 \
#  --pi-delta 0.5028068 \

  
mkdir -p fig
./gelman-diag.sh
# R CMD BATCH compareParms.r 
# R CMD BATCH phi-alp-del.r 
# R CMD BATCH c.r
# R CMD BATCH phi.r
# R CMD BATCH alp.r
# R CMD BATCH del.r
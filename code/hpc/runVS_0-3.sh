#!/bin/bash
# Script for running 6 chains of models 0-3 for variable selection

echo "Starting lc_runMod script"
#FILES=./lc/data/strat_15pct/*.Rdump
FILES=./lc/data/3km/*.Rdump
d=($FILES)
dfull=($(basename -a $FILES))
echo "loaded data:"
echo $FILES
echo " "
for f in {0..3}
do
  df="${dfull[$f]}"
  dname="${df%.Rdump}"
  out="./lc/out_3km/vs_"
  out+=$dname
  out+="_"
  for i in {1..6}
  do
    echo " "
    echo "||||---------- Starting" $f $dname "chain" $i
    echo " "
    outi=$out
    outi+=$i
  
    lc/lc_VS sample \
        algorithm=hmc \
          engine=nuts \
            metric=diag_e \
        num_samples=5000 \
        num_warmup=5000 \
      data file="${d[$f]}" \
      init=0 \
      random seed=`expr $i + 4337` \
      output file=$outi.csv \
        refresh=100 &

    sleep 10
    fm=`free | awk 'FNR == 2 {print int($4)}'`
    initTime=10
    echo "||||---------- Waiting..."
    while [ $fm -lt 100000000 ]
    do
      sleep 10
      fm=`free | awk 'FNR == 2 {print int($4)}'`
      initTime=`expr $initTime + 10`
    done
    echo "||||---------- Init time for" $f $dname $i "=" $initTime "s"
  done
done

echo " "
echo "||||---------- All gradients evaluated"
echo "||||---------- Finishing sampling"
echo " "
wait
	

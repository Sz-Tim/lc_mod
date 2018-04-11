#!/bin/bash
# Script for running 24 chains of the climate model

echo "Starting run_full script"
FILES=./lc/data/all/Clim.Rdump
d=($FILES)
dfull=($(basename -a $FILES))
echo "loaded data:"
echo $FILES
echo " "
f=0
  df="${dfull[$f]}"
  dname="${df%.Rdump}"
  out="./lc/out/full_"
  out+=$dname
  out+="_"
  for i in {1..24}
  do
    echo " "
    echo "||||---------- Starting" $f $dname "chain" $i
    echo " "
    outi=$out
    outi+=$i
  
    lc/lc_full sample \
        algorithm=hmc \
          engine=nuts \
            metric=diag_e \
        num_samples=2 \
        num_warmup=2 \
      data file="${d[$f]}" \
      init=0 \
      random seed=`expr $i + 4337` \
      output file=$outi.csv \
        refresh=1 &

    sleep 30
    fm=`free | awk 'FNR == 2 {print int($4)}'`
    initTime=30
    echo "||||---------- Waiting..."
    while [ $fm -lt 100000000 ]
    do
      sleep 10
      fm=`free | awk 'FNR == 2 {print int($4)}'`
      initTime=`expr $initTime + 10`
    done
    echo "||||---------- Init time for" $f $dname $i "=" $initTime "s"
  done

echo " "
echo "||||---------- All gradients evaluated"
echo "||||---------- Finishing sampling"
echo " "
wait
	

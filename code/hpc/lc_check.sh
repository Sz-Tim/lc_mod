#!/bin/bash

# This script monitors the progress of chains on the cluster
# It assumes the output is passed in as the first argument


echo "Chains at 5%:" `grep -o ' 5%' $1 | wc -l`
echo "Chains at 10%:" `grep -o '10%' $1 | wc -l`
echo "Chains at 15%:" `grep -o '15%' $1 | wc -l`
echo "Chains at 20%:" `grep -o '20%' $1 | wc -l`
echo "Chains at 25%:" `grep -o '25%' $1 | wc -l`
echo "Chains at 30%:" `grep -o '30%' $1 | wc -l`
echo "Chains at 35%:" `grep -o '35%' $1 | wc -l`
echo "Chains at 40%:" `grep -o '40%' $1 | wc -l`
echo "Chains at 45%:" `grep -o '45%' $1 | wc -l`
echo "Chains at 50%:" `grep -o '50%' $1 | wc -l`
echo "Chains at 55%:" `grep -o '55%' $1 | wc -l`
echo "Chains at 60%:" `grep -o '60%' $1 | wc -l`
echo "Chains at 65%:" `grep -o '65%' $1 | wc -l`
echo "Chains at 70%:" `grep -o '70%' $1 | wc -l`
echo "Chains at 75%:" `grep -o '75%' $1 | wc -l`
echo "Chains at 80%:" `grep -o '80%' $1 | wc -l`
echo "Chains at 85%:" `grep -o '85%' $1 | wc -l`
echo "Chains at 90%:" `grep -o '90%' $1 | wc -l`
echo "Chains at 95%:" `grep -o '95%' $1 | wc -l`
echo "Chains at 100%:" `grep -o '100%' $1 | wc -l`

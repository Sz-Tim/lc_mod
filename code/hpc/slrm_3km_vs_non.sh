#!/bin/bash

## Note - Slurm script comments require two hash symbols (##).  A single
## hash symbol immediately followed by SBATCH indicates an SBATCH
## directive.  "##SBATCH" indicates the SBATCH command is commented
## out and is inactive.

## For jobs running on a single node using multiple threads, the number of 
## tasks should be 1.  This reflects how many processes are running (1), and
## not how many threads that process will use.
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=21

## If it's likely your job will use more than 128GB of RAM, be sure
## to specify a minimum above this to ensure you are allocated a node
## with 512GB of RAM. Note: this value is specified in megabytes.
##SBATCH --mem=512000

## Normal Slurm options
## SBATCH -p shared
#SBATCH --job-name="vs_non"
#SBATCH --output=logs/vs_non.output

## Load the appropriate modules first.  Linuxbrew/colsa contains most
## programs, though some are contained within the anaconda/colsa
## module.  Refer to http://premise.sr.unh.edu for more info.
module purge
module load linuxbrew/colsa

## Instruct your program to make use of the number of desired threads.
## As your job will be allocated an entire node, this should normally
## be 24.
srun Rscript R/3km_vs_non.R

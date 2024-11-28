#!/bin/bash
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=80
#SBATCH --time=00:16:00
#SBATCH --job-name E3n4ht2
#SBATCH --output=/recurrent_project/E1/sim%j.txt
#SBATCH --mail-type=FAIL

module load intel/2019u4 gcc/8.3.0 r/4.1.2

Rscript Simulation.R --adjusted 0 --n 400 --time 2 --tv 0 --experiment 1 --BaseEvent 1 --BetaDeath 0 --BetaEvent 0 --reps 1000 --out "/recurrent_project/E1/sim/"

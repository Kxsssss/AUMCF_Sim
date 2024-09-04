#!/bin/bash
fin=~/Documents/AUMCF_Sim/SimConfig.txt

sed 1d ${fin} | while read line
do
	# Read parameter configuration.
	ad=$(echo ${line} | awk "{print \$1}")
	n=$(echo ${line} | awk "{print \$2}")
	t=$(echo ${line} | awk "{print \$3}")
	c=$(echo ${line} | awk "{print \$4}")
	f=$(echo ${line} | awk "{print \$5}")
	d0=$(echo ${line} | awk "{print \$6}")
	d1=$(echo ${line} | awk "{print \$7}")
	e0=$(echo ${line} | awk "{print \$8}")
	e1=$(echo ${line} | awk "{print \$9}")
	bd=$(echo ${line} | awk "{print \$10}")
	be=$(echo ${line} | awk "{print \$11}")
	

	# Run simulation.
	# Specify the number of simulation replicates 'reps', the number of bootstrap samples per replicate 'boot', 
	# and the output directory 'out'.
	Rscript ~/Documents/AUMCF_Sim/Simulation.R  --adjusted ${ad} --n ${n} --time ${t} --censor ${c} --BaseDeath0 ${d0} --BaseDeath1 ${d1} --BaseEvent0 ${d0} --BaseEvent1 ${d1} --BetaDeath ${bd} --BetaEvent ${be} --reps 1 --boot 2 --out "~/Documents/AUMCF_Sim/Sim/";
done
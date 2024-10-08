#!/bin/bash
fin=~/Documents/AUMCF_Sim/SimConfig.txt

sed 1d ${fin} | while read line
do
	# Read parameter configuration.
	ad=$(echo ${line} | awk "{print \$1}")
	n=$(echo ${line} | awk "{print \$2}")
	t=$(echo ${line} | awk "{print \$3}")
	f=$(echo ${line} | awk "{print \$4}")
	e1=$(echo ${line} | awk "{print \$5}")
	bd=$(echo ${line} | awk "{print \$6}")
	be=$(echo ${line} | awk "{print \$7}")
	tv=$(echo ${line} | awk "{print \$8}")
	ei=$(echo ${line} | awk "{print \$9}")
	

	# Run simulation.
	# Specify the number of simulation replicates 'reps', the number of bootstrap samples per replicate 'boot', 
	# and the output directory 'out'.
	Rscript ~/Documents/AUMCF_Sim/Simulation.R  --tv ${tv} --ei ${ei} --adjusted ${ad} --n ${n} --time ${t} --frailtyVar ${f} --BaseEvent1 ${d1} --BetaDeath ${bd} --BetaEvent ${be} --reps 1 --out "~/Documents/AUMCF_Sim/test/";
done

#!/bin/bash
fin=~/Desktop/recurrent_events/SimConfig.txt

sed 1d ${fin} | while read line
do
	# Read parameter configuration.
	n=$(echo ${line} | awk "{print \$1}")
	t=$(echo ${line} | awk "{print \$2}")
	c=$(echo ${line} | awk "{print \$3}")
	f=$(echo ${line} | awk "{print \$4}")
	d0=$(echo ${line} | awk "{print \$5}")
	d1=$(echo ${line} | awk "{print \$6}")
	e0=$(echo ${line} | awk "{print \$7}")
	e1=$(echo ${line} | awk "{print \$8}")
	bd=$(echo ${line} | awk "{print \$9}")
	be=$(echo ${line} | awk "{print \$10}")

	# Run simulation.
	# Specify the number of simulation replicates 'reps', the number of bootstrap samples per replicate 'boot', 
	# and the output directory 'out'.
	Rscript ~/Desktop/recurrent_events/Simulation.R --n ${n} --time ${t} --censor ${c} --BaseDeath0 ${d0} --BaseDeath1 ${d1} --BaseEvent0 ${d0} --BaseEvent1 ${d1} --BetaDeath ${bd} --BetaEvent ${be} --reps 1 --boot 2 --out "~/Desktop/recurrent_events/Sim/";
done

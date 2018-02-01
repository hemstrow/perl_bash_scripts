#!/bin/bash -l

if [ "$#" -ne 10 ]; then
	echo -e "Runs sim_mut2.pl multiple times with mulitple parameter sets by calling run_sim_mut2_ms.sh.\n\nRequired Arguments:\n\tnum_runs\n\tnum_gene_copies\n\tpop_size\n\tmutation_rate\n\tparameter_to_change\n\tchange_amount\n\tending_parameter_value\n\toutput_style\n\tsbatch_time\n\tTS_aTMRCA_reps\n\n\nparameter_to_change options:\n\t1: Number of gene copies.\n\t2: Population size\n\noutput_style options, divisible by:\n\t2: ms style\n\t3: raw gene copy info\n\t5: Tree Statistics"
	exit
fi

nruns=$1
ncop=$2
npop=$3
mut=$4
change=$5
by=$6
end=$7
os=$8
stime=$9
tsreps=$10

#set which variable to change from set to set.
if [ $change = 1 ]; then
	cur=$ncop
	echo -e "Changing number of gene copies by ${by} from ${cur} to ${end}.\n"
elif [ $change = 2 ]; then
	cur=$npop
	echo -e "Changing pop size by ${by} from ${cur} to ${end}.\n"

else
	exit
fi

echo -e "Submitting:"

#run sims for all of the variable ranges.
while [ $cur -le $end ] 
do
	echo -e "\n\t${cur}."

	#run sims
	if [ $change = 1 ]; then
	       sbatch -t ${stime} ~/scripts/run_sim_mut2.sh $nruns $cur $npop $mut $os sm2_out $tsreps
	elif [ $change = 2 ]; then
	       sbatch -t ${stime} ~/scripts/run_sim_mut2.sh $nruns $ncop $cur $mut $os sm2_out $tsreps
	else
	        exit
	fi

	#increase variable of interest. 
	cur=$(( $cur + $by ))
done


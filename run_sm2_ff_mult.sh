#!/bin/bash -l

if [ "$#" -lt 4 ]; then
	echo -e "Runs sim_mut2_ploidy.pl multiple times with mulitple parameter sets by calling run_sim_mut2p.sh.\n\nRequired Arguments:"
	echo -e "\tParameter File"
	echo -e "\tsbatch time argument"
	echo -e "\tOutput Filename"
	echo -e "\tOutput Style"
	echo -e "\tNumber of pairwise times to most recent common ancestor to calculate. Required if tree stats requested in output style."
	echo -e "\nParameter File contains the following information, delimited by a single space. Multiple lines can be given for multiple parameter sets to run."
	echo -e "\tnumber of gene copies"
	echo -e "\tpopulation size"
	echo -e "\tmutation rate"
	echo -e "\tploidy"
	echo -e "\tselfing, yes or no"
	echo -e "\tnumber of runs"
	echo -e "\nOutput Styles, number divisible by:"
	echo -e "\t2: ms"
	echo -e "\t3: gene copy info"
	echo -e "\t5: tree statistics"
	exit
fi

file=$1
stime=$2
output=$3
os=$4
nts=$5

echo -e "\n\nInput file: ${file}.\nRuntime: ${stime}.\nOutput: ${output}.\nos: ${os}.\nnts: ${nts}.\n\n"


if (( $os % 5 == 0 )); then
	re='^[0-9]+$'
	if [[ $nts =~ $re ]]; then
		if [[ "$nts" -lt 100 ]]; then
			if [[ "$nts" -eq 0 ]];then
				echo "ts requested, nts integer must be greater than 0!"
				exit
			else
				echo "WARNING: ts requested, nts should be greater than 100!."
			fi
		fi
	else
		echo "ts requested, nts must be an integer greater than 0."
		exit
	fi
fi


end=`wc -l ${file}`
end=${end% *}

echo -e "Submitting:"

x=1
while [ $x -le $end ]
do
	#get variables from current line of input file
        line=`sed "${x}q;d" ${file}`
        IFS=$' ' read -a params <<< $line
        n="${params[0]}"
        N="${params[1]}"
	m="${params[2]}"
	p="${params[3]}"
	s="${params[4]}"
	nruns="${params[5]}"

	echo -e "Run: ${x}. Parameters:\n\tgene copies: ${n}.\n\tpop size: ${N}.\n\tmutation rate: ${m}.\n\tploidy: ${p}.\n\tSelfing: ${s}.\n\tIterations: ${nruns}." 

	#run sims
	sbatch -t $stime --mem=16G ~/scripts/run_sim_mut2p.sh $nruns $n $N $m $p $s $os $output $nts

	echo -e "\n"

	#increase variable of interest. 
        x=$(( $x + 1 ))
done


#!/bin/bash -l

if [ "$#" -lt 3 ]; then
	echo -e "Runs ms multiple times with mulitple parameter sets.\n\nRequired Arguments:"
	echo -e "\tParameter File"
	echo -e "\tsbatch time argument"
	echo -e "\tOutput Filename"
	echo -e "\nParameter File contains the following information, delimited by a single space. Multiple lines can be given for multiple parameter sets to run."
	echo -e "\tnumber of gene copies"
	echo -e "\tpopulation size"
	echo -e "\tmutation rate"
	echo -e "\tnumber of runs"
	exit
fi

file=$1
stime=$2
output=$3

echo -e "\n\nInput file: ${file}.\nRuntime: ${stime}.\nOutput: ${output}.\n\n\n"


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
	nruns="${params[3]}"

	theta="4*${N}*${m}/2"
        theta=`perl -E "say ${theta}"`

	echo -e "Run: ${x}. Parameters:\n\tgene copies: ${n}.\n\tpop size: ${N}.\n\tmutation rate: ${m}.\n\tIterations: ${nruns}.\n\tTheta: ${theta}.\n" 

	#run sims
	mkdir ${n}_${N}_${m}
	echo "#!/bin/bash -l
        ~/ms/msdir/ms ${n} ${nruns} -t ${theta} > ${n}_${N}_${m}/${output}.ms" > ${n}_${N}_${m}/ms_run.sh
	sbatch -t $stime --mem=16G ${n}_${N}_${m}/ms_run.sh

	echo -e "\n"

        x=$(( $x + 1 ))
done


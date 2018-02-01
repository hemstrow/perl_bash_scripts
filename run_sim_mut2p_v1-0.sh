#!/bin/bash -l

nruns=$1
ncop=$2
npop=$3
mut=$4
ploid=$5
self=$6
os=$7
output=$8
tsreps=$9


outdir="${ncop}_${npop}_${mut}_${ploid}"
mkdir $outdir
cd $outdir

if (( $os % 2 == 0 )); then
	theta="4*${npop}*${mut}/2"
	theta=`perl -E "say ${theta}"`
	echo -e "ms ${ncop} ${nruns} -t ${theta}\n\n" > ${output}.ms
fi

x=1
while [ $x -le $nruns ] 
do
	
	c1=$x
	~/scripts/sim_mut2_ploidy.pl $ncop $npop $mut $ploid $self $os $tsreps > ${c1}.raw
	if (( $os % 2 == 0 )); then	
		grep -v "TreeStats" ${c1}.raw | grep -v "GC" > ${c1}.ms
		cat ${output}.ms ${c1}.ms > ${output}.ms.temp
		rm ${output}.ms
		rm ${c1}.ms
		mv ${output}.ms.temp ${output}.ms
	fi
	if (( $os % 3 == 0 )); then
		grep "GC" ${c1}.raw > ${c1}.gc
		echo -e "//" >> ${output}.gc
                cat ${output}.gc ${c1}.gc > ${output}.gc.temp
		rm ${output}.gc
                rm ${c1}.gc
                mv ${output}.gc.temp ${output}.gc
	fi
	if (( $os % 5 == 0 )); then
		grep "TreeStats" ${c1}.raw > ${c1}.ts
                cat ${output}.ts ${c1}.ts > ${output}.ts.temp
                rm ${output}.ts
                rm ${c1}.ts
                mv ${output}.ts.temp ${output}.ts
	fi
	rm ${c1}.raw
	x=$(( $x + 1 ))

done


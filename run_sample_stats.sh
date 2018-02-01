#!/bin/bash -l

list=$1
output=$2

wc=$(wc -l ${list} | awk '{print $1}')

x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p ${list}" 
	str=$($string)

	var=$(echo $str | awk -F"\t" '{print $1}')   
	set -- $var
	c1=$1
	echo "$c1"
	cat ${c1} | ~/ms/msdir/sample_stats > ${c1}.ss
	cat $output ${c1}.ss > ${output}.temp
	rm $output
	rm ${c1}.ss
	mv ${output}.temp $output
	x=$(( $x + 1 ))

done



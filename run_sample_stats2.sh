#!/bin/bash -l
#version of run_sample_stats that incorporates first three or four bits of information split by underscore in file names to columns in final output

if [ "$#" -lt 2 ]; then
        echo -e "Gets sample stats from ms outputs and adds columns containing the first three or four parts of information contained in the file handles.\n\nRequired Arguments:"
        echo -e "\tlist of ms files"
        echo -e "\toutput name"
        exit
fi




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
	reps=`wc -l ${c1}.ss`
	reps=${reps% *}
	~/scripts/coal_info_grabber.pl ${c1} $reps
	paste ${c1}.ss infofile.txt >> $output
	rm ${c1}.ss
	rm infofile.txt
	x=$(( $x + 1 ))

done



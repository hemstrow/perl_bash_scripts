#!/bin/bash -l
space=$1
name=$2

re='^[0-9]+$'
if [[ $name =~ $re ]]; then
	if [[ "$name" -le 5 ]]; then
		echo "integer, less than or equal to 5"
	else
		echo "integer, greater than 5"
	fi
else
	echo "not an integer"
fi

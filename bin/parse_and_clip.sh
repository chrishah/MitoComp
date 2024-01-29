#!/bin/bash

f=$1
log=$2
readlength=$3
outprefix=$4
## could allow to choose different criteria
mode=$5

clip=$(cat $log | grep "clip points" | awk '{print $9"\t"$12"\t"$22}' | sed 's/;//' | sed "s/'//" | sort -n -k 2 -k 1)

if [[ "$(echo "$clip" | wc -c)" -eq 1 ]]
then
	echo "no candidate found"
	touch $outprefix.fasta.missing
else
	total=$(cat $log | grep "^Processing" | cut -d " " -f 5)
	if [[ "$(echo "$clip" | wc -l)" -eq 1 ]]
	then
		echo "found single candidate (see $log for details):"
		echo "$clip"
		froto=$(echo "$clip" | cut -f 3)
		length=$(echo "$clip" | cut -f 1)
	else
		clip=$(echo "$clip" | tail -n 1)
		echo "found multiple candidates (see $log for details) - trying best:"
		echo "$clip"
		froto=$(echo "$clip" | cut -f 3)
		length=$(echo "$clip" | cut -f 1)
	fi

	if [[ $mode == "strict" ]]
	then
		echo "check if candidate fullfills our criteria"

		diff=$(( total - length ))
	
		if [[ $diff -lt $(( readlength + readlength + readlength )) ]]
		then
			echo "Candidate fullfills our criteria"
			circules.py -f $f -c $froto --prefix $outprefix
			sed -i "s/^>.*/>${outprefix}_circular_${length}/" $outprefix.circular.${length}.fasta
			echo "result in file: $(pwd)/$outprefix.fasta"
			ln -s $outprefix.circular.${length}.fasta $outprefix.fasta
		else
			echo "Candidate doesn't fullfill our criteria - we recommend manual curation"
			sed "s/^>.*/>${outprefix}_circular_${length}/" $f > $outprefix.fasta.needs_attention
			echo "see file: $(pwd)/$outprefix.fasta.needs_attention"
		fi
	elif [[ $mode == "best" ]]
	then
		echo "Extracting candidate of length $length bp"
		circules.py -f $f -c $froto --prefix $outprefix
		sed -i "s/^>.*/>${outprefix}_circular_${length}/" $outprefix.circular.${length}.fasta
		echo "result in file: $(pwd)/$outprefix.fasta"
		ln -s $outprefix.circular.${length}.fasta $outprefix.fasta
	fi
fi


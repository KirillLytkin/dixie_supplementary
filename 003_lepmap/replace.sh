#!/bin/bash

for file in *.vcf;
do

	cd "${file%%.*}"
	
	sed 's/1 2/12/' genotype.txt > variable.txt
	sed 's/1 2/12/' variable.txt > genotype.txt
	
	if grep -q 'LG = 2' order.txt
		then grep -B1000 '#*** LG = 2' order.txt > variable.txt
		sed '$d' variable.txt > order.txt
	fi

	rm variable.txt

	cd ..

done

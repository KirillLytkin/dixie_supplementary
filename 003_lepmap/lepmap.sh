#!/bin/bash
leppath=/home/user/Program/lepmap3
for file in *.vcf;
do
	mkdir "${file%%.*}"
	cd "${file%%.*}"
	java -cp $leppath/bin ParentCall2 data=../pedigree.txt vcfFile=../$file > data.call
	java -cp $leppath/bin Filtering2 data=data.call removeNonInformative=1 dataTolerance=0.001 > data_fil.call
	java -cp $leppath/bin SeparateChromosomes2 data=data_fil.call lodLimit=5 > map.txt
	sort map.txt|uniq -c|sort -n
	java -cp $leppath/bin OrderMarkers2 data=data_fil.call map=map.txt scale=2M/N > order.txt
	cut -f1,2 data.call > cut_call.txt
	awk -vfullData=1 -f ../map2genotypes.awk order.txt > genotype.txt
	cd ..
done

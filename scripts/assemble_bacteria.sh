#!/bin/bash

prefix=`echo $1 | rev | cut -d '/' -f 1 | rev`
echo $prefix

##Determine genome length
if [[ $prefix == *"Ecoli"* ]] ; then
    ##k-12 size
    gsize=4.6m
elif [[ $prefix == *"KLPN"* ]] ; then
    ##klpn size used for all others
    gsize=5.3m
elif [[ $prefix == *"AB"* ]] ; then
    ##Acinetobacter baumannii
    gsize=4m
elif [[ $prefix == *"cloacae"* ]]; then
    ##E. cloacae
    gsize=5.3m
elif [[ $prefix == *"Citrobacter"* ]]; then
    gsize=5.2m
elif [[ $prefix == *"Pantoea"* ]]; then
    gsize=3.9m 
elif [[ $prefix == *"amalon"* ]]; then
    gsize=5.6m
elif [[ $prefix == *"radio"* ]]; then
    gsize=5.1m
elif [[ $prefix == *"aero"* ]]; then
    gsize=4.8m
fi


##Assemble if it's not already done
if [ -f $1/fastqs/$prefix.fq ] ; then
    canu \
	-p $prefix -d $1/canu_assembly \
	-gridOptions="--time=22:00:00 --account=mschatz1" \
	genomeSize=$gsize \
	corMemory=30 \
	-nanopore-raw $1/fastqs/$prefix.fq
fi

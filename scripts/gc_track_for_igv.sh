#!/bin/bash

# Split the input into filename and its path
fbase=$(dirname $1)
fname=$(basename $1)

# generate a fasta index (including sizes in second column)
samtools faidx ${fbase}/${fname}
 
# extract sequence names (chromosomes) and their lengths to a new file
cut -f 1,2 ${fbase}/${fname}.fai > ${fbase}/${fname}_sizes

# create a folder to store all results
basedir=${fbase}/gcpercent
mkdir -p ${basedir}
 
# if no other input was given make 1kb wide bins
defvalue=1000
width=${2:-$defvalue}

# create a BED file with windows of wished width
bedtools makewindows \
	-g ${fbase}/${fname}_sizes \
	-w ${width} \
	> ${basedir}/${fname}_${width}bps.bed

#create a bedtools nuc file
bedtools nuc \
	-fi ${fbase}/${fname} \
	-bed ${basedir}/${fname}_${width}bps.bed  \
	> ${basedir}/${fname}_nuc.txt

# width was defined earlier (default 1000), hence the title column filled with 'GCpc_1000bps'
gawk -v w=${width} 'BEGIN{FS="\t"; OFS="\t"}
{
if (FNR>1) {print $1,$2,$3,"GCpc_"w"bps",$5}
}' ${basedir}/${fname}_nuc.txt > ${basedir}/${fname}_${width}bps.igv

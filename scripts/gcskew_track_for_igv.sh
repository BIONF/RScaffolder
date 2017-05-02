#!/bin/bash

# Written by Bardya Djahanschiri, based on suggestions on http://wiki.bits.vib.be/index.php/Create_a_GC_content_track
# This script takes a (genomic) DNA sequence and uses samtools, bedtools and gawk to generate a file for the IGV
# It first collects the sequences and the sizes within the (multisequence -) fastafile (ARG1) and then uses bedtools to
# generate non-overlapping windows of defined length (ARG2). Lastly it computes a 


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

# width was defined earlier (default 1000), hence the title column filled with 'GCpc_widthbps'
gawk -v w=${width} 'BEGIN{FS="\t"; OFS="\t"}
{
if (FNR>1) {print $1,$2,$3,"GC_skew_"w"bps", (($7-$8)/100)/$5}
}' ${basedir}/${fname}_nuc.txt > ${basedir}/${fname}_${width}bps_gcskew.igv

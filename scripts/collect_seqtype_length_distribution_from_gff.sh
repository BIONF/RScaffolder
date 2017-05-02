#!/bin/bash

# Split the input into filename and its path
fbase=$(dirname $2)
fname=$(basename $2)

# if no other input was given select gene as seqtype
deftype='gene'
seqtype=${1:-$deftype}

defout=${fbase}'/'$1'_length_distribution.txt'
outpath=${3:-$defout}

defin=${fbase}'/mcoords_mapped_contigs.gff3'
inpath=${2:-$defin}

awk -F'\t' -v st=$(echo "$seqtype") '{ if ( $3 == st) print; }' $inpath | cut -s -f 4,5 | perl -ne '@v = split(/\t/); if ($v[1] > $v[0]) { printf("%d\n", $v[1] - $v[0] + 1) } else { printf("%d\n", $v[0] - $v[1] + 1)}' > $outpath ;

sort $outpath -o $outpath ;

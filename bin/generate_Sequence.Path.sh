#!/usr/bin/env zsh

genome=$1

if [[ -z "$genome" ]]; then
   echo ""
   echo "Usage:"
   echo "./generate_Sequence.Path.sh [genome folder]"
   echo ""
   echo "Example:"
   echo "./generate_Sequence.Path.sh /desktop/h19" 
   echo ""
   exit
fi

BASEDIR=$(readlink -f $genome)
cat <(echo "Chromosome name""\t""path") > path_header.tmp
ls $genome | awk '{print $1}' | sed 's/\.fa//g' > Sequence.Path.tmp1
ls $genome | awk '{print aa "/"$1}' aa=$(echo $BASEDIR) > Sequence.Path.tmp2
paste Sequence.Path.tmp1 Sequence.Path.tmp2 | tr ' ' \\t > Sequence.Path.tmp3
cat path_header.tmp Sequence.Path.tmp3 > Sequence.Path
rm -r -f path_header.tmp
rm -r -f Sequence.Path.tmp1
rm -r -f Sequence.Path.tmp2
rm -r -f Sequence.Path.tmp3





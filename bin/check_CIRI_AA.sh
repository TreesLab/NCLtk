#!/usr/bin/env bash
CIRI_col5=$1
ambiguous_align=$2
echo -n > "ambigous_align_"$CIRI_col5
head=$(echo "circRNA_ID strand gene_id	#junction_reads	junction_reads_ID #ambiguous_reads ambiguous_reads_ID")
cat <(echo $head) | tr ' ' \\t > "ambigous_align_"$CIRI_col5
num=$(cat $CIRI_col5 | wc -l) 
for i in $(seq 1 $num)
do
   takeOne=$(cat $CIRI_col5 | head -n $i | tail -n 1)
   AA_count=$(join <(echo $takeOne| awk '{print $5}' | tr ',' \\n | sed '/^$/d' | sort -k1,1) <(cat $ambiguous_align| sort -k1,1) | sort | uniq  | wc -l)  

   if test $AA_count -gt 0 ; then
      AA_read=$(join <(echo $takeOne | awk '{print $5}' | tr ',' \\n | sed '/^$/d' | sort -k1,1) <(cat $ambiguous_align | sort -k1,1) | sort | uniq  | tr '\n' '\,')
   else
      AA_read=$(echo "na")
   fi
   
   paste <(echo $takeOne) <(echo $AA_count)  <(echo $AA_read) | tr ' ' \\t >> "ambigous_align_"$CIRI_col5
done



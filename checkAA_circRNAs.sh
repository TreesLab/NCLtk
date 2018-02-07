#!/usr/bin/env zsh
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -circRNAs) 
     circRNAs=$(readlink -f $2)
     shift
     ;; 
     -thread) 
     thread=$2
     shift
     ;; 
     -genome|) 
     genome=$(readlink -f $2)
     shift
     ;; 
     -genome_chr|) 
     genome_chr=$(readlink -f $2)
     shift
     ;; 
     -others) 
     others=$(readlink -f $2)
     shift
     ;; 
     -blat) 
     blat_dir=$(readlink -f $2)
     shift
     ;;
     -tools)
     tools_bin=$(readlink -f $2)
     shift
     ;;
     
     *)

esac
shift
done

if [[ -z "$circRNAs" ]]; then
   echo ""
   echo "Usage:"
   echo "./checkAA_circRNAs.sh -circRNAs [circRNAs.txt] -thread [number of thread] -genome [genome.fa] -genome_chr [chromosome folder] -others [transcriptome folder] -blat [blat link] -tools [bin path]"
   echo ""
   exit
fi

nameout=$(echo $circRNAs | sed 's/\//\n/g' | tail -n 1)
echo "Step1 : fetch flanking +- 100 read"
cat $circRNAs | sed 's/\r//g' > circ.tmp1
cat circ.tmp1 | awk '{print $1 "\t" $2}' | sed 's/:/\t/g' | sed 's/|/\t/g' > circ.tmp2 
cat circ.tmp2 | tr ' ' \\t | grep '+' | awk '$2>$3 {print $0}' > circ_plus.tmp1 
cat circ.tmp2 | tr ' ' \\t | grep '+' | awk '$2<$3 {print $1 "\t" $3 "\t" $2 "\t" $4}' > circ_plus.tmp2
cat circ.tmp2 | tr ' ' \\t | grep '-' | awk '$2<$3 {print $0}' > circ_minus.tmp1
cat circ.tmp2 | tr ' ' \\t | grep '-' | awk '$2>$3 {print $1 "\t" $3 "\t" $2 "\t" $4}' > circ_minus.tmp2
cat circ_plus.tmp1 circ_plus.tmp2 circ_minus.tmp1 circ_minus.tmp2 | awk '{print $0}' > circ.tmp3
cat circ.tmp3 | grep '+' | awk '{print $0 "\t" $1 "\t" $2-99 "\t" $2}' > circ_plus_100.tmp1
cat circ.tmp3 | grep '+' | awk '{print $0 "\t" $1 "\t" $3 "\t" $3+99}' > circ_plus_100.tmp2
cat circ.tmp3 | grep '-' | awk '{print $0 "\t" $1 "\t" $2 "\t" $2+99}' > circ_minus_100.tmp1
cat circ.tmp3 | grep '-' | awk '{print $0 "\t" $1 "\t" $3-99 "\t" $3}' > circ_minus_100.tmp2
cat circ_plus_100.tmp1 circ_minus_100.tmp1 | awk '{print $1":"$2":"$4":"$3".1" "\t" $5 "\t" $6 "\t" $7 "\t" $4}' > circ_100.tmp1
cat circ_plus_100.tmp2 circ_minus_100.tmp2 | awk '{print $1":"$2":"$4":"$3".2" "\t" $5 "\t" $6 "\t" $7 "\t" $4}' > circ_100.tmp2
cat circ_100.tmp1 circ_100.tmp2 | sort -k1,1  > circ_200_intervel.txt

$tools_bin/generate_Sequence.Path.sh $genome_chr
$tools_bin/Getseq circ_200_intervel.txt Sequence.Path circ_200.fa
$tools_bin/merge_paired_sequences.py circ_200.fa circ_200.merged.fa

echo ""
echo  "Step2: runs blat"
mkdir BLAT_tmp
cat $genome | $tools_bin/SeqOut chrM.list 1 > chrM.fa
cat <(echo ">RepChrM") <(cat chrM.fa | sed -e '1d' ) <(cat chrM.fa | sed -e '1d') > RepChrM.fa  
cat $others/*.fa RepChrM.fa > others.fa

$tools_bin/mp_blat.py $genome circ_200.merged.fa circ.rG.1.psl -p $thread --blat_bin $blat_dir --tmp_path BLAT_tmp
$tools_bin/mp_blat.py $genome circ_200.merged.fa circ.rG.2.psl -p $thread --blat_bin $blat_dir --blat_opt "-tileSize=9 -stepSize=9 -repMatch=32768" --tmp_path BLAT_tmp
$tools_bin/mp_blat.py others.fa circ_200.merged.fa circ.rO.1.psl -p $thread --blat_bin $blat_dir --tmp_path BLAT_tmp
$tools_bin/mp_blat.py others.fa circ_200.merged.fa circ.rO.2.psl -p $thread --blat_bin $blat_dir --blat_opt "-tileSize=9 -stepSize=9 -repMatch=32768" --tmp_path BLAT_tmp

echo ""
echo "Step3: evaluate blat results"
cat circ.rG.1.psl | sed -e '1,5d' | awk 'substr($14,0,2)=="GL"' | awk '($1+$3)>=100' | awk '{print $10}' | sort  | uniq  > GL.rG.1.list
cat circ.rG.1.psl | sed -e '1,5d' | awk '($13-$12)/$11 > 0.9 {print $10}' | sort  | uniq  > colinear.rG.1.list
cat circ.rG.1.psl | sed -e '1,5d' | $tools_bin/PslChimeraFilter 30 5 > circ.chi0.1.bed
cat circ.rG.1.psl | sed -e '1,5d' | $tools_bin/RemoveInList 10 circ.chi0.1.bed 4  | $tools_bin/RemoveInList 10 colinear.rG.1.list 1 | $tools_bin/RemoveInList 10 GL.rG.1.list 1 > circ.rG.1_ambigous.psl 
cat circ.rG.1_ambigous.psl | awk '{print $10}' | sort  | uniq -c | awk '$1=="1" {print $2}' > badBlat.rG.1.list
cat circ.rG.1_ambigous.psl | $tools_bin/RemoveInList 10 badBlat.rG.1.list 1 | awk '{print $10}' | sort  | uniq > multipleHit.rG.1.list

cat circ.rG.2.psl | sed -e '1,5d' | awk 'substr($14,0,2)=="GL"' | awk '($1+$3)>=100' | awk '{print $10}' | sort  | uniq  > GL.rG.2.list
cat circ.rG.2.psl | sed -e '1,5d' | awk '($13-$12)/$11 > 0.9 {print $10}' | sort  | uniq > colinear.rG.2.list
cat circ.rG.2.psl | sed -e '1,5d' | $tools_bin/PslChimeraFilter 30 5 > circ.chi0.2.bed
cat circ.rG.2.psl | sed -e '1,5d' | $tools_bin/RemoveInList 10 circ.chi0.2.bed 4  | $tools_bin/RemoveInList 10 colinear.rG.2.list 1 | $tools_bin/RemoveInList 10 GL.rG.2.list 1 > circ.rG.2_ambigous.psl 
cat circ.rG.2_ambigous.psl | awk '{print $10}' | sort  | uniq -c | awk '$1=="1" {print $2}' > badBlat.rG.2.list
cat circ.rG.2_ambigous.psl | $tools_bin/RemoveInList 10 badBlat.rG.2.list 1 | awk '{print $10}' | sort  | uniq > multipleHit.rG.2.list

cat circ.rO.1.psl | sed -e '1,5d' | awk '($13-$12)/$11 > 0.9 {print $10}' | awk '{print $10}' | sort  | uniq  > colinear.rO.1.list
cat circ.rO.2.psl | sed -e '1,5d' | awk '($13-$12)/$11 > 0.9 {print $10}' | sort  | uniq  > colinear.rO.2.list
cat GL.rG.1.list GL.rG.2.list | sort | uniq > GL.rG.list
cat colinear.rG.1.list colinear.rG.2.list colinear.rO.1.list colinear.rO.2.list | sort  | uniq | sed '/^$/d'  > colinear.rGO.list
cat multipleHit.rG.1.list multipleHit.rG.2.list | sort  | uniq  > multipleHit.rG.list
cat colinear.rGO.list | sed '/^$/d'> colinear.clear.list
cat multipleHit.rG.list GL.rG.list | sort  | uniq > multipleHit.rG.GL.list 
join multipleHit.rG.GL.list colinear.clear.list -v1 | sed '/^$/d' > multipleHit.clear.list
cat colinear.clear.list multipleHit.clear.list | sort | uniq > ambiguous_align_circRNAs.tmp
cat ambiguous_align_circRNAs.tmp | grep '+' | sed 's/:/\t/g' | awk '{print $1 ":" $4 "|" $2 "\t" $3}' > ambiguous_align_circRNAs_sense.tmp
cat ambiguous_align_circRNAs.tmp | grep '-' | sed 's/:/\t/g' | awk '{print $1 ":" $2 "|" $4 "\t" $3}' > ambiguous_align_circRNAs_antisense.tmp
cat ambiguous_align_circRNAs_sense.tmp ambiguous_align_circRNAs_antisense.tmp | sort -k1,1 > ambiguous_align_$nameout

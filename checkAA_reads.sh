#!/usr/bin/env zsh
while [[ $# -gt 1 ]]
do
key="$1"

case $key in
     -circReads) 
     circReads=$(readlink -f $2)
     shift
     ;; 
     -thread) 
     Thread=$2
     shift
     ;; 
     -genome|) 
     genome=$(readlink -f $2)
     shift
     ;; 
     -others) 
     others=$(readlink -f $2)
     shift
     ;; 
     -read1) 
     read1=$(readlink -f $2)
     shift
     ;; 
     -read2) 
     read2=$(readlink -f $2)
     shift
     ;;
     -blat) 
     blat_dir=$(readlink -f $2)
     shift
     ;;
     -seqtk)
     seqtk_dir=$(readlink -f $2)
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

if [[ -z "$circReads" ]]; then
   echo ""
   echo "Usage:"
   echo "./checkAA_reads.sh -circReads [circReads.txt] -thread [number of thread] -genome [genome.fa] -others [transcriptome folder] -read1 [read1 fastq.gz file] -read2 [read2 fastq.gz file]  -blat [blat link] -seqtk [seqtk link] -tools [bin path]"
   echo ""
   exit
fi

echo ""
echo "Step1: extract backspliced junction reads"
cat $circReads | awk '{print $5}' | sed -e '1d' | tr ',' \\n | sed '/^$/d' | sort | uniq > readID.txt


check_readID=$(cat $read1 | grep '/1')
if [[ -z "$check_readID" ]]; then 
   cat readID.txt > readID_1.tmp
   cat readID.txt > readID_2.tmp
fi

cat readID.txt | awk '{print $0"/1 "}' > readID_1.tmp
cat readID.txt | awk '{print $0"/2 "}' > readID_2.tmp
$seqtk_dir subseq  $read1 readID_1.tmp > readID_1.fastq
$seqtk_dir subseq  $read2 readID_2.tmp > readID_2.fastq
$tools_bin/AssembleFastq readID_1.fastq readID_2.fastq > readID.fa

echo ""
echo "Step2: runs blat"
mkdir BLAT_tmp
cat $genome | $tools_bin/SeqOut chrM.list 1 > chrM.fa
cat <(echo ">RepChrM") <(cat chrM.fa | sed -e '1d' ) <(cat chrM.fa | sed -e '1d') > RepChrM.fa  
cat $others/*.fa RepChrM.fa > others.fa

$tools_bin/mp_blat.py $genome readID.fa circ.rG.1.psl -p $Thread --blat_bin $blat_dir --tmp_path BLAT_tmp
$tools_bin/mp_blat.py $genome readID.fa circ.rG.2.psl -p $Thread --blat_bin $blat_dir --blat_opt "-tileSize=9 -stepSize=9 -repMatch=32768" --tmp_path BLAT_tmp
$tools_bin/mp_blat.py others.fa readID.fa circ.rO.1.psl -p $Thread --blat_bin $blat_dir --tmp_path BLAT_tmp
$tools_bin/mp_blat.py others.fa readID.fa circ.rO.2.psl -p $thread --blat_bin $blat_dir --blat_opt "-tileSize=9 -stepSize=9 -repMatch=32768" --tmp_path BLAT_tmp

echo ""
echo "Step3: evaluate blat results"
cat circ.rG.1.psl | sed -e '1,5d' | awk 'substr($14,0,2)=="GL"' | awk '($1+$3)>=100' | awk '{print $10}' | sort  | uniq  > GL.rG.1.list
cat circ.rG.1.psl | sed -e '1,5d' | awk '($13-$12)/$11 > 0.8 {print $10}' | sort  | uniq  > colinear.rG.1.list
cat circ.rG.1.psl | sed -e '1,5d' | $tools_bin/PslChimeraFilter 30 5 > circ.chi0.1.bed
cat circ.rG.1.psl | sed -e '1,5d' | $tools_bin/RemoveInList 10 circ.chi0.1.bed 4 | $tools_bin/RemoveInList 10 colinear.rG.1.list 1 | $tools_bin/RemoveInList 10 GL.rG.1.list 1 > circ.rG.1_ambigous.psl 
cat circ.rG.1_ambigous.psl | awk '{print $10}' | sort  | uniq -c | awk '$1=="1" {print $2}' > badBlat.rG.1.list
cat circ.rG.1_ambigous.psl | $tools_bin/RemoveInList 10 badBlat.rG.1.list 1 | awk '{print $10}' | sort  | uniq > multipleHit.rG.1.list

cat circ.rG.2.psl | sed -e '1,5d' | awk 'substr($14,0,2)=="GL"' | awk '($1+$3)>=100' | awk '{print $10}' | sort  | uniq  > GL.rG.2.list
cat circ.rG.2.psl | sed -e '1,5d' | awk '($13-$12)/$11 > 0.8 {print $10}' | sort  | uniq > colinear.rG.2.list
cat circ.rG.2.psl | sed -e '1,5d' | $tools_bin/PslChimeraFilter 30 5 > circ.chi0.2.bed
cat circ.rG.2.psl | sed -e '1,5d' | $tools_bin/RemoveInList 10 circ.chi0.2.bed 4  | $tools_bin/RemoveInList 10 colinear.rG.2.list 1 | $tools_bin/RemoveInList 10 GL.rG.2.list 1 > circ.rG.2_ambigous.psl 
cat circ.rG.2_ambigous.psl | awk '{print $10}' | sort  | uniq -c | awk '$1=="1" {print $2}' > badBlat.rG.2.list
cat circ.rG.2_ambigous.psl | $tools_bin/RemoveInList 10 badBlat.rG.2.list 1 | awk '{print $10}' | sort  | uniq > multipleHit.rG.2.list

cat circ.rO.1.psl | sed -e '1,5d' | awk '($13-$12)/$11 > 0.8 {print $10}' | awk '{print $10}' | sort  | uniq  > colinear.rO.1.list
cat circ.rO.2.psl | sed -e '1,5d' | awk '($13-$12)/$11 > 0.8 {print $10}' | sort  | uniq  > colinear.rO.2.list
cat GL.rG.1.list GL.rG.2.list | sort | uniq > GL.rG.list
cat colinear.rG.1.list colinear.rG.2.list colinear.rO.1.list colinear.rO.2.list | sort  | uniq | sed '/^$/d'  > colinear.rGO.list
cat multipleHit.rG.1.list multipleHit.rG.2.list | sort  | uniq  > multipleHit.rG.list
cat colinear.rGO.list | sed '/^$/d' > colinear.clear.list
cat multipleHit.rG.list GL.rG.list | sort  | uniq > multipleHit.rG.GL.list 
join multipleHit.rG.GL.list colinear.rGO.clear.list -v1 | sed '/^$/d' > multipleHit.clear.list
cat colinear.clear.list multipleHit.clear.list  > ambiguous_align.list

echo ""
echo "Final output: ambiguous_align_InputFileName "
$tools_bin/check_CIRI_AA.sh $circReads ambiguous_align.list 


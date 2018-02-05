# NCLtk
#### a toolkit to examine ambiguous alignment in non-co-linear trascripts (circular RNAs and trans-splicing)
--------------------------------------------------
### two modules
##### 1. checking ambiguous alignment in the genomic positions of circRNAs (checkAA_circRNAs.sh):
##### to concatenate the flanking the circRNAs (within -100 nucleotides to +100 nucleotides for each circRNA's junction) against the reference genome and NCBI Refseq-/GENCODE-/Ensembl- identified mRNAs and ncRNAs using BLAT.
##### its criteria of ambiguous alignment: 
  + alternative co-linear explanation
      ##### to detect alternative co-linear explanation ovelapped more than 90% of the concatenated seq
  + multiple hits
     ##### to mapped to multiple positions with similiar BLAT mapping scores (difference of mapping scores < 5)
##### 2. checking ambiguous alignemnt in the supporting backedspliced junctions' RNA-seqs (checkAA_read.sh): 
##### to concatenate the supporting backspliced junctions' paired-end reads (the concatenated read may be shorter than length of the paired-end reads if the paired-end reads have the overlapped part) against the reference genome and  NCBI Refseq-/GENCODE-/Ensembl- identified mRNAs and ncRNAs using BLAT.
##### its criteria of ambiguous alignment: 
  + alternative co-linear explanation
      ##### to detect alternative co-linear explanation ovelapped more than 80% of the concatenated seq
  + multiple hits
     ##### to mapped to multiple positions with similiar BLAT mapping scores (difference of mapping scores < 5)

#### 1. System requirements 
The two shell scripts (check_circRNAs.sh and check_reads.sh) runs under Linux-like environment (i.e. Bio-Linux, also see http://environmentalomics.org/bio-linux/) with at least 30 GB RAM. 

#### 2. Installation
```sh
$ tar zxvf NCLtk.tar.gz
$ cd NCLtk
$ chmod +x *.sh
$ cd bin
$ make 
$ chmod +x *
```

#### 3. Installation external tools
(1) BLAT (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/)
(2) seqtk (https://github.com/lh3/seqtk)

#### 4. Preparation
(1) Chromosome fasta files 
```
$ mkdir hg19
$ cd hg19
$ wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/*
$ gunzip *.gz
```
(2) Genome fasta file 
```
$ wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz
$ gunzip GRCh37.p13.genome.fa.gz
```
(3) Transcritome files
```
$ mkdir others
$ cd others
$ wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
$ wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz
$ wget ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.lncRNA_transcripts.fa.gz
$ wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/mrna.fa.gz
$ wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/refMrna.fa.gz
$ wget ftp://ftp.ensembl.org/pub/grch37/release-91/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh37.cdna.all.fa.gz
$ wget ftp://ftp.ensembl.org/pub/grch37/release-91/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh37.ncrna.fa.gz
$ gunzip *.gz
```
    
#### 5. Execution of checkAA_circRNAs
Usage : 
```sh
$ ./checkAA_circRNAs.sh -circRNAs [circRNAs.txt] -thread [number of threads] -genome [genome] -genome_chr [chromosome folder] -others [transcriptome folder] -blat [blat link] -tools [bib path]
```
An example:
```sh
$ ./checkAA_circRNAs.sh -circRNAs examples/circ1046.txt -thread 6 -genome /path/to/GRCh37.p13.genome.fa -genome_chr /path/to/hg19/ -others /path/to/others/ -blat /path/to/blat -tools /path/to/bin/
```
Input format (e.g. circ1046.txt)

|circID                 | strand|
|-----------------------|-------|
|chr1:93297917\|93307481 | + |
|chr1:52954581\|53019130 | - |

Output file: ambiguous_align_circ1046.txt
#### 6. Execution of checkAA_reads 
Usage : 
```sh
$ ./checkAA_reads.sh -circReads [circReads.txt] -genome [genome.fa] -others [transcriptome folder] -read1 [read1.fastq or fastq.gz] -read2 [read2.fastq or fastq.gz] -blat [blat link] -seqtk [seqtk link ] -tools [bin path]
```
An example:
```sh
$ ./checkAA_reads.sh -circReads examples/CIRI_positive.txt -genome /path/to/GRCh37.p13.genome.fa -others /path/to/others/ -read1 /path/to/pos_1.fastq.gz -read2 /path/to/pos_2.fastq.gz -blat /path/to/blat -seqtk /path/to/seqtk -tools /path/to/bin/
```
Input format (e.g. CIRI_positive.txt)

|circRNA_ID | strand | gene_id | #junction_reads | junction_reads_ID |
|-------|--------|-----------|------------------------------------|---------------------------------|
|chr1:8928047\|8932045|		-|	ENSG00000074800.9,|	3|	simulate:21841,simulate:21844,simulate:21846,|

Output format (e.g. ambiguous_align_CIRI_positive.txt)

|circRNA_ID|	strand|	gene_id|	#junction_reads|	junction_reads_ID|	#ambiguous_reads| ambiguous_reads_ID|
|--------| --------|--------|-----|---------|---------|--------| 
|chr1:865535\|866469|	+	|ENSG00000187634.6,|	6|	simulate:69646,simulate:69647,simulate:69648,simulate:69650,simulate:69653,simulate:69657,|	0|	na|
|chr1:871152\|874509|	+	|ENSG00000187634.6,|	5|	simulate:69678,simulate:69679,simulate:69682,simulate:69683,simulate:69686,|	1|	simulate:69679,|





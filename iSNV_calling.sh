#!/bin/bash


echo 
echo "********************************************************"
echo "*       Intrahost SNV calling for amplicon-seq         *"
echo "********************************************************"
echo "nim, 201709 v2"
echo

if [ $# == 0 ]; then
	customStep="1,2,3,4,5"
else
	customStep=$1
	echo
	echo ~~~~~~~ custom step selection = $1 ~~~~~~~
	echo 
fi

fastqPath="./sample_reads"
outPath="./clean_fastq"
#sickleBinPath="~/bin/sickle-master"
#spadesBinPath="~/bin/SPAdes-3.7.1-Linux/"
noThread=16
if [[ $customStep =~ "1" ]];then
bash step1_fastqQC.sh $fastqPath $outPath  $noThread
fi

cleanPath="./clean_fastq"
samtoolsThread=8
bowtieThread=16
reffasta="./sample_ref/Ebola_genome1.fa"
bowtie2indexpath="./sample_ref/Ebola_genome"
outputpath="./mpileup_and_ntfreq"
if [[ $customStep =~ "2" ]];then
bash step2_generate_mpileup.sh $cleanPath $samtoolsThread $bowtieThread $reffasta $bowtie2indexpath $outputpath
fi

noThread=8
filepath="./mpileup_and_ntfreq"
if [[ $customStep =~ "3" ]];then
bash step3_mpileup2ntfreq.sh $noThread $filepath
fi

#Usage: perl f_generate_summary_tables.pl -in <ntfreq dir> -ref <ref fasta> [-out output_path ]
tablePath="big_tables"
if [[ $customStep =~ "4" ]];then
perl step4_generate_summary_tables.pl  -in $filepath -ref $reffasta  -out $tablePath
fi

#Usage:	perl filtering_bigtables.pl -in <summary table dir> [options]
#Options:
#	-in <dir>          the path where *ntfreq files are.
# Thresholds:
#	-d <int>           the minimum depth of a genome position. Default = 100
#	-c <int>           the minimum depth of minor allele. Default = 5
#	-s <int>           the minimum number of valid (pass -d and -c) genome positions of a sample. Dafault = 3000
#	-f <fload|0-1>     the minimum minor allele frequency to identify a intrahost single nucleotide variations (iSNVs).If f > 1 - cut, output as SNP. Default = 0.02
#	-st-t <int>        the 4 types to measure stranded bias. 0: total. 1: max allele. 2: minor allele. 3: both max and minor alleles. 4: max or minor alleles. 
#					   Only type 2 available now. Default = 2
#	-st-c <fload|0-1>  the stranded bias ratio cutoff. Default = 0.1
#	-pf <file>         the file of genome region to exclude in iSNV identification, such as the primer regions. Format: stard-end\n... . Default = NULL
#	-pf-d <int>        the upstream and downstream adjcent region of -pf are also excluded, usually for -pf is primer regions.Default = 0
if [[ $customStep =~ "5" ]];then
perl step5_filtering_bigtables.pl -in $tablePath -d 100 -c 5 -s 3000 -f 0.02 -st-t 2 -st-c 0.1 
fi

echo ================ done ====================


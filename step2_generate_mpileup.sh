#!/bin/bash

cleanPath=$1 #"./clean_fastq"
samtoolsThread=$2 #8
bowtieThread=$3 #16
reffasta=$4 #"./ref/Ebola_genome1.fa"
bowtie2indexpath=$5 #"./ref/Ebola_genome"
outputpath=$6 #"./mpileup_and_ntfreq"

#bowtie2binpath=""
#samtoolsbinpath=""

echo 
echo ============== step 2: generate mpileup "(samtools 1.3 required)" ==============
echo

if [ "${cleanpath:0-1}"x = "/"x ]; then
	cleanpath=${cleanpath%?}
fi

echo "clean fastq path: "$cleanPath
echo "reference fasta file: "$reffasta
echo "bowtie2 index: "$bowtie2indexpath
if [ ! -e $bowtie2indexpath.1.bt2 ];then
	echo "!!! cannot find bowtie2 index file"
	exit
fi
echo "output path: "$outputpath
if [ ! -d $outputpath ];then
	echo "  !"$outputpath does not exsit, mkdir
	mkdir $outputpath
fi
echo "no. of threads for bowtie2 (-p): "$bowtieThread
echo "no. of threads for samtools (-@ for sort, no. of task for mpileup): "$samtoolsThread
echo

# ----------------------------------
echo -------------- step 2.1: bowtie2-align and sort bam -------------
for file in $cleanPath/*R1_*fastq
do
	samp=${file%%R1_*}
	samp=${samp##*\/}
	echo ~~~ SampId = $samp ~~~~
#	if [ ! -e $outputpath/$samp.sam ];then
	if [ ! -e $outputpath/$samp.sorted.bam ];then
		bowtie2 -x $bowtie2indexpath -1 $file -2 ${file/R1/R2} -p $bowtieThread -S $outputpath/$samp.sam
		samtools view -bS $outputpath/$samp.sam | samtools sort - -@ $samtoolsThread -o $outputpath/$samp.sorted.bam 
	else
		echo $outputpath/$samp.sorted.bam already exists, skip bowtie2, samtools view and samtools sort. To perform samtools mpileup.
	fi
	samtools index $outputpath/$samp.sorted.bam
done

echo
echo --------------- step 2.2: generate mpileup file -----------------
for file in $outputpath/*sorted.bam
do
	echo $file
	while [ $(ps -Af|grep "samtools mpileup -f"|wc -l) -gt $samtoolsThread ]
	do
		sleep 2
	done
	if [ ! -e ${file/sorted.bam/mpileup} ];then
		samtools mpileup -f $reffasta -d 1000000 $file  > ${file/sorted.bam/mpileup} &
	else
		echo ${file/sorted.bam/mpileup} already exists, skip samtools mpileup.
	fi
done
	
while [ $(ps -Af|grep "samtools mpileup -f"|wc -l) -gt 1 ]; do
	sleep 2
done

echo
if [ -e $outputpath/*sam ]; then
	echo delete sam files
	rm $outputpath/*sam
fi

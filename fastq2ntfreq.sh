#!/bin/bash

# ---------------------------------------------------------------------------
# Load clean data of viral genome amplicon-seq, align them to reference genome,
# generate mpileup file, and obtain the nucleotide frequencies and other
# results (by mpileup2ntFreq.pl). 
# ---------------------------------------------------------------------------

# usage: fastq2ntfreq <read 1 fastq file in sampl*_1.fastq format>

ref=/public1/home/nim/EBOV_iSNV_OurSamp_20150401/ebov_ref/KJ660346.2
reffasta=/public1/home/nim/EBOV_iSNV_OurSamp_20150401/ebov_ref/KJ660346.2.fasta
binpath=/public1/home/nim/bin
outputpath=/public1/home/nim/EBOV_iSNV_OurSamp_20150401/20150526_NewQC_results/bam


samp=${file%%\_1*}
samp=${samp##*\/}

echo ~~~ SampId = $samp ~~~~

$binpath/bowtie2-2.2.5/bowtie2 -x $ref -p 16 -1 $file -2 ${file/_1/_2} -S $outputpath/$samp.sam
$binpath/samtools-1.2/samtools view -bS $outputpath/$samp.sam | $binpath/samtools-1.2/samtools sort - $outputpath/$samp.sorted;
$binpath/samtools-1.2/samtools index $outputpath/$samp.sorted.bam;
$binpath/samtools-1.2/samtools mpileup -f $reffasta -d 1000000 $outputpath/$samp.sorted.bam  > $outputpath/$samp.mpileup;

cat $outputpath/$samp.mpileup | awk '{print $2"\t"$4}' > $outputpath/$samp.pos_dep # just for visualization

perl $binpath/iSNV_script/run_mpileup2ntFreq_addStranded.pl $outputpath/$samp.mpileup > ${outputpath/bam/ntfreq}/$samp.mpileup.ntfreq

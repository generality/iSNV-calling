#!/bin/bash

#amplicon_qc.sh s1 /public1/home/yefq/data/test/rawdata /public1/home/yefq/data/test/sickleoutput /public1/home/yefq/data/test/errcorroutput /public1/home/yefq/data/test/assemblyoutput

sampleid=$1
sourcedir=$2
trimdatadir=$3
sickleoutputdir=$4
errcorroutputdir=$5
assemblyoutputdir=$6

forward=$sampleid"_1.fastq"
reverse=$sampleid"_2.fastq"

echo "Sample ID is "$sampleid

mkdir -p $trimdatadir/$sampleid
mkdir -p $sickleoutputdir/$sampleid
mkdir -p $errcorroutputdir/$sampleid
mkdir -p $assemblyoutputdir/$sampleid

#Step0 trim the 10 base at 5'
perl trim_fastq_first_N_base.pl $sourcedir/$forward $trimdatadir/$sampleid
perl trim_fastq_first_N_base.pl $sourcedir/$reverse $trimdatadir/$sampleid


#Step1 run sickle for adaptive triming
sickle pe -f $trimdatadir/$sampleid/$forward -r $trimdatadir/$sampleid/$reverse -o $sickleoutputdir/$sampleid/$forward -p $sickleoutputdir/$sampleid/$reverse -s $sickleoutputdir/$sampleid/trimmed_single.fastq -t sanger -q 20 -l 100

#Step2 run Bayeshammer for error correction
spades.py --only-error-correction --disable-gzip-output -o $errcorroutputdir/$sampleid --pe1-1 $sickleoutputdir/$sampleid/$forward --pe1-2 $sickleoutputdir/$sampleid/$reverse --pe1-s $sickleoutputdir/$sampleid/trimmed_single.fastq -t 16 -m 20


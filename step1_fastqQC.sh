#!/bin/bash

# ----------------------------------- fastq QC: the pipeline for iSNV -----------------------------------------------------------------------------------------------
fastqPath=$1 #"./"
outPath=$2 #"./clean_fastq"
noThread=$3 #16


echo 
echo ============== step 1: fastq qc ==============
echo

if [ "${fastqPath:0-1}"x = "/"x ]; then
	fastqPath=${fastqPath%?}
fi

if [ "${outPath:0-1}"x = "/"x ]; then
	outPath=${outPath%?}
fi

if [ "${sickleBinPath:0-1}"x = "/"x ]; then
	sickleBinPath=${sickleBinPath%?}
fi

echo "Path to fastq file(s): "$fastqPath
echo "Output fastq after QC: "$outPath
if [ ! -d $outPath ]; then
	echo '  !' output path does not exist, mkdir
	mkdir $outPath
fi
echo "number of threads: "$noThread

# trim frist 10 base for all fastq, both *fastq and *fastq.gz cannot be recognized
echo
echo ------------------ step 1.1: trim the first 10 bases ------------
count=0
for file in $fastqPath/*fastq
do
	while [ $(ps -Af|grep "awk NR%4==2"|wc -l) -gt $noThread ]
	do
		sleep 2
	done
	if [ -e $file -a ! -d $file ]; then
		let n=n+1
		samp=${file##*\/}
		samp=${samp%%\.*}
		echo "  "$file, sampId = $samp
		cat  $file | awk 'NR%4==2||NR%4==0{print substr($0,11);} NR%4==1||NR%4==3{print}' > $outPath"/"Atrim10.$samp.fastq &
	fi
done

for file in $fastqPath/*fastq.gz
do
	while [ $(ps -Af|grep "awk NR%4==2"|wc -l) -gt $noThread ]
	do
#		echo  $(ps -Af|grep "awk NR%4==2"|wc -l)
#		ps -Af|grep "awk NR%4==2"
		sleep 2
	done
	if [ -e $file ]; then
		let n=n+1
		samp=${file##*\/}
		samp=${samp%%\.*}
		echo "  "$file, sampId = $samp
		gzip -dc $file | awk 'NR%4==2||NR%4==0{print substr($0,11);} NR%4==1||NR%4==3{print}' > $outPath"/"Atrim10.$samp.fastq &
	fi
done

while [ $(ps -Af|grep "awk NR%4==2"|wc -l) -gt 1 ]
do
	sleep 5
done


echo
echo ----------------- step 1.2: run sickle -------------------------- 
echo "  " ! discriminate read 1 and read 2 by \"R1_\"
cd $outPath


for file in Atrim10*R1_*fastq
do
	while [[ $(ps -Af|grep "sickle pe"|wc -l) -gt $noThread ]]
	do
		sleep 2
	done
#	echo $file
	sickle pe -f $file -r ${file/R1/R2} -o Bsickle.$file -p Bsickle.${file/R1/R2} -s Bsickle.${file/R1/se}  -t sanger -q 20 -l 100 | tee sickleStdout.txt &
done

while [[ $(ps -Af|grep "sickle pe"|wc -l) -gt 1 ]]
do
	sleep 5
done

echo ---------------- step 1.3: run spades ----------------------

if [ -e spadesStdout.txt ]; then
	rm spadesStdout.txt
fi

for file in Bsickle*R1_*fastq
do 
	echo "  "$file, ${file/R1/R2}
	echo spades.py  --only-error-correction --disable-gzip-output -o ./  --pe1-1 $file  --pe1-2 ${file/R1/R2}  --pe1-s ${file/R1/se}  -t $noThread -m 20
	time spades.py  --only-error-correction --disable-gzip-output -o ./  --pe1-1 $file  --pe1-2 ${file/R1/R2}  --pe1-s ${file/R1/se}  -t $noThread -m 20 >>spadesStdout.txt
	for cfile in corrected/*cor.fastq
	do
		mv $cfile ./${cfile##*Atrim10.}
	done
done

echo
echo clean output tmp files
rm Atrim10*fastq
rm Bsickle*fastq
rm -r corrected
rm -r tmp
rm input_dataset.yaml
rm spades.log
rm params.txt

cd -

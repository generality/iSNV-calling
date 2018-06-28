#!/bin/bash

noThread=$1 #8
filepath=$2 #./mpileup_and_ntfreq

echo
echo ============== step 3: process mpileup file ===============
echo

if [ "${filepath:0-1}"x = "/"x ]; then
	filepath=${filepath%?}
fi

for file in $filepath/*mpileup
do
	echo $file
	while [ $(ps -Af|grep "perl f_mpileup2ntFreq_step3.pl"|wc -l) -gt $noThread ]
	do
		sleep 2
	done
	echo perl f_mpileup2ntFreq_step3.pl -in $file
	echo
	time perl f_mpileup2ntFreq_step3.pl -in $file &
done


while [ $(ps -Af|grep "perl f_mpileup2ntFreq_step3.pl"|wc -l) -gt 1 ]; do
	sleep 2
done

# The batch file to mapping libaries of YFV amplicon-seq fastq data.

use warnings;
use strict;

my @id = qw/
MiLib052_S1
MiLib053_S2
MiLib054_S3
MiLib055_S4
MiLib056_S5
MiLib057_S6
MiLib058_S7
MiLib059_S8
MiLib060_S9
MiLib076_S1
MiLib077_S2
MiLib078_S3
MiLib079_S4
MiLib080_S5
MiLib081_S6
MiLib082_S7
MiLib083_S8
MiLib084_S9
MiLib085_S10
MiLib086_S11
MiLib087_S12
MiLib092_S5
/;

my $suf1 = '_L001_R1_001.fastq.gz';
my $suf2 = '_L001_R2_001.fastq.gz';
my $ref = '/home/nim/serial_iSNV_phasing_YFV/yf1';
foreach (@id){
#	my $fq1 = '/data/ebola/big_tables/EBOV_135samp_Clean_Reads/'.$_.$suf1;
#	my $fq2 = '/data/ebola/big_tables/EBOV_135samp_Clean_Reads/'.$_.$suf2;
	my $fq1 = '/data/ebola/exchange/YFV/0628/'.$_.$suf1;
	my $fq2 = '/data/ebola/exchange/YFV/0628/'.$_.$suf2;

	print "bowtie2 -x $ref -p 4 -1 $fq1 -2 $fq2 -S $_.sam\n";
	system "bowtie2 -x $ref -p 4 -1 $fq1 -2 $fq2 -S $_.sam\n";
}

#!/usr/bin/perl -w
use warnings;
use strict;

#-------------------------------------
# Trim the first N base of a read
#-------------------------------------

#perl trim_fastq_first_N_base.pl raw_data/shi1_1.fastq raw_data_trimmed

my $in=$ARGV[0];
my $dir=$ARGV[1];

my @eles=split/\//,$in;
my $file_name=$eles[$#eles];
my $out="$dir/$file_name";

open FP1, "< $in";
open FP2, "> $out";
my $count=0;
while(<FP1>){
	chomp;
	$count++;
	if($count % 4 ==2 || $count % 4 ==0 ){
		my $line=substr($_,10);
		print FP2 "$line\n";
	}
	else {
		print FP2 "$_\n";
	}
}
close FP1;
close FP2;


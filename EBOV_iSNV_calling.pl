use warnings;
use strict;

# -------------------------------------------------------
# The iSNV calling for our amplicon-seq protocol of cDNA 
# extracted from whole blood samples of ebola patients in
# West Africa, 2014.
# -------------------------------------------------------

# =============== The thresholds ============
my $THRES_siteDepCut  = 100; #  >=
my $THRES_sampSiteCut = 3000; # >=
my $THRES_minorAlleDepCut = 5; # >=
my $THRES_iSNVCut = 0.05; # >=
my $THRES_strandedRatioCut = 0.1;

#my $path = '/data/ebola/big_tables/Merged/';
#my $path = '/data/ebola/big_tables/Merged_addNewReads/';
#my $path = '/data/ebola/big_tables/MiSeq_135samp/';
my $path = './'; #'/data/ebola/big_tables/bigtables_miseqAfterQC_withNewReads/';

my $fileDep       = $path."out_cov_noIndel_QC.bigtables.txt";
my $fileMinorDep  = $path."out_minor_cov.bigtables.txt";
my $fileMinorFreq = $path."out_minor_freq.bigtables.txt";
my $fileMaxStrandRatio = $path."out_max_stranded_posRatio.bigtables.txt";
my $fileSecStrandRatio = $path."out_second_stranded_posRatio.bigtables.txt";
#my $fileDate       = $path."238.date.info.sort";
my $fileDate       = $path."192.date.info.sort";

print "# read big tables ...\n";
print "$fileDep\n";
my $rh_siteDep         = read_tables($fileDep);
print "fileMinorDep\n";
my $rh_siteMinorDep    = read_tables($fileMinorDep);
print "fileMinorFreq\n";
my $rh_siteMinorFreq   = read_tables($fileMinorFreq);
print "fileMaxStrandRatio\n";
my $rh_siteMaxStranded = read_tables($fileMaxStrandRatio);
print "fileSecStrandRatio\n";
my $rh_siteSecStranded = read_tables($fileSecStrandRatio);
print "read file date from $fileDate\n";
my ($rh_sample2mon, $rh_sample2date) = read_sample_date($fileDate);
print "load ampliconseq primer regions\n";
my $rh_pos2isprimer_20 = load_primer_region_20p();
my $rh_pos2isprimer_39 = load_primer_region_39p();
print "load sample ids that using 39 pairs of ampliconseq primers\n";
my $rh_sample39p = load_samples_39primers();
print "\n";


my @samples = sort keys %{$rh_siteDep};
my @positions = sort {$a<=>$b} keys %{$rh_siteDep->{$samples[0]}};

my $rh_for_iSNV;
my $rh_iSNV;

my $cSampValid = 0;
open OUT, ">samples.statistics.txt" or die "$!: samples.staitstics.txt\n";
my $rh_samplesToRemove;
printf "# filtering and output ...\n\n";
print "#sampleId(formatted)\tsampleId(originId)\tvalidPosNo(Dep>=$THRES_siteDepCut|MAF>=$THRES_minorAlleDepCut)\tiSNV_No\tmeanSiteDep(noIndelCov&onlyCountValidpos)\tmonth\tdate\tIsSampleIncluded\tExcluded_iSNV_no_within_primerRegion\tprimerSet(20pair|39pair)\n";
print OUT "#sampleId(formatted)\tsampleId(originId)\tvalidPosNo(Dep>=$THRES_siteDepCut|MAF>=$THRES_minorAlleDepCut)\tiSNV_No\tmeanSiteDep(noIndelCov&onlyCountValidpos)\tmonth\tdate\tIsSampleIncluded\tExcluded_iSNV_no_within_primerRegion\tprimerSet(20pair|39pair)\n";
my $cOurSample = 0;
my $cDownloadSample = 0;
for (my $i = 0; $i < @samples; $i ++){
	my $cValidPos   = 0;
	my $cValid_iSNV = 0;
	my $siteMeanDep = 0;
	my $ciSNVtoExcludeInPrimer = 0;
	my $formatsample = format_sample_name($samples[$i]);
	foreach my $p (@positions){
		if ( ($rh_siteDep->{$samples[$i]}->{$p} ne 'na' && $rh_siteDep->{$samples[$i]}->{$p} >= $THRES_siteDepCut)  || ($rh_siteMinorDep->{$samples[$i]}->{$p} ne 'na' && $rh_siteMinorDep->{$samples[$i]}->{$p} >= $THRES_minorAlleDepCut) ) {
			$cValidPos ++;
			$siteMeanDep += $rh_siteDep->{$samples[$i]}->{$p};
			$rh_for_iSNV->{$samples[$i]}->{$p} = int($rh_siteMinorFreq->{$samples[$i]}->{$p}*10000)/10000;
 			my $flag_iSNVpassStrandedFilter =  $rh_siteSecStranded->{$samples[$i]}->{$p} ne 'na' &&  $rh_siteSecStranded->{$samples[$i]}->{$p} >= $THRES_strandedRatioCut && $rh_siteSecStranded->{$samples[$i]}->{$p} <= 1-$THRES_strandedRatioCut;
			if ($flag_iSNVpassStrandedFilter == 1 && $rh_siteMinorFreq->{$samples[$i]}->{$p} >= $THRES_iSNVCut){
				if (defined $rh_sample39p->{$formatsample}){
					if (! defined $rh_pos2isprimer_39->{$p}){
						$cValid_iSNV ++;
						$rh_iSNV->{$samples[$i]}->{$p} = int($rh_siteMinorFreq->{$samples[$i]}->{$p}*10000)/10000;				
					} else {
						$ciSNVtoExcludeInPrimer ++;
					}
				} else {
					if (! defined $rh_pos2isprimer_20->{$p}){
						$cValid_iSNV ++;
						$rh_iSNV->{$samples[$i]}->{$p} = int($rh_siteMinorFreq->{$samples[$i]}->{$p}*10000)/10000;				
					} else {
						$ciSNVtoExcludeInPrimer ++;
					}
				}
#				$cValid_iSNV ++;
#				$rh_iSNV->{$samples[$i]}->{$p} = int($rh_siteMinorFreq->{$samples[$i]}->{$p}*10000)/10000;
				#print "$rh_siteSecStranded->{$samples[$i]}->{$p}\t$flag_iSNVpassStrandedFilter\n";
			} 
 			#my $flag_iSNVpassStrandedFilter =  $rh_siteSecStranded->{$samples[$i]}->{$p} ne 'na' &&  $rh_siteSecStranded->{$samples[$i]}->{$p} >= $THRES_strandedRatioCut && $rh_siteSecStranded->{$samples[$i]}->{$p} <= 1-$THRES_strandedRatioCut) &&  ($rh_siteMaxStranded->{$samples[$i]}->{$p} ne 'na' && $rh_siteMaxStranded->{$samples[$i]}->{$p} >= $THRES_strandedRatioCut && $rh_siteMaxStranded->{$samples[$i]}->{$p} <= 1-$THRES_strandedRatioCut) ) 
		}
	}
	if ($cValidPos != 0){
		$siteMeanDep = int(($siteMeanDep/$cValidPos)*10)/10;
	} else {
		$siteMeanDep = 0;
	}
	print "$formatsample\t$samples[$i]\t$cValidPos\t$cValid_iSNV\t$siteMeanDep\t";
	print OUT "$formatsample\t$samples[$i]\t$cValidPos\t$cValid_iSNV\t$siteMeanDep\t";
	if (defined $rh_sample2mon->{$formatsample} && defined $rh_sample2date->{$formatsample}){
		print $rh_sample2mon->{$formatsample}."\t".$rh_sample2date->{$formatsample}."\t";
		print OUT $rh_sample2mon->{$formatsample}."\t".$rh_sample2date->{$formatsample}."\t";
	} else {
		print "-\t-\t";
		print OUT "-\t-\t";
	}
	if ($cValidPos >= $THRES_sampSiteCut){
		print "Yes\t";
		print OUT "Yes\t";
		if ($formatsample =~ /^J/){
			$cOurSample ++;
		} else {
			$cDownloadSample ++;
		}
		$cSampValid ++;
	} else {
		$rh_samplesToRemove->{$samples[$i]} = 1;
		print "No\t";
		print OUT "No\t";
	}
	print "$ciSNVtoExcludeInPrimer\t";
	print OUT "$ciSNVtoExcludeInPrimer\t";
	if (defined $rh_sample39p->{$formatsample}){
		print "39-pair\n";
		print OUT "39-pair\n";
	} else {
		print "20-pair\n";
		print OUT "20-pair\n";
	}
}

print "# Total $cSampValid samples are valid with cuts as: valid_site#>$THRES_sampSiteCut (dep>=$THRES_siteDepCut || minor_dep>=$THRES_minorAlleDepCut); iSNV cut as minor freq > $THRES_iSNVCut && (sec allele positive stranded $THRES_strandedRatioCut to ",1-$THRES_strandedRatioCut,")\n\tOur valid sample No: $cOurSample\n\tDownload vadlid sample No: $cDownloadSample\n";
print OUT "# Total $cSampValid samples are valid with cuts as: valid_site#>$THRES_sampSiteCut (dep>=$THRES_siteDepCut || minor_dep>=$THRES_minorAlleDepCut); iSNV cut as minor freq > $THRES_iSNVCut && (sec allele positive stranded $THRES_strandedRatioCut to ",1-$THRES_strandedRatioCut,")\n\tOur valid sample No: $cOurSample\n\tDownload vadlid sample No: $cDownloadSample\n";
#print OUT "# Total $cSampValid samples are valid with cuts as: valid_site#>$THRES_sampSiteCut (dep>=$THRES_siteDepCut || minor_dep>=$THRES_minorAlleDepCut); iSNV cut as minor freq > $THRES_iSNVCut && (sec allele positive stranded $THRES_strandedRatioCut to ",1-$THRES_strandedRatioCut,")\n";
close OUT;

while (my ($s,$v) = each %$rh_for_iSNV){
	delete $rh_for_iSNV->{$s} if (defined $rh_samplesToRemove->{$s});
	delete $rh_iSNV->{$s} if (defined $rh_samplesToRemove->{$s});
}

print "\n#output datasets for display, minorFreq.all.txt/formatlab ...\n";
open OUT, "> minorFreq.all.txt" or die "$!: minorFreq.all.txt\n";
open OUT2, "> minorFreq.all.formatlab" or die "$!: minorFreq.all.formatlab\n";
#open OUT3, "> iSNV.all.txt" or die "$!: iSNV.all.txt\n";
#open OUT4, "> iSNV.all.formatlab" or die "$!: iSNV.all.formatlab\n";
print OUT "samples/pos\t";
@samples = keys %$rh_for_iSNV;
for (my $i = 0; $i < @positions; $i ++){
	print OUT "$positions[$i]\t";
#	print OUT3 "$positions[$i]\t";
} 
print OUT "\n";
#print OUT3 "\n";

for (my $i = 0; $i < @samples; $i ++){
	my $formatsample = format_sample_name($samples[$i]);
	print OUT "$formatsample\t";
#	print OUT3 "$formatsample\t";
	for (my $j = 0; $j < @positions; $j ++){
		my $p = $positions[$j];
		if (defined $rh_for_iSNV->{$samples[$i]}->{$p} ){
			print OUT $rh_for_iSNV->{$samples[$i]}->{$p}."\t";
			print OUT2 $rh_for_iSNV->{$samples[$i]}->{$p}."\t";
		} else {
			print OUT "-1.0\t";
			print OUT2 "-1.0\t";
		}
		
#		if (defined $rh_iSNV->{$samples[$i]}->{$p} ){
#			print OUT3 $rh_iSNV->{$samples[$i]}->{$p}."\t";
#			print OUT4 $rh_iSNV->{$samples[$i]}->{$p}."\t";
#		} else {
#			print OUT3 "-1.0\t";
#			print OUT4 "-1.0\t";
#		}		
	} 
	print OUT "\n";
	print OUT2 "\n";
#	print OUT3 "\n";
#	print OUT4 "\n";
}
close OUT;
close OUT2;
#close OUT3;
#close OUT4;



print "\n#output only iSNV dataset for display, iSNV.all.txt/formatlab ...\n";
open OUT, "> iSNV.all.txt" or die "$!: iSNV.all.txt\n";
open OUT2, "> iSNV.all.formatlab" or die "$!: iSNV.all.formatlab\n";
print OUT "samples/pos\t";
@samples = keys %$rh_iSNV;
for (my $i = 0; $i < @positions; $i ++){
	print OUT "$positions[$i]\t";
} print OUT "\n";

for (my $i = 0; $i < @samples; $i ++){
	my $formatsample = format_sample_name($samples[$i]);
	print OUT "$formatsample\t";
	for (my $j = 0; $j < @positions; $j ++){
		my $p = $positions[$j];
		if (defined $rh_iSNV->{$samples[$i]}->{$p} ){
			print OUT $rh_iSNV->{$samples[$i]}->{$p}."\t";
			print OUT2 $rh_iSNV->{$samples[$i]}->{$p}."\t";
		} else {
			print OUT "-1.0\t";
			print OUT2 "-1.0\t";
		}		
	} 
	print OUT "\n";
	print OUT2 "\n";
}
close OUT;
close OUT2;


#-----------------------------------------------------------------

sub read_tables{
	open FP, $_[0] or die "$!: $_[0]\n";
	my $rh_dat;
	my @samples = ();
	while (<FP>){
		chomp;
		if(/^#/){
			@samples = split /\s+/,$_;
		} elsif (/^\d/) {
			my @row = split /\s+/, $_;
			for(my $i = 1; $i < @row; $i++){
				$rh_dat->{$samples[$i]}->{$row[0]} = $row[$i]; 
			}
		}
	}
	close FP;
	return $rh_dat;
}


sub format_sample_name{
	my $s = $_[0];
	if ($s =~ /^\d/){
		$s =~ s/i//g;
		$s = '1646' if ($s eq '1046');
		if (length($s) == 1 ){
			return 'J000'.$s;
		} elsif (length($s) == 2) {
			return 'J00'.$s;
		} elsif (length($s) == 3) {
			return 'J0'.$s;
		} elsif (length($s) == 4) {
			return 'J'.$s;
		}
	} else {
		return $s;
	}
}

sub read_sample_date{
	open FP, $_[0] or die "$!: $_[0]\n";
	my $rh_month;
	my $rh_date;
	while (<FP>){
		chomp;
		next if (/^Sample/ || length($_) < 5);
		my @row = split /\t/, $_;
		my $samp = $row[0];
		my $mon  = $row[4];
		my $date = lc($row[6]);
		$date =~ s/^\s+//g;
		$date =~ s/\s+$//g;
		$rh_month->{$samp} = $mon;
		if ($date =~ /^no$/ || ($date =~ /^\d+$/ && length($date) == 8)){
			$rh_date->{$samp} = $date;
		} elsif ($date=~ /(\d+)\-(\w+)\-14/){
			my $d = $1;
			my $m = $2;
			my $y = '2014';
			$d = '0'.$d if (length($d) == 1);
			if ($m eq 'jun') {
				$m = '06';
			} elsif ($m eq 'jul'){
				$m = '07';
			} elsif ($m eq 'aug'){
				$m = '08';
			} elsif ($m eq 'sep'){
				$m = '09';
			} elsif ($m eq 'oct'){
				$m = '10';
			} elsif ($m eq 'nov'){
				$m = '11';
			} elsif ($m eq 'dec'){
				$m = '12';
			} else {
				$m = 'unknown';
				#die "cannot reg month format in $date\n$_\n";
			}
			$date = $y.$m.$d;
			$rh_date->{$samp} = $date;
		} elsif ($date =~ /2014\-(\d+)\-([\w\d]+)/){
			$date = "2014".$1.$2;
			$rh_date->{$samp} = $date;
		} elsif ($date =~ /2014\/(\d+)\/(\d+)/){
			my $m = $1;
			my $d = $2;
			$m = '0'.$m if (length($m)==1);
			$d = '0'.$d if (length($d)==1);
			$date = "2014".$1.$2;
			$rh_date->{$samp} = $date;
		}else { 
			$rh_date->{$samp} = 'no';
			#die "cannot reg date format in $date\n$_\n";
		}
	}	
	return ($rh_month, $rh_date);
}

# sample id list from XuPeisong, WeiXin Pic. 
# * 146 are next to 145, in the PCR box
# * All in the first batch. The 2nd batch, all were 20-pair primers
sub load_samples_39primers{
	my @samples = qw/J2065
J0146
J2011
J2040
J2007
J0100
J2021
J0120
J2053
J2048
J2062
J0140
J2054
J0034
J2014
J0037
J0066
J0038
J2019
J2008
J2052
J0111
J2038
J0163
J2029
J2012
J0129
J2064
J2051
J2045
J2009
J2034
	/;
	my $rh_samp;
	foreach (@samples){
		$rh_samp->{$_} = 1;
	}
	return $rh_samp;
}

sub load_primer_region_20p{
	my $rh_pos2isprimer;
	#[nim@localhost test]$ cat ampliseq_primerSet_20p.toRef.sam | awk '!/^@/{cigar=$6;gsub("M","",cigar);  print $4"-"$4+cigar-1}'
	my @regions = qw(130-151 1407-1426 1231-1252 2449-2469 2255-2274 3431-3452 3213-3233 4239-4260 4078-4098 5161-5180 4958-4977 6011-6031 5911-5930 6965-6984 6823-6842 7997-8019 7750-7771 9028-9046 8807-8828 9759-9782 9576-9596 10670-10692 10470-10490 11607-11626 11290-11309 12273-12292 12205-12225 12961-12980 12772-12796 13876-13898 13714-13733 14821-14842 14653-14671 15823-15842 15721-15741 16867-16887 16738-16756 17846-17865 17634-17653 18851-18869);
	foreach (@regions) {
		if (/(\d+)\-(\d+)/){
			for (my $i = $1; $i <= $2; $i++){
				$rh_pos2isprimer->{$i} = 1;
			}
		}
	}
	return $rh_pos2isprimer;
}

sub load_primer_region_39p{
	my $rh_pos2isprimer;
	#[nim@localhost test]$ cat ampliseq_primerSet_39p.toRef.sam | awk '!/^@/{cigar=$6;gsub("M","",cigar);  print $4"-"$4+cigar-1}'
	my @regions = qw/41-59
586-606
468-487
1107-1126
949-968
1595-1613
1551-1573
2246-2265
2107-2130
2742-2763
2595-2615
3119-3137
2975-2995
3525-3546
3445-3466
4124-4145
3989-4010
4562-4583
4495-4516
5115-5134
4976-4994
5744-5765
5579-5603
6250-6269
6198-6217
6782-6801
6662-6683
7269-7290
7190-7212
7772-7794
7664-7682
8181-8200
8030-8048
8672-8690
8552-8573
9189-9210
9081-9099
9602-9620
9544-9566
10136-10155
10031-10050
10565-10585
10517-10538
11161-11180
11042-11060
11658-11676
11589-11610
12200-12221
12122-12141
12777-12796
12653-12670
13244-13263
13102-13124
13692-13713
13620-13641
14304-14326
14161-14183
14729-14750
14660-14679
15140-15159
15068-15088
15609-15631
15473-15493
16080-16099
15954-15975
16518-16537
16357-16376
16867-16888
16736-16755
17313-17332
17233-17254
17840-17860
17629-17650
18189-18211
18039-18058
18491-18508
18290-18310
18859-18879/;
	foreach (@regions) {
		if (/(\d+)\-(\d+)/){
			for (my $i = $1; $i <= $2; $i++){
				$rh_pos2isprimer->{$i} = 1;
			}
		}
	}
	return $rh_pos2isprimer;
}

use warnings;
use strict;

# Load summary tables and filter iSNVs and SNPs by various characteristics 
# Output iSNV summary table both in txt and matlab data formats, and 
# also output all characteristics extraced from the summary tables for
# the iSNVs and SNPs.
# SNP is defined as MuAF > freq cutoff, and passing depth and stranded filters
# niming, 20170906

print "\n============ step5: iSNV and SNP calling =============\n\n";
#my $THRES_siteDepCut = 100; #  >=
#my $THRES_sampSiteCut = 3000; # >=
#my $THRES_minorAlleDepCut = 5; # >=
#my $THRES_iSNVCut = 0.02; # >=
#my $THRES_strandedRatioCut = 0.1;
#$path = './'; #'/data/ebola/big_tables/bigtables_miseqAfterQC_withNewReads/';

my ($path, $THRES_siteDepCut, $THRES_minorAlleDepCut, $THRES_sampSiteCut,
		$THRES_iSNVCut,$TYPE_strandedRatio, $THRES_strandedRatioCut,
		$FILE_excludeRegion, $INT_excludeDist) = load_arg(\@ARGV);

if ($path eq 'NULL'){
	help_info();
	die "[error] input folder error.\n";
}

my @bigtables = <$path/*bigtables.txt>;

my $rh_bigtables;
foreach my $file (@bigtables){
	if ($file =~ /out_([\w\d\_]+)\.big/){
		print "Load summary table from $file.\n";
		$rh_bigtables->{$1} = read_tables($file);
	}
}
my $rh_siteDep         = $rh_bigtables->{cov_noIndel_QC};
my $rh_secDep          = $rh_bigtables->{second_cov}; # second means second, minor means sum of second, third (if any) ...
my $rh_siteMinorDep    = $rh_bigtables->{minor_cov};
my $rh_siteMinorFreq   = $rh_bigtables->{minor_freq};
my $rh_siteMaxStranded = $rh_bigtables->{max_stranded_posRatio};
my $rh_siteSecStranded = $rh_bigtables->{second_stranded_posRatio};


my $rh_pos2isprimer;
if ($FILE_excludeRegion ne 'NULL'){
	my ($rh_pos2isprimer,$excludeSize) = load_exclude_region($FILE_excludeRegion, $INT_excludeDist);
	print "\nLoad excluding region from $FILE_excludeRegion, with up/downstream $INT_excludeDist bp.\n  Total region size to exclude: $excludeSize\n";
}

print "\n";


my @samples = sort keys %{$rh_siteDep};
my @positions = sort {$a<=>$b} keys %{$rh_siteDep->{$samples[0]}};

my $rh_iSNV;
my $rh_SNP;
my $rh_validPos;
my $rh_posHasiSNV;

my $cSampValid = 0;
my $rh_samplesToRemove;
open OUT, ">samples.statistics.txt" or die "$!: samples.staitstics.txt\n";
printf "# filtering and output ...\n\n";
printf("%-15s\t%-10s\t%-10s\t%-10s\t%-10s\t%-10s\t", "sampleId","validPosNo","#iSNV(final)","#SNP","meanSiteDep","IsSampleValid?");
print "#iSNV_in_excluding_region" if ($FILE_excludeRegion ne 'NULL');
print "\n";

for (my $i = 0; $i < @samples; $i ++){
	my $cValidPos   = 0;
	my $cValid_iSNV = 0;
	my $cValid_SNP  = 0;
	my $siteMeanDep = 0;
	my $ciSNVtoExcludeInPrimer = 0;
	foreach my $p (@positions){
		if ( ($rh_siteDep->{$samples[$i]}->{$p} ne 'na' && $rh_siteDep->{$samples[$i]}->{$p} >= $THRES_siteDepCut)  ||
		($rh_siteMinorDep->{$samples[$i]}->{$p} ne 'na' && $rh_siteMinorDep->{$samples[$i]}->{$p} >= $THRES_minorAlleDepCut) ) {
			$cValidPos ++;
			$rh_validPos->{$samples[$i]}->{$p} = 1;
			$siteMeanDep += $rh_siteDep->{$samples[$i]}->{$p};
 			my $flag_iSNVpassStrandedFilter =  $rh_siteSecStranded->{$samples[$i]}->{$p} ne 'na' &&
				                               $rh_siteSecStranded->{$samples[$i]}->{$p} >= $THRES_strandedRatioCut &&
											   $rh_siteSecStranded->{$samples[$i]}->{$p} <= 1-$THRES_strandedRatioCut;
			if ($flag_iSNVpassStrandedFilter == 1 && $rh_siteMinorFreq->{$samples[$i]}->{$p} >= $THRES_iSNVCut){
					if (! defined $rh_pos2isprimer->{$p}){
						$cValid_iSNV ++;
						if (defined $rh_posHasiSNV->{$p}){
							$rh_posHasiSNV->{$p} ++;
						} else {
							$rh_posHasiSNV->{$p} = 1;
						}
						$rh_iSNV->{$samples[$i]}->{$p} = int($rh_secDep->{$samples[$i]}->{$p}/$rh_siteDep->{$samples[$i]}->{$p}*10000)/10000; 
														#int($rh_siteMinorFreq->{$samples[$i]}->{$p}*10000)/10000;
					} else {
						$ciSNVtoExcludeInPrimer ++;
					}
			}  elsif ($rh_siteMinorFreq->{$samples[$i]}->{$p} < $THRES_iSNVCut) {
				my @patt = split //, $rh_bigtables->{ref_max_sec}->{$samples[$i]}->{$p};
				if ($patt[0] ne $patt[2] && 
				    $rh_siteMaxStranded->{$samples[$i]}->{$p} >= $THRES_strandedRatioCut &&
                    $rh_siteMaxStranded->{$samples[$i]}->{$p} <= 1-$THRES_strandedRatioCut &&
					! defined $rh_pos2isprimer->{$p}){
						#$rh_SNP->{$samples[$i]}->{$p} = int((1 - $rh_siteMinorFreq->{$samples[$i]}->{$p})*10000)/10000;
						$rh_SNP->{$samples[$i]}->{$p} = int(($rh_siteMinorFreq->{$samples[$i]}->{$p})*10000)/10000;
						$cValid_SNP ++;
				}
			}
		}
	}
	if ($cValidPos != 0){
		$siteMeanDep = int(($siteMeanDep/$cValidPos)*10)/10;
	} else {
		$siteMeanDep = 0;
	}
	printf("%-15s\t%-10d\t%-10d\t%-10d\t%-10d\t",$samples[$i], $cValidPos, $cValid_iSNV, $cValid_SNP, $siteMeanDep);
	if ($cValidPos >= $THRES_sampSiteCut){
		print "Yes       \t";
		$cSampValid ++;
	} else {
		$rh_samplesToRemove->{$samples[$i]} = 1;
		while (my ($tmpP, $tmp) = each %{$rh_iSNV->{$samples[$i]}}) {
			$rh_posHasiSNV->{$tmpP} --;
			delete $rh_posHasiSNV->{$tmpP} if ($rh_posHasiSNV->{$tmpP} == 0);
		}
		print "No        \t";
	}
	print "$ciSNVtoExcludeInPrimer" if ($FILE_excludeRegion ne "NULL");
	print "\n";
}

print "\nTotally, $cSampValid samples are valid.\n";
print "\tvalid_site: dep >= $THRES_siteDepCut (ignoring indels) or minor_dep>=$THRES_minorAlleDepCut\n";
print "\t#valid_site >= $THRES_sampSiteCut\n";
print "iSNV identification cutoff:\n";
print "\tminor freq >= $THRES_iSNVCut\n";
print "\tsecond allele positive stranded belongs to [$THRES_strandedRatioCut to ",1-$THRES_strandedRatioCut,"]\n";
print "SNP identification cutoff:\n";
print "\tminor freq < $THRES_iSNVCut\n";
print "\tmajor allele positive stranded belongs to [$THRES_strandedRatioCut to ",1-$THRES_strandedRatioCut,"]\n";

close OUT;

while (my ($s,$v) = each %$rh_iSNV){
	delete $rh_iSNV->{$s} if (defined $rh_samplesToRemove->{$s});
}



print "\n#output summary of iSNVs to iSNV.all.txt and iSNV.all.formatlab ...\n";
open OUT, "> iSNV.all.txt" or die "$!: iSNV.all.txt\n";
open OUT2, "> iSNV.all.formatlab" or die "$!: iSNV.all.formatlab\n";

print OUT "#NOTE: NO = no isnv or freq < cutoff, denoted as -1.0 in iSNV.all.formatlab; NA = this pos is invalid for iSNV calling, denoted as -2.0 in iSNV.all.formatlab\n";
print OUT "pos\tiSNV#\t";
@samples = sort keys %$rh_iSNV;
my @iSNVPos = sort {$a<=>$b} keys %$rh_posHasiSNV;
for (my $i = 0; $i < @samples; $i ++){
	print OUT "$samples[$i]\t";
} print OUT "\n";


for (my $j = 0; $j < @iSNVPos; $j ++){
	my $p = $iSNVPos[$j];
	print OUT "$p\t$rh_posHasiSNV->{$p}\t";
	for (my $i = 0; $i < @samples; $i ++){
		if (defined $rh_iSNV->{$samples[$i]}->{$p} ){
			print OUT $rh_iSNV->{$samples[$i]}->{$p}."\t";
			print OUT2 $rh_iSNV->{$samples[$i]}->{$p}."\t";
		} elsif (defined $rh_validPos->{$samples[$i]}->{$p}) {
			print OUT "NO\t";
			print OUT2 "-1.0\t";
		} else {
			print OUT "NA\t";
			print OUT2 "-2.0\t";
		}
	} 
	print OUT "\n";
	print OUT2 "\n";
}
close OUT;
close OUT2;


print "\nOutput all iSNV info to iSNV_info.txt\n";
print "Output all SNP info to SNP_info.txt\n\n";
list_all_info($rh_iSNV, $rh_bigtables, 'iSNV_info.txt');
list_all_info($rh_SNP, $rh_bigtables, 'SNP_info.txt');

# ==========================================================================================
sub load_arg {
	my @args   = @{$_[0]};
	my $input  = 'NULL';
	my $d = 100;
	my $c = 5;
	my $s = 3000;
	my $f = 0.02;
	my $st_t = 2;
	my $st_c = 0.1;
	my $pf = 'NULL';
	my $pf_d = 0;
	if (@args < 2){
		help_info();
		die "[error] invalid arg.\n";
	}
	for(my $i = 0; $i < @args; $i+=2){
		if ($args[$i] eq '-in'){
			if ( !(-d $args[$i+1])){
				help_info();
				die "\n[error] the folder \"$args[$i+1]\" does not exist.\n";
			} else {
				$input  = $args[$i+1];
				if ($input =~ /\/$/){
					$input =~ s/\/$//;
				}
			}
		} elsif ($args[$i] eq '-d'){
			$d = $args[$i+1];
			if (! $d =~ /^\d+$/){
				help_info();
				die "[error] -d $d: -d requires <int>\n";
			}
		} elsif ($args[$i] eq '-c'){
			$c = $args[$i+1];
			if (! $c =~ /^\d+$/){
				help_info();
				die "[error] -c $c: -c requires <int>\n";
			}
		} elsif ($args[$i] eq '-s'){
			$s = $args[$i+1];
			if (! $s =~ /^\d+$/){
				help_info();
				die "[error] -s $s: -s requires <int>\n";
			}
		} elsif ($args[$i] eq '-f'){
			$f = $args[$i+1];
			if (! $f =~ /^[\d\.]+$/ || $f < 0 || $f > 1){
				help_info();
				die "[error] -f $args[$i+1]: -d requires <float> belonging to [0, 1]\n";
			}
		} elsif ($args[$i] eq '-st-t'){
			$st_t = $args[$i+1];
			if (! $st_t =~ /^\d+$/ || !($st_t == 0 || $st_t == 1 || $st_t == 2 || $st_t == 3 || $st_t == 4) ) {
				help_info();
				die "[error] -st-t $st_t: -d requires <int> = 0, 1, 2, 3, or 4.\n";
			}
		} elsif ($args[$i] eq '-st-c'){
			$st_c = $args[$i+1];
			if (! $st_c =~ /^[\d\.]+$/ || $st_c < 0 || $st_c > 1){
				help_info();
				die "[error] -st_c $args[$i+1]: -d requires <float> belonging to [0, 1]\n";
			}
		} elsif ($args[$i] eq '-pf'){
			$pf = $args[$i+1];
			if (! -e $pf){
				help_info();
				die "[error] -pf $pf: $pf does not exist\n";
			}
		} elsif ($args[$i] eq '-pf-d'){
			$pf_d = $args[$i+1];
			if (! $pf_d =~ /^\d+$/){
				help_info();
				die "[error] -pf-d $pf_d: -pf-d requires <int>\n";
			}
		} else {
			help_info();
			die "[error]: cannot recognize option $args[$i].\n";
		}
	}
	return ($input, $d, $c, $s, $f, $st_t, $st_c, $pf, $pf_d);
}


sub help_info{
	print "\n*** filter summary table to identify iSNVs ***\n";
	print "\nUsage:\tperl filtering_bigtables.pl -in <summary table dir> [options]\n";#-d pos_dep_cut -s valid_pos_no_cut -c minor_allele_dep_cut -f  iSNV_cut -st strand_bias_cut ]\n";
	print "\nOptions:\n";
	print "\t-in <dir>          the path where *ntfreq files are.\n";
	print "\n  Thresholds:\n";
	print "\t-d <int>           the minimum depth of a genome position. Default = 100\n";
	print "\t-c <int>           the minimum depth of minor allele. Default = 5\n";
	print "\t-s <int>           the minimum number of valid (pass -d and -c) genome positions of a sample. Dafault = 3000\n";
	print "\t-f <fload|0-1>     the minimum minor allele frequency to identify a intrahost single nucleotide variations (iSNVs).If f > 1 - cut, output as SNP. Default = 0.02\n";
	print "\t-st-t <int>        the 4 types to measure stranded bias. 0: total. 1: max allele. 2: minor allele. 3: both max and minor alleles. 4: max or minor alleles. Default = 2\n";
	print "\t-st-c <fload|0-1>  the stranded bias ratio cutoff. Default = 0.1\n";
	print "\t-pf <file>         the file of genome region to exclude in iSNV identification, such as the primer regions. Format: stard-end\\n... . Default = NULL\n";
	print "\t-pf-d <int>        the upstream and downstream adjcent region of -pf are also excluded, usually for -pf is primer regions.Default = 0\n";
	print "\n";
}


# ----------------------------------------------------------
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
	my $s = shift;
	return $s;
}


# sample id list from XuPeisong, WeiXin Pic. 
# * 146 are next to 145, in the PCR box
# * All in the first batch. The 2nd batch, all were 20-pair primers

sub load_exclude_region{
	my $file = shift;
	my $d    = shift;
	my $rh_pos2isprime;
	my $size = 0;
	open FP, "$file" or die "[error] -pf $file: $!\n";
	while (<FP>){
		chomp;
		if(/(\d+)\-(\d+)/){
			for (my $i = $1-$d; $i <= $2+$d; $i++){
				$rh_pos2isprimer->{$i} = 1;
				$size ++;
			}
		}
	}
	close FP;
	return ($rh_pos2isprimer, $size);
}

sub list_all_info{
	my ($rh_isnv,$rh_bigtables, $file) = @_;
	my @terms = qw/
	ref_max_sec
	cov_from_mpileup
	cov_noIndel
	cov_noIndel_QC
	ranked_alle
	minor_cov
	minor_freq
	second_cov
	second_freq
	second_stranded
	second_stranded_posRatio
	max_stranded
	max_stranded_posRatio
	entropy
	nt_pattern/;
	open OUT, "> $file" or die "[error] $file: $!\n"; 
	keys %$rh_bigtables;
	print OUT "#samp\tpos\tisnvFreq\tMuAF\t";
	for (my $i = 0; $i< @terms; $i++){
		print OUT "$terms[$i]\t";
	} print OUT "\n";
	while (my ($samp, $rh_dat) = each %$rh_isnv){
		my @p = sort {$a<=>$b} keys %$rh_dat;
		for(my $i = 0; $i < @p; $i ++){
			print OUT "$samp\t$p[$i]\t$rh_dat->{$p[$i]}\t";
			my @patt = split //, $rh_bigtables->{ref_max_sec}->{$samp}->{$p[$i]};

			if ($patt[0] eq $patt[2]){
				print OUT "$rh_dat->{$p[$i]}\t";
			} else {
				print OUT 1-$rh_dat->{$p[$i]},"\t";
			}
			for (my $j = 0; $j< @terms; $j++){
#				print '\'' if ($rh_bigtables->$terms[$j]->{$samp}->{$p[$i]} =~ /\d+\/\d+/);
				print OUT "$rh_bigtables->{$terms[$j]}->{$samp}->{$p[$i]}\t";
			}
			print OUT "\n";
		}
		print OUT "\n";
	}
	close OUT;
}

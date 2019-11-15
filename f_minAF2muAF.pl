# This script transfer the minor allele frequencies in iSNV.all.txt to mutated allele frequencies.
# The mutated allele frequencies are loaded from both iSNV_info.txt and SNP_info.txt.
# Note that the SNP_Info.txt will not contain the loci whose minor allele failed QC such as strand-bias filter.

# Output is the same format as iSNV_info.txt but with iSNV and SNP MuAF. Output to stdout.

# nim, 20191114

use warnings;
use strict;

my $isnvInfo = "iSNV_info.txt";
my $snpInfo  = "SNP_info.txt";
my $isnvall  = "iSNV.all.txt";

my $rh_isnv_muaf = load_info($isnvInfo);
my $rh_snp_muaf  = load_info($snpInfo);

open FP, "$isnvall" or die "$!:$isnvall\n";
my @header;
while(<FP>){
    chomp;
    if (/^\#/){
	print $_."\n";
	next;
    }
    if (/^pos/){
	print $_."\n";
	@header = split /\s+/, $_;
	next;
    }
    my @row = split /\s+/, $_;
    my $pos = $row[0];
    print "$row[0]\t$row[1]\t";
    for(my $i=2; $i<@row; $i++){
	my $v = $row[$i];
	if (defined $rh_isnv_muaf->{$header[$i]}->{$pos}){
	    print $rh_isnv_muaf->{$header[$i]}->{$pos}."\t";
	} elsif (defined $rh_snp_muaf->{$header[$i]}->{$pos}){
	    print $rh_snp_muaf->{$header[$i]}->{$pos}."\t";
	} else {
		print "$v\t"; # should be NA (not available) no NO (no isnv or snp)
        }
    }
    print "\n";
}


sub load_info{
    open FP, shift or die "iSNV/SNP_info.txt: $!\n";
    my @headers;
    my $rh_muaf;
    while(<FP>){
	chomp;
	next if (/^$/);
	if(/^\#/){ # header line
	    @headers = split /\s+/, $_;
	} else {
	    my @row  = split /\s+/, $_;
	    my $samp = $row[0];
	    my $pos  = $row[1];
	    my $minf = $row[2];
	    my $muaf = $row[3];
	    $rh_muaf->{$samp}->{$pos} = $muaf;
	}
    }
    close FP;
    return $rh_muaf;
}

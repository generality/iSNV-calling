use warnings;
use strict;

# ------------------------------------------------------------------------
# Generate site nucleotide frequencies and other values from mpileup files 
# ------------------------------------------------------------------------

die "usage: perl mpileup2ntFreq.pl <mpileup file>\n" if (@ARGV!=1); 
my $file = $ARGV[0];
die "usage: perl mpileup2ntFreq.pl <mpileup file>\nerror: file = $file" if (!($file =~ /mpileup/)); 

my $QTHRES = 20;
my $DEBUG = 0;

open FP, "$file" or die "$!: $file\n";
# output the headers
print "#chr\tpos\tref\tdep(mpileup)\tdep(-indel)\tntfreqPatt\tA\tG\tC\tT\ttot(Q$QTHRES)\tFa\tFg\tFc\tFt\t"."A5\'/3\'\tG5/3\tC5/3\tT5/3\tA5\'(na_if_dep<5)\tG5\tC5\tT5\tentropy\n";
while (<FP>){
	chomp;
	my @row = split /\s+/, $_;
	my $des = $row[0];
	my $ref = $row[2];
	my $pos = $row[1];
	my $dep = $row[3];
	my $seqPatt = $row[4];
	my $qval= $row[5];
	my ($formatSeqPatt, $formatSeqQval,$depNoIndel);
	if (defined $seqPatt){
		($formatSeqPatt,$formatSeqQval) = mpileup_seq_tran($seqPatt, $qval, $dep);
		$depNoIndel = length($formatSeqPatt);
		my ($filSeqPatt, $filqval) = filterNtByQ($formatSeqPatt,$formatSeqQval,$QTHRES);
		my $freStat = fre_stat($filSeqPatt,$ref);
		print "$des\t$pos\t$ref\t$dep\t$depNoIndel\t$freStat\n";#\t$formatSeqPatt\n";
      		if ($DEBUG == 1){
			print "1:original\n$seqPatt\n"; # for debug
			print "2:formatted\n$formatSeqPatt\n"; # for debug
			print "3:filtered\n$filSeqPatt\n\n"; # for debug
       		}
	}
}
close FP;


#delete chars that representing indel
sub mpileup_seq_tran{
	my $pat = $_[0];
	my $qva = $_[1];
	my $dep = $_[2];
	my $originalPat = $pat;
	$pat =~ s/\^.//g;  
	$pat =~ s/\$//g;
	if (length($pat) == length($qva)){
		return ($pat,$qva);
	} else {
		my ($newpat,$newqva) = del_insert_or_deletion($pat,$qva);
#		$pat =~ s/\+\d+[ACGTNacgtn]+//g;
#		$pat =~ s/\-\d+[ACGTNacgtn]+//g;
		if (length($newpat) == length($newqva)){
			return ($newpat,$newqva);
		} else {
			print "[error] dep=",$dep,"\tformatPatLen=",length($newpat),"\tQvalLen=",length($newqva),"\n";
			print "[error] ",$originalPat."\n".$newpat."\n$newqva\n\n";
			die;
		}
	}
	
}

sub del_insert_or_deletion{
	my @sp = split /[\+\-]/,$_[0];
	my @qvals = split //, $_[1];
	my $seq = '';
	my @occupys = (); # the occupying pos for indel to be removed
	push @occupys, length($sp[0])-1; # occupys before the +/-
	for(my $i = 0; $i < @sp; $i++){
		if ($sp[$i] =~ /^(\d+)/){
			my $num = $1;
			$sp[$i] =~ s/^\d+//;
        		my $spnew = substr($sp[$i],$num,1000000); # the this position has a q-value
			$sp[$i] = $spnew;
			my $lastp = @occupys[@occupys-1];
			push @occupys, $lastp+length($spnew)-1;
		}
	}
	for(my $i = 0; $i < @sp; $i++){
		$seq .= $sp[$i];
	}
	my @seqs = split //, $seq;
	if (@seqs == @qvals){
		for (my $i = 0; $i < @occupys-1; $i++){
			$seqs[$occupys[$i]] = '' if (defined $seqs[$occupys[$i]]); # te if is for the last fragment has no occupys, e.g.: ..,+1g,+1g
			$qvals[$occupys[$i]] = '' if (defined $seqs[$occupys[$i]]);
		}
	} else {
		die "in sub del_insert_or_deletion(), length of pvals and seqs are not equal.\n$seq\n$_[1]\n";
	}
	my $newseq = join('',@seqs);
	my $newqval = join('',@qvals);
	return ($newseq, $newqval);
}


sub fre_stat{
	my $seq = $_[0];
	my $ref = $_[1];
	# --------- count postitive and negative strand -----------
	my $cRefPos = $seq =~ tr/\.//;
	my $cRefNeg = $seq =~ tr/\,//;
	my $cAPos = $seq =~ tr/A//;
	my $cANeg = $seq =~ tr/a//;
	my $cGPos = $seq =~ tr/G//;
	my $cGNeg = $seq =~ tr/g//;
	my $cCPos = $seq =~ tr/C//;
	my $cCNeg = $seq =~ tr/c//;
	my $cTPos = $seq =~ tr/T//;
	my $cTNeg = $seq =~ tr/t//;
	my $strandedInfo = '';
	if ($ref eq 'A'){
		$cAPos = $cRefPos; $cANeg = $cRefNeg;
	} elsif ($ref eq 'G'){
		$cGPos = $cRefPos; $cGNeg = $cRefNeg;
	} elsif ($ref eq 'C'){
		$cCPos = $cRefPos; $cCNeg = $cRefNeg;
	} elsif ($ref eq 'T'){
		$cTPos = $cRefPos; $cTNeg = $cRefNeg;
	}
	$strandedInfo .= $cAPos.'/'.$cANeg."\t";
	$strandedInfo .= $cGPos.'/'.$cGNeg."\t";
	$strandedInfo .= $cCPos.'/'.$cCNeg."\t";
	$strandedInfo .= $cTPos.'/'.$cTNeg;
    my $APosRatio = 'na';
    my $GPosRatio = 'na';
    my $CPosRatio = 'na';
    my $TPosRatio = 'na';
    my $demc = 1000;
    $APosRatio = int($cAPos/($cAPos+$cANeg)*$demc)/$demc if ($cAPos + $cANeg >= 5);
    $GPosRatio = int($cGPos/($cGPos+$cGNeg)*$demc)/$demc if ($cGPos + $cGNeg >= 5);
    $CPosRatio = int($cCPos/($cCPos+$cCNeg)*$demc)/$demc if ($cCPos + $cCNeg >= 5);
    $TPosRatio = int($cTPos/($cTPos+$cTNeg)*$demc)/$demc if ($cTPos + $cTNeg >= 5);
    $strandedInfo .= "\t$APosRatio\t$GPosRatio\t$CPosRatio\t$TPosRatio";

	$seq =~ s/\,/\./g;
	$seq =~ s/\./$ref/g;
	$seq = uc($seq);
	my @seq = split //, $seq;
	my $rh_cou;
	foreach (@seq){
		if (defined $rh_cou->{$_}){
			$rh_cou->{$_} ++;
		} else {
			$rh_cou->{$_} = 1;
		}
	}
	my $returnString = '';

	my $tot=0;
	if (defined $rh_cou->{'A'}){
		$returnString .= 'A:'.$rh_cou->{'A'}.';';
		$tot += $rh_cou->{'A'};
	} else {
		$returnString .= 'A:0;';
	}
	if (defined $rh_cou->{'G'}){
		$returnString .= 'G:'.$rh_cou->{'G'}.';';
		$tot += $rh_cou->{'G'};
	} else {
		$returnString .= 'G:0;';
	} 
	if (defined $rh_cou->{'C'}){
		$returnString .= 'C:'.$rh_cou->{'C'}.';';
		$tot += $rh_cou->{'C'};
	} else {
		$returnString .= 'C:0;';
	} 
	if (defined $rh_cou->{'T'}){
		$returnString .= 'T:'.$rh_cou->{'T'}.';';
		$tot += $rh_cou->{'T'};
	} else {
		$returnString .= 'T:0;';
	}
	$returnString .= 'tot:'.$tot;
	$returnString .= "\t$rh_cou->{'A'}" if ( defined $rh_cou->{'A'});
	$returnString .= "\t0"              if (!defined $rh_cou->{'A'});
	$returnString .= "\t$rh_cou->{'G'}" if ( defined $rh_cou->{'G'});
	$returnString .= "\t0"              if (!defined $rh_cou->{'G'});
	$returnString .= "\t$rh_cou->{'C'}" if ( defined $rh_cou->{'C'});
	$returnString .= "\t0"              if (!defined $rh_cou->{'C'});
	$returnString .= "\t$rh_cou->{'T'}" if ( defined $rh_cou->{'T'});
	$returnString .= "\t0"       if (!defined $rh_cou->{'T'});
	$returnString .= "\t$tot";
	my $fa = 0;
	my $fg = 0;
	my $fc = 0;
	my $ft = 0;
	my $entropy = 0;
    $demc = 10000;
	if (defined $rh_cou->{'A'}){
		$fa = int($rh_cou->{'A'}/$tot*$demc)/$demc;
		$entropy += $fa*log($fa) if ($fa > 0);
	}
	if (defined $rh_cou->{'G'}){
		$fg = int($rh_cou->{'G'}/$tot*$demc)/$demc;
		$entropy += $fg*log($fg) if ($fg > 0);
	} 
	if (defined $rh_cou->{'C'}){
		$fc = int($rh_cou->{'C'}/$tot*$demc)/$demc;
		$entropy += $fc*log($fc) if ($fc > 0);
	} 
	if (defined $rh_cou->{'T'}){
		$ft = int($rh_cou->{'T'}/$tot*$demc)/$demc;
		$entropy += $ft*log($ft) if ($ft > 0);
	}
	$entropy *= -1;
	$returnString .= "\t$fa\t$fg\t$fc\t$ft";
	return $returnString."\t".$strandedInfo."\t".$entropy;
}

sub log2{
	return log($_[0])/log(2);
}

sub filterNtByQ{
	my ($seq,$q,$QThres) = @_;
	my @seq = split //,$seq;
	my @q = split //,$q;
	my $filterSeq = '';
	my $filterQ = '';
	return "NULL" if(@seq != @q);
	for (my $i = 0; $i < @seq; $i++){
		if ((ord($q[$i])-33) >= $QThres){
			$filterSeq .= $seq[$i];
			$filterQ  .= $q[$i];
		}
	}
	return ($filterSeq, $filterQ);
}

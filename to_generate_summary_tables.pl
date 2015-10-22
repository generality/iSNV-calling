use warnings;
use strict;

# --------------------------------------------------------------------------
# Integrate the results (frequencies, depth, strand biased, etc produced by 
# fastq2ntfreq.sh) for a batch of samples for the following iSNV calling.
# ---------------------------------------------------------------------------


my @files =<the-path/*ntfreq>;


##print join("\n",@files),"\n";

my @sites = 1..18959; # length of the reference EBOV genome

# $rh_all->{site}->{samp}->{'cov'}
# $rh_all->{site}->{samp}->{'cov'}
my $rh_all;
my $rh_samp;
my $rh_sample2pos;
my @metrics;
for (my $i = 0; $i < @files; $i++){
	print "#",$i+1,": $files[$i] loading...\n";
	if ($files[$i] =~ /([\w\d\-]+)\.mpileup/ || $files[$i] =~ /([\w\d\-]+)\.0504\.xls/) {
	    my $samp = $1;
	    my $rh_siteinfo = read_ntfreq($files[$i]);
	    my @pos = keys %$rh_siteinfo;
	    if (defined $rh_samp->{$samp}){
		print "\t[idRepeatWarning] sample ID repeated: $samp\n";#$files[$i]\n$rh_samp->{$samp}\n";
		if (defined $pos[0] && @pos > $rh_sample2pos->{$samp}){
			print "\t[idRepeatWarning] select $files[$i] (posNum=".@pos.")\n\t[idRepeatWarning] overwrite values from $rh_samp->{$samp} (posNum=".$rh_sample2pos->{$samp}.")\n";
	    		$rh_samp->{$samp}        = $files[$i];
			$rh_sample2pos->{$samp}  = @pos;
		} else {
			print "\t[idRepeatWarning] skip $files[$i] (posNum=".@pos.")\n\t[idRepeatWarning] keep values from $rh_samp->{$samp} (posNum=".$rh_sample2pos->{$samp}.")\n";
			next;
		}
	    } else {
	    	$rh_samp->{$samp} = $files[$i];
		$rh_sample2pos->{$samp} = @pos;
            }
		
	    @metrics = sort keys %{$rh_siteinfo->{$pos[0]}} if (! defined $metrics[1]);
	    for (my $p = 0; $p < @sites; $p++){
       		if (defined $rh_siteinfo->{$sites[$p]}){
     		       $rh_all->{$sites[$p]}->{$samp} = $rh_siteinfo->{$sites[$p]};
  	        } else {
          	       #$rh_all->{$sites[$p]}->{$samp} = "na";
  	     	}
    	    }
	    print "\t$samp\tpos num\t".@pos."\n";
	} else {
		die "Cannot extract sample id from file name: $files[$i]\n";
	}
}

my @samps = sort keys %$rh_samp;
for (my $i = 0; $i < @metrics; $i ++){
	print "# output metric $metrics[$i] to $metrics[$i].bigtables.txt\n";
	open OUT, "> $metrics[$i].bigtables.txt" or die "$!: $metrics[$i] output error\n";
	print OUT "#$metrics[$i]\t";
	for (my $j = 0; $j < @samps; $j++){
		print OUT "$samps[$j]\t";
	} print OUT "\n";
	for (my $p = 0; $p < @sites; $p++){
		print OUT $sites[$p]."\t";
		for (my $j = 0; $j < @samps; $j++){
			if (defined $rh_all->{$sites[$p]}->{$samps[$j]}->{$metrics[$i]}) {
				print OUT $rh_all->{$sites[$p]}->{$samps[$j]}->{$metrics[$i]}."\t";
			} else {
				print OUT "na\t";
			}
		} print OUT "\n";
	}
	close OUT;
}



# =======================================================
# read a given mpipeup2ntfreq file
sub read_ntfreq{
	open FP, "$_[0]" or die "sub read_ntfreq, $!: $_[0]\n";
	my $rh_siteinfo;
	my @headers;
	my $DEBUG = 0; 
#chr	pos	ref	dep(mpileup)	dep(-indel)	ntfreqPatt	A	G	C	T	tot(Q20)	Fa	Fg	Fc	Ft	A5'/3'	G5/3	C5/3	T5/3	A5'(na_if_dep<5)	G5	C5	T5	entropy
	while (<FP>){
        	chomp;
		if (/^#/){
			s/^#//g;
			@headers = split /\s+/,$_;
			for (my $i = 0; $i < @headers; $i++){
				$headers[$i] =~ s/[\(\)\-\<\>\/\']/_/g;		
			}
		} elsif (/^[a-zA-Z]/) {
			my @row = split /\s+/, $_;
			my $pos = $row[1];
#			for (my $i = 2; $i < @headers; $i++){ # pos->otherHeader = values
#				$rh_siteinfo->{$row[1]}->{$headers[$i]} = $row[$i];
#			}
#
			my @ntCov = ($row[6],$row[7],$row[8],$row[9]);
			my @ntCovSort = sort {$b <=> $a} @ntCov;
			my ($ra_rank,$ra_rankAlle) = get_des_rank(\@ntCov);#@rank[ sort { $ntCov[$b] <=> $ntCov[$a] } 0..$#ntCov] = 0..$#ntCov;		
			my @rank = @$ra_rank;
			print "\nAGCT cov for debug\nntcov:@ntCov\nrank:@rank\nrankalle:"."@{$ra_rankAlle}"."\n" if($DEBUG==1);
			my @bases = qw/A G C T/;
			if ($row[10] > 0){
				$rh_siteinfo->{$pos}->{'out_cov_from_mpileup'}= $row[3];
				$rh_siteinfo->{$pos}->{'out_cov_noIndel'}     = $row[4];
				$rh_siteinfo->{$pos}->{'out_cov_noIndel_QC'}  = $row[10];
				$rh_siteinfo->{$pos}->{'out_entropy'}         = $row[@row-1];
				$rh_siteinfo->{$pos}->{'out_minor_cov'}       = $ntCovSort[1] + $ntCovSort[2] + $ntCovSort[3]; 
				$rh_siteinfo->{$pos}->{'out_second_cov'}      = $ntCovSort[1]; 
  		      		$rh_siteinfo->{$pos}->{'out_minor_freq'}      = $rh_siteinfo->{$pos}->{'out_minor_cov'}/$row[10];
        			$rh_siteinfo->{$pos}->{'out_second_freq'}     = $rh_siteinfo->{$pos}->{'out_second_cov'}/$row[10];
				$rh_siteinfo->{$pos}->{'out_nt_pattern'}      = $row[5];
				$rh_siteinfo->{$pos}->{'out_max_stranded'}    = $row[15+$rank[0]] ; 
				$rh_siteinfo->{$pos}->{'out_second_stranded'} = $row[15+$rank[1]]; 
				$rh_siteinfo->{$pos}->{'out_max_stranded_posRatio'}    = $row[19+$rank[0]] ; 
				$rh_siteinfo->{$pos}->{'out_second_stranded_posRatio'} = $row[19+$rank[1]]; 
				$rh_siteinfo->{$pos}->{'out_ref_max_sec'}     = "$row[2]:$bases[$rank[0]]-$bases[$rank[1]]"; 
				$rh_siteinfo->{$pos}->{'out_ranked_alle'}     = join("-",@$ra_rankAlle);
			}
			if ($DEBUG == 1){
				print join("\t",@row),"\n";
				while (my ($h,$v) = each %{$rh_siteinfo->{$pos}}){
					print "$h\t$v\n";
				}
			}
       		}
    	}
   	close FP;
	my @poss = keys %$rh_siteinfo;
	print "# total key no. = ".@poss."\n" if ($DEBUG == 2 );
   	return $rh_siteinfo;
}

sub min_freq{
    my @cou = sort {$a<=>$b} @_;
    my $tot = 0;
    foreach (@cou){
        $tot += $_;
    }
    return ("noRead","noRead") if ($tot <= 0); 
    return (1-$cou[@cou-1]/$tot, $cou[@cou-2]/$tot);
}

sub get_des_rank{
	my @arr = @{$_[0]};
	my $DEBUG = 0;
	print "get_des_rank();@arr\n" if ($DEBUG == 1);
	my @bases = qw/A G C T/;
	my $rh_oldrank;
	for(my $i = 0; $i < @bases; $i ++){
		$rh_oldrank->{$bases[$i]} = $i;
	}
	my $rh_val2key;
	for(my $i = 0; $i < @arr; $i++){
		if (defined $rh_val2key->{"$arr[$i]"}){
			$rh_val2key->{"$arr[$i]"} .= "/".$bases[$i];
		} else {
			$rh_val2key->{"$arr[$i]"} = $bases[$i];
		}
	}
	my @vals = sort {$b<=>$a} keys %$rh_val2key;
	my @newrank = ();
	my @newrankAlle = ();

	print "get_des_rank():@vals\n" if ($DEBUG == 1);
	for (my $i = 0; $i < @vals; $i++ ) {
		my $alle = $rh_val2key->{$vals[$i]};
		push @newrankAlle, $alle;
		# push the index
		print "get_des_rank():$alle\n" if ($DEBUG == 1);
		if (length($alle) == 1){
			push @newrank,$rh_oldrank->{$alle};
		} elsif (length($alle) > 1) {
			my @bs = split /\//,$alle;
		print "get_des_rank():-->@bs\n" if ($DEBUG == 1);
			my $degranks = '';
			foreach my $ba (@bs){
				die "get_des_rank(): cannot find rank for base $ba\n" if (! defined $rh_oldrank->{$ba});
				$degranks .= '/'.$rh_oldrank->{$ba};
			}
			$degranks =~ s/^\///;
			$degranks = $bs[0]; # <---------- for simple in its parent sub, just get the first
			push @newrank, $rh_oldrank->{$degranks};
		}
	} 
	print "get_des_rank():@newrank\n" if ($DEBUG == 1);
	print "get_des_rank():@newrankAlle\n" if ($DEBUG == 1);
	return (\@newrank, \@newrankAlle);
}

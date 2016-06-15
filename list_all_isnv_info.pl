use warnings;
use strict;

my $isnvFile = "iSNV.all.txt";
my $rh_isnv = read_isnv_file($isnvFile);

my $bigtablePath = ".";
my @bigtableFiles = <$bigtablePath/*bigtables.txt>;
my ($rh_isnv_bigtables, $ra_terms) = read_isnv_from_bigtable_files(\@bigtableFiles,$rh_isnv);


print_isnv_hash($rh_isnv, $rh_isnv_bigtables, $ra_terms);




sub read_isnv_file{
	open FP, "$_[0]" or die "read_isnv_file: $!\n";
	my @pos;
	my $rh_isnv;
	while (<FP>){
		chomp;
		if (/^samples\/pos/){
			@pos = split /\s+/,$_;
		} elsif (/^[\w\d\_\-]/) {
			my @row = split /\s+/, $_;
			for(my $i = 1; $i < @row; $i ++){
				if ($row[$i]>-0.1){
					$rh_isnv->{$row[0]}->{$pos[$i]} = $row[$i];
				}
			}
		}
	}
	close FP;
	return $rh_isnv;
}

sub read_isnv_from_bigtable_files{
	my ($ra_file, $rh_isnv) = @_;
	my @headers;
	my @terms;
	my $rh_isnv_table;
	foreach my $file (@$ra_file){
		open FP, "$file" or die "$file\n";
		while (<FP>){
			chomp;
			if (/^\#/){
				@headers = split /\s+/, $_;
				push @terms,$headers[0];
			} elsif (/^\d+/) {	
				my @row = split /\s+/, $_;
				for(my $i = 1; $i < @row; $i ++){
					if (defined $rh_isnv->{$headers[$i]}->{$row[0]}){
						$rh_isnv_table->{$headers[$i]}->{$row[0]}->{$headers[0]} = $row[$i];
						#print "$headers[$i]\t$row[0]\t$headers[0]\t$row[$i]\n";
					}
				}
			}
		}
		close FP;
	}
	return ($rh_isnv_table, \@terms);
}

sub print_isnv_hash{
	my ($rh_isnv,$rh_bigtables, $ra_terms) = @_;
	print "#samp\tpos\tisnvFreq\t";
	for (my $i = 0; $i< @{$ra_terms}; $i++){
#		$$ra_terms[$i] =~ s/\#//;	
		print "$$ra_terms[$i]\t";
	} print "\n";
	while (my ($samp, $rh_dat) = each %$rh_isnv){
		my @p = sort {$a<=>$b} keys %$rh_dat;
		print "\/\/ sample=$samp:\n";
		for(my $i = 0; $i < @p; $i ++){
			print "$samp\t$p[$i]\t$rh_dat->{$p[$i]}\t";
			for (my $j = 0; $j< @{$ra_terms}; $j++){
				print '\'' if ($rh_bigtables->{$samp}->{$p[$i]}->{$$ra_terms[$j]} =~ /\d+\/\d+/);
				print "$rh_bigtables->{$samp}->{$p[$i]}->{$$ra_terms[$j]}\t";
			}
			print "\n";
		}
		print "\n";
	}
}

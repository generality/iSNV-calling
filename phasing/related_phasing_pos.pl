use warnings;
use strict;

#my @files = <./serial_isnv_phasing/*phasing>;

my @files = <./*phasing>;

@files = sort @files;
print "#pos\tref>mut\tfreq\tphasedRatioOfMutant\tphasedMutX\tunphasedMutX\twildtypeX\tphasedMuXF\tunphasedMutF\twildtypeF\n";
for(my $i=0; $i<@files; $i++){
    print "\n$files[$i]\n";
    pos_serial_corr($files[$i]);
}

sub pos_serial_corr{
    open FP,"$_[0]" or die "$!: $_[0]\n";
    my @pos;
    my ($rh_ref,$rh_mut,$rh_f);
    my ($r1,$r2,$m1,$m2);
    my ($rh_corr, $rh_uncorr, $rh_wildtype);
    while (<FP>){
        chomp;
        if(/^\d/){
            @pos = split /\s+/,$_;
            for(my $i=0;$i<@pos-1;$i++){
                $rh_corr->{$i} = 0;
                $rh_uncorr->{$i} = 0;
                $rh_wildtype->{$i} = 0;
            }
        } elsif (/\w\>\w\|/){
            ($rh_ref,$rh_mut,$rh_f) = handle_line2($_);
        } elsif (/^[\w\.]/) {
            my @bases = split /\s+/,$_;
            
            for(my $i=0;$i<@bases-1;$i++){
                $r1 = $rh_ref->{$i};
                $r2 = $rh_ref->{$i+1};
                $m1 = $rh_mut->{$i};
                $m2 = $rh_mut->{$i+1};
                if ($bases[$i] eq $m1 && $bases[$i+1] eq $m2){
                    $rh_corr->{$i} ++;
                } elsif ($bases[$i] ne "." && $bases[$i+1] ne ".") {
                    if ($bases[$i] eq $m1 || $bases[$i+1] eq $m2){
                        $rh_uncorr->{$i} ++;
                    } else {
                        $rh_wildtype->{$i} ++;
                    }
                }
            }
        }
    }
    close FP;
    for(my $i=0;$i<@pos;$i++){
        print "$pos[$i]\t$rh_ref->{$i}>$rh_mut->{$i}\t$rh_f->{$i}\t";
        if ($i != @pos-1){
            my $tot = $rh_corr->{$i} + $rh_uncorr->{$i} + $rh_wildtype->{$i};
	    if ($tot>0){
                print int($rh_corr->{$i}/($rh_corr->{$i}+$rh_uncorr->{$i})*1000)/1000,"\t";
	    } else {
		print "-\t";
	    }
            print "$rh_corr->{$i}\t";
            print "$rh_uncorr->{$i}\t";
            print "$rh_wildtype->{$i}\t";
            if ($tot>0){
                print int($rh_corr->{$i}/$tot*1000)/1000,"\t";
                print int($rh_uncorr->{$i}/$tot*1000)/1000,"\t";
                print int($rh_wildtype->{$i}/$tot*1000)/1000,"\n";
            } else {
                print "-\t-\t-\n";
            }
        } else {print "\n";}
    }
}

sub handle_line2{
    my ($rh_ref,$rh_mut,$rh_f);
    my @row = split /\s+/,$_[0];
    for(my $i=0;$i<@row;$i++){
        if ($row[$i] =~ /(\w)\>(\w)\|([\d\.]+)/){
            $rh_ref->{$i} = $1;
            $rh_mut->{$i} = $2;
            $rh_f->{$i}   = $3;
        } else {
            die "[error] cannot reg $_[0]\n"
        }
    }
    return ($rh_ref, $rh_mut, $rh_f);
}


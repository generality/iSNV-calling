use warnings;
use strict;

my $rh_patient2group2pos = load_patGrouPos();

#find_read(sam_file, rh_position_group);
while (my ($patient, $rh_group) = each %$rh_patient2group2pos){
    my @groups = sort {$a<=>$b} keys %$rh_group;
    for (my $i=0;$i<@groups;$i++){
	my $outfile =  $patient."_".$groups[$i].".sam.phasing";
        print "\n\npatient=$patient, group=$groups[$i], output to $outfile\n";

	find_read($patient.".sam", $rh_group->{$groups[$i]},$outfile);
    }
}

#$rh_test->{1241} = 1;
#$rh_test->{1243} = 1;
#find_read('tt.sam',$rh_test);

sub find_read{
    open FP,"$_[0]" or die "$!:$_[0]\n";
    my $rh_pos = $_[1];
    my $outfile = $_[2];
    open OUT, ">$outfile" or die "$!: $outfile\n";
    my @pos =  sort {$a<=>$b} keys %$rh_pos;
    my $min = $pos[0];
    my $max = $pos[@pos-1];
    print OUT join(" ",@pos),"\n";
    for (my $i=0;$i<@pos;$i++){
    	print OUT $rh_pos->{$pos[$i]}," ";
    } print OUT "\n";
    #print "@pos\n$min&$max\n";
    while (<FP>){
        chomp;
        next if(/^\@/ || /^\:/ || /^ana/);
        my @row = split /\s+/,$_;
        if(!/^\@/ && $row[2] ne "*"){
            my $start = $row[3];
            my $cigar = cigarTransfor($row[5]);
            my @bases = split //,$row[9];
            #print "\n$row[0]\n@pos\n$row[3]\t$cigar\n$row[9]\n";
            next if ($row[5] eq '*' || $start+@bases+10 < $min || $start > $max );

            my @cigarPos = split /[A-Z]+/, $cigar; 
            my @cigarLab = split /\d+/, $cigar;
            @cigarLab = @cigarLab[1..@cigarLab-1];
            # ------------------------------------------------------------------
            # -- align the alignment of a reads, with insert del and del='-' ---
            my $align='';
            my $walk = 0;
            for(my $i=0;$i<@cigarPos;$i++){
                if ($cigarLab[$i] eq 'M'){
                    die "[error] $row[5]\t$cigar\twalk=$walk\t$cigarPos[$i]\treadlen=".@bases."\n$_\n" if ($walk+$cigarPos[$i]>@bases);
#                    print @bases[$walk.$walk+$cigarPos[$i]-1]."\n";
                    $align .= join("",@bases[$walk..$walk+$cigarPos[$i]-1]);
                    $walk += $cigarPos[$i];
                } elsif ($cigarLab[$i] eq 'I'){
                    $walk += $cigarPos[$i];
                } elsif ($cigarLab[$i] eq 'D'){
                    for(my $j=0;$j<$cigarPos[$i];$j++){$align .= "-";};
                } else {
                    die "[error] extended cigar: $row[5], $cigar, $cigarLab[$i]\n|@cigarPos|\n|@cigarLab|\n$_\n";
                }
            }
            # ------------------------------------------------------------------
            # ---- the read covers >=2 position -----
            my @align = split //, $align;
            #print $align."\n";
            for (my $j=0;$j<@pos;$j++){
                if ($pos[$j]-$start>=0 && $pos[$j]-$start<@align){
                    print OUT $align[$pos[$j]-$start]." ";
                } else {
                    print OUT ". ";
                }
            }
            print OUT "\n";
        }
    }
    close FP;
    close OUT;
}


# ===== eg. 5MD10M = 5M1D10M ========
sub cigarTransfor{
    my $cigar = $_[0];
    my @cigars = split //,$cigar;
    my $newCigar = $cigars[0]; 
    for(my $i=1;$i<@cigars;$i++){
        if ($cigars[$i-1]=~ /[A-Z]/ && $cigars[$i]=~ /[A-Z]/){
            $newCigar .= "1".$cigars[$i];
        } else {
            $newCigar .= $cigars[$i];
        }
    }
    return $newCigar;
}

sub load_patGrouPos{
    #my $file = "patientId_SerialGroupId_pos.txt";
    my $file = "patientId_SerialGroupId_pos_20161217.txt";
    #my $file = "tmp.txt";
    open FP,"$file" or die "load_patGrouPos, $!: $file\n";
    my $rh;
    while (<FP>){
        chomp;
        if (/(\d+)\s+(\d+)\s+(\d+)\s+([\w\d\.\-\>\|]+)/){
            $rh->{$1}->{$2}->{$3} = $4;
        }
    }
    close FP;
    return $rh;
}

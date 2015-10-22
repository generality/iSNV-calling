#!/usr/bin/perl -w

# -------------------------------------------------------
#	Select the properly paired reads from the clean data
# ------------------------------------------------------

my $inputFile1 = './test_data/1.fastq';
my $inputFile2 = './test_data/2.fastq';
my $outfile1='./test_data/1.fastq.pair';
my $outfile2='./test_data/2.fastq.pair';
my $outfile3='./test_data/1.fastq.unpair';
my $outfile4='./test_data/2.fastq.unpair';

my $temp_time=localtime(time);								# start time
print "running check_pair.pl at:$temp_time\n";
#----------------------------------------main------------------------------------------
$inputFile1=$ARGV[0];
$inputFile2=$ARGV[1];
$outfile1=$ARGV[2];
$outfile2=$ARGV[3];

my %list1=set_read($inputFile1);
my %list2=set_read($inputFile2);
my $old1=keys %list1;
my $old2=keys %list2;
my $intersection=0;
my $single1=0;
my $single2=0;
if($old1>0&&$old2>0)#not empty
{
	#out_hash(\%list1,$outfile1);
	#out_hash(\%list1,$outfile2);
	open MY_OUT1,"> $outfile1" or die "Can't open $outfile1: $!";#open output file
	open MY_OUT2,"> $outfile2" or die "Can't open $outfile2: $!";#open output file
	open MY_OUT3,"> $outfile1" . ".single" or die "Can't open $outfile1 .single: $!";#open output file
	open MY_OUT4,"> $outfile2" . ".single" or die "Can't open $outfile2 .single: $!";#open output file
	#my @res1;
	#my @res2;
	#my $res1_count=0;
	#my $res2_count=0;

	print "test1\n";
	foreach $key (keys %list1)								#get the intersection
	{
		if(exists $list2{$key})								#intersection
		{
		#	$res1[$res1_count++]=$list1{$key};
		#	$res2[$res2_count++]=$list2{$key};
			print MY_OUT1 $list1{$key};
			print MY_OUT2 $list2{$key};
			$intersection++;
		}
		else
		{
			print MY_OUT3 $list1{$key};
			$single1++;
		}
	}
	print "test2\n";
	foreach $key (keys %list2)								#
	{
		if(exists $list1{$key})								#skip
		{
		#	$res1[$res1_count++]=$list1{$key};
		#	$res2[$res2_count++]=$list2{$key};
		}
		else
		{
			print MY_OUT4 $list2{$key};
			$single2++;
		}
	}

#@res1=sort { $a cmp $b } @res1;
#@res2=sort { $b cmp $a } @res2;
#out_array(\@res1,$outfile1);
#out_array(\@res2,$outfile2);
	close MY_OUT1;
	close MY_OUT2;
	close MY_OUT3;
	close MY_OUT4;

	print "old1\t$old1(",flt_to_pct($intersection/$old1),")\t";
	print "old2\t$old2(",flt_to_pct($intersection/$old2),")\t";
	print "intersection\t$intersection\t";
	print "single1\t$single1(",flt_to_pct($single1/$old1),")\t";
	print "single2\t$single2(",flt_to_pct($single2/$old2),")\n";
}
add_ipc();#
$temp_time=localtime(time);#
print "Exit successfully at:$temp_time\n";
#----------------------------------------sub proc------------------------------------------
sub set_read
{
	my $inputfile=$_[0];
	my %temp_result;
	open MY_IN, "$inputfile" or die "Can't open $inputfile: $!\n";#
	my $count=0;
	while (my $line=<MY_IN>)#
	{
		$line=~s/[\n\r]//g;
		if($line ne "")
		{
			my @temp0=split(/\|/,$line);
			my @temp=split(/[\ _\/]/,$temp0[$#temp0]);
			my $line2=<MY_IN>;
			my $line3=<MY_IN>;
			my $line4=<MY_IN>;
			$line2=~s/[\n\r]//g;
			$line3=~s/[\n\r]//g;
			$line4=~s/[\n\r]//g;
			$temp_result{$temp[0]}=join "\n",$line,"$line2","$line3","$line4\n";
		}
	}
	close MY_IN;
	return %temp_result;
}
sub out_hash												#print hash list
{
	my (%list)=%{$_[0]};									#input file
	my $outputfile=$_[1];									#output file

	open MY_OUT,"> $outputfile" or die "Can't open $outputfile: $!";#open output file
	foreach $key (keys %list)#
	{
		print MY_OUT $key,"\t",$list{$key},"\n";
	}
	close MY_OUT;
}
sub out_array												#print array list
{
	my @list=@{$_[0]};										#input file
	my $outputfile=$_[1];									#output file

	open MY_OUT,"> $outputfile" or die "Can't open $outputfile: $!";#open output file
	for(my $i=0;$i<=$#list;$i++)#
	{
		print MY_OUT $list[$i];
#		print $list[$i],"\n";
	}
	close MY_OUT;
}
sub flt_to_pct
{
    sprintf( "%.4f", shift ) * 100 . '%';
}
sub add_ipc
{
	my $pid=$$;
	my $ipc_path = './tmp_ipc';
	system("echo add >$ipc_path/$pid");
}

#!/usr/bin/perl
if($#ARGV == 2){
	$infile = $ARGV[0];
	$discards = $ARGV[1];
	$outfile = $ARGV[2];
}
else{
	die("Removes sequences from a migrate-n infile according to specified indexes.\n\nInputs:\n\tmigrate-n infile\n\tFile containing indices to discard\n\toutfile\n\nFile of indices to discard should have one line, with indices tab seperated.\n");
}


use Data::Dumper qw(Dumper);

open(INFILE, "<", $infile) or die("Can't find infile $infile\n");
open(DISCARDS, "<", $discards) or die("Can't find discards $discards\n");

#get discards
my $d = <DISCARDS>;
chomp($d);
my @da = split("\t", $d);
@da = sort(@da);

#correct for array index starting at zero.
my $i = 0;
while($i <= $#da){
	$da[$i]--;
	$i++;
}
my %dels = map { $_ => 1 } @da; #map to hash.


close(DISCARDS);

open(OUTFILE, ">", $outfile);

my $start = 1;
my $npop = 0;
OUTER: while (<INFILE>){
	chomp($_);
	#deal with header info and n per locus stuff
	if($start == 1){
		$start++;
		@header = split(" ", $_);
		#change the number of loci to account for those removed.
		$header[1] -= ($#da + 1);
		#print Dumper \@header;
		print OUTFILE ("   $header[0] $header[1]  $header[2]\n");
		next OUTER;
	}
	if($start == 2){
		#print the lengths after shortening the array by the number to be removed.
		#print("Doing lengths.\n");
		@lengths = split(' ', $_);
		my $j = 0;
		while($j <= $#da){
			pop(@lengths);
			$j++;
		}
		print OUTFILE (join(" ", @lengths), "\n");
		$start++;
		next OUTER;
	}
	my $line = $_;
	my @la = split(' ', $line);

	#if the line is a list of number of inds sequenced per loci...
	if($#la > 0){
		#print("New pop.\n");
		$npop++;
		$tloc = 0; #tracker for current loci
		$pop = pop(@la);
		#print("Current counts per loci:");
		#print Dumper \@la;
		@locn = @la;
		my $j = 0;
		#print("\n");
		#print this out other than any loci to delete.
		while($j <= $#la){
			if(!exists($dels{$j})){
				print OUTFILE ("$la[$j]");
				print OUTFILE ("   ");
			}
			$j++;
		}
		print OUTFILE ("$pop\n");
		$i = 1;
		next OUTER;
	}

	#check if we have changed loci
	if($i > $locn[$tloc]){
		#print("We changed loci. i = $i\n");
		#print("Resetting i.\n");
		$i = 1;
		$tloc ++;
	}

	#print("This Locus: $tloc. This sample: $i.\n");
	#print it out as long as it isn't a loci to remove
	if(!exists($dels{$tloc})){
		#print("Printing the line: $line\n");
		print OUTFILE ("$line\n");
		$i++;
	}
	else{
		#print("Not printing the line.\n");
		$i++;
	}
}

close OUTFILE;
close INFILE;
	 		

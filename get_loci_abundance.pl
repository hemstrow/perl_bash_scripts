#!/usr/bin/perl
#get inputs if 3 inputs, else die
if($#ARGV == 0 || $#ARGV == 2){
	$blist = $ARGV[0]; #get the list of input BAM files
	$exclude = $ARGV[1]; #file of positions which indicate reads to exclude.
	$rlength = $ARGV[2]; #read length, for exclusion only
}
else{
	die("Generates a table containing the number of individuals that each loci is genotyped in.\nUseage: get_rand_reads.pl input_bamlist optional_excluded_positions optional_read_length\nArguments 2 and 3 are optional.");
}

use Data::Dumper qw(Dumper);

#add exclusion data if given
if($exclude){
	print("Generating exclusion hash.\n");
	%e; #intialize exclusion hash	
	open(EXCLUDE, "<", $exclude) or die("Can't locate exclusion hash\n"); #open the file of exclude positions if it is given
	while(<EXCLUDE>){
		chomp($_);
		my @tsnp = split(/\s/, $_); #split the line
		push(@{$e{$tsnp[0]}}, $tsnp[1]);
		#print("$tsnp[0]: ");
		#foreach $element (@{$e{$tsnp[0]}}) {print("$element\t");}
		#print("\n");
	}
	close(EXCLUDE);
	print("Done.\n");
}

print("Counting loci...\n");

#loop through the bamlist, count number of individuals sequenced for each loci
%trc;
open(BLIST, "<", $blist) or die("Can't find bamlist.\n");
my $ind_prog = 1;
OUTER: while(<BLIST>){
	my @dat = split(/\n/, `samtools view $_`); #convert and split the bam file into an array, one read per element
	my %rc;
	if($ind_prog % 5){
		print("\tIndividual: $ind_prog\n");
	}
	MID: foreach $read (@dat){ #for each read in this individual
		chomp($read); 
		my @rarray = split(/\t/, $read); 
		my $tchr = $rarray[2]; #chr
		my $tpos = $rarray[3]; #position
		my $tseq = $rarray[9]; #sequence
		my $rid = $tchr . "__" . $tpos;
		if($rc{$rid} == -1){
			next MID; #if already caught by exclusion hash, skip
		}
		
		#check against exclusion list
		if($exclude){
			if(exists($e{$tchr})){
	                	my $epos = $tpos + $rlength - 1; #get ending postion
        	                #print("This read: $tpos to $epos.\n");
        	                foreach my $vpos (@{$e{$tchr}}){ #for each excluding position on this lg/chr
                	                #print("Checking position $vpos against read...\n");
                        	        if ($tpos <= $vpos && $epos >= $vpos){ #if excluding postion is on this read...
						$rc{$rid} = -1; #set output to -1 for exclusion
						next MID; #skip to the next read
					}
				}
			}
		}
				
		#Increase count for this locus in this individual. Could use this later to get number of unique reads per ind or something.
		$rc{$rid}++;
	}
	
	#read %rc, store info for final output
	foreach my $locus (sort keys %rc){
		if ($rc{$locus} == -1){
			$trc{$locus} = -1;
		}
		elsif ($rc{$locus} == 1){
			$trc{$locus}++;
		}
	}
	$ind_prog++;
}

##########################################################################################################################
#print output
open(RINFO, ">", "locusinfo.txt");

foreach my $rid (sort keys %trc){
	my @locus = split(/__/, $rid);
	my $tchr = $locus[0];
	my $tpos = $locus[1];
	print RINFO ("$tchr\t$tpos\t$trc{$rid}\n");
}

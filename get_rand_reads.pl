#!/usr/bin/perl
#get inputs if 3 inputs, else die
if($#ARGV >= 3){
	$nrds = shift(@ARGV); #get the number of reads to grab
	$blist = shift(@ARGV); #get the list of input BAM files
	$rlength = shift(@ARGV); #read length
	$outfile = shift(@ARGV); #output filename
}
else{
	die("\nUseage for get_rand_reads.pl:\n\nRequried Arguments, in order:\n\tNumber of reads per individual to select.\n\tList of bam files containing aligned sequences. Different pops can be seperated with a line containing the name for the following pop in the tab seperated format POP:<popname>.\n\tLength of reads.\n\tName of outfile.\n\nOther Arguments:\n\t-e: file with reads to exclude in tab delim format <linkage_group><position_on_read>\n\t-r: file with candidate reads to use. If specified, will select reads ONLY from this set. Tab seperated format <linkage_group><starting_position>.\n\t-m: Required if -r not set. Number of individuals in which a locus must be sequenced to be kept.\n\t-maf: minimum minor allele frequency to accept.\n\n");
}

use Data::Dumper qw(Dumper);

use List::Util qw(max);

print("Checking arguments...\n");
######################################################################################################################################################################
#check remaining arguments and prep exclusion hash and candreads if given
if($#ARGV >= 0){
	my @aa = @ARGV;
	foreach $element (@aa){
		$element =~ s/\-//;
	}
	my %args = @aa;
	
	if(exists($args{m})){
		$max_miss = $args{m};
	}

	if(exists($args{maf})){
		$maf = $args{maf};
	}

	#check if exclusion info given
	if(exists($args{e})){
		print("Generating exclusion hash.\n");
		%e; #intialize exclusion hash
		open(EXCLUDE, "<", $args{e}) or die("Can't locate exclusion file.\n"); #open the file of exclude positions if it is given
		while(<EXCLUDE>){
			chomp($_);
			my @tsnp = split(/\s/, $_); #split the line
			push(@{$e{$tsnp[0]}}, $tsnp[1]);
		}
		close(EXCLUDE);
		print("Done.\n");
	}


	#check if candidate read list given.
	if(exists($args{r})){
		print("Generating hash of candidate loci from provided list...\n");
		#get array of read identifiers to use
		open (CANDREADS, "<", $args{r}) or die("Can't find input $args{r}\n");
		while(<CANDREADS>){
			chomp($_);
			my @read = split(/\t/, $_);
			my $chr = $read[0];
			my $pos = $read[1];
			my $coord = $chr . "__" . $pos;
			#check against exclusion hash if exists
			if(exists($e{$chr})){
					my $epos = $pos + $rlength - 1; #get ending position
					#print("This read: $tpos to $epos.\n");
					foreach my $vpos (@{$e{$chr}}){ #for each excluding position on this lg/chr
						if ($pos <= $vpos && $epos >= $vpos){ #if excluding position is on this read...
						$rc{$rid} = -1; #set output to -1 for exclusion
					}
				}
			}
			else{
				push(@candreads, $coord);
			}
		}
		close (CANDREADS);
	}
}


#initialize a few things, open blist and array it, since we'll be looping through it a bunch and it's going to be small memory wise.
my @barray;
open (BLIST, "<", $blist) or die("Can't locate bam list $blist\n");
while(<BLIST>){
	push(@barray, $_);#add this line to the array.
}

close(BLIST); #close the now unused connections


######################################################################################################################################################################
#loop through each file and compile list of candidate loci if candreads not given.
if(!@candreads){
	if(!$max_miss){
		die("If no list of candidate loci to use is provided (-r), minimum number of individuals sequenced needed to keep loci must be set (-m).\n");
	}
	print("Generating list of candidate loci from data, minimum individuals = $max_miss.\n");
	if($maf){
		print("Filtering by minor allele frequency: $maf.\n");
		%lmaf; #intialize maf info if necissary.
	}
	OUTER: foreach my $ind (@barray){
		if ($ind =~ /^POP\:/){ #skip pop names for now.
			next OUTER;
		}
		my @dat = split(/\n/, `samtools view $ind`); #convert and split the bam file into an array, one read per element
		my %rc;
		MID: foreach $read (@dat){ #for each read in this individual
			chomp($read); 
			my @rarray = split(/\t/, $read); 
			my $chr = $rarray[2]; #chr
			my $pos = $rarray[3]; #position
			my $rid = $chr . "__" . $pos;
			if($rc{$rid} == -1){
				next MID; #if already caught by exclusion hash, skip
			}

			#check against exclusion list if given
			if(exists($e{$chr})){
				my $epos = $pos + $rlength - 1; #get ending postion
				foreach my $vpos (@{$e{$chr}}){ #for each excluding position on this lg/chr
					#print("Checking position $vpos against read...\n");
					if ($pos <= $vpos && $epos >= $vpos){ #if excluding position is on this read...
						$rc{$rid} = -1; #set output to -1 for exclusion
						next MID; #skip to the next read
					}
				}
			}
			
			#Increase count for this locus in this individual. Could use this later to get number of unique reads per ind or something.
			#crucially, could use this filter out loci with more than X number of unique reads in an individual
			$rc{$rid}++;

			if($maf){
				#Increase the allele count for this allele.
				$lmaf{$rid}{$rarray[9]}++;
			}
		}
		
		#read %rc, store info for final set of loci to use.
		foreach my $locus (sort keys %rc){
			if ($rc{$locus} == -1){
				$trc{$locus} = -1;
			}
			elsif($rc{$locus} > 2){
				$trc{$locus} = -1;
			}
			else{
				$trc{$locus}++;
			}
		}
	}
	
	#get candreads
	foreach my $rid(sort keys %trc){
		if($trc{$rid} >= $max_miss){ #check number sequenced in
			if($maf){ #if maf exists, check minor allele frequency.
				my %thash = %{$lmaf{$rid}};
				my $tac = 0;
				my $mac = 0;
				my $maj = 0;
				#calculate minor allele frequency (maf)
				foreach my $af (keys(%thash)){
					if($thash{$af} > $maj){
						$mac += $maj;
						$maj = $thash{$af};
						$tac += $maj;
					}
					else{
						$mac += $thash{$af};
						$tac += $thash{$af};
					}
				}
				$tmaf = $mac / $tac;
				
				if($tmaf >= $maf){
					push(@candreads, $rid);
				}
			}
			#otherwise just add the read, since it already passed minimum sequencing filter.
			else{
				push(@candreads, $rid);
			}
		} 
	}
	
	#if candreads is empty, print out locus info and die.
	print("$#candreads candidate reads found.\n");
	if(!@candreads || $#candreads < $nrds){
		open(RINFO, ">", "locusinfo.txt");
		foreach my $rid (sort keys %trc){
			my @locus = split(/__/, $rid);
			my $tchr = $locus[0];
			my $tpos = $locus[1];
			print RINFO ("$tchr\t$tpos\t$trc{$rid}\n");
		}
		if(!@candreads){
			die("No loci sequenced at at least $max_miss individuals found. Printing table of number of sequenced individuals per loci to locusinfo.txt\n");
		}
		else{
			print("Not enough loci sequenced at at least $max_miss individuals found. Printing table of number of sequenced individuals per loci to locusinfo.txt\n");
		}
	}
}

print("Getting reads...\n");



######################################################################################################################################################################
#Pick random read from each individual for nrds number of loci from candreads.
#while enough reads haven't been found, select a random read, loop through output, grab matches while one is found, skip to another random read if any reads don't find a match
my $nsaved = 0;
my %kh; #initialize hash to of reads to keep. Structure: kh(hash), pop(hash), locus(array), <ind><read>(concatenated variable).
my %ci; #keep number of individuals sequenced per read per pop for output. Structure: ci(hash), pop(hash), locus (array), number sequenced (variable).
OUTER:while($nsaved < $nrds && @candreads){
	my $indnum = 1; #initialize individual number per pop tracker for output.
	my $cpop = "popxx"; #set default pop name in case there is only one pop.

	my $rnum = int(rand($#candreads + 1));
	my $rinf = $candreads[$rnum];
	my @rid = split(/__/, $rinf);
	my $chr = $rid[0];
	my $pos = $rid[1];
	splice(@candreads, $rnum, 1); #remove the random read from the array of possible reads, don't want to check it again	
	print("$chr\t$pos\n");
	#loop through bam files and grab a copy of the read if it exists.
	MID:foreach my $line (@barray){
		chomp($line);
		#if row is a pop name, set the current pop to this.
		if ($line =~ /^POP\:/){
			my @temp = split(/\s/, $line);
			$indnum = 1;
			$cpop = $temp[1]; #grab pop name
			if(length($cpop) > 5){ #shorten the pop name if it is longer than 7 characters.
				$cpop = substr($cpop, 0, 5);
			}
			elsif(length($cpop) < 5){
				$cpop = sprintf($cpop, %xd5);
			}
			next MID;
		}
		
		my $schr = `samtools view $line | grep $chr | grep $pos`; #get the lines from this converted bam file are in the same chr and have a match for pos somewhere 

		#if something matched, store data. otherwise do nothing (just less individuals for this locus)
		if ($schr) {
			my @schra = split(/\n/, $schr);
			#for each match, make sure that they actually match by checking the position of the pos match, then pick a random one to save.
			my @trcand; #initialize true candidates array
			#foreach $element (@schra){print("$element\n")}
			INNER:foreach my $pcand (@schra){
				chomp($pcand);
				my @thiscand = split(/\t/,$pcand); #split the candidate
				#check if the target position is the same as the starting BUT NOT the paired end tag position (to get rid of merged tags).
				if($thiscand[3] eq $pos && $thiscand[3] ne $thiscand[6]){
					#print("True match: $tpos, $tchr, $thiscand[3], $thiscand[2]\n");
					push(@trcand, $thiscand[9]); #if the candidate is actually the same tag (starting pos is the same), add to the true candidates array
				}
			}

			#if there is at least one match, save one random result
			if($#trcand >= 0){
				my $keep = int(rand($#trcand + 1));
				my $adj_inum = sprintf("%05d", $indnum); #add leading zeros if necessary to ensure the individual name is the correct length
				my $ind_id = $cpop . $adj_inum . "\t" . $trcand[$keep]; #get the info to store. Individual identifier, tab, then read.
				push(@{$kh{$cpop}{$rinf}}, $ind_id); #save one random result
				$ci{$cpop}{$rinf}++;
				undef(@trcand); #forget trcand array
			}
			else{
				undef(@trcand); #forget trcand array
				#note in %ci if a value isn't already there.
				if (!$ci{$cpop}{$rinf}){
					$ci{$cpop}{$rinf} = 0;
				}
			}
		}
		$indnum++;
		#note in %ci if a value isn't already there.
		if (!$ci{$cpop}{$rinf}){
			$ci{$cpop}{$rinf} = 0;
		}

	}
	$nsaved = $nsaved + 1; #add one to nsaved, signifiying that we have a complete match
}


######################################################################################################################################################################
#print the output information in the correct migrate-n format
open(OUTFILE, ">", $outfile); #open the output connection
#print header info
$npop = scalar keys %kh;
print OUTFILE ("   $npop $nrds  $blist\n", join(" ", ($rlength) x $nrds), "\n");

#print info for each pop.
foreach my $pop (sort keys %kh){
	my %tpop = %{$kh{$pop}}; #get the sub-hash for this pop
	
	#print the number of individual at each locus
	foreach my $locus (sort keys %{$ci{$pop}}){
		print OUTFILE ("$ci{$pop}{$locus}   ");
	}	
	print OUTFILE ("$pop\n");
	foreach my $locus (sort keys %tpop){
		my @dat = @{$tpop{$locus}}; #get the array of individuals with reads
		foreach my $ind (@dat){
			@idat = split(/\t/, $ind); #split individual into id and data
			print OUTFILE (join("",@idat), "\n"); #print the individual followed by the data.
		}
	}
}




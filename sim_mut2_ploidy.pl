#!/usr/bin/perl
$n = $ARGV[0]; #n: gene copies
$N = $ARGV[1]; #N: pop size, individuals (total pop gene copies is 2*N for diploids, N for haploids, for example)
$m = $ARGV[2]; #m: mutation rate
$p = $ARGV[3]; #ploidy
$s = $ARGV[4]; #selfing, yes or no
$os = $ARGV[5]; #os: output style. Divisible by 2 for ms style output, by 3 for raw gene copy info, and/or by 5 for raw tree statistics.
$nts = $ARGV[6]; #$nts: number of pairwise TMRCAs to do for tree stats


use Data::Dumper qw(Dumper);
use List::Util qw(shuffle);


#get output styles
if($os % 2 == 0){$osms = 1}
if($os % 3 == 0){$osgc = 1}
if($os % 5 == 0){$osts = 1}



my @popsize = ($N) x 1; #population size for each gen, gens larger than array assume last array element size

my @starts = (0) x $n; #start gen of each branch
my @stops = () x $n; #stop gen of each branch
my @decends = ( 0 .. $n ); #decendent gcs of each branch, original branches only have original gcs as decendents
my @anc = (0) x $n; #stores random ancestor each gen, -1 when gc has already coalesced
my @ploidy = ($p) x 1; #ploidy size for each gen
my $theta = 4*($popsize[0]/2)*$m; #save this for output
my $init_copies = $n; #save this for output.


if($s eq "yes"){#set ploidy to haploid in each gen if selfing is true
	my $i = 0;
	while ($i <= $#ploidy){
		$ploidy[$i] = 1;
	$i++;
	}
}

#sort intial gene copies into individuals
#intialize array for randomly pulling gene copies from different individuals correctly.
my @sa = (1 .. $ploidy[0]*2);
#load starting gene copies into individuals
my $ns_inds = int(($n/$ploidy[0])+0.5); #get number of individuals represented by copies. Assumes that individuals are maximally grouped, $p copies per ind. Could write this to take individual-copy identity from file or something.
my $i = 0;
while ($i < $n){
	my $k = int($i/$ploidy[0]); #ind to push to
	push(@{$is{$k}}, $i); #add the gene copy ID to the correct individual.
	$i++;
}

my $g = 1;
while ($n > 1) {
	
	#intialize is_added tracker for this gen
	my @is_added = (0) x $#anc;
	
	#print("\n===================================================================================================================\nGen: $g\n");
	#print("Individuals:\n");
	#print Dumper \%is;
	undef $tploidy;
	my $lg = $g - 1;
	if($g > $#ploidy){ #use last ploidy as this ploidy if not specified for this gen.
		$tploidy = $ploidy[$#ploidy];
	}
	elsif($ploidy[$g] != $ploidy[$tg]){ #if ploidy changed and is now poly, remake @sa. Elsif here since ploidy can't have changed if it isn't specified 
		$tploidy = $ploidy[$g];
		if($tploidy != 1){
			@sa = (1 .. $tploidy*2);
		}
	}
	else{
		$tploidy = $ploidy[$g];
	}

	
	######################################################################################################################################
	#pull ancestral gene copies for each current gene copy
	#for haploids
	if($tploidy == 1){
		my $x = 0; #number of gene copies with a drawn ancestor
		while ($x <= $#anc) { #while all copies don't have an ancestor...
			if ($anc[$x] != -1) { #if this copy hasn't already coalesced (-1)
				if ($g <= $#popsize) { #if a popsize is specificied for this gen
					$anc[$x] = int(rand($popsize[$g])); #get a random draw between zero and the pop size (origin individual)
				} else { #otherwise
					$anc[$x] = int(rand($popsize[$#popsize])); #get a random draw between zero and the final pop size
				}
			}
			$x++;
		}
	}
	
	#for polyploids:
	else{ #for polyploids.
		#print("--Polyploid this gen.\n");
		foreach my $ind (keys %is){ #for each individual
			#print ("------------------\nInd: $ind\n");
			#get two parents for each individual
			undef $p1;
			undef $p2;
			if($g <= $#popsize){ #pull a given pop size
				$p1 = int(rand($popsize[$g])); #parent 1
				$p2 = int(rand($popsize[$g])); #parent 2
				while ($p1 == $p2){ #make sure they aren't the same
					$p2 = int(rand($popsize[$g]));
				}
			}
			else{ #no popsize specified, use final specified
			$p1 = int(rand($popsize[$#popsize])); #parent 1
				$p2 = int(rand($popsize[$#popsize])); #parent 2
				while ($p1 == $p2){ #make sure they aren't the same
					$p2 = int(rand($popsize[$#popsize]));
				}
			}
			#print("\tParent 1: $p1\n\tParent 2: $p2\n");
			#get ancestors for each copy in those individuals
			my @gcs = @{$is{$ind}}; #get the gene copies to check for coalescence in this individual
			#print("\tGene copies:\n");
			#print Dumper \@gcs;
			my $tgc = 0;
			#need to see which parent each came from. For ALL copies in the end (including those we aren't interested in), half need to come from each parent
			@sa = shuffle @sa;
			#print("\tShuffle array:\n");
			#print Dumper \@sa;
			#undef $sgc;
			while ($tgc <= $#gcs){ #for each copy in this individual...
				if($anc[$gcs[$tgc]] != -1){ #if this copy hasn't already coalesced...
					#print("\tCopy $gcs[$tgc], index $tgc. Source: ");
					if($sa[$tgc] <= $tploidy){ #check if the copy comes from parent one or two, parent one if the value is less than ploidy
						$sgc = $p1 . "_" . $sa[$tgc]; #source is parent one, chr copy the suffled position.
					}
					else{
						$sgc = $p2 . "_" . ($sa[$tgc] - $tploidy); #source is parent two, chr copy is the shuffled result (if 5 and ploiy is 3, it would be position two, ind two).
					}
					#print("$sgc\n");
					$anc[$gcs[$tgc]] = $sgc; #set the source individual/gene copy
				}
				$tgc++;
			}
		}
	}
		
	my $x = 0;
	my $stop = $#anc;
	undef %is;
	#print("Ancestors:\n");
	#print Dumper \@anc;
	
	#################################################################################################################
	#check for coalescence
	while ($x <= $stop) {
		if ($anc[$x] != -1 && $anc[$x] != -2) { #if hasn't coalesced or is new
			#print("\n=======\nChecking coalescence on copy $x\n");
			my @col = ($x) x 1;
			
			#if polyploid, grab ancestral individual and copy for this copy for comparison. If haploid, no reason to (it's just a number that we access) 
			undef @inf; #explicitly undef any global variables created below.
			undef $xind;
			undef $xcp;
			if($tploidy != 1){ #if polyploid
				@xinf = split(/_/, $anc[$x]); #get info on copy
				$xind = $xinf[0]; #individual source
				$xcp = $xinf[1]; #copy source 
			}
		
			#check to see if the ancestor is the same in each other copy
			my $y = $x+1; #start with the copy 1+ the current
			while ($y <= $#anc) { #for each remaining comparison involving this copy
				
				#if haploid, see if random number is the same:
				if($tploidy == 1){ #if haploid
					if ($anc[$x] == $anc[$y]) { #if they have the same ancestor
						#print("\tCoalesced with copy $y.\n");
						push(@col, $y); #add the id of the copy with the same ancestor
						$anc[$y] = -1; $stops[$y] = $g; #mark that the other copy coalesced, put the gen that it stopped
					}
				}
				
				#if polyploid, see if ancestral individual and copy are the same
				else{ #if polyploid, same deal but split the source ind and source copy and check each.
					my @yinf = split(/_/, $anc[$y]);
					my $yind = $yinf[0];
					if($yind == $xind){ #if y is sourced from the same individual as x
						my $ycp = $yinf[1];
						if($xcp == $ycp){ #coalesced, save info. as above.							
							#print("\tCoalesced with copy $y.\n");
							push(@col, $y);
							$anc[$y] = -1; $stops[$y] = $g;
						}
						else{ #otherwise note that another copy, ID y, came from the same individual (didn't coalesce though). This is used in the next gen.
							#print("\tSame source individual as copy $y.\n");
							if ($is_added[$y] >= 0){ #if this hasn't been previously added to an individual
								push(@{$is{$x}}, $y); #add y to an individual ID'd by loop iteration.
								@is_added[$y] = -1;
							}
						}
					}
				}	
				$y++;
			}
			
			
			
			################################################################################################################
			#do book-keeping for next generation by noting coalescences, creating a "new" copy for the ancestral source for each coalescent event
			if ($#col > 0) { #if this copy coalesced
				$anc[$x] = -1; $stops[$x] = $g; #note that it coalesced, add gen when it happened
				$n = $n - $#col; #change the number of gene copies to the current minus the number which coalesced
				push(@anc,-2); $starts[$#anc] = $g; #add info for the "new" copy/node, set the gen where it started
				push(@{$is{$x}}, $#anc); #add the new copy to an individual ID'd by loop iteration.
				$stops[$#anc] = ""; @{$decends[$#anc]} = @col; #add a filler spot for when the new copy/node ends, store descendant info for all the copies that coalesced into new copy.
			}
			if ($anc[$x] != -1 && $anc[$x] != -2 && $is_added[$x] >= 0){ #if hasn't already coalesced or is new (including this gen)
				$is_added[$x] = -1;
				unshift(@{$is{$x}}, $x); #add/put this copy into an individual, ID'd by loop iteration.
			}
		}
		$x++;
	}
	$g++;
}

#######################################################################################################################################################################
#print outputs


#ms and gc
if($osms == 1 | $osgc == 1){
	@mutations = (0) x ($#anc+1); #intialize mutation array

	#add mutations to line, print out results
	$x = 0;
	while ($x <= $#mutations) { #for each copy
		$y = $starts[$x]; #y is the starting gen
		while ($y < $stops[$x]) { #for each gen that the branch existed...
			if (rand(1) <= $m) { #check to see if a mutation occured
				$mutations[$x]++; #add one to mutaiton counter if so
				@init_array = (0) x ($#mutations + 1); #initialize mutation array
				$init_array[$x] = 1; #mark which branch the mutation occured in
				$tpos = rand(1); #get the mutation position
				while(exists $muthash{$tpos}){ #ensure no duplicate positions
					$tpos = rand(1);
				}
				push(@{$muthash{$tpos}}, @init_array); #add this array to the muthash.
			}
			$y++;
		}
		#print results: branch ID, start gen, end gen, number of mutations, which branches are the ancestor
		
		#gc	
		if($osgc == 1){
			print "GC$x\t$starts[$x]\t$stops[$x]\t$mutations[$x]\t" . join(',',@{$decends[$x]}) . "\n";
		}
		$x++;
	}
	#print Dumper \%muthash;

	#ms
	if($osms == 1){
		foreach $mut (sort keys %muthash){
			@tmut = @{$muthash{$mut}};
			$x = $#tmut;
			#print("\n\nThis mut: $mut\n");
			while ($x >= 0){
				if ($tmut[$x] == 1){
					#print("\tCopy $x decendats:\t");
					foreach $d (@{$decends[$x]}){
							#print("$d\t");
							$tmut[$d] = 1;
					}
					#print("\n");
				}
				$x = $x - 1;
			}
			#print Dumper \@tmut;
			$i = 0;
			while($i <= $#tmut){
					push(@{$oAoA[$i]}, $tmut[$i]);
					$i++;
			}
		}
		#print("===========================================\nEnding AoA:\n");
		#print Dumper \@oAoA;

		#print out ms style output
		$segsites = scalar keys %muthash;
		#print("ms $init_copies 1 -t $theta\n\n\n"); #don't want to print this if doing many sims, only needed at the top.
		print("//" . "\n");
		print("segsites: $segsites\n");
		if($segsites > 0){
			print("\npositions: ");
			if(%muthash){
				foreach $mut (sort keys %muthash){
					print("$mut ");
				}
			}
			print("\n");
			$i = 0;
			while ($i <= ($init_copies - 1)){
				@tgc = @{$oAoA[$i]};
				print(join("",@tgc) . "\n");
				$i++;
			}
		}
		print("\n\n");
	}
}

#ts
if($osts == 1){
	my $i = $init_copies; #i starts at the first non-terminal branch
	while($i <= $#anc){
		my @tds = (0) x ($i + 1); #set the decendant array, contains only the possible decendants and this branch
		#print("\n-------------------------------------------------------------------\nCopy $i decendants:\n");
		my $x = $i; #set current position in array, move down to figure out decends
		$tds[$#tds] = 1;
		#print Dumper \@tds;
		while ($x >= ($init_copies - 1)){ #for each non-terminal branch less than the starting place
			if ($tds[$x] == 1){
				#print("\tCopy $x decendats:\t");
				foreach my $d (@{$decends[$x]}){
					#print("$d\t");
					$tds[$d] = 1;
				}
				#print("\n");
			}
			$x = $x - 1;
		}
		#print Dumper \@tds;
		push(@{$full_decends[$i]}, @tds);
		$i++;
	}
	#print Dumper \@full_decends;

	#calculate mean and var for time to ancestor for each pair of terminal gene copies
	my $x = 0;
	my $tot_tmrca = 0;
	my $n_pwc = 0;
	while ($x < $init_copies && $n_pwc <= $nts){
		#print("\n================\ntrmcas for $x :\n");
		my $y = $x + 1; #set comparison pair
		COMPS: while ($y < $init_copies && $n_pwc <= $nts){
			#print("gc $y :\n");
			#find the least ancestral branch with a @full_decends element with both copies listed
			my $z = $init_copies; #starting check element, first non-terminal branch
			while($z <= $#anc){
				#print("\n$z");
				my @tbrch = @{$full_decends[$z]}; #get the array for this branch
				#print Dumper \@tbrch;
				if($tbrch[$x] == 1 & $tbrch[$y] == 1){ #if both are listed as decendants
					#print("\tB: $z\t tmrca: $starts[$z]\n");
					push(@tmrca, $starts[$z]); #add tmrca to array
					$tot_tmrca += $starts[$z]; #add tmrca to total
					$y = $y + 1;
					$n_pwc++;
					next COMPS;
				}
				$z++;
			}
			print("\nERROR: no decendant detected for gcs $x and $y , run $tr of $r for input line $this_run_set !\n");
			#print("z: $z out of $#anc . $#full_decends tbrch : \n");
			#print Dumper \@{$full_decends[$init_copies]};
			die;
		}
		$x++;
	}
	my $ave_tmrca = $tot_tmrca/($#tmrca + 1);
	my $SSE = 0;
	foreach $gc (@tmrca){ 
		$SSE += ($gc - $ave_tmrca)**2; #squared error for this copy
	}
	my $var_tmrca = $SSE/($#tmrca + 1); #variance of tmrca
	
	#get total tree length, note that time to complete coalescence is $g - 1.
	$i = 0;
	my $TTL = 0;
	while($i <= $#stops - 1){ #last branch has no stop time.
		$TTL += $stops[$i] - $starts[$i];
		$i++;
	}
	$g -= 1;
	#print results
	print ("TreeStats:\t$TTL\t$g\t$ave_tmrca\t$var_tmrca\n");
}

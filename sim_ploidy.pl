#!/usr/bin/perl
#conducts coalescent simulations with no assumption of low n/N ratio or haploidy. Output is in MS format if no output format
#arugment is given, can be quickly analyzed via the sample_stats utility packaged with ms (output | /path/to/ms/directory/sample_stats)


if($#ARGV == 4 | $#ARGV == 5){
	$n = $ARGV[0]; #n: gene copies
	$N = $ARGV[1]; #N: pop size
	$m = $ARGV[2]; #m: mutation rate
	$p = $ARGV[3]; #p: ploidy
	$s = $ARGV[4]; #s: selfing, 0 = no, 1 = yes.
	if($#ARGV = 5){$os = $ARGV[5];} #os: output style. Empty for ms, 1 for only times and number of mutations, any other number for both.
}
else{die("Usage: sim_ploidy.pl	n_gene_copies	popsize	mutation_rate	ploidy	selfing	optional_output_styles\n");}
use List::Util qw(shuffle);
use Data::Dumper qw(Dumper);

#intialize stuff
@popsize = ($N) x 1; #population size for each gen, gens larger than array assume last array element size
@starts = (0) x $n; #start gen of each branch
@stops = () x $n; #stop gen of each branch
@decends = ( 0 .. $n ); #decendent gcs of each branch, original branches only have original gcs as decendents
@anc = (0) x $n; #stores random ancestor each gen, -1 when gc has already coalesced
@ploidy = ($p) x 1; #ploidy size for each gen
$init_copies = $n; #save this for output.
$theta = 4*$popsize[0]*$m; #save this for output

if($s == 1){#set ploidy to haploid in each gen if selfing is true
	$i = 0;
	while ($i <= $#ploidy){
		$ploidy[$i] = 1;
	$i++;
	}
}

#sort intial gene copies into individuals
#intialize array for randomly pulling gene copies from different individuals correctly.
@sa = (1 .. $ploidy[0]*2);
#load starting gene copies into individuals
$ns_inds = int(($n/$ploidy[0])+0.5); #get number of individuals represented by copies. Assumes that individuals are maximally grouped, $p copies per ind. Could write this to take individual-copy identity from file or something.
$i = 0;
while ($i < $n){
	$k = int($i/$ploidy[0]); #ind to push to
	push(@{$is{$k}}, $i); #add the gene copy ID to the correct individual.
	$i++;
}
#print("Starting individuals:\n");
#print Dumper \%is;




#run coalescent simulation
$g = 1;
while ($n > 1) {
	#print("-----------------------------------------------------------\nGen $g:\n\tPloidy: ");
	
	#intialize is_added tracker for this gen
	undef @is_added;
	@is_added = (0) x $#anc;

	#get ploidy for this gen
	$lg = $g - 1;
	if($g > $#ploidy){ #use last ploidy as this ploidy if not specified for this gen.
		$tploidy = $ploidy[$#ploidy];
        }
        elsif($ploidy[$g] != $ploidy[$tg]){ #if ploidy changed and is now poly, remake @sa. Else if here since ploidy can't have changed if it isn't specified $
		$tploidy = $ploidy[$g];
                if($tploidy != 1){
			@sa = (1 .. $tploidy*2);
		}
        }
        else{
                $tploidy = $ploidy[$g];
	}


	#print("$tploidy\n\tPop Size: ");
	#if($g <= $#popsize){
	#	print("$popsize[$g]\n");
	#}
	#else{
	#	print("$popsize[$#popsize]\n");
	#}

	#for haploids:	
	if($tploidy == 1){
		#print("--Haploid this gen.\n");
		$x = 0; #number of gene copies with a drawn ancestor
		while ($x <= $#anc) { #while all copies don't have an ancestor...
			if ($anc[$x] >= 0) { #if this copy hasn't already coalesced
				if ($g <= $#popsize) { #if a popsize is specificied for this gen
					$anc[$x] = int(rand($popsize[$g])); #get a random draw between zero and the pop size (origin individual)
				} else { #otherwise
					$anc[$x] = int(rand($popsize[$#popsize])); #get a random draw between zero and the final pop size
				}				
			}
			$x++;
		}
	}

	else{ #for polyploids.
		#print("--Polyploid this gen.\n");
		foreach $ind (keys %is){ #for each individual
			#print("-----\nInd: $ind\n");
			#get two parents for each individual
			if($g <= $#popsize){
				$p1 = int(rand($popsize[$g])); #parent 1
				$p2 = int(rand($popsize[$g])); #parent 2
				while ($p1 == $p2){ #make sure they aren't the same
					$p2 = int(rand($popsize[$g]));
				}
			}
			else{
			        $p1 = int(rand($popsize[$#popsize])); #parent 1
                                $p2 = int(rand($popsize[$#popsize])); #parent 2
                                while ($p1 == $p2){ #make sure they aren't the same
                                        $p2 = int(rand($popsize[$#popsize]));
                                }
                        }
			#print("\tParent 1: $p1\n\tParent 2: $p2\n");
			#get ancestors for each copy in those individuals
			@gcs = @{$is{$ind}}; #get the gene copies to check for coalescence in this individual
			#print("\tGene copies:\n");
			#print Dumper \@gcs;
			$tgc = 0;
			#need to see which parent each came from. For ALL copies in the end (including those we aren't interested in), half need to come from each parent
			@sa = shuffle @sa;
			#print("\tShuffle array:\n");
			#print Dumper \@sa;
			while ($tgc <= $#gcs){ #for each copy in this individual...
				if($anc[$gcs[$tgc]] >= 0){ #if this copy hasn't already coalesced...
					#print("\tCopy $tgc source: ");
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
	
	#print("Anc:\n");
	#print Dumper \@anc;
	$x = 0;
        $stop = $#anc;
	undef %is; #remove the old individual hash, need to make a new one for the next gen back to use
        #compare each gene copy to every other gene copy to see if coalesced
        while ($x <= $stop) {
		@col = ($x) x 1;
		if ($anc[$x] >= 0) { #if hasn't coalesced before this gen
			#print("=======\nChecking coalescence on copy $x\n");
			if($tploidy != 1){ #if polyploid
				@xinf = split(/_/, $anc[$x]); #get info on copy
				$xind = $xinf[0]; #individual source
				$xcp = $xinf[1]; #copy source 
			}
                	$y = $x+1; #start with the copy 1+ the current
                        while ($y <= $#anc) { #for each remaining comparison involving this copy
                        	if($tploidy == 1){ #if haploid
					if ($anc[$x] == $anc[$y]) { #if they have the same ancestor
						#print("\tCoalesced with copy $y.\n");
                        		       	push(@col, $y); #add the id of the copy with the same ancestor
                        		        $anc[$y] = -1; $stops[$y] = $g; #mark that the other copy coalesced, put the gen that it stopped
                        		}
				}
				else{ #if polyploid, same deal but split the source ind and source copy and check each.
					@yinf = split(/_/, $anc[$y]);
					$yind = $yinf[0];
					if($yind == $xind){ #if y is sourced from the same individual as x
						$ycp = $yinf[1];
						if($xcp == $ycp){ #coalesced, save info. as above.							
							#print("\tCoalesced with copy $y.\n");
							push(@col, $y);
							$anc[$y] = -1; $stops[$y] = $g;
						}
						else{ #otherwise note that another copy, ID y, came from the same individual (didn't coalesce though).
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
                }
		#print("n col $#col\n");
                if ($#col > 0) { #if this copy coalesced this gen
			$anc[$x] = -1; $stops[$x] = $g; #note that it coalesced, add gen when it happened
                        $n = $n - $#col; #change the number of gene copies to the current minus the number which coalesced
                        push(@anc,0); $starts[$#anc] = $g; #add info for the "new" copy/node, set the gen where it started
                        push(@{$is{$x}}, $#anc); #add the new copy to an individual ID'd by loop iteration.
			$stops[$#anc] = ""; @{$decends[$#anc]} = @col; #add a filler spot for when the new copy/node ends, store descendant info f$
		}
		if ($anc[$x] >= 0 && $is_added[$x] >= 0){ #if hasn't already coalesed (including this gen)
			$is_added[$x] = -1;
			unshift(@{$is{$x}}, $x); #add/put this copy into an individual, ID'd by loop iteration.
		}
		$x++;
	}

	#print("Individuals for next gen:\n");
	#print Dumper \%is;
	#print(" Number of copies remaining: $n\n");
	$g++;
}

@mutations = (0) x ($#anc+1); #intialize mutation array

#add mutations, store in muthash for propagation later.
$x = 0;
while ($x <= $#mutations) { #for each copy
	$y = $starts[$x]; #y is the starting gen
	while ($y < $stops[$x]) { #for each gen that the branch existed...
		if (rand(1) <= $m) { #check to see if a mutation occured
			@init_array = (0) x ($#mutations + 1); #initialize mutation array
			$init_array[$x] = 1; #mark which branch the mutation occured in
			$mutations[$x]++; #add one to mutaiton counter if so. Can remove this later.
			$tpos = rand(1); #get the mutation position
			while(exists $muthash{$tpos}){ #ensure no duplicate positions
				$tpos = rand(1);
			}
			push(@{$muthash{$tpos}}, @init_array); #add this array to the muthash.
		}
		$y++;
	}
	#print results: branch ID, start gen, end gen, number of mutations, which branches are the ancestor
        if($os){print "$x\t$starts[$x]\t$stops[$x]\t$mutations[$x]\t" . join(',',@{$decends[$x]}) . "\n";}
	$x++;
}


#propagate mutations, create ms style output. Pseudo:
#foreach key in muthash
#	loop through each gene copy, starting with the most ancestral.
#	if 1 (mutation present), set to 1 in each decends copy.
#by stepping down from most ancestral to most derived, this will ensure proper propagation. If it seems to produce errors, can loop through until the array stops
#changing.

#print("===========================================\nPropagating mutations.\n");
#print Dumper \%muthash;

if($os != 1){
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
	print("ms $init_copies 1 -t $theta\n\n\n");
	print("//" . "\n");
	print("segsites: $segsites");
	if(%muthash){
		print ("\npositions: ");
		foreach $mut (sort keys %muthash){
			print("$mut ");
		}
	}
	print("\n");
	$i = 0;
	if(%muthash){
		while ($i <= ($init_copies - 1)){
			@tgc = @{$oAoA[$i]};
			print(join("",@tgc) . "\n");
			$i++;
		}
	}
}

#!/usr/bin/perl
#conducts coalescent simulations with no assumption of low n/N ratio or haploidy. Output can be given in MS format,
#which can be quickly analyzed via the sample_stats utility packaged with ms (output | /path/to/ms/directory/sample_stats).
#Can also give gene copy information (branch start and end generations and number of mutations) or tree stats
#(TTL, TMRCA for all copies, average TMRCA for each pairwise comparison of copies, and variance of the latter).


#run info file format: <gene_copes><pop_size_file><mutation_rate><ploidy_file><selfing><n_runs_with_these_settings>
#tokens delimited by whitespace

if($#ARGV == 2 | $#ARGV == 3){
	$runinf = $ARGV[0]; #runinf: file containing run info
	$output_prefix = $ARGV[1];
	$os = $ARGV[2]; #os: output style. Empty for ms, 1 for only times and number of mutations, 2 for tree stats.
	$nts = $ARGV[3]; #$nts: number of pairwise TMRCAs to do for tree stats
}
else{die("Usage: sim_ploidy_multirun.pl	run_info_file output_prefix output_style opt_num_pairwise_TMRCAs\nrun info file format, tokens whitespace delimited:\n<gene_copes><pop_size_file><mutation_rate><ploidy_file><selfing><n_runs_with_these_settings>\nOutput Styles: 2: ms. 3: gene copy. 5: tree stats. Multiple of these returns both/all.\n");}
use List::Util qw(shuffle);
use Data::Dumper qw(Dumper);

if($os % 2 == 0){$osms = 1}
if($os % 3 == 0){$osgc = 1}
if($os % 5 == 0){$osts = 1}

if(not($osms) & not($osgc) & not($osts)){
	die("Output style not valid. Must be divisible by 2, 3, and/or 5.\n");
}

#loop through info file, do runs listed
open (RUNINFO, "<", $runinf) or die("Can't find run info file $runinf\n");
$this_run_set = 1;
while (<RUNINFO>){
	chomp($_);
	

	#print("-" x 60, "\n");
	#{
	# no strict 'refs';
	#
	#	foreach my $entry ( sort keys %main:: )
	#	{
	#		print "$entry\n";
	#	}
	#}


	my @meta = split(/\s/,$_);
	#foreach $element(@meta){print("$element\n");}
	if(not($#meta == 5 | $#meta == 6)){die("Run info format incorrect at line $this_run_set.\n");}
	my $n = $meta[0]; #n: gene copies
	my $N = $meta[1]; #N: pop size
	my $m = $meta[2]; #m: mutation rate
	my $p = $meta[3]; #p: ploidy
	my $s = $meta[4]; #s: selfing
	my $r = $meta[5]; #r: number of reps


	my $init_copies = $n; #save this for output.
	my @popsize = ($N) x 1; #population size for each gen, gens larger than array assume last array element size
	my @ploidy = ($p) x 1; #ploidy size for each gen
	my $theta = 4*($popsize[0]/2)*$m; #save this for output, pop size in gene copies as a diploid per ms.

	my $tout = $output_prefix . "_" . $this_run_set; 
	
	if($osms){
		my $msnm = "ms_" . $tout;
		open (MSOUT, ">", $msnm);
		print MSOUT ("ms $init_copies $r -t $theta\n\n\n");
	}

	if($osgc){
		my $gcnm = "gc_" . $tout;
		open (GCOUT, ">", $gcnm);
	}

	if($osts){
		my $tsnm = "ts_" . $tout;
		open (TSOUT, ">", $tsnm);
	}

	my $tr = 1;

	while ($tr <= $r){
		my $n = $init_copies;
	

		#intialize stuff
		my @starts = (0) x $n; #start gen of each branch
		my @stops = () x $n; #stop gen of each branch
		my @decends = ( 0 .. $n ); #decendent gcs of each branch, original branches only have original gcs as decendents
		my @anc = (0) x $n; #stores random ancestor each gen, -1 when gc has already coalesced


		if($s == 1){#set ploidy to haploid in each gen if selfing is true
			$i = 0;
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
		$i = 0;
		my %is;
		while ($i < $n){
			my $k = int($i/$ploidy[0]); #ind to push to
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
			my @is_added = (0) x ($#anc + 1);
			#print("$#anc\n");
			#print Dumper \@is_added;

			#get ploidy for this gen
			my $lg = $g - 1;
			if($g > $#ploidy){ #use last ploidy as this ploidy if not specified for this gen.
				$tploidy = $ploidy[$#ploidy];
				}
			elsif($ploidy[$g] != $ploidy[$g - 1]){ #if ploidy changed and is now poly, remake @sa. Else if here since ploidy can't have changed if it isn't specified $
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
				foreach my  $ind (keys %is){ #for each individual
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
					my @gcs = @{$is{$ind}}; #get the gene copies to check for coalescence in this individual
					#print("\tGene copies:\n");
					#print Dumper \@gcs;
					my $tgc = 0;
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
			my $stop = $#anc;
			undef %is; #remove the old individual hash, need to make a new one for the next gen back to use
			#compare each gene copy to every other gene copy to see if coalesced
			while ($x <= $stop) {
				my @col = ($x) x 1;
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
										$is_added[$y] = -1;
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


		
		if($osms | $osgc){
			my @mutations = (0) x ($#anc+1); #intialize mutation array
			my %muthash;
			#add mutations, store in muthash for propagation later.
			$x = 0;

			while ($x <= $#mutations) { #for each copy
				my $y = $starts[$x]; #y is the starting gen
				while ($y < $stops[$x]) { #for each gen that the branch existed...
					if (rand(1) <= $m) { #check to see if a mutation occured
						my @init_array = (0) x ($#mutations + 1); #initialize mutation array
						$init_array[$x] = 1; #mark which branch the mutation occured in
						$mutations[$x]++; #add one to mutaiton counter if so. Can remove this later.
						my $tpos = rand(1); #get the mutation position
						while(exists $muthash{$tpos}){ #ensure no duplicate positions
							$tpos = rand(1);
						}
						push(@{$muthash{$tpos}}, @init_array); #add this array to the muthash.
					}
					$y++;
				}
				#print results: branch ID, start gen, end gen, number of mutations, which branches are the ancestor
				if($osgc){print GCOUT ("$x\t$starts[$x]\t$stops[$x]\t$mutations[$x]\t" . join(',',@{$decends[$x]}) . "\n");}
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

			if($osms){
				my @oAoA;
				foreach my $mut (sort keys %muthash){
					my @tmut = @{$muthash{$mut}};
					$x = $#tmut;
					#print("\n\nThis mut: $mut\n");
					while ($x >= 0){
						if ($tmut[$x] == 1){
							#print("\tCopy $x decendats:\t");
							foreach my  $d (@{$decends[$x]}){
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
				my $segsites = scalar keys %muthash;
				print MSOUT ("//" . "\n");
				print MSOUT ("segsites: $segsites");
				if(%muthash){
					print MSOUT ("\npositions: ");
					foreach my $mut (sort keys %muthash){
						print MSOUT ("$mut ");
					}
				}
				print MSOUT ("\n");
				if(scalar keys %muthash >= 1){
					$i = 0;
					while ($i <= ($init_copies - 1)){
						my @tgc = @{$oAoA[$i]};
						print MSOUT (join("",@tgc) . "\n");
						$i++;
					}
				}
			}
		} #ends doing mutation stuff if not $osms or $osgc
		


		if($osts){ #$osts, do tree stats.
			undef @full_decends;
			undef @tmrca;
			my @full_decends;
			my @tmrca;
			#for debugging:		
			#$x = 0;
			#while($x <= $#anc){
			#	print ("$x\t$starts[$x]\t$stops[$x]\t" . join(',',@{$decends[$x]}) . "\n");
			#	$x++;
			#}
			#print ("\n");
			

				
			#first need to propagete full decends to each branch
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
			$x = 0;
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
			print TSOUT ("$TTL\t$g\t$ave_tmrca\t$var_tmrca\n");
		}	
		
		if($osms){print MSOUT ("\n\n");}
		if($osgc){print GCOUT ("\n\n");}
		$tr++;
	} #end stuff to do each run set
	$this_run_set++; #agument run set counter
} 

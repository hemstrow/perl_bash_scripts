#!/usr/bin/perl
$n = $ARGV[0]; #n: gene copies
$N = $ARGV[1]; #N: pop size, haploid
$m = $ARGV[2]; #m: mutation rate
$os = $ARGV[3]; #os: output style. Divisible by 2 for ms style output, by 3 for raw gene copy info, and/or by 5 for raw tree statistics.
$nts = $ARGV[4]; #$nts: number of pairwise TMRCAs to do for tree stats

#get output styles
if($os % 2 == 0){$osms = 1}
if($os % 3 == 0){$osgc = 1}
if($os % 5 == 0){$osts = 1}



@popsize = ($N) x 1; #population size for each gen, gens larger than array assume last array element size

@starts = (0) x $n; #start gen of each branch
@stops = () x $n; #stop gen of each branch
@decends = ( 0 .. $n ); #decendent gcs of each branch, original branches only have original gcs as decendents
@anc = (0) x $n; #stores random ancestor each gen, -1 when gc has already coalesced
$theta = 4*($popsize[0]/2)*$m; #save this for output
$init_copies = $n; #save this for output.
use Data::Dumper qw(Dumper);

$g = 1;
while ($n > 1) {
	$x = 0; #number of gene copies with a drawn ancestor
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

	$x = 0;
	$stop = $#anc;
	#compare each gene copy to every other gene copy to see if coalesce
	while ($x <= $stop) {
		if ($anc[$x] >= 0) { #if hasn't coalesced or is new
			$y = $x+1; #start with the copy 1+ the current
			@col = ($x) x 1;
			my @storage = ($anc[$x]);
			while ($y <= $#anc) { #for each remaining comparison involving this copy
				if ($anc[$x] == $anc[$y]) { #if they have the same ancestor
					push(@col, $y); #add the id of the copy with the same ancestor
					$anc[$y] = -1; $stops[$y] = $g; #mark that the other copy coalesced, put the gen that it stopped
				}
				$y++;
			}
			if ($#col > 0) { #if this copy coalesced
				$anc[$x] = -1; $stops[$x] = $g; #note that it coalesced, add gen when it happened
				$n = $n - $#col; #change the number of gene copies to the current minus the number which coalesced
				push(@anc,-2); $starts[$#anc] = $g; #add info for the "new" copy/node, set the gen where it started 
				$stops[$#anc] = ""; @{$decends[$#anc]} = @col; #add a filler spot for when the new copy/node ends, store descendant info for all the copies that coalesced into new copy.
			}
		}
		$x++;
	}
	$g++;
}

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
		
		if($osgc == 1){
			print "GC$x\t$starts[$x]\t$stops[$x]\t$mutations[$x]\t" . join(',',@{$decends[$x]}) . "\n";
		}
		$x++;
	}
	#print Dumper \%muthash;


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

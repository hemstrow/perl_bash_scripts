#!/usr/bin/perl
$n = $ARGV[0]; #n: gene copies
$N = $ARGV[1]; #N: pop size, haploid
$m = $ARGV[2]; #m: mutation rate

@popsize = ($N) x 1; #population size for each gen, gens larger than array assume last array element size

@starts = (0) x $n; #start gen of each branch
@stops = () x $n; #stop gen of each branch
@decends = ( 0 .. $n ); #decendent gcs of each branch, original branches only have original gcs as decendents
@anc = (0) x $n; #stores random ancestor each gen, -1 when gc has already coalesced

$g = 1;
while ($n > 1) {
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

	$x = 0;
	$stop = $#anc;
	#compare each gene copy to every other gene copy to see if coalesce
	while ($x <= $stop) {
		if ($anc[$x] >= 0) { #if hasn't coalesced
			$y = $x+1; #start with the copy 1+ the current
			@col = ($x) x 1;
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
				push(@anc,0); $starts[$#anc] = $g; #add info for the "new" copy/node, set the gen where it started 
				$stops[$#anc] = ""; @{$decends[$#anc]} = @col; #add a filler spot for when the new copy/node ends, store descendant info for all the copies that coalesced into new copy.
			}
		}
		$x++;
	}
	$g++;
}

@mutations = (0) x ($#anc+1); #intialize mutation array

#add mutations to line, print out results
$x = 0;
while ($x <= $#mutations) { #for each copy
	$y = $starts[$x]; #y is the starting gen
	while ($y < $stops[$x]) { #for each gen that the branch existed...
		if (rand(1) <= $m) { #check to see if a mutation occured
			$mutations[$x]++; #add one to mutaiton counter if so
		}
		$y++;
	}
	#print results: branch ID, start gen, end gen, number of mutations, which branches are the ancestor
	print "$x\t$starts[$x]\t$stops[$x]\t$mutations[$x]\t" . join(',',@{$decends[$x]}) . "\n";
	$x++;
}




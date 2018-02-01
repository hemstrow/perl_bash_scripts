#!/usr/bin/perl
#get inputs if 3 inputs, else die
use Data::Dumper qw(Dumper);
use List::Util qw/shuffle/;

if($#ARGV == 3){
        $ssam = $ARGV[0];  #get merged SAM file
	$outfile = $ARGV[1]; #output filename
	$sigma = $ARGV[2]; #sigma accuracy of pvalue
	$max_3plus = $ARGV[3]; #maximun number of individuals with more than three copies detected.
	#$Npermutations = $ARGV[2]; #number of permutations
}
else{die("Usage: ./HWE_from_sam.pl merged.sam outfile p-value_accuracy_in_sigma number_inds_with_3plus_alleles_before_vio\n");}

open(SAM, "<", $ssam) or die("Couldn't find merged sam file.\n");



%h;
%loc;
%vloci;
%tph;
#loop though sam, store reads per loci per ind in hash. Also store #seqs per loci (more than two = reject).
print("Looping through input...\n");
while (<SAM>){
	$line = $_;
	chomp($line);
	@read = split(/\t/, $line); #all of the info from the line
	@dat = split(/;/, $read[0]); #just the sample ID stuff
	#foreach $element(@dat){print("$element\t");}
	#print("\n");
	$sampID = $dat[0]; #sample name
	$chr = $read[2]; #chr/lg
	$pos = $read[3]; #position
	$locus = $chr . "__" . $pos; #get succinct sequence locaiton identifier
	$seq = $read[9]; #sequence
	#print("$locus\t$seq\n");
	$h{$sampID}{$locus}{$seq} = 0; #set seq as in sample
	if(exists($loc{$locus}{$seq}) == 0){
		$loc{$locus}{$seq} = 0; #add sequence to locus
		$tph{$locus} = 0; #add locus to hash which keeps track of number of times an individual is discovered to have more than three alleles.
	}
	#foreach $element (@{$loc{$locus}}){print("$element\n");}
	#print("\n\n\n");
}
close (SAM);

open (OUTFILE, ">", $outfile) or die("Can't open outfile.\n");;

print("Getting genotypes for each individual at each locus.\n");
#get genotypes for each locus

#print Dumper \%h;
OUTER:foreach $ind (sort keys %h){
	#print("-------------------------------------------------------------------------------------------\n$ind\n");
	%ahash = %{$h{$ind}};
	#print Dumper \%ahash;
	MID:foreach $locus (sort keys %ahash){
		#print("\n$locus\n============\n");
		if (exists($vloci{$locus})){ #if this locus has already been found to violate, just skip it.
			#print("Already added to violating locus list...\n");
			next MID;
		}
		%bhash = %{$ahash{$locus}};
		#print Dumper \%bhash;
		$nalleles = scalar keys %bhash; #the number of different reads observed at this locus
		#print("number of alleles: $nalleles\n");
		if($nalleles > 2){ #if there were more than two alleles
			#print("More than two alleles in an individual at $locus, adding to violating list.\n");
			#foreach $element (sort keys %{$h{$ind{$locus}}}){print("$h{$ind}{$locus}{$element}\t");}
			#print("\n");
			$tph{$locus}++;
			if($tph{$locus} >= $max_3plus){
                		@l = split(/__/, $locus); #recover the locus chr and position
				print OUTFILE ("$l[0]\t$l[1]\t-$tph{$locus}\n"); #print the locus chr and position to the output
                        	$vloci{$locus} = 1; #add the locus to the violation loci list to save computation if we find this locus again during read in.
                        	delete $loc{$locus}; #remove this locus from the locus list
                        }
			
			#skip this individual for HWE calcs.
			next MID;
                }
		@a;
		foreach $seq (sort keys %bhash){
			#print("$seq\n");
			if ($nalleles == 1){ #if just one allele, homozygote for that
				#print("1 allele: $seq.\n");
				$loc{$locus}{$seq}++; #add one to this allele's homozygote count
				next MID;
			}
			push(@a, $seq);	#add this allele to the list
		}
		#print("2 alleles:");
		$het = $a[0] . "_" . $a[1]; #concat the strings for saving as het
		#print("$het\n");
		$loc{$locus}{$het}++; #add one to this heterozygote count
	}
}

#print Dumper \%loc;

if($sigma == 0){
	print("Sigma is 0: skipping HWE calculations.\n");
	close OUTFILE;
	exit;
}

$Npermutations = (2.576/(2*$sigma))**2;
$Npermutations = int($Npermutations + 0.5);
print("Doing $Npermutations permutations to achieve $sigma accuracy.\n");

#calc HWE divergence for each.
print("Checking HWE for each locus...\n");
$nloci = scalar keys %loc;
$prog_count = 1;
$vio_count = scalar keys %vloci;
$n_tot_loci = $nloci + $vio_count;
foreach $locus (sort keys %loc){
	#print("-------------------------------------------------------------------------------------------\nLocus: $locus\n");
	if($prog_count % 100 == 0){print("$prog_count out of $nloci loci complete.\n");}
	$prog_count++;
	undef %geno;
	undef %allele;
	undef @Ai;
	%allele;
	@Ai;
	%geno;
	
	
	#redefine obs genotype counts for convieneince.
	%geno = %{$loc{$locus}};
	#print Dumper \%geno;
	$nalleles = scalar keys %geno;
	if($nalleles == 1){
		#print("Only one allele detected, skipping HWE calculation.\n");
		next;
	}	
	
	#get obs allele counts and total
	$n = 0;
	foreach $genotype (sort keys %geno){
		@a = split(/_/, $genotype);
		chomp(@a);
		$n += $geno{$genotype};
		foreach $a1(@a){
			#print("Observed allele: $a1\t");
			if($#a == 0){ #homozygotes
				$allele{$a1} = $allele{$a1} + (2*$geno{$genotype});
			}
			else{ #heterozygotes
				$allele{$a1} = $allele{$a1} + $geno{$genotype}; 
			}
		}
		#print("\n");
	}
	
	#get array of i,j for downstream processing
	foreach $a1 (sort keys %allele){
		$i = 1;
		while($i <= $allele{$a1}){ #for the number of observed instances of $a1
			push(@Ai, $a1); #add that many $a1 alleles to @Ai
			$i++;
		}
	}
	#print("Ai: \t");
	#print Dumper \@Ai;
	#print("\n");
	#print("Observed genotype counts:\n");
	#foreach $element (sort keys %geno){print("\n$element : $geno{$element}\n");}
	#print("Observed allele counts:\n");
	#print Dumper \%allele;

	#check HWE
	
	#function to get factorials:
	sub fact{
		my $input = $_[0];
		$i = 1;
		$factorial = 1;
		while ($i < $input){
			$i++;
			$factorial *= $i;
		}
		return ($factorial)
	}

	#get n factorial, used in all P calcs
	$nfact = fact($n);
	$twon = 2*$n;
	$twonfact = fact($twon);
	
	#get the top factorial for P(g), which doesn't change
	$topfact = 1;
        foreach $a1 (sort keys %allele){
                $topfact *= fact($allele{$a1});
        }
        $topfact *= $nfact;

	#print("\n Stats:\n\t2n:$twon\n\t2n!:$twonfact\n\tTop Factorial Product: $topfact\n");

	#check using monte carlo method, as described in Guo and Thompson (1992). 
	#P(g) is adapted directly from Levene (1949), since it is, I believe, incorrect in Guo and Thompson (1992).
	#This version agrees with results produced by Wigginton et al 2005 for biallelic systems whereas that in Guo produces probabilites above one.
	#function to calculate P(g)
	#takes $topfact, $twonfact, %geno from the global env, redone for whatever configuration is being tested.
	sub Pg{
		#get the bottom and right portions (same loop)
		$botfact = $twonfact;
		$right_exp = 0;
		foreach $genotype (sort keys %geno){
			$botfact *= fact($geno{$genotype});
			@a = split(/_/, $genotype);
			chomp(@a);
			if ($#a == 1){ #if a heterozygote
				$right_exp += $geno{$genotype};
			}
		}
		
		#put it together
		$right = 2**$right_exp;
		$prob = ($topfact / $botfact)*$right;
		return($prob);
	}

	#get P(f)
	$Pf = Pg();
	#print("\tBottom Factorial: $botfact\n\t Right Multiplier: $right\n\tP(f):$Pf\n");	

	#perform permutations, get p value
	$perm = 1;
	$K = 0; #counter variable
	while($perm <= $Npermutations){
		#print("\tPermutation: $perm\n");
		undef %geno;
		@deck = shuffle(0..$#Ai); #get an array of random indices
		@Ai = @Ai[@deck]; #re-order the Ai array by the new order
		
		#print("Ai: \t");
        	#print Dumper \@Ai;
        	#print("\n");
		#print("This Ai:\t");
		#foreach $element (@Ai){print("$element\t");}
		#print("\n");
		#chop up @Ai to get new genotypes.
		$j = 0;
		while ($j <= ($#Ai - 1)){
			$k = $j + 1;
			undef @ta;
			push(@ta, $Ai[$j]); #get the first allele
			push(@ta, $Ai[$k]); #get the second allele
			#print("Genotype from $j and $j+1:\n");
			#print("\t\t$ta[0]\n\t\t$ta[1]\n");
			if ($ta[0] eq $ta[1]){ #if a homozygote
				#print("Homozygote.\n");
				$genotype = $ta[0];
				$geno{$genotype}++;
			}
			else{ #if a heterozygote
				@ta = sort(@ta); #sort them
				$genotype = $ta[0] . "_" . $ta[1]; #get the genotype
				$geno{$genotype}++; #increase the genotype count.
			}
			$j = $j + 2;
		}
		#foreach $element (keys %geno){print("\t\tGenotype:$element\t$geno{$element}\n");}
		$tPg = Pg();
		#print("\t\tP(g) = $tPg\n");
		if ($tPg <= $Pf){ #if this is less probable, add one to $K.
			$K++;
		}
		$perm++;
	}
	$pval = $K / $Npermutations; #get the P-value for this locus
	#print("\tP-value: $pval\n");

	#add to violating loci if lower that 0.05.
	if($pval <= 0.05){
		@l = split(/__/, $locus);
                print OUTFILE ("$l[0]\t$l[1]\t$pval\n"); #print the locus chr and position to the output
		$vio_count++;
        }
}
print("Completed. $vio_count out of $n_tot_loci loci found in violation of HWE or with more than two alleles in one sample.\n");




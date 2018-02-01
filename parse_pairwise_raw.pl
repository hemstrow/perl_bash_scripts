#!/usr/bin/perl
$infile = $ARGV[0];
$outfile = $ARGV[1];


%h = ();
open(IN, $infile);
foreach $line (<IN>)
	{
	chomp($line);
	@w_array = split(' ', $line);
	$i=0;
	#foreach $e (@w_array){print("$e is element $i\t"); $i++;}
	#print("\n");
	if ($w_array [0] =~ /Locus/)
		{
		#print("$w_array[0] is a match for Locus.\n");
		($snp) = $w_array[1] =~ /SNP_(\w+)/;
		#print("$snp\n");
		}
	elsif ($w_array[0] =~ /2/)
		{
		$fst = $w_array[1];
		$h{$snp} = $fst;
		}
	else
		{
		next;
		}
	}
close (IN);

open (OUT, '>', $outfile);
$snp = ();
$fst = ();
foreach $snp (sort {$a <=> $b} keys %h)
	{
	print OUT ("$snp\t", "$h{$snp}", "\n");
	}

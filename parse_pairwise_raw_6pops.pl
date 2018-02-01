#!/usr/bin/perl
$infile = $ARGV[0];
$outfile = $ARGV[1];


%hoh = ();
open(IN, $infile);
foreach $line (<IN>)
	{
	chomp($line);
	@w_array = split(/\s+/, $line);
	#$i=0;
	#foreach $e (@w_array){print("$e is element $i\t"); $i++;}
	#print("\n");
	if ($w_array [1] =~ /Locus/)
		{
		#print("$w_array[1] is a match for Locus.\n");
		($snp) = $w_array[2] =~ /SNP_(\w+)/;
		#print("$snp\n");
		}
	elsif ($w_array[0] =~ /2/)
		{
		$fst1_2 = $w_array[1];
		$hoh{$snp}{"1_2"} = $fst1_2;
		}
	elsif ($w_array[0] =~ /3/)
		{
		$fst1_3 = $w_array[1];
		$fst2_3 = $w_array[2];
		$hoh{$snp}{"1_3"} = $fst1_3;
		$hoh{$snp}{"2_3"} = $fst2_3;
		}
	elsif ($w_array[0] =~ /4/)
		{
		$fst1_4 = $w_array[1];
                $fst2_4 = $w_array[2];
                $fst3_4 = $w_array[3];
                $hoh{$snp}{"1_4"} = $fst1_4;
                $hoh{$snp}{"2_4"} = $fst2_4;
		$hoh{$snp}{"3_4"} = $fst3_4;
		}
	elsif ($w_array[0] =~ /5/)
		{
		$fst1_5 = $w_array[1];
                $fst2_5 = $w_array[2];
                $fst3_5 = $w_array[3];
		$fst4_5 = $w_array[4];
                $hoh{$snp}{"1_5"} = $fst1_5;
                $hoh{$snp}{"2_5"} = $fst2_5;
                $hoh{$snp}{"3_5"} = $fst3_5;
		$hoh{$snp}{"4_5"} = $fst4_5;
		}
	elsif($w_array[0] =~ /6/)
		{
		$fst1_6 = $w_array[1];
                $fst2_6 = $w_array[2];
                $fst3_6 = $w_array[3];
                $fst4_6 = $w_array[4];
                $fst5_6 = $w_array[5];
                $hoh{$snp}{"1_6"} = $fst1_6;
                $hoh{$snp}{"2_6"} = $fst2_6;
                $hoh{$snp}{"3_6"} = $fst3_6;
                $hoh{$snp}{"4_6"} = $fst4_6;
		$hoh{$snp}{"5_6"} = $fst5_6;
		}
	else
		{
		next;
		}
	}
close (IN);

open (OUT, '>', $outfile);
print OUT ("snp\t1_2\t1_3\t1_4\t1_5\t1_6\t2_3\t2_4\t2_5\t2_6\t3_4\t3_5\t3_6\t4_5\t4_6\t5_6\n");
$snp = ();
foreach $snp (sort {$a <=> $b} keys %hoh)
	{
	print OUT ("$snp\t");
	foreach $comp (sort keys %{$hoh{$snp}})
		{
		print OUT ("$hoh{$snp}{$comp}", "\t");
		}
	print OUT ("\n");
	}

#!/usr/bin/perl

#Same deal as Hash Seqs, but takes a list of fastqs and hashes all of them together. For doing HWE on reads and such. Does not split by pop.
$file = $ARGV[0]; #input list of fastqs
$min = $ARGV[1]; #min number of reads to keep a seq
$max = $ARGV[2]; #max number of reads to keep a seq
$lib = $ARGV[3]; #identifier for output


open(FILE, "<", $file) 
	or die("Can't find list of fastq files\n");

while(<FILE>){
	$tfile = $_;
	chomp($tfile);
	#print("$tfile\n");
	open(TFILE, "<", $tfile);
	while (<TFILE>) {
		$seq_line = <TFILE>;
		chomp($seq_line);
		$tags{$seq_line}++;

		$_ = <TFILE>;
		$_ = <TFILE>;
	}
	close TFILE;
}

$z = 1;
foreach $key (sort { $tags{$b} <=> $tags{$a} } keys %tags) {

	if ($tags{$key} >= $min && $tags{$key} <= $max) {	
		print ">$lib;$z;$tags{$key}\n";
		print "$key\n";
	}

	$z++;
}



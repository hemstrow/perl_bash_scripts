#!/usr/bin/perl

$file = $ARGV[0];
$lib = $ARGV[1];
$min = 10;
$max = 100;

open(FILE, "<$file")
	or die;

while (<FILE>) {
		
	$seq_line = <FILE>;
	chomp($seq_line);
	$tags{$seq_line}++;

	$_ = <FILE>;
	$_ = <FILE>;
}
close FILE;


$z = 1;
foreach $key (sort { $tags{$b} <=> $tags{$a} } keys %tags) {

	if ($tags{$key} >= $min && $tags{$key} <= $max) {	
		print ">$lib;$z;$tags{$key}\n";
		print "$key\n";
	}

	$z++;
}



#!/usr/bin/perl

############################################################################################
$max_count = 200;
############################################################################################

$file = $ARGV[0];

open(FILE, "<$file")
	or die;

while (<FILE>) {
	$line = $_; chomp($line);
	$_ = <FILE>;

	($lib_ID, $seq_ID, $seq_count) = split(/\;/,$line);

	if ($seq_count > $max_count) { $seq_count = $max_count + 1; }
	
	$array[$seq_count]++;
}
close FILE;

$x = 1;
while ($x <= $max_count + 1) {
	print "$x\t$array[$x]\n";
	$x++;
}




#!/usr/bin/perl

my $info = $ARGV[0];

if($#ARGV == 1){
	$reps = $ARGV[1];
}
else{
	$reps = `wc -l $info`;
}

my @params = split(/_/, $info);
my $gcs = $params[0];
my $N = $params[1];

if ($#params > 2){
	$mut = $params[2];
	@ploida = split(/\//, $params[3]);
	$ploid = $ploida[0];
}
else{
	my @muta = split(/\//, $params[2]);
	$mut = $muta[0];
}

open(INFOFILE, ">", "infofile.txt");

$i = 1;
while($i <= $reps){
	if ($ploid){
		 print INFOFILE ("$gcs\t$N\t$mut\t$ploid\n");
	}
	else{
		print INFOFILE ("$gcs\t$N\t$mut\n");
	}
	$i++;
}	

close(INFOFILE);	

#!/usr/bin/perl
open(FILE, "$ARGV[0]") or die ("Can't open $ARGV[0]!\n"); 
open(FILE2, "$ARGV[1]") or die ("Can't open $ARGV[1]!\n");
while(<FILE2>) {
	chomp($_); 
	$_ =~ s/ /_/g;
       	@array2 = split(/\t/, $_);
	#$n = $array2[0];
	#$n =~ s/\w\d+$//g;
	#$hash{$array2[0]} = "$n"."_$array2[3]";
       	$hash{$array2[0]} = $array2[3];
}
while(<FILE>) {
	chomp($_);
       	@array = split(/\t/, $_);
       	for ($i = 0; $i < @array; $i++) { 
		if ($hash{$array[$i]} ne "") {
			print ("$hash{$array[$i]}\t");
		} 
		else {
			print("$array[$i]\t");
		}
	}
}

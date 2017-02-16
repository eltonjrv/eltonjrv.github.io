#!/usr/bin/perl
#Programmer: Elton Vasconcelos (16/Feb/2017)
#Usage: perl get_num_reads_per_otu.pl [QIIME *_otus.txt file]

open(FILE, "$ARGV[0]") or die ("Can't open infile $ARGV[0]\n");
$line=<FILE>; 
chomp($line); 
while($line ne "") { 
	@array = split(/\t/, $line); 
	$num_reads = @array -1;
       	print("$num_reads\t$line\n");
       	$line=<FILE>;
       	chomp($line);
}

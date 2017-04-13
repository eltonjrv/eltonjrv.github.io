#!/usr/bin/perl
# Programmer: Elton Vasconcelos (11/Apr/2017)
# Usage: perl reads_perOTU_perSample.pl [seq_otus.txt2] [sample_IDs.txt] >outfile.tsv
#####################################################################################
# ATTENTION: The "_\d+" at the end of the reads' names must be removed on the seqs_otus.txt qiime output
# Do the following to remove:
# sed 's/_[0-9]*\t/\t/g' seqs_otus.txt >seqs_otus.txt2
# perl -pi -e 's/_\d+$//g' seqs_otus.txt2

my (@array, %hash, $line2, $cluster);

open(FILE,"$ARGV[0]") or die ("Can't open $ARGV[0]'\n");
my $line=<FILE>;
chomp($line);
open(FILE2,"$ARGV[1]") or die ("Can't open $ARGV[1]'\n");
my @array2=<FILE2>;
chomp(@array2);
my $header = join("\t", @array2);
print("\t$header\n");

while($line ne ""){
	for ($i = 0; $i < @array2; $i++) {
		$hash{$array2[$i]} = 0;
	} 
	@array = split(/\t/,$line);
	$cluster = shift(@array);
	for ($j = 0; $j < @array; $j++) {
		$hash{$array[$j]}++;
	}
	print("$cluster");
	for ($i = 0; $i < @array2; $i++) {
		print("\t$hash{$array2[$i]}");
	} 
	print("\n");
	$line=<FILE>;
	chomp($line);
}

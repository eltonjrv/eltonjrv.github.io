#!/usr/bin/perl -w
#Programmer: ELton Vasconcelos (27/Jan/2017)

foreach $file (<*p.fastq>) { 
	system("cat $file | paste - - - - >$file.tab"); 
	open(FILE, "$file.tab") or die ("Can't open $file.tab\n");
	$line = <FILE>;
	chomp($line);
	$outfile = $file;
	$outfile =~ s/\.fastq/\.lengths/g;
	open(OUT, ">$outfile"); 
	while($line ne "") {
		@array = split(/\t/, $line);
		$len = length($array[1]);
		print OUT ("$len\n");
		$line = <FILE>;
		chomp($line);
	}
	`rm $file.tab`; 
}

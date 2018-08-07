#!/usr/bin/perl
# Programmer: Elton Vasconcelos (12/Jul/2017)
# Usage: perl prepBCmapFiles4Qiime2.pl [iNext-barcodes.tab] [samples-map.tab]
# Script first run within the following dir:

my(@array, @array2, %hash, $samp, @iF, @iR, @samples, $iFid, $iFbc, $iRid, $iRbc, $outfile);
######## Working on samples-map.tab ###################
#######################################################
open (FILE2, "$ARGV[1]") or die ("Can't open file $ARGV[1]!\n");
my $line2 = <FILE2>;
chomp($line2);
while ($line2 ne "") {
   @array2 = split(/\t/, $line2);
   $array2[1] =~ s/_/\-/g;
   $samp = "$array2[1]"."-"."$array2[0]";
   push(@samples, $samp);
   $hash{$samp} = "$array2[2]"."$array2[3]";
   $line2 = <FILE2>;
   chomp($line2);
}

######## Working on iNext-barcodes.tab ################
########################################################
open (FILE, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
@array = <FILE>;
chomp(@array);
for ($i = 0; $i < @array; $i++) {
	if ($array[$i] =~ m/^519/) {
		push(@iF, $array[$i]);	
	}
	elsif ($array[$i] =~ m/^926/) {
		push(@iR, $array[$i]);	
	}
}

for ($a = 0; $a < @samples; $a++) {
	$outfile = "barcodes_"."$samples[$a]".".tab";
	open(OUT, ">$outfile") or die ("Can't open $outfile\n");
	print OUT ("#SampleID\tBarcodeSequence\tLinkerPrimerSequence\tDescription\n");        #CAGCMGCCGCGGTAATWC  CCGTCAATTCCTTTRAGGTT 
	for ($i = 0; $i < @iF; $i++) {
		$iFid = $iF[$i];
		$iFid =~ s/\t.*//g;
		$iFbc = uc($iF[$i]);
		$iFbc =~ s/.*\t//g;
		for ($j = 0; $j < @iR; $j++) {
			$iRid = $iR[$j];
			$iRid =~ s/\t.*//g;
			$iRbc = uc($iR[$j]);
			$iRbc =~ s/.*\t//g;
			if ($hash{$samples[$a]} eq "$iFid$iRid") {
				print OUT ("$samples[$a]\t$iFbc$iRbc\t\t\n");
			}
			#else {
			#	print OUT ("BARCODE\t$iFbc\t$iRbc\tignore\n");
			#}
		}
	}
}


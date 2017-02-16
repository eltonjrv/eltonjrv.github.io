#!/usr/bin/perl
## Programmer: Elton Vasconcelos (08/Feb/2017)
## Script that takes the demultiplex_fasta.py qiime output as input (file.fasta) and catches the sequences from samples informed on the metadata.tab
## Usage: perl catching_samples_after_demultiplexing-QIIMEpipe.pl [file.fasta]] [metadata.tab] >outfile.fasta
################################################################################################################
## ATTENTION:
## The fasta file name must start with its respective metaSampleID followed by a dash or underscore (e.g.: M10-demultiplexed.fna)
## The metadata.tab must contain four colunms like the following:
#metaSample	study_label	iNext519F	iNext926R
#M1	B_flea	519_A	926_1
#M10	NC-014	519_A	926_10
#M11	FOF-012	519_A	926_11
#M12	DS-006	519_A	926_12
#M2	X_flea	519_A	926_2
#M3	F_flea	519_A	926_3
#################################################################
# The headers from the fasta file must have the following structure
# >iNextR-iNextF-metaSampleID 
# See the example below:
#>926_3_1 Reversed: 519_A_1 orig_bc=CCACGTCC new_bc=CCACGTCC bc_diffs=0
#>926_3_2 Reversed: 519_B_2 orig_bc=CCACGTCC new_bc=CCACGTCC bc_diffs=0
#>926_3_3 Reversed: 519_B_3 orig_bc=CCACGTCC new_bc=CCACGTCC bc_diffs=0
################################################################

my($line, $line2, @array, @array2, %hash, $metaSample, $metaID, $iNextF, $iNextR, $header);

########## Working on the metadata file ############
####################################################
open(FILE2,"$ARGV[1]") or die ("Can't open file $ARGV[1]!\n");
$line2 = <FILE2>;
chomp($line2);

while ($line2 ne "") {
	@array2 = split(/\t/,$line2);
	$metaSample = $array2[3]."-".$array2[2]."-".$array2[0];	#Concatenating iNextR to iNextF to metaSampleID
	$hash{$metaSample} = $array2[1];	#Associating the iNext_pairs and metaSampleID to the study label
	$line2 = <FILE2>;
	chomp($line2);
}
##############################################################


######### Working on the demultiplexed fasta file ############
##############################################################
open (FILE, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
$line = <FILE>;
chomp($line);
$metaID = $ARGV[0];
$metaID =~ s/[_\-].*$//g;	#Grabbing the metasampleID from the beginning of the fasta file name
$metaSample = "";
my $c = 1;

while ($line ne "") {
	if ($line =~ m/^>/) {	
		$line =~ s/^>//g;
		@array = split(/\s/, $line);
		$iNextR = $array[0];
		$iNextR =~ s/_\d+$//g;
		$iNextF = $array[2];
		$iNextF =~ s/_\d+$//g;
		$metaSample = $iNextR."-".$iNextF."-".$metaID;
		if ($hash{$metaSample} ne "") {
			$header = ">$hash{$metaSample}\-$metaID";
			$header =~ s/_/\-/g; #ATTENTION: This replacement is very important, because Sample IDs must not contain "_" for the plotting charts qiime scripts. The "_" must only be placed at the end of the sequence ID, followed by any digit(s), like I am doing at the print function below ("$header"."_$c ...)
			#print(">$hash{$metaSample}\-$metaSample\_$c $line\n");
			print("$header"."_$c $line\n");
			$line = <FILE>;
			chomp($line);
			$c++;
			until ($line =~ m/^>/ || $line eq "") {
				print("$line\n");
				$line = <FILE>;
				chomp($line);
			}
		}
		else {
			 $line=<FILE>;
			chomp($line);
		}
	}
	else {
		$line=<FILE>;
		chomp($line);	
	} 
}

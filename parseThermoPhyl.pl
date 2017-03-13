#!/usr/bin/perl
# Programmer: Elton Vasconcelos (10/Mar/2017)
# Usage: perl parseThermoPhyl.pl [target_species_names.txt] [raw_search_results.txt] [sorted_search_results.txt] [nonTargets_cutoff]
################################################################################################################
# The target_species_names.txt file must be a list of species names (Genus species) one per line.
# The other two inputs are the default outputs from ThermoPhyl.
# The last argument is an integer corresponding to the maximum number of non-target sequences the user wants to see.
# The output will be named best_candidates.tsv, which will contain the primer_pairs that amplify all the target species and less than [nonTargets_cutoff] non-target sequences.

use List::MoreUtils qw(uniq);

if (@ARGV != 4) {
	die ("*** You must specify 4 arguments ***\nUsage: perl parseThermoPhyl.pl [target_species_names.txt] [raw_search_results.txt] [sorted_search_results.txt] [nonTargets_cutoff]\nPlease, read script's header for more instructions.\n");
}

my $old_file = "best_candidates.tsv";
if (-f $old_file) {	#If there is already a best_candidates.tsv file within the current directory, it will be removed!
	`rm best_candidates.tsv`;
}

open (FILE1, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
my @species = <FILE1>;
chomp(@species);

open (FILE2, "$ARGV[1]") or die ("Can't open file $ARGV[1]!\n");
my $line = <FILE2>;
chomp($line);

my(@array, $c, $primerID, @caught_species, @uniq_caught_species);

`echo -e "Primer_pair\tTarget_matches\tNon_target_matches\tF\tF_pos\tR_motif\tR_pos\tamp_length" >best_candidates.tsv`;

while ($line ne "") {
	if ($line =~ m/^[YN]\t/) {
		$c = 0;
		@caught_species = ();
		@array = split(/\t/, $line);
		$primerID = $array[1];
		until ($array[1] != $primerID || $line eq "") {	# Walking on a primerID at a time
			if ($array[0] eq "N") {		# Counting the number of non-targets
				$c++;
			}
			elsif ($array[0] eq "Y") {					
	   			for ($i = 0; $i < @species; $i++) {	
					if ($line =~ m/$species[$i]/) {	# Finding a match for the target species' name within the raw_search_results.txt
						push(@caught_species, $species[$i]);
					}
				}
			}
			$line = <FILE2>;
			chomp($line);
			@array = split(/\t/, $line);
		}
		@uniq_caught_species = uniq(@caught_species);
		if (@uniq_caught_species == @species && $c <= $ARGV[3]) {   # All target species were caught together with less than 1000 non-target species
			print ("Primer_pair $primerID caught all target species and less than $ARGV[3] non-target sequences\n");
			`grep -P \'^$primerID\t\' $ARGV[2] >>best_candidates.tsv`;
		}

	}	
	else {
		$line = <FILE2>;
		chomp($line);
	}
}

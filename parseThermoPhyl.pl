#!/usr/bin/perl
# Programmer: Elton Vasconcelos (10/Mar/2017)
# Usage: perl parseThermoPhyl.pl [target_species_names.txt] [raw_search_results.txt] [sorted_search_results.txt]
################################################################################################################
# The target_species_names.txt file must be a list of species names (Genus species) one per line.
# The other two inputs are the default outputs from ThermoPhyl.
# The output will be named best_candidates.tab, which will contain the primer_pairs that amplify all the target species and less than 1000 non-target sequences.

use List::MoreUtils qw(uniq);

open (FILE1, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
my @species = <FILE1>;
chomp(@species);

open (FILE2, "$ARGV[1]") or die ("Can't open file $ARGV[1]!\n");
my $line = <FILE2>;
chomp($line);

my(@array, $c, $primerID, @caught_species, @uniq_caught_species);

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
   			for ($i = 0; $i < @species; $i++) {	
				if ($line =~ m/$species[$i]/) {	# Finding a match for the target species' name within the raw_search_results.txt
					push(@caught_species, $species[$i]);
				}
			}
			$line = <FILE2>;
			chomp($line);
			@array = split(/\t/, $line);
		}
		@uniq_caught_species = uniq(@caught_species);
		if (@uniq_caught_species == @species && $c <= 1000) {	# All target species were caught together with less than 1000 non-target species
			print ("Primer_pair $primerID caught all target species and less than 1000 non-target ones\n");
			`grep -P \'^$primerID\t\' $ARGV[2] >>best_candidates.tab`;
			print ("--\n");
		}

	}	
	else {
		$line = <FILE2>;
		chomp($line);
	}
}

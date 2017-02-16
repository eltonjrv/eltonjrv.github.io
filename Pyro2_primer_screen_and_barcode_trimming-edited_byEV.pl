#!/usr/bin/perl
# Usage: perl script.pl [input.fasta] [hashtable_file-barcodes.tab]	# EV
use Bio::SeqIO;
use Bio::Tools::IUPAC;

#############################################################################################
### Working first on the hash table to associate each file with their respective barcodes ###
open(FILE2,"$ARGV[1]") or die ("Can't open $ARGV[1]!\n");
$line2=<FILE2>;
chomp($line2);
while($line2 ne ""){
	@array2 = split(/\t/,$line2);
	$hash_i5{$array2[0]} = $array2[2];
	$hash_i7{$array2[0]} = $array2[3];
	$line2=<FILE2>;
	chomp($line2);
}
############################################################################################# EV

#=================================
#Script to screen for primers
#Will trim off sequences external to primer(s)
#briandotoakleyatarsdotusdadotgov
#=================================

#---define input file----
my $usage = "\nUsage:  $0 filename\n\n";
$input_file = $ARGV[0] or die $usage;
chomp $input_file;
my $prefix = $input_file;	# EV
$prefix =~ s/\-*n*n*\.[fastanq]+$//g; # EV

#----open output file----
#Define output file for forward primer found
$output_file1 = $prefix.'-forward_primer_found.fas';	# EV
$out_f = Bio::SeqIO->new(-file => ">$output_file1",-format => 'fasta');
#Define output file for reverse primer found
$output_file2 = $prefix.'-reverse_primer_found.fas';	# EV
$out_r = Bio::SeqIO->new(-file => ">$output_file2",-format => 'fasta');
#Define output file for both primers found
$output_file3 = $prefix.'-both_primers_found.fas';	# EV
$out_both = Bio::SeqIO->new(-file => ">$output_file3",-format => 'fasta');

print "\nScript is running...\n\n";
	
#---Primer sequence definitions---  MODIFY TO SUIT	
#tag=AACCCACTTT
#$forward_motif=GTGTAGCAGTTACCGCA;#  - dsr forward primer
#$reverse_motif=CCCTGKGTATGRACGATGTTG;

#my 16S primers
#$forward_motif=GAGTTTGATCNTGGCTCAG; #27F primer
#$forward_motif = GRGTTTGATCNTGGCTCAG; # 27F primer more degenerate
#$forward_motif=ATMTTTAYGGTGSTGCWGCT; #porA_103F
#$forward_motif=CTGARGTAAAGCGTTTGCG;#pgm
#$forward_motif=TGYACACACCGCCCGT; #1406F
#$reverse_motif=TACCTTGTTACGACTT; #1492R
#$reverse_motif=GTTACGACTT; #1492R truncated
#$reverse_motif=AAGGAGGTGWTCCARCC;#1541R
#$forward_motif=CGAAATTCCTTGTCRGKTAA; #1804F(23S)
#$forward_motif=CAGCMGCCGCGGTAATWC; #519F
#$forward_motif=CNGCNGCC; #519F - shortened
$forward_motif = $hash_i5{$prefix};	# EV

#$forward_motif=GGCGVACGGGTGAGTAA; #104F primer
#$forward_motif=AGCGCACTGCTCAGTAA; #104F primer Archaeaa-like
#$reverse_motif=CCGCNGCNGCTGGCAC; # 530R primer
#$reverse_motif=GTGCCAGCNGCCGC; # 515F primer shortened Turner et al. 1999
#$reverse_motif=ACCGCCCCAGTCAAACTRC; #2109R(23S)
#$reverse_motif=GGGTTBCCCCATTCRG; #23Sr
#$reverse_motif=GAARCAGTTCGYkTWGGTG; #CMP1178
#$forward_motif=GGATGACACTTTTCGGAG; #410F
#$reverse_motif=AATCCATCTGCCTCTCC; #661R
#$reverse_motif=CCYTATCCWCARCTTTTAAYCAA; #pgm
#$reverse_motif=CACCWAMRCGAACTGYTT; #porA_549R
#$reverse_motif=CCGTCAATTCCTTTRAGGTT; #926R
#$reverse_motif=CCGTCAATT; #926R - shortened
$reverse_motif = $hash_i7{$prefix};	# EV

$F_offset=1;   #Will trim anything before primer (e.g. pyrosequencing tag), leaves primer sequence
$R_offset=length($reverse_motif); #Will trim anything after primer, leaves primer sequence

#---Allow for IUPAC ambiguity codes---	
 my $F_ambiseq = new Bio::Seq (-seq => $forward_motif, -alphabet => 'dna');
 my $F_stream  = new Bio::Tools::IUPAC(-seq => $F_ambiseq);
 print "\nSearching for forward primers:\n\n";
	while ($unique_seq = $F_stream->next_seq() ) {
	$forward_motif=$unique_seq->seq();
	push (@forward_motif, $forward_motif);
	print "$forward_motif\n";
	}
	
 my $R_ambiseq = new Bio::Seq (-seq => $reverse_motif, -alphabet => 'dna');
 my $R_stream  = new Bio::Tools::IUPAC(-seq => $R_ambiseq);
 print "\nSearching for reverse motif:\n\n";
	while ($unique_seq = $R_stream->next_seq() ) {
	$unique_seq=$unique_seq->revcom();  #reverse complement primers if desired
	$reverse_motif=$unique_seq->seq();
	push (@reverse_motif, $reverse_motif);
	print "$reverse_motif\n";
	}

print "\nSearching for ",scalar(@forward_motif)," forward primers, and ",scalar(@reverse_motif)," reverse primers\n\n";

#---heavy lifting--------
$forward_counter=0;
$reverse_counter=0;	
my $seq_in  = Bio::SeqIO->new( -format => 'fasta',-file => $input_file);
while ($seq_obj = $seq_in->next_seq()) {
 	$seq_string = $seq_obj->seq();
	$length=$seq_obj->length();
	$F_pos=-1;
 	$R_pos=-1;
	
	foreach $forward_motif(@forward_motif) {
	if ($seq_string =~ /$forward_motif/) {
	$forward_counter+=1;
	$F_pos = index($seq_string, $forward_motif);
	}
	}
	
	foreach $reverse_motif(@reverse_motif) {
	if ($seq_string =~ /$reverse_motif/) {
	$reverse_counter+=1;
	$R_pos = index($seq_string, $reverse_motif);
	}
	}
	
	if (($F_pos>=0)){
	$out_f->write_seq($seq_obj->trunc($F_pos+$F_offset,$length));
	}

	if (($R_pos>=0)){
	$out_r->write_seq($seq_obj->trunc(1,$R_pos+$R_offset));
	}
	
	if (($F_pos>=0)&&($R_pos>=0)){
#	$out_both->write_seq($seq_obj);
	$out_both->write_seq($seq_obj->trunc($F_pos+$F_offset,$R_pos+$R_offset));
	$both_primers_found_counter+=1;
	}
}

#---Print results--------

print "\n", $forward_counter," forward primers found\n";
print $reverse_counter," reverse primers found\n";
print "\n",$both_primers_found_counter," sequences had both primers present\n";


print "\nResults written to:\n\n$output_file1\n$output_file2\n$output_file3\n\n";



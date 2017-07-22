#!/usr/bin/perl
# Programmer: Elton Vasconcelos (20/Jul/2017)
# Usage: perl autoMergeSamples_qiime2.pl [qzaTables-ls.txt] [qzaRepSeqs-ls.txt]

open (FILE1, "$ARGV[0]") or die ("Can't open file $ARGV[0]!\n");
my @tables = <FILE1>;
chomp(@tables);
open (FILE2, "$ARGV[1]") or die ("Can't open file $ARGV[1]!\n");
my @repseqs = <FILE2>;
chomp(@repseqs);

`source activate qiime2-2017.6`;

###########################################################################
#### Working on merging feature-table.biom from each individual sample ####
###########################################################################
`qiime feature-table merge --i-table1 $tables[0] --i-table2 $tables[1] --o-merged-table mergedTable-1.qza`;

for ($i = 2; $i < @tables; $i++) {
	$a = $i - 1;
	`qiime feature-table merge --i-table1 mergedTable-$a.qza --i-table2 $tables[$i] --o-merged-table mergedTable-$i.qza`;
}
$a = $i - 1;
`mv mergedTable-$a.qza allSamples-mergedTable.qza`;
`rm mergedTable*qza`;

###########################################################################
### Working on merging representative seqs from each individual sample ####
###########################################################################
`qiime feature-table merge-seq-data --i-data1 $repseqs[0] --i-data2 $repseqs[1] --o-merged-data mergedRepSeqs-1.qza`;

for ($i = 2; $i < @repseqs; $i++) {
	$a = $i - 1;
	`qiime feature-table merge-seq-data --i-data1 mergedRepSeqs-$a.qza --i-data2 $repseqs[$i] --o-merged-data mergedRepSeqs-$i.qza`;
}
$a = $i - 1;
`mv mergedRepSeqs-$a.qza allSamples-mergedRepSeqs.qza`;
`rm mergedRepSeqs*qza`;

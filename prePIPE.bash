#!/usr/bin/bash
### This is a pilot pipeline using both ad-hoc PERL scripts (placed in the bin branch) and QIIME tools (http://qiime.org/scripts/index.html)

### One must have QIIME installed and run the following command before starting the pipeline
source activate qiime1

### Quality control (trimming) and demultiplexing
bash qiimePipe-demultiplexing.bsh

### Catching our samples and relabelling each sequence (adding Diniz's study labels)
for i in `ls *fna`; do perl catching_samples_after_demultiplexing-QIIMEpipe.pl $i metaSampleID-studyLabel-barCodes-prep.tab >`echo $i | sed 's/\.fna/\-DinizLabels.fasta/g'`; done

### Preparing a summary of the sequences content in each sample
wc -l *joined/fastqjoin.join.fastq | sed -r 's/^ +//g' | sed 's/ /\t/g' | awk '{ print $2 "\t" $1 / 4 }' >numSeqs-joined.txt
grep -c '>' *fastaqual/fastqjoin.join.fna | sed 's/\:/\t/g' >numSeqs-fastaqual.txt 
grep -c '>' *F_splitOUT/seqs.fna | sed 's/\:/\t/g' >numSeqs-iNextF-demultiplexed.txt 
grep -c '>' *R_splitOUT/seqs.fna | sed 's/\:/\t/g' >numSeqs-iNextR-demultiplexed.txt 
grep -c '>' *Labels.fasta | sed 's/\:/\t/g' >numSeqs-DinizLabels.txt
paste numSeqs-joined.txt numSeqs-fastaqual.txt numSeqs-iNextF-demultiplexed.txt numSeqs-iNextR-demultiplexed.txt numSeqs-DinizLabels.txt >qiimePipe-Summary.tsv

### Screening for chimeras using usearch aligner and PATRIC 16S db as reference
#$ for i in `ls *Labels.fasta`; do identify_chimeric_seqs.py -i $i -r ../../16S-DBs/PATRIC/16s_RNA_PATRIC.frn -m usearch61 -o `echo $i | sed 's/\-.*//g'`-chimScreenOUT; done
for i in `ls *Labels.fasta`; do identify_chimeric_seqs.py -i $i -r ../../16S-DBs/ARB-SILVA/SSURef_NR99_128_SILVA_07_09_16_opt-typeStrains-UbyT.fasta -m usearch61 --non_chimeras_retention intersection -o `echo $i | sed 's/\-.*//g'`-chimScreenOUT; done
for i in `ls -d *chimScreenOUT`; do  perl /home/elton/bioinformatics-tools/perl-scripts/seqs1.pl -outfmt fasta -incl $i/non_chimeras.txt -seq `echo $i | sed 's/\-chimScreenOUT//g'`-demultiplexed-DinizLabels.fasta  >`echo $i | sed 's/\-chimScreenOUT//g'`-demultiplexed-DinizLabels-nonChim.fasta; done
## Adding the number of non-chimera seqs in each sample to the Summary table created on line 18
wc -l *chimScreenOUT/non_chimeras.txt | head -26 | sed -r 's/^ +//g' | sed 's/ /\t/g' | awk '{ print $2 "\t" $1 }' >numSeqs-non_chimeras.txt
paste numSeqs-joined.txt numSeqs-fastaqual.txt numSeqs-iNextF-demultiplexed.txt numSeqs-iNextR-demultiplexed.txt numSeqs-DinizLabels.txt numSeqs-non_chimeras.txt >qiimePipe-Summary.tsv 

### Picking OTUs
#$ for i in `ls *nonChim.fasta`; do  parallel_pick_otus_usearch61_ref.py -i $i -r ../../16S-DBs/PATRIC/16s_RNA_PATRIC.frn -o `echo $i | sed 's/\-.*$//g'`-OTUs; done
for i in `ls *nonChim.fasta`; do  parallel_pick_otus_usearch61_ref.py -i $i -r ../../16S-DBs/ARB-SILVA/SSURef_NR99_128_SILVA_07_09_16_opt-typeStrains-UbyT.fasta -o `echo $i | sed 's/\-.*$//g'`-OTUs; done
wc -l *OTUs/*im_otus.txt | head -26 | sed -r 's/^ +//g' | sed 's/ /\t/g' | awk '{ print $2 "\t" $1 }' >numOTUs-DinizLabels.txt
# 1160  for i in `ls *OTUs/*nonChim_otus.txt`; do cut -f 1 $i | xargs -i grep '{} ' ../../16S-DBs/PATRIC/16s_RNA_PATRIC.frn | sort >`echo $i | sed 's/\/.*/\-DinizLabels\.txt/g'`; done
for i in `ls *OTUs/*nonChim_otus.txt`; do cut -f 1 $i | xargs -i grep '{} ' ../../16S-DBs/ARB-SILVA/SSURef_NR99_128_SILVA_07_09_16_opt-typeStrains-UbyT.fasta | sort >`echo $i | sed 's/\/.*/\-DinizLabels\.txt/g'`; done

### Generating OTU tables
for i in `ls *nonChim.fasta`; do pick_closed_reference_otus.py -i $i -r ../../16S-DBs/ARB-SILVA/SSURef_NR99_128_SILVA_07_09_16_opt-typeStrains-UbyT.fasta -s -p parameter.txt -f -o `echo $i | sed 's/\-.*$//g'`-closed_ref_OTUs; done 
## Adding the number of reads per OTU as the first column on the the *_otus.txt output from pick_closed_reference_otus.py script
for i in `ls *ref_OTUs/usearch61_ref_picked_otus/*_otus.txt`; do perl get_num_reads_per_otu.pl $i >`echo $i | sed 's/\.txt$/\.numberReads\.txt/g'`; done;

## Plotting summary charts for all the samples
merge_otu_tables.py -i `ls *ref_OTUs/otu_table.biom | xargs | sed 's/ /,/g'` -o otu_table-All.biom 
summarize_taxa_through_plots.py -i otu_table-All.biom -m newSampleIDs.tab -p parameter.txt -f -o OTUsPlots-silva

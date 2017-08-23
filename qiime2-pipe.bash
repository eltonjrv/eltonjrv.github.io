#!/usr/bin/bash
# Programmer: Elton Vasconcelos, DVM, PhD
# July, 2017
### This is a pilot pipeline using both ad-hoc PERL scripts (placed in the "bin" branch) and QIIME2 tools (https://qiime2.org).
### Accesories files used in some commands are placed within "accFiles" branch
### If you use this whole tool or part of it, please cite this github page acknowledging the author (Vasconcelos, EJ) as well as QIIME2.
################################################################################################################
#### Extracting barcodes for metasamples coming from a Nextera 2-steps library prep (paired-end sequencing) ####
## This initial step is performed with QIIME1 extract_barcodes.py script
source activate qiime1
mkdir 00-bc-extraction
cd 00-bc-extraction/
## Put all your fastq files into this current directory
for i in `ls | sed 's/R[12].*/R/g' | uniq`; do extract_barcodes.py -c barcode_paired_end -f `echo "$i""1_001.fastq.gz"` -r `echo "$i""2_001.fastq.gz"` --bc1_len 8 --bc2_len 8 -o `echo $i | sed 's/_.*/\-barcodes/g'`; done
for i in `ls -d *barcodes`; do gzip $i/*fastq; mv $i/reads1.fastq.gz $i/forward.fastq.gz; mv $i/reads2.fastq.gz $i/reverse.fastq.gz; done
source deactivate
cd ../

######################################
#### Demultiplexing and denoising ####
######################################
mkdir 01-demux
cd 01-demux/
## Prepapring the barcodes mapFile for each sample with an ad-hoc PERL script (prepBCmapFiles4Qiime2.pl)
perl prepBCmapFiles4Qiime2.pl iNext-barcodes.tab samples-map.tab 
## Activating QIIME2
source activate qiime2-2017.6
source tab-qiime
## Importing sequencing data (from 00-bc-extraction dir) to QIIME2
for i in `ls -d ../00-bc-extraction/M*barcodes`; do qiime tools import --type EMPPairedEndSequences --input-path $i --output-path `echo $i | sed 's/\.\..*M/M/g' | sed 's/\-barcodes/\-input4demux/g'`; done
## Demultiplexing 
for i in `ls barcodes_*tab`; do Mbase=`echo $i | sed s'/barcodes_.*M/M/g' | sed 's/\.tab//g'`; qiime demux emp-paired --m-barcodes-file $i --m-barcodes-category BarcodeSequence --i-seqs `echo $Mbase`-input4demux.qza --o-per-sample-sequences `echo $Mbase`-demuxOUT; done
for i in `ls *OUT.qza`; do qiime demux summarize --i-data $i --o-visualization `echo $i | sed 's/qza$/qzv/g'`; done
## Denoising (truncation parameters were decided after visualizing the sequences quality boxplot chart, loading *demuxOUT.qzv files onto view.qiime2.org)
for i in `ls *demuxOUT.qza`; do qiime dada2 denoise-paired --i-demultiplexed-seqs $i --p-trunc-len-f 260 --p-trunc-len-r 220 --o-representative-sequences `echo $i | sed 's/\.qza/\-denoised_repseqs.qza/g'` --o-table `echo $i | sed 's/\.qza/\-denoised_table.qza/g'` --p-n-threads 6; done
cd ../

################################################
#### Merging samples in a single biom table ####
################################################
mkdir 02-merge
cd 02-merge/
ln -s ../01-demux/*repseqs.qza .
ln -s ../01-demux/*table.qza .
qiime feature-table merge --i-table1 M1-demuxOUT-denoised_table.qza --i-table2 M2-demuxOUT-denoised_table.qza --o-merged-table test-M1M2_mergedTable.qza
qiime feature-table merge-seq-data --i-data1 M1-demuxOUT-denoised_repseqs.qza --i-data2 M2-demuxOUT-denoised_repseqs.qza --o-merged-data test-M1M2_mergedRepSeqs.qza
 
## Automating the merging process with an ad-hoc PERL script (autoMergeSamples_qiime2.pl)
ls *table.qza >qzaTables-ls.txt
ls *repseqs.qza >qzaRepSeqs-ls.txt
perl autoMergeSamples_qiime2.pl qzaTables-ls.txt qzaRepSeqs-ls.txt 
qiime feature-table summarize --i-table allSamples-mergedTable.qza --o-visualization allSamples-mergedTable.qzv --m-sample-metadata-file ../sample-metadata.tsv 

### Philogenetic analysis for grouping repseqs and then performing diversity metrics
qiime alignment mafft --i-sequences allSamples-mergedRepSeqs.qza --o-alignment allSamples-alignedRepSeqs.qza
qiime alignment mask --i-alignment allSamples-alignedRepSeqs.qza --o-masked-alignment allSamples-maskedAlignedRepSeqs.qza
qiime phylogeny fasttree --i-alignment allSamples-maskedAlignedRepSeqs.qza --o-tree allSamples-unrooted-tree.qza
qiime phylogeny midpoint-root --i-tree allSamples-unrooted-tree.qza --o-rooted-tree allSamples-rooted-tree.qza
cd ../

############################
#### Diversity Analysis ####
############################
mkdir 03-diversity
cd 03-diversity/
ln -s ../02-merge/allSamples-mergedTable.qza 
ln -s ../02-merge/allSamples-rooted-tree.qza 
qiime diversity core-metrics --i-phylogeny allSamples-rooted-tree.qza --i-table allSamples-mergedTable.qza --p-sampling-depth 47 --output-dir core-metrics-results
# -> 47 is the lowest number of sequences in a sample. So the diversity core metrics are going to use 47 seqs for all samples. Check that by visualizing allSamples-mergedTable.qzv on view.qiime2.org

## Alpha-diversity
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/faith_pd_vector.qza --m-metadata-file ../sample-metadata.tsv --o-visualization core-metrics-results/faith_pd-group-significance.qzv 
qiime diversity alpha-group-significance --i-alpha-diversity core-metrics-results/evenness_vector.qza --m-metadata-file ../sample-metadata.tsv --o-visualization core-metrics-results/evenness-group-significance.qzv 

## Beta-diversity
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza --m-metadata-file ../sample-metadata.tsv --m-metadata-category Description --o-visualization core-metrics-results/unweighted_unifrac-group-significance.qzv --p-pairwise 
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza --m-metadata-file ../sample-metadata.tsv --m-metadata-category Description --o-visualization core-metrics-results/weighted_unifrac-group-significance.qzv --p-pairwise 
qiime diversity beta-group-significance --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza --m-metadata-file ../sample-metadata.tsv --m-metadata-category Description --o-visualization core-metrics-results/bray_curtis-group-significance.qzv --p-pairwise 

## PCoA Emperor plots for the three most used bDiv metrics
qiime emperor plot --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza --m-metadata-file ../sample-metadata.tsv  --o-visualization core-metrics-results/unweighted_unifrac_pcoa_emperor.qzv 
qiime emperor plot --i-pcoa core-metrics-results/weighted_unifrac_pcoa_results.qza --m-metadata-file ../sample-metadata.tsv  --o-visualization core-metrics-results/weighted_unifrac_pcoa_emperor.qzv 
qiime emperor plot --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza --m-metadata-file ../sample-metadata.tsv  --o-visualization core-metrics-results/bray_curtis_pcoa_emperor.qzv
cd ../

###################################
#### Taxonomic classification #####
###################################
mkdir 04-taxClass
cd 04-taxClass/
ln -s ../02-merge/allSamples-mergedTable.qza
ln -s ../02-merge/allSamples-mergedRepSeqs.qza 
wget -O "gg-13-8-99-515-806-nb-classifier.qza" "https://data.qiime2.org/2017.6/common/gg-13-8-99-515-806-nb-classifier.qza"
qiime feature-classifier classify-sklearn --i-classifier gg-13-8-99-515-806-nb-classifier.qza --i-reads allSamples-mergedRepSeqs.qza --o-classification taxonomy.qza
qiime taxa tabulate --i-data taxonomy.qza --o-visualization taxonomy.qzv
qiime taxa barplot --i-table allSamples-mergedTable.qza --i-taxonomy taxonomy.qza --m-metadata-file ../sample-metadata.tsv --o-visualization taxa-bar-plots.qzv
cd ../


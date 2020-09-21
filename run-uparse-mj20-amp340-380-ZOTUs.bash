#!/bin/bash
# Robert Edgar pipeline presented at STAMPS 2017 course (Aug/2017)
# Adjusted by Elton Vasconcelos for own purposes (Jan/2018)
# Please refer to https://www.drive5.com/usearch/download.html in order to download USEARCH tools
# Script that takes paired-end Illumina fastq files as input and performs sequencing quality filtering, ZOTUs/ESVs assembly, taxonomic classification, and some diversity analyses.
# Usage: bash run-uparse-mj20-amp340-380-ZOTUs.bash [path-to-input-fastq-files]
###########################################################################################################################
# NOTE-1: One may tune parameters on each command below according to his/her needs.
# NOTE-2: On the first uparse command "usearch -fastq_mergepairs" (line 35), we have set a minimum of 20 bp for merging R1 and R2 mates (accepting a maximum difference of 5 bases within the overlapped region, as set by default), as well as a minimum and maximum merged sequence length of 340 and 380, respectively. This is because our V4-V5 target amplicon region is ~ 400 bp long and, after primers are trimmed, we get a ~ 360 bp-long full amplicon to be joined, allowing an arbitrary  +/- 20 bp range. 
# NOTE-3: This script uses a customized RDP refDB, which we have added Ehrlichia_canis, Ehrlichia_chafeensis, Anaplasma_platys, Anaplasma_phagocytophilum, Mycoplasma_haemocanis,Mycoplasma_haematoparvum 16S rRNA sequences.

if [ xusearch == x ] ; then
	echo Must set \usearch >> /dev/stderr
	exit 1
fi

version=`usearch -version | sed "-es/usearch //" | sed "-es/v11.*/v11/"`

if [ x$version != xv11 ] ; then
	echo "usearch version too old, need v11" >> /dev/stderr
	exit 1
fi

if [ ! -d $1 ] ; then
	echo "Directory $1 not found." >> /dev/stderr
	exit 1
fi

out=outputs

rm -rf $out
mkdir -p $out

# Assemble paired reads, put sample names into read labels
usearch -fastq_mergepairs $1/*_R1_*.fq -fastq_minmergelen 250 -fastq_maxmergelen 450 -relabel @ -fastqout $out/merged.fq

# Discard reads which probably have errors (quality filtering)
usearch -fastq_filter $out/merged.fq -fastq_maxee 1.0 -relabel Filt -fastaout $out/filtered.fa

# Find unique sequences and abundances (dereplication)
usearch -fastx_uniques $out/filtered.fa -sizeout -relabel Uniq -fastaout $out/uniques.fa

# Create ZOTUs by denoising (error-correction)
usearch -unoise3 $out/uniques.fa -zotus $out/zotus.fa

# Create a ZOTUs table (ZOTUs are 100% identical OTUs)
perl -pi -e 's/Zotu/Otu/g' $out/zotus.fa
usearch -otutab $out/merged.fq -zotus $out/zotus.fa -strand plus -otutabout $out/zotus_table_uparse.tsv -mapout $out/zotus_map.txt
perl -pi -e 's/Otu/Zotu/g' $out/zotu*

### Taxonomic classification ###
usearch -makeudb_sintax rdp_16s_extra_seqs.fa -output rdp_16s_extra_seqs.udb 	#formatting refDB
usearch -sintax $out/zotus.fa -db rdp_16s_extra_seqs.udb -tabbedout $out/zotus.sintax -strand both -sintax_cutoff 0.8
perl -pi -e 's/\t$/\td:Bacteria/g' $out/zotus.sintax
usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotus_table_uparse.tsv -output $out/phylum_summary_zotus.txt -rank p
usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotus_table_uparse.tsv -output $out/family_summary_zotus.txt -rank f
usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotus_table_uparse.tsv -output $out/genus_summary_zotus.txt -rank g

### Alpha and Beta Diversity ###
usearch -alpha_div $out/zotus_table_uparse.tsv -output $out/zotus-tab_alpha.txt
mkdir betaDiv-zotus
usearch -beta_div $out/zotus_table_uparse.tsv -filename_prefix betaDiv-zotus/
mv betaDiv-zotus/ $out

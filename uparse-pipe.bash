#!/bin/bash
# Robert Edgar pipeline presented at STAMPS 2017 course (Aug/2017)
# Adjusted by Elton Vasconcelos
# usage: run-uparse.bash [path-to-input-fastq-files]

if [ xusearch == x ] ; then
	echo Must set \usearch >> /dev/stderr
	exit 1
fi

version=`usearch -version | sed "-es/usearch //" | sed "-es/v10.*/v10/"`

if [ x$version != xv10 ] ; then
	echo "usearch version too old, need v10" >> /dev/stderr
	exit 1
fi

if [ ! -d $1 ] ; then
	echo "Directory $1 not found." >> /dev/stderr
	exit 1
fi

out=outputs

rm -rf $out
mkdir -p $out

#cd $out

# Assemble paired reads, put sample names into read labels
usearch -fastq_mergepairs $1/*_R1_*.fastq -relabel @ -fastqout $out/merged.fq

# Discard reads which probably have errors (quality filtering)
usearch -fastq_filter $out/merged.fq -fastq_maxee 1.0 -relabel Filt -fastaout $out/filtered.fa

# Find unique sequences and abundances
usearch -fastx_uniques $out/filtered.fa -sizeout -relabel Uniq -fastaout $out/uniques.fa

# Create 97% OTUs
usearch -cluster_otus $out/uniques.fa -relabel Otu -otus $out/otus.fa

# Create ZOTUs by denoising (error-correction)
usearch -unoise3 $out/uniques.fa -zotus $out/zotus.fa

# Create OTU table for 97% OTUs
usearch -otutab $out/merged.fq -otus $out/otus.fa -strand plus -otutabout $out/otu_table_uparse.tsv -mapout $out/otu_map.txt
# Create OTU table for ZOTUs	(ZOTUs are 100% identical OTUs)
perl -pi -e 's/Zotu/Otu/g' $out/zotus.fa	#According to work-around for bug 10.1 (https://www.drive5.com/usearch/manual/bugs.html)
usearch -otutab $out/merged.fq -zotus $out/zotus.fa -strand plus -otutabout $out/zotu_table_uparse.tsv -mapout $out/zotu_map.txt

#########################################################################################################

### Taxonomic classification (for OTUs 97%)
#./setup_sintax.bash	# Downloading RDP-16S
#usearch -makeudb_sintax sintax_16s.fa -output rdp_16s.udb	#formatting the database
usearch -sintax $out/otus.fa -db rdp_16s.udb -tabbedout $out/otus.sintax -strand both -sintax_cutoff 0.8
perl -pi -e 's/\t$/\td:Bacteria/g' $out/otus.sintax
usearch -sintax_summary $out/otus.sintax -otutabin $out/otu_table_uparse.tsv -output $out/phylum_summary_otu97.txt -rank p
usearch -sintax_summary $out/otus.sintax -otutabin $out/otu_table_uparse.tsv -output $out/family_summary_otu97.txt -rank f
usearch -sintax_summary $out/otus.sintax -otutabin $out/otu_table_uparse.tsv -output $out/genus_summary_otu97.txt -rank g
# Alpha and Beta Diversity 
usearch -alpha_div $out/otu_table_uparse.tsv -output $out/otutab_alpha.txt
#usearch -otutab_norm $out/otutab.txt -sample_size 5000 -output $out/otutab_norm.txt
#usearch -alpha_div otutab_norm.txt -output otutab_norm_alpha.txt
mkdir bDiv-otu97
usearch -beta_div $out/otu_table_uparse.tsv -filename_prefix bDiv-otu97/

## Taxonomic classification (for ZOTUs)
#usearch -sintax $out/zotus.fa -db rdp_16s.udb -tabbedout $out/zotus.sintax -strand both -sintax_cutoff 0.8
#perl -pi -e 's/\t$/\td:Bacteria/g' $out/zotus.sintax
#usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotu_table_uparse.tsv -output $out/phylum_summary_zotu.txt -rank p
#usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotu_table_uparse.tsv -output $out/family_summary_zotu.txt -rank f
#usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotu_table_uparse.tsv -output $out/genus_summary_zotu.txt -rank g
# Alpha and Beta Diversity
#usearch -alpha_div $out/zotu_table_uparse.tsv -output $out/zotutab_alpha.txt
#usearch -otutab_norm $out/otutab.txt -sample_size 5000 -output $out/otutab_norm.txt
#usearch -alpha_div otutab_norm.txt -output otutab_norm_alpha.txt
#mkdir bDiv-zotu
#usearch -beta_div $out/zotu_table_uparse.tsv -filename_prefix bDiv-zotu/

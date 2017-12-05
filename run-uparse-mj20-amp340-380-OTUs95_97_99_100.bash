#!/bin/bash
# Robert Edgar pipeline presented at STAMPS 2017 course (Aug/2017)
# Adjusted by Elton Vasconcelos
# usage: run-uparse-4plates.bash [path-to-input-fastq-files]
# The only change on this script is the use of a Mycoplasma_haemocanis-containing RDP refDB (see line 57)

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
usearch -fastq_mergepairs $1/*_R1_*.fastq -fastq_minovlen 20 -fastq_minmergelen 340 -fastq_maxmergelen 380 -relabel @ -fastqout $out/merged.fq

# Discard reads which probably have errors (quality filtering)
usearch -fastq_filter $out/merged.fq -fastq_maxee 1.0 -relabel Filt -fastaout $out/filtered.fa

# Find unique sequences and abundances (dereplication)
usearch -fastx_uniques $out/filtered.fa -sizeout -relabel Uniq -fastaout $out/uniques.fa

# Create 97% OTUs
usearch -cluster_otus $out/uniques.fa -relabel Otu -otus $out/otusDef.fa

# Create OTU table for OTUs default (97% id from the reads)
usearch -otutab $out/merged.fq -otus $out/otusDef.fa -strand plus -otutabout $out/otusDef_table_uparse.tsv -mapout $out/otusDef_map.txt
############################################################
# ZOTUs-clustered OTUs (following Robert Edgar's recommendation at https://www.drive5.com/usearch/manual/uparse_otu_radius.html)
# Create ZOTUs by denoising (error-correction)
usearch -unoise3 $out/uniques.fa -zotus $out/zotus.fa
perl -pi -e 's/Zotu/Otu/g' $out/zotus.fa
java -jar /home/elton/bioinformatics-tools/readseq.jar -inform fasta -f fasta -o $out/zotus.fa2  $out/zotus.fa
sed -r 's/ [0-9]+ bp/\;size\=&/g' $out/zotus.fa2 | sed 's/\= /\=/g' | sed 's/ bp//g' >$out/zotus.fa3
usearch -sortbysize $out/zotus.fa3 -fastaout $out/zotus_sorted.fa
rm $out/zotus.fa2 $out/zotus.fa3

# Create 99% OTUs 
usearch -cluster_smallmem $out/zotus_sorted.fa -id 0.99 -centroids $out/otus99.fa

# Create 97% OTUs
usearch -cluster_smallmem $out/zotus_sorted.fa -id 0.97 -centroids $out/otus97.fa

# Create 95% OTUs
usearch -cluster_smallmem $out/zotus_sorted.fa -id 0.95 -centroids $out/otus95.fa

# Create OTU table for 95% OTUs
usearch -otutab $out/merged.fq -otus $out/otus95.fa -strand plus -otutabout $out/otus95_table_uparse.tsv -mapout $out/otus95_map.txt

# Create OTU table for 97% OTUs
usearch -otutab $out/merged.fq -otus $out/otus97.fa -strand plus -otutabout $out/otus97_table_uparse.tsv -mapout $out/otus97_map.txt

# Create OTU table for 99% OTUs
usearch -otutab $out/merged.fq -otus $out/otus99.fa -strand plus -otutabout $out/otus99_table_uparse.tsv -mapout $out/otus99_map.txt

# Create OTU table for ZOTUs	(ZOTUs are 100% identical OTUs)
#perl -pi -e 's/Zotu/Otu/g' $out/zotus.fa
usearch -otutab $out/merged.fq -zotus $out/zotus.fa -strand plus -otutabout $out/zotus_table_uparse.tsv -mapout $out/zotus_map.txt
perl -pi -e 's/Otu/Zotu/g' $out/zotu*
#########################################################################################################

### Taxonomic classification (for default OTUs: 97%)
#./setup_sintax.bash	# Downloading RDP-16S
#usearch -makeudb_sintax sintax_16s.fa -output rdp_16s.udb	#formatting the database
usearch -sintax $out/otusDef.fa -db rdp_16s-wMhaemocanis.udb -tabbedout $out/otusDef.sintax -strand both -sintax_cutoff 0.8
perl -pi -e 's/\t$/\td:Bacteria/g' $out/otusDef.sintax
usearch -sintax_summary $out/otusDef.sintax -otutabin $out/otusDef_table_uparse.tsv -output $out/phylum_summary_otusDef.txt -rank p
usearch -sintax_summary $out/otusDef.sintax -otutabin $out/otusDef_table_uparse.tsv -output $out/family_summary_otusDef.txt -rank f
usearch -sintax_summary $out/otusDef.sintax -otutabin $out/otusDef_table_uparse.tsv -output $out/genus_summary_otusDef.txt -rank g
# Alpha and Beta Diversity 
usearch -alpha_div $out/otusDef_table_uparse.tsv -output $out/otusDef-tab_alpha.txt
#usearch -otutab_norm $out/otutab.txt -sample_size 5000 -output $out/otutab_norm.txt
#usearch -alpha_div otutab_norm.txt -output otutab_norm_alpha.txt
mkdir bDiv-otusDef
usearch -beta_div $out/otusDef_table_uparse.tsv -filename_prefix bDiv-otusDef/
##########################################################################################################

### Taxonomic classification (for ZOTUs-clustered OTUs: 95%)
#./setup_sintax.bash	# Downloading RDP-16S
#usearch -makeudb_sintax sintax_16s.fa -output rdp_16s.udb	#formatting the database
usearch -sintax $out/otus95.fa -db rdp_16s-wMhaemocanis.udb -tabbedout $out/otus95.sintax -strand both -sintax_cutoff 0.8
perl -pi -e 's/\t$/\td:Bacteria/g' $out/otus95.sintax
usearch -sintax_summary $out/otus95.sintax -otutabin $out/otus95_table_uparse.tsv -output $out/phylum_summary_otus95.txt -rank p
usearch -sintax_summary $out/otus95.sintax -otutabin $out/otus95_table_uparse.tsv -output $out/family_summary_otus95.txt -rank f
usearch -sintax_summary $out/otus95.sintax -otutabin $out/otus95_table_uparse.tsv -output $out/genus_summary_otus95.txt -rank g
# Alpha and Beta Diversity 
usearch -alpha_div $out/otus95_table_uparse.tsv -output $out/otus95-tab_alpha.txt
#usearch -otutab_norm $out/otutab.txt -sample_size 5000 -output $out/otutab_norm.txt
#usearch -alpha_div otutab_norm.txt -output otutab_norm_alpha.txt
mkdir bDiv-otus95
usearch -beta_div $out/otus95_table_uparse.tsv -filename_prefix bDiv-otus95/

### Taxonomic classification (for ZOTUs-clustered OTUs: 97%)
#./setup_sintax.bash	# Downloading RDP-16S
#usearch -makeudb_sintax sintax_16s.fa -output rdp_16s.udb	#formatting the database
usearch -sintax $out/otus97.fa -db rdp_16s-wMhaemocanis.udb -tabbedout $out/otus97.sintax -strand both -sintax_cutoff 0.8
perl -pi -e 's/\t$/\td:Bacteria/g' $out/otus97.sintax
usearch -sintax_summary $out/otus97.sintax -otutabin $out/otus97_table_uparse.tsv -output $out/phylum_summary_otus97.txt -rank p
usearch -sintax_summary $out/otus97.sintax -otutabin $out/otus97_table_uparse.tsv -output $out/family_summary_otus97.txt -rank f
usearch -sintax_summary $out/otus97.sintax -otutabin $out/otus97_table_uparse.tsv -output $out/genus_summary_otus97.txt -rank g
# Alpha and Beta Diversity 
usearch -alpha_div $out/otus97_table_uparse.tsv -output $out/otus97-tab_alpha.txt
#usearch -otutab_norm $out/otutab.txt -sample_size 5000 -output $out/otutab_norm.txt
#usearch -alpha_div otutab_norm.txt -output otutab_norm_alpha.txt
mkdir bDiv-otus97
usearch -beta_div $out/otus97_table_uparse.tsv -filename_prefix bDiv-otus97/

### Taxonomic classification (for ZOTUs-clustered OTUs: 99%)
#./setup_sintax.bash	# Downloading RDP-16S
#usearch -makeudb_sintax sintax_16s.fa -output rdp_16s.udb	#formatting the database
usearch -sintax $out/otus99.fa -db rdp_16s-wMhaemocanis.udb -tabbedout $out/otus99.sintax -strand both -sintax_cutoff 0.8
perl -pi -e 's/\t$/\td:Bacteria/g' $out/otus99.sintax
usearch -sintax_summary $out/otus99.sintax -otutabin $out/otus99_table_uparse.tsv -output $out/phylum_summary_otus99.txt -rank p
usearch -sintax_summary $out/otus99.sintax -otutabin $out/otus99_table_uparse.tsv -output $out/family_summary_otus99.txt -rank f
usearch -sintax_summary $out/otus99.sintax -otutabin $out/otus99_table_uparse.tsv -output $out/genus_summary_otus99.txt -rank g
# Alpha and Beta Diversity 
usearch -alpha_div $out/otus99_table_uparse.tsv -output $out/otus99-tab_alpha.txt
#usearch -otutab_norm $out/otutab.txt -sample_size 5000 -output $out/otutab_norm.txt
#usearch -alpha_div otutab_norm.txt -output otutab_norm_alpha.txt
mkdir bDiv-otus99
usearch -beta_div $out/otus99_table_uparse.tsv -filename_prefix bDiv-otus99/

## Taxonomic classification (for ZOTUs)
usearch -sintax $out/zotus.fa -db rdp_16s-wMhaemocanis.udb -tabbedout $out/zotus.sintax -strand both -sintax_cutoff 0.8
perl -pi -e 's/\t$/\td:Bacteria/g' $out/zotus.sintax
usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotus_table_uparse.tsv -output $out/phylum_summary_zotus.txt -rank p
usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotus_table_uparse.tsv -output $out/family_summary_zotus.txt -rank f
usearch -sintax_summary $out/zotus.sintax -otutabin $out/zotus_table_uparse.tsv -output $out/genus_summary_zotus.txt -rank g
#A lpha and Beta Diversity
usearch -alpha_div $out/zotus_table_uparse.tsv -output $out/zotus-tab_alpha.txt
#usearch -otutab_norm $out/otutab.txt -sample_size 5000 -output $out/otutab_norm.txt
#usearch -alpha_div otutab_norm.txt -output otutab_norm_alpha.txt
mkdir bDiv-zotus
usearch -beta_div $out/zotus_table_uparse.tsv -filename_prefix bDiv-zotus/

mkdir aDiv-all
mv $out/*tab_alpha.txt aDiv-all/
mv aDiv-all/ $out
mv bDiv* $out

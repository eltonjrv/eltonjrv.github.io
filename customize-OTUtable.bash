#!/usr/bin/bash
# Programmer: Elton Vasconcelos (01/Aug/2018)
# Usage: bash customize-OTUtable.bash [zotus_table_uparse.tsv] [zotus.sintax] [your_sample-metadata_full_path]
##############################################################################################
echo "Taxa" >x
cut -f 1 $1 | xargs -i grep -P '{}\t' $2 | cut -f 4 >y
cat x y >z
paste zotus_table_uparse.tsv z  >zotus_table_uparse-wTaxa.tsv
rm x y z
head -1 zotus_table_uparse-wTaxa.tsv >x
perl sampleID-to-sampleDescription.pl x $3 >y
perl -pi -e "s/\t$/\n/g" y
perl -pi -e "s/\#OTU ID/OTU_ID/g" y
tail -n $((`wc -l zotus_table_uparse-wTaxa.tsv | cut -d ' ' -f 1` - 1))	zotus_table_uparse-wTaxa.tsv >z
cat y z >zotus_table_uparse-wTaxa.tsv2		#same file with a new header (sample_Descriptions, instead of sample_IDs)
rm x y z
Rscript customize-OTUtable.R zotus_table_uparse-wTaxa.tsv2
grep -v -P '^OTU_ID' zotus_table_uparse-customized.tsv >x
rm zotus_table_uparse-customized.tsv
mv x zotus_table_uparse-customized.tsv

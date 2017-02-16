#!/usr/bin/sh

count=0;

#Find files of type .fastq, assign each to a file

for file in *.fastq.gz
do
count=$((count+1))
if [ "$count" -eq 1 ]; then
R1=$file
fi
if [ "$count" -eq 2 ]; then
R2=$file
#fi

#Make new variable to hold base file name
#R3=${R1%_S*_L001_R*_001.fastq}
R3=$(echo $file | sed -e 's/_S[0-9]*_L001_R2_001\.fastq\.gz//g')

#Run trimmomatic with mate pairs - trims adapters, quality trim, leaves paired reads
java -jar /usr/local/bin/trimmomatic-0.36.jar PE -threads 3 -trimlog logfile.log -phred33 $R1 $R2 $R3-R1p.fastq $R3-R1u.fastq $R3-R2p.fastq $R3-R2u.fastq ILLUMINACLIP:/home/elton/bioinformatics-tools/Trimmomatic-0.36/adapters/adapters.fasta:2:30:10 LEADING:6 TRAILING:6 MINLEN:30

#Run flash with mate pairs
#/opt/local/bin/flash/./flash -M 100 $R1 $R2 #for ARS laptop install
flash -M 200 -t 3 $R3-R1p.fastq $R3-R2p.fastq 

#Do some quality trimming on merged file and convert to fasta format
file=out.extendedFrags.fastq
cp out.extendedFrags.fastq $R3.flashOUT.fastq
fastq_quality_trimmer -i $file -t 20 -l 250 -Q 33 | fastq_to_fasta -Q 33 > $R3.fas


#Rename sequences sequentially
#my_123_fasta_renamer.pl $R3.fas > ${R3}_123.fas
#mv ${R3}_123.fas $R3.fas


count=0


fi
done

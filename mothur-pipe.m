# Mothur pipeline learned from Pat Schloss during "Mothur Workshop", April 2017.
# Adjustments by Elton Vasconcelos for my particular data type.
################################################################################
# Before running this pipeline with $ mothur mothur-pipe.m, I had to run make.contigs individually for each fastq pair, with the following bash command:
# $ for o in `ls oligos-prep/oligos_*`; do f=M`echo $o | sed 's/.*M//g' | sed 's/\.tab//g'`_S`echo $o | sed 's/.*M//g' | sed 's/\.tab//g'`_L001_R1_001.fastq; r=M`echo $o | sed 's/.*M//g' | sed 's/\.tab//g'`_S`echo $o | sed 's/.*M//g' | sed 's/\.tab//g'`_L001_R2_001.fastq; mothur "#make.contigs(ffastq=$f, rfastq=$r, oligos=$o, bdiffs=1, pdiffs=2, checkorient=t,  processors=6)"; done
# Merging all the outputs from the individual make.contigs runs
# mkdir mkctg-individual-outputs
# mv *.contigs.* mkctg-individual-outputs/
# mv mothur.149435* mkctg-individual-outputs/
# mv mothur.150285* mkctg-individual-outputs/
# cat mkctg-individual-outputs/*.groups >fleas.contigs.groups
# cat mkctg-individual-outputs/*.report >fleas.contigs.report
# cat mkctg-individual-outputs/*.scrap.contigs.fasta >fleas.scrap.contigs.fasta
# cat mkctg-individual-outputs/*.scrap.contigs.qual >fleas.scrap.contigs.qual
# cat mkctg-individual-outputs/*.trim.contigs.qual >fleas.trim.contigs.qual
# cat mkctg-individual-outputs/*.trim.contigs.fasta >fleas.trim.contigs.fasta
###############################################################################
### Running mothur whole pipeline
### Removing undesireable sequences
summary.seqs(fasta=fleas.trim.contigs.fasta)
get.current()
screen.seqs(fasta=fleas.trim.contigs.fasta, group=fleas.contigs.groups, maxambig=0, maxlength=552, processors=6)
get.current()
summary.seqs()
unique.seqs(fasta=current)
get.current()
count.seqs(name=current, group=current)
summary.seqs(count=current)
align.seqs(fasta=current, reference=silva.bacteria.fasta)
summary.seqs(fasta=current, count=current)
get.current()
screen.seqs(fasta=current, count=current, summary=current, start=13875, end=27654, maxhomop=8)
summary.seqs(fasta=current, count=current)
filter.seqs(fasta=current, vertical=T, trump=.)
summary.seqs(fasta=current, count=current)
unique.seqs(fasta=current, count=current)
summary.seqs(fasta=current, count=current)
pre.cluster(fasta=current, count=current, diffs=5)
summary.seqs(fasta=current, count=current)
chimera.vsearch(fasta=current, count=current, dereplicate=t)
remove.seqs(fasta=current, accnos=current)
summary.seqs(fasta=current, count=current)

### Classifying sequences with RDP
classify.seqs(fasta=current, count=current, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax, cutoff=80)
remove.lineage(fasta=current, count=current, taxonomy=current, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.tax(taxonomy=current, count=current)
get.current()

### Assigning error rates (if one has a Mock community sample)
#get.groups(count=current, fasta=current, groups=Mock)
#seq.error(fasta=current, count=current, reference=HMP_MOCK.v35.fasta, aligned=F)
#remove.groups(count=current, fasta=current, taxonomy=current, groups=Mock)

### Preparing for analysis (catching OTUs)
#get.current()
dist.seqs(fasta=current, cutoff=0.20)
get.current()
cluster(column=current, count=current)
cluster.split(fasta=current, count=current, taxonomy=current, splitmethod=classify, taxlevel=4, cutoff=0.03)
get.current()
make.shared(list=current, count=current)   # OTU table generation (the *.shared outfile is an otu table)
classify.otu(list=current, count=current, taxonomy=current)        # OTUs classification (two outputs: *cons.tax.summary and *cons.taxonomy)

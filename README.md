# Customized RDP reference database
The original fasta file was downloaded from http://drive5.com/sintax/rdp_16s_v16_sp.fa.gz

Then, after uncompressing with "gunzip" command, we added, at the end of the file, the following 16S rRNA sequences from SILVA-DB: 
>*Mycoplasma haemocanis*  (SILVA ID: H0HHaemo)

>*Ehrlichia canis* (SILVA ID: I8UCani3)

>*Anaplasma platys* (SILVA ID: IE3Plat5)

Those are vector-borne pathogens of interest in Vet Med that were missing at the original RDP file from Robert Edgar (Uparse developer).

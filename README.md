# Customized RDP reference database
The original fasta file was downloaded from http://drive5.com/sintax/rdp_16s_v16_sp.fa.gz

Then, after uncompressing with "gunzip" command, we appended, at the end of the file, 16S rRNA sequences from the following species: 

>*Ehrlichia canis* and *Ehrlichia chafeensis*

>*Anaplasma platys* and *Anaplasma phagocytophilum*

>*Mycoplasma haemocanis* and *Mycoplasma haematoparvum*

Those are vector-borne pathogens of interest in Vet Med that were missing at the original RDP file from Robert Edgar (Uparse developer). Those sequences are exactly the same as inserted in our positive controls plasmids for the multi-infection simulations. Please see Methods section "Positive controls and standard curves" subsection for details on GenBank accesssion numbers (Vasconcelos et al., 2020 - in preparation). M. haematoparvum sequence is derived from GenBank AY383241.1.

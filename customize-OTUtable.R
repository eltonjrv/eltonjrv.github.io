# Programmer: Elton Vasconcelos (17/Oct/2017)
# Usage: Rscript customize-OTUtable.R [zotus_table_uparse-wTaxa.tsv2]
#####################################################################
args <- commandArgs(TRUE)
x = read.table(args[1], header=T, check.names=F)
out = gsub("wTaxa.tsv2", "customized.tsv", args[1])
for(j in 2:dim(x)[2]-1) { 
	for(i in 1:length(x[,j][which(x[,j] > 0)])) { 
		y = matrix(c(colnames(x)[j], as.character(x[,1][which(x[,j] > 0)][i]), as.character(x[,dim(x)[2]][which(x[,j] > 0)][i]), x[,j][which(x[,j] > 0)][i]), nrow=1, ncol=4); 
		write.table(y, file = out, sep = "\t", append = T, quote=F, row.names=F, col.names=F);
	}
}


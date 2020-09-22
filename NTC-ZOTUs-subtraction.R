### Programmer: Chayan Roy (Jul/2020)
### Usage: Rscript NTC-ZOTUs-subtraction.R zotus_table_uparse-customized.tsv

args <- commandArgs(TRUE)
# Take input (change your input file name)
My.Data <- read.table(args[1], header= F, sep = "\t", check.names = F)

# process
control<-My.Data[grep("^H2O*", My.Data$V1), ]
sample<-My.Data[!grepl("^H2O*",My.Data$V1),]
library(data.table)
setDT(sample)[control, V4 := V4 - i.V4, on = .(V2)]
outDT<-rbind(sample, setDT(control)[!sample, on = .(V2)])
outDT_positives <- outDT[outDT$V4 >= 0, ]

# Output (change your output file name)
write.table(outDT_positives, "output.tsv")

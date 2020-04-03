#!/usr/bin/Rscript

fam=as.character(commandArgs(TRUE)[1])
pop=as.character(commandArgs(TRUE)[2])
output=as.character(commandArgs(TRUE)[3])

fam_file=read.table(fam,header = F,stringsAsFactors = F)
colnames(fam_file)[2]="samples"
pop_file=read.table(pop,header = F,stringsAsFactors = F)
colnames(pop_file)[1]="samples"

library(dplyr)
merged_data=left_join(fam_file[,"samples",drop=F],pop_file,by="samples")
merged_data[is.na(merged_data[,3]),3]="-"
merged_data[is.na(merged_data[,2]),2]="-"

write.table(data.frame(merged_data$V3),output,quote = F,row.names = F, col.names = F)

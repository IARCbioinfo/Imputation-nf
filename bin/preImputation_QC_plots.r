#!/usr/bin/Rscript

library(data.table)
library(ggplot2)
library(gridExtra)
library(parallel)
library(hexbin)

pop=as.character(commandArgs(TRUE)[1]) #"CHB_JPT"  #"CEU" "YRI"
pop_ref=as.character(commandArgs(TRUE)[2]) #"EAS" #"EUR" "AFR"

pop="CHB_JPT"  #"CEU" "YRI"
pop_ref="EAS" #"EUR" "AFR"


af_diff=0.2
af_fc=5

# Read in the frequency files
chip <- fread(paste0("target_freq_",pop,".frq"), header = TRUE)
panel <- fread("ref_freq_withHeader.txt", header = TRUE)

custom_f=function(id,allele){
  ref=unlist(strsplit(id,":"))[3]
  alt=unlist(strsplit(id,":"))[4]
  return(get(allele))
}
panel[, ref := mcmapply(custom_f,id,"ref",mc.cores=5L)]
panel[, alt := mcmapply(custom_f,id,"alt",mc.cores=5L)]
head(panel)
summary(panel$ref==panel$a0)

# read SNPs that have ben filtered out
filtered_SNPs=read.table("filtered_snps.txt",header = F)
chip=chip[-which(chip$SNP %in% filtered_SNPs$V1 ),]

# Read affy annotation to retrieve SNPs positions and   -------------------
probes_IDs=fread(paste0("ID-target_",pop,"-1000G.txt"),header = F) #ancestry_run_04022020/ID-all_cohorts_freq-1000G.txt = ID-target5-1000G.txt
head(probes_IDs)
colnames(probes_IDs)=c("affy_IDs","ref_IDs")

flip_als=c("A","T","C","G")
names(flip_als)=c("T","A","G","C")
chip_updated=merge(chip,probes_IDs,by.x="SNP",by.y="affy_IDs",all.x=T)
chip_updated[A1 == "A", A1_flip := "T"]
chip_updated[A1 == "T", A1_flip := "A"]
chip_updated[A1 == "G", A1_flip := "C"]
chip_updated[A1 == "C", A1_flip := "G"]

chip_updated[A2 == "A", A2_flip := "T"]
chip_updated[A2 == "T", A2_flip := "A"]
chip_updated[A2 == "G", A2_flip := "C"]
chip_updated[A2 == "C", A2_flip := "G"]

# check if snp is ambiguous
custom_f=function(a1,a2){
  strand_case=ifelse(paste(sort(c(a1,a2)),collapse = ":") %in% c("C:G","A:T"),T,F )
  if(strand_case ){return(T)
  }else{return(F)}
}
chip_updated[, strand_pb := mcmapply(custom_f,A1,A2,mc.cores=5L)]
chip_updated=chip_updated[!which(strand_pb),]

for (i in 1:length(chip_updated$ref_IDs)){
  if(is.na(chip_updated$ref_IDs[i])){
    chip_updated$ref_IDs[i]<-paste0(chip_updated$SNP,':',chip_updated$CHR,':',chip_updated$A1,':',chip_updated$A2)
  }
}
  

# Flip MAF ref don't match the ref of the reference
custom_f=function(id,a1,a1_flip,maf){
  ref=unlist(strsplit(id,":"))[4]
  if(a1!=ref & a1_flip!=ref){
    return(1-maf)
  }else{
    return(maf)}
}
chip_updated[, flip_MAF := mcmapply(custom_f,ref_IDs,A1,A1_flip,MAF,mc.cores=2L)]

# Take an intersection of the panel and chip data
isec <- merge(chip_updated[,.SD,.SDcols=c(7,11)], panel[,.SD,.SDcols=c("id",pop_ref)], by.x="ref_IDs",by.y="id",all.x=T)
colnames(isec)
colnames(isec)[c(3,2)] <- c("AF_PANEL", "AF_CHIP")

isec$AF_diff <- abs(isec$AF_PANEL - isec$AF_CHIP)
isec$AF_ratio <- isec$AF_CHIP/isec$AF_PANEL
isec$AF_FC <- log2(isec$AF_ratio)

# Check if AFs are within the pp and fold change ranges
af_ok <- (isec$AF_diff < af_diff)
exclude <- !af_ok
# id_test=isec$ref_IDs[which(af_ok==F & isec$AF_PANEL<0.2 & grepl("A:G",isec$ref_IDs))][4]
# panel[which(panel$id==id_test),]
# chip_updated[which(chip_updated$ref_IDs==id_test),]

# Generate an exclusion list for variants not in the panel
nonpanel <- chip[!(chip$SNP) %in% (isec$SNP)]
# Generate AF list for discordant variants
af_discrepant <- isec[exclude]

# Save the plots
pdf(paste0("AF_beforeImputation_",pop_ref,".pdf"),w=10,h=5)

p1 = ggplot(isec,aes(x=AF_PANEL,y=AF_CHIP)) + stat_binhex(bins = 500)+
  labs(title=paste0("Target AF vs reference AF (",pop_ref," ancestry)"),x="AF in ref data", y = "AF in target data")+
  theme(legend.position = c(0.8, 0.2))

p2 = ggplot(isec[!exclude], aes(x=AF_CHIP)) +
  geom_histogram(color="black", fill="white",bins=100)+
  labs(title=paste0("AF in target data (",pop_ref," ancestry)"),x="AF in target data", y = "Count")

grid.arrange(p1, p2, nrow = 1)

dev.off()

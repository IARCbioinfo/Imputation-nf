#!/usr/bin/Rscript
library(BSgenome.Hsapiens.UCSC.hg38)
seqlengths(Hsapiens)
options("scipen"=100, "digits"=4)
chr=as.integer(commandArgs(TRUE)[1])
print(seqlengths(Hsapiens)[chr])
chunks=data.frame(strat=vector(),end=vector())
positions=seq(1,seqlengths(Hsapiens)[chr],20000000) #19999999

zz=gzfile(paste0("isec_chr_",chr,".vcf.gz"))
chr_map=read.table(zz,header = F)
i=1
inc=1
while(i<=length(positions)){
  if(i==1){
    chr_map_subset=chr_map[which(chr_map$V2>positions[i] & chr_map$V2< positions[i+inc]-1),]
    dim(chr_map_subset)
    if(nrow(chr_map_subset)>1000){
      start=positions[i]
      end=positions[i+inc]-1
      i=i+inc
      inc=1
      chunks=rbind(chunks,data.frame(strat=start,end=end))
    }else{inc=inc+1}
  }else if(i!=length(positions)){
    chr_map_subset=chr_map[which(chr_map$V2>positions[i] & chr_map$V2< positions[i+inc]-1),]
    dim(chr_map_subset)
    if(nrow(chr_map_subset)>1000){
      start=positions[i]
      end=positions[i+inc]-1
      i=i+inc
      inc=1
      chunks=rbind(chunks,data.frame(strat=start,end=end))
    }else{inc=inc+1}
  }else{
    chr_map_subset=chr_map[which(chr_map$V2>positions[i] & chr_map$V2< seqlengths(Hsapiens)[chr]),]
    dim(chr_map_subset)
    if(nrow(chr_map_subset)>1000){
      start=positions[i]
      end=seqlengths(Hsapiens)[chr]
      i=i+inc
      inc=1
      chunks=rbind(chunks,data.frame(strat=start,end=end))
    }else{
      chunks[nrow(chunks),]=c(positions[i-1]+1,seqlengths(Hsapiens)[chr])
      i=i+1
    }
  }
}
write.table(chunks,file=paste0("chunk_split_chr",chr,".txt"),quote=F,row.names = F,col.names = F)

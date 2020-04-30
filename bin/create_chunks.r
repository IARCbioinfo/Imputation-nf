#!/usr/bin/Rscript
library(BSgenome.Hsapiens.UCSC.hg19)
seqlengths(Hsapiens)

chr=as.integer(commandArgs(TRUE)[1])

print(seqlengths(Hsapiens)[chr])
chunks=data.frame(strat=vector(),end=vector())
positions=seq(1,seqlengths(Hsapiens)[chr],20000000)
for(i in 1:length(positions)){
  if(i==1){chunks=rbind(chunks,data.frame(strat=positions[i],end=positions[i+1]))
  }else if(i!=length(positions)){chunks=rbind(chunks,data.frame(strat=positions[i]+1,end=positions[i+1]))
  }else{chunks=rbind(chunks,data.frame(strat=positions[i]+1,end=seqlengths(Hsapiens)[chr]))
  }
}

write.table(chunks,file=paste0("chunk_split_chr",chr,".txt"),quote=F,row.names = F,col.names = F)

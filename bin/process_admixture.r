#!/usr/bin/Rscript

K=3
admixture_Q=read.table(paste0("merge.",K,".Q"),header = F)
n_groups=ncol(admixture_Q)

target_name = as.character(commandArgs(TRUE)[1])
target_name='NA'
ref_samples=read.table("ref_samples.txt",stringsAsFactors=F,header=F)

infered_ancestry_samples=read.table("merge.fam",stringsAsFactors=F,header=F)
admixture_Q$samples_cohort=infered_ancestry_samples[,1]
admixture_Q$samples_ids=infered_ancestry_samples[,2]
# admixture_Q$participant=sapply(infered_ancestry_samples[,2],function(i) paste(unlist(strsplit(i,"-"))[1:3],collapse = ".") )
# retrieve admixture groups
admixture_Q$group_admixture=sapply(1:nrow(admixture_Q),function(i) which.max(admixture_Q[i,1:n_groups]))
admixture_Q$certainty=sapply(1:nrow(admixture_Q), function(i) admixture_Q[i,order(admixture_Q[i,1:n_groups],decreasing = T)[1]]/admixture_Q[i,order(admixture_Q[i,1:n_groups],decreasing = T)[2]] >=1.5)
table(admixture_Q$group_admixture)
table(admixture_Q$group_admixture[admixture_Q$certainty])

# save admixture groups
ref_res=admixture_Q[which(admixture_Q$samples_ids %in% ref_samples$V1 ),]
# read B ancestry file
rel_file=read.table("relationships_w_pops_121708_2.txt",header = T,stringsAsFactors = F)
ref_res=merge(ref_res,rel_file[,c("IID","population")],by.x = "samples_ids",by.y = "IID")
table(ref_res$group_admixture)
table(ref_res$population)

list_origins=list()
for(i in 1:K){
  print(unique(ref_res$population[which(ref_res$group_admixture==i)]))
  list_origins[[i]]=unique(ref_res$population[which(ref_res$group_admixture==i)])
}

target_res=admixture_Q[which(grepl(target_name ,admixture_Q$samples_ids)==FALSE),]
group_summary=table(target_res$group_admixture)
table(target_res$group_admixture)
table(target_res$group_admixture[target_res$certainty])

target_res$admixture_ancestry=NA
for(i in 1:K){
  print(summary(target_res[which(target_res$group_admixture==i),paste0("V",i)]))
  target_res$admixture_ancestry[target_res$group_admixture==i]=paste(list_origins[[i]],collapse = "_")
}
table(target_res$admixture_ancestry)

colnames(target_res)[1:K]=sapply(1:K,function(i) paste0("p",i))
write.table(target_res,"admixture_results_withGroups.txt",row.names = F,col.names = T,quote = F)

# select samples IDs in the biggest group
biggest_group=names(group_summary)[which.max(group_summary)]
write.table(target_res[which(target_res$group_admixture==as.numeric(biggest_group)),c("samples_cohort","samples_ids")],file="samples_principalGroup.txt",row.names = F,col.names = F,quote = F)

# select samples IDs in the biggest group
write.table(target_res[which(target_res$group_admixture==as.numeric(biggest_group) & target_res$certainty==T),c("samples_cohort","samples_ids")],file="samples_principalGroup_withCertainty.txt",row.names = F,col.names = F,quote = F)

system(command = "mkdir out_pop_admixture")
for(pop in unique(target_res$admixture_ancestry)){
  subset_pop=target_res[which(target_res$admixture_ancestry==pop & target_res$certainty==T),c("samples_cohort","samples_ids")]
  print(dim(subset_pop))
  system(command = paste0("mkdir out_pop_admixture/1000G_checking_",pop,"/"))
  write.table(subset_pop,file=paste0("out_pop_admixture/1000G_checking_",pop,"/","samples_",pop,"_withCertainty.txt"),row.names = F,col.names = F,quote = F)

  subset_pop=target_res[which(target_res$admixture_ancestry==pop ),c("samples_cohort","samples_ids")]
  write.table(subset_pop,file=paste0("out_pop_admixture/1000G_checking_",pop,"/","samples_",pop,".txt"),row.names = F,col.names = F,quote = F)
}

# return all samples
subset_pop=target_res[which(target_res$certainty==T),c("samples_cohort","samples_ids")]
print(dim(subset_pop))
system(command = paste0("mkdir out_pop_admixture/1000G_checking_ALL/"))
write.table(subset_pop,file=paste0("out_pop_admixture/1000G_checking_ALL/","samples_ALL_withCertainty.txt"),row.names = F,col.names = F,quote = F)

subset_pop=target_res[c("samples_cohort","samples_ids")]
write.table(subset_pop,file=paste0("out_pop_admixture/1000G_checking_ALL/","samples_ALL.txt"),row.names = F,col.names = F,quote = F)

#!/usr/bin/Rscript
rm(list=ls())

# Libraries ---------------------------------------------------------------
library(readxl)
library(viridis)
library(grid)
library(ggplot2)
library(gridExtra)

# Parameters --------------------------------------------------------------
#path=as.character(commandArgs(TRUE)[1]);

# Read previous infos on the selected samples -----------------------------
selected_samples=read.table("admixture_results_withGroups.txt",header = F,sep=" ",dec=",")

samples_to_annotate=c()
list_points=list()
list_hist=list()
selected_samples$sex_checking_F=NA
selected_samples$SNP_SEX=NA
for(pop in c("ALL")){
  #qc_path=paste0("out_pop_admixture/1000G_checking_",pop,"/")

  # Read sex checking results -----------------------------------------------
  sex_checking=read.table(paste0("target_sexCheck_",pop,".sexcheck"),header = T,stringsAsFactors = F)
  sex_checking_problems=sex_checking[which(sex_checking$STATUS=="PROBLEM"),]
  print("PEDSEX distributions in the detected problems")
  print(table(sex_checking_problems$PEDSEX))
  print("SNPSEX distributions in the detected problems")
  print(table(sex_checking_problems$SNPSEX))
  print(table(sex_checking_problems$SNPSEX[which(sex_checking_problems$PEDSEX!=0)]))
  samples_sex_problems=sex_checking_problems$IID

  sex_checking_merged=merge(sex_checking,selected_samples,by.x = "IID",by.y = "V5",all.x = T)
  # col_sex=rep(adjustcolor("blue",alpha.f = 0.7),nrow(sex_checking_merged))
  # col_sex[which(sex_checking_merged$PEDSEX==1)]=adjustcolor("red",alpha.f = 0.7)
  # col_sex[which(sex_checking_merged$PEDSEX==2)]=adjustcolor("black",alpha.f = 0.7)
  # plot(1:nrow(sex_checking_merged),sex_checking_merged$F,col=col_sex,pch=16,main=pop,xlab="samples",ylab="Fstat")
  # legend('topright', title = 'Reported sex', legend = c("Missing","Male","Female"),col=c(adjustcolor("blue",alpha.f = 0.7),adjustcolor("red",alpha.f = 0.7),adjustcolor("black",alpha.f = 0.7)), pch=21)
  # abline(h=c(0.2,0.8))

  legend_size=2
  text_size=8
  sex_checking_merged$PEDSEX=factor(sex_checking_merged$PEDSEX, levels = c(0,1,2), labels = c("Missing", "Male", "Female"))
  list_points[[paste0(pop,"_sex")]]=grid.arrange(ggplot(sex_checking_merged,aes(x=1:nrow(sex_checking_merged),y=F,col=PEDSEX)) + geom_point() +
    labs(title=pop,x="Samples", y = "F") +
    scale_colour_viridis_d(alpha = 0.5) +
    theme(legend.title = element_text(size = text_size),
          legend.text  = element_text(size = text_size),
          legend.key.size = unit(0.5, "lines"),legend.position = c(0.8, 0.7)),nrow=1 )

  # col_sex=rep(adjustcolor("blue",alpha.f = 0.7),nrow(sex_checking_merged))
  # col_sex[which(sex_checking_merged$SNPSEX==1)]=adjustcolor("red",alpha.f = 0.7)
  # col_sex[which(sex_checking_merged$SNPSEX==2)]=adjustcolor("black",alpha.f = 0.7)
  # plot(1:nrow(sex_checking_merged),sex_checking_merged$F,col=col_sex,pch=16,main=pop,xlab="samples",ylab="Fstat")
  # legend('topright', title = 'SNPs sex', legend = c("Missing","Male","Female"),col=c(adjustcolor("blue",alpha.f = 0.7),adjustcolor("red",alpha.f = 0.7),adjustcolor("black",alpha.f = 0.7)), pch=21)
  # abline(h=c(0.2,0.8))

  list_hist[[paste0(pop,"_sex")]]=grid.arrange(ggplot(sex_checking_merged, aes(x=F)) +
    geom_histogram(color="black", fill="white",bins=50)+
    labs(title=pop,x="F", y = "Frequency"),nrow=1)

  sex_checking_problems_merged=merge(sex_checking_problems,selected_samples,by.x = "IID",by.y = "V5",all.x = T)
  table(sex_checking_problems_merged$samples_cohort)
  table(sex_checking_problems_merged$samples_cohort)/table(selected_samples$samples_cohort)

  samples_to_annotate=c(samples_to_annotate,sex_checking_problems$IID)

  selected_samples$sex_checking_F[which(selected_samples$V5 %in% sex_checking$IID)]=sapply(selected_samples$V5[which(selected_samples$V5 %in% sex_checking$IID)],function(i) sex_checking[which(sex_checking$IID == i),"F"])
  selected_samples$SNP_SEX[which(selected_samples$V5 %in% sex_checking$IID)]=sapply(selected_samples$V5[which(selected_samples$V5 %in% sex_checking$IID)],function(i) sex_checking[which(sex_checking$IID == i),"SNPSEX"])
}

selected_samples$sex_checking_problem=ifelse(selected_samples$V5 %in% samples_to_annotate,yes = T,no = F)
selected_samples$sex_checking_problem[which(selected_samples$apt_QC_removed==T | selected_samples$callRate_removed==T)]=NA
table(selected_samples$sex_checking_problem,useNA = "always")

# Read relatedness checking results ---------------------------------------
samples_to_annotate_rel=c()
summary_rel=c()
selected_samples$relative=NA
selected_samples$PI_HAT=NA
for(pop in c("ALL")){
  #qc_path=paste0("out_pop_admixture/1000G_checking_",pop,"/")

  relatedness=read.table(paste0("target_rel_",pop,".genome"),header = T,stringsAsFactors = F)
  if(nrow(relatedness)!=0){
    all_related_samples=unique(c(relatedness$IID1,relatedness$IID2))
    summary_rel=rbind(summary_rel,relatedness)
    samples_to_annotate_rel=c(samples_to_annotate_rel,all_related_samples)
    for(i in 1:nrow(relatedness)){
      selected_samples$relative[which(selected_samples$V5 == relatedness$IID1[i])]=relatedness$IID2[i]
      selected_samples$relative[which(selected_samples$V5 == relatedness$IID2[i])]=relatedness$IID1[i]
      selected_samples$relative[which(selected_samples$V5 == relatedness$IID1[i])]=relatedness$IID2[i]
      selected_samples$PI_HAT[which(selected_samples$V5 %in% c(relatedness$IID1[i],relatedness$IID2[i]))]=relatedness$PI_HAT[i]
    }
  }
}

selected_samples$relatives_problem=ifelse(selected_samples$V5 %in% samples_to_annotate_rel,yes = T,no = F)
selected_samples$relatives_problem[which(selected_samples$apt_QC_removed==T | selected_samples$callRate_removed==T)]=NA
table(selected_samples$relatives_problem,useNA = "always")


# Heterozygosity and missingness ------------------------------------------
selected_samples$prop_HET=NA
selected_samples$F_MISS=NA
list_het=list()
for(pop in c("ALL")){
  #qc_path=paste0("out_pop_admixture/1000G_checking_",pop,"/")

  het_imiss=read.table(paste0("het_",pop,".imiss.txt"),header = T,stringsAsFactors = F)
  #colnames(het_imiss)=c('FID','IID','F_MISS','prop_HET')
  print(summary(het_imiss$F_MISS))

  het=read.table(paste0("het_",pop,".het"),header = T,stringsAsFactors = F)
  print(summary(het$F))
  print(summary(het$F>0.2))

  het_imiss=merge(het_imiss,het[,c("IID","F")],by="IID")
  colnames(het_imiss)[which(colnames(het_imiss)=="F")]="F_stat"

  selected_samples$prop_HET[which(selected_samples$V5 %in% het_imiss$IID)]=sapply(selected_samples$V5[which(selected_samples$V5 %in% het_imiss$IID)],function(i) het_imiss[which(het_imiss$IID == i),"prop_HET"])
  selected_samples$F_MISS[which(selected_samples$V5 %in% het_imiss$IID)]=sapply(selected_samples$V5[which(selected_samples$V5 %in% het_imiss$IID)],function(i) het_imiss[which(het_imiss$IID == i),"F_MISS"])

  list_het[[paste0(pop,"_het")]] = grid.arrange(ggplot(het_imiss,aes(x=F_MISS,y=prop_HET,col=F_stat)) + geom_point() + ylim(0,0.5)+
    labs(title=paste0(pop),x="Samples missingness", y = "Heterozigosity rate")+
    theme(legend.title = element_text(size = text_size),
          legend.text  = element_text(size = text_size),
          legend.key.size = unit(0.5, "lines"),legend.position = c(0.85, 0.2)),nrow=1)

}

write.table(selected_samples,file="selected_samples_afterSamplesQCsChecking.txt",row.names = F,col.names = T,quote = F,sep="\t",dec=",")

pdf("postGenotyping_samples_QCs.pdf",w=10,h=9)
grid.arrange(grobs=c(list_hist,list_points,list_het),nrow=1)
dev.off()

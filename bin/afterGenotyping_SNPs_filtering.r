#!/usr/bin/Rscript
library(data.table)
library(limma)

# Parameters --------------------------------------------------------------
cluster_path="/home/gabriela/remoteDir"
# cluster_path=""
project_path=paste0(cluster_path,"/data/gep/MR_Signatures/work/gabriela/GIT_Rstudio_project/")
output_path_figures=paste0(project_path,"results/figures/imputation_QCs/")

# Check_bim ---------------------------------------------------------------

pop_QC_path=paste0(cluster_path,"/data/gep/MR_Signatures/work/gabriela/imputation/January_2020/plink_files/Genotyping_R2/combined/")

summary=data.frame(pop=vector(),n_exclude_SNPs=vector(),n_Chr_SNPs=vector(),n_StrandFlip_SNPs=vector(),n_Pos_SNPs=vector())
summary_withoutFreqFilt=data.frame(pop=vector(),n_exclude_SNPs=vector(),n_Chr_SNPs=vector(),n_StrandFlip_SNPs=vector(),n_Pos_SNPs=vector())

pop_dependent_files=c("Exclude-all_cohorts_subpop_split-1000G.txt")
non_pop_dependent_files=c("Chromosome-all_cohorts_subpop_split-1000G.txt","Strand-Flip-all_cohorts_subpop_split-1000G.txt","Position-all_cohorts_subpop_split-1000G.txt")

exclude_snps_allPop=c()
for(pop in c("CEU","CHB_JPT","YRI")){
  freq_filter="withFreqFiltering/"
  nrow_check_files=c(pop)
  for(file in pop_dependent_files){
    nrow_check_files=c(nrow_check_files,length(readLines(paste0(pop_QC_path,"1000G_",pop,"checking/",freq_filter,file))) )
    assign(paste0("exclude_snps_",pop),read.table(paste0(pop_QC_path,"1000G_",pop,"checking/",freq_filter,file),stringsAsFactors = F)[,1])
    exclude_snps_allPop=c(exclude_snps_allPop,get(paste0("exclude_snps_",pop)))
  }
  for(file in non_pop_dependent_files){
    nrow_check_files=c(nrow_check_files,length(readLines(paste0(pop_QC_path,"1000G_",pop,"checking/",freq_filter,file))) )
  }
  names(nrow_check_files)=c("pop","n_exclude_SNPs","n_Chr_SNPs","n_StrandFlip_SNPs","n_Pos_SNPs")
  summary=rbind(summary,data.frame(t(nrow_check_files),stringsAsFactors = F))

  freq_filter="withoutFreqFiltering/"
  nrow_check_files=c(pop)
  for(file in pop_dependent_files){
    nrow_check_files=c(nrow_check_files,length(readLines(paste0(pop_QC_path,"1000G_",pop,"checking/",freq_filter,file))) )
  }
  for(file in non_pop_dependent_files){
    nrow_check_files=c(nrow_check_files,length(readLines(paste0(pop_QC_path,"1000G_",pop,"checking/",freq_filter,file))) )
  }
  names(nrow_check_files)=c("pop","n_exclude_SNPs","n_Chr_SNPs","n_StrandFlip_SNPs","n_Pos_SNPs")
  summary_withoutFreqFilt=rbind(summary_withoutFreqFilt,data.frame(t(nrow_check_files),stringsAsFactors = F))

}
as.numeric(summary$n_exclude_SNPs) - as.numeric(summary_withoutFreqFilt$n_exclude_SNPs)

sum_all_probes=data.frame(all_snps=unique(exclude_snps_allPop))
sum_all_probes$CEU=ifelse(sum_all_probes$all_snps %in% exclude_snps_CEU,T,F)
sum_all_probes$CHB_JPT=ifelse(sum_all_probes$all_snps %in% exclude_snps_CHB_JPT,T,F)
sum_all_probes$YRI=ifelse(sum_all_probes$all_snps %in% exclude_snps_YRI,T,F)

pdf(paste0(output_path_figures,"Probes_filtering_PerAncestryCheckBim.pdf"))
venn_obj=vennCounts(data.frame(sum_all_probes$CEU,sum_all_probes$CHB_JPT,sum_all_probes$YRI))
vennDiagram( venn_obj,names = c("CEU" , "CHB_JPT" , "YRI"))
dev.off()

# Hwe filter -----------------------------------------------------
hwe_filtered_probes=c()
for(pop in c("CEU","CHB_JPT","YRI")){
  before_filtering=fread(paste0(pop_QC_path,"1000G_",pop,"checking/all_cohorts_subpop_split.bim"),header = F,stringsAsFactors = F)
  probes_before_filtering=setDF(before_filtering[,2])[,1]
  after_filtering=fread(paste0(pop_QC_path,"1000G_",pop,"checking/all_cohorts_subpop_hwe.bim"),header = F,stringsAsFactors = F)
  probes_after_filtering=setDF(after_filtering[,2])[,1]

  assign(paste0("hwe_filtered_probes_",pop),probes_before_filtering[ which( !(probes_before_filtering %in% probes_after_filtering)) ])
  hwe_filtered_probes=c(hwe_filtered_probes,get(paste0("hwe_filtered_probes_",pop)))
}

sum_all_probes_hwe_filtered=data.frame(all_snps=unique(hwe_filtered_probes),stringsAsFactors = F)
sum_all_probes_hwe_filtered$CEU=ifelse(sum_all_probes_hwe_filtered$all_snps %in% hwe_filtered_probes_CEU,T,F)
sum_all_probes_hwe_filtered$CHB_JPT=ifelse(sum_all_probes_hwe_filtered$all_snps %in% hwe_filtered_probes_CHB_JPT,T,F)
sum_all_probes_hwe_filtered$YRI=ifelse(sum_all_probes_hwe_filtered$all_snps %in% hwe_filtered_probes_YRI,T,F)

pdf(paste0(output_path_figures,"Probes_filtering_PerAncestryGenoHWE.pdf"))
venn_obj=vennCounts(data.frame(sum_all_probes_hwe_filtered$CEU,sum_all_probes_hwe_filtered$CHB_JPT,sum_all_probes_hwe_filtered$YRI))
vennDiagram( venn_obj,names = c("CEU" , "CHB_JPT" , "YRI"))
dev.off()

summary(sum_all_probes_hwe_filtered$all_snps %in% sum_all_probes$all_snps)
summary(sum_all_probes$all_snps %in% sum_all_probes_hwe_filtered$all_snps)

length(unique(c(as.character(sum_all_probes_hwe_filtered$all_snps) , as.character(sum_all_probes$all_snps))))


filtered_snps=unique(c(as.character(sum_all_probes_hwe_filtered$all_snps) , as.character(sum_all_probes$all_snps)))
length(filtered_snps)

write.table(data.frame(filtered_snps),file=paste0(pop_QC_path,"filtered_snps.txt"),quote = F,col.names = F,row.names = F)

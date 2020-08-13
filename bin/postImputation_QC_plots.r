#!/bin/env Rscript --no-save

# Required libraries
library(data.table) 
library(ggplot2)         
library(gridExtra)
library(viridis)
library(doParallel)
library(ggplotify)

# Input variables
args <- commandArgs(TRUE)
pop <- args[1]       
INFO_thr <- as.numeric(args[2])

ref_AF=c("EUR_AF","EUR_AF","AFR_AF","EAS_AF")
names(ref_AF)=c("ALL","CEU","YRI","CHB_JPT")

ref_AF_col=ref_AF[pop]
if(pop=="ALL"){pop=""}else{pop=paste0("_",pop)}

cl = makeCluster(22)
registerDoParallel(cl)

par_res <-  foreach(chr=c(1:22),.packages = c("data.table","ggplot2","gridExtra","viridis","ggplotify")) %dopar%{ 
  panel.frq=paste0("1000GP_chr",chr,"_imputation_all.frq")
  infile=paste0("INFO_group_chr",chr,pop,".txt")
  
  panel <- fread(panel.frq, header = TRUE)
  imp_vars <- fread((infile), header = TRUE)
  colnames(imp_vars)[which(colnames(imp_vars) == "AF")] = "AF_target"
  colnames(panel)[which(colnames(panel) == ref_AF_col)] = "AF_ref"
  
  # Create random sample for INFO-value plot
  rand1 <- sample(1:nrow(imp_vars), 100000, replace = FALSE)
  temp1 <- imp_vars[rand1,]
  temp1$AF_GROUP <- factor(temp1$AF_GROUP, levels = c("1", "2","3"), labels = c("MAF > 5%", "MAF 0.5-5%", "MAF < 0.5%"))

  # Merge by common variants ------------------------------------------------
  # SNP column in format CHR_POS_REF_ALT
  isec <- merge(panel[,.SD,.SDcols=c("SNP","AF_ref")], imp_vars, by="SNP")
  
  # Plot absolute AF difference as imputed AF - panel AF
  isec$POS <- as.numeric(as.character(data.frame(do.call('rbind', strsplit(as.character(isec$SNP),'_',fixed=TRUE)))[,2]))
  
  isec=isec[isec$INFO >= INFO_thr,]
  # Order the variants by position
  sisec <- isec[order(isec$POS),]
  
  text_size=8
  legend_size=2
  # Plot INFO-value distributions
  p1 = ggplot(temp1, aes(x=as.numeric(temp1$INFO), color=AF_GROUP, fill=AF_GROUP)) + 
    geom_histogram(aes(y=..density..), alpha=0.4,position="identity",bins = 100) + 
    geom_density(alpha=.2) +
    labs(title=paste0("Distribution of R2: chr ",chr),x="R2", y = "Density") +
    guides(shape = guide_legend(override.aes = list(size = legend_size)),
           color = guide_legend(override.aes = list(size = legend_size))) +
    scale_color_viridis(discrete=TRUE) +
    scale_fill_viridis(discrete=TRUE) +
    theme(legend.title = element_text(size = text_size), 
          legend.text  = element_text(size = text_size),
          legend.key.size = unit(0.3, "lines"),legend.position = c(0.5, 0.5))
  
  # Plot AF of panel vs. imputed variants
  p2 = ggplot(isec,aes(x=as.numeric(AF_ref),y=as.numeric(AF_target))) + stat_binhex(bins = 300) +
    labs(title=paste0("Imputed AF vs. reference panel AF: chr ",chr),x="AF in ref data", y = "AF in target data") +
    guides(shape = guide_legend(override.aes = list(size = legend_size)),
           color = guide_legend(override.aes = list(size = legend_size))) +
    scale_fill_viridis() +
    theme(legend.title = element_text(size = text_size), 
          legend.text  = element_text(size = text_size),
          legend.key.size = unit(0.5, "lines"),legend.position = c(0.85, 0.3))
  
  # Imputed data AF histogram for intersecting variants
  p3 = ggplot(isec, aes(x=as.numeric(isec$AF_target))) + 
    geom_histogram(color="black", fill="white",bins=100)+
    labs(title=paste0("AF distribution of imputed variants: chr ",chr),x="Imputed AF", y = "Count")
 
  p4 = ggplot(sisec,aes(x=as.numeric(POS),y=as.numeric(sisec$AF_target)-as.numeric(sisec$AF_ref))) + stat_binhex(bins = 300) +
    labs(title=paste0("Absolute AF differences along the chromosome"),x="Chromosome position", y = "AF difference (imputed - panel)") +
    guides(shape = guide_legend(override.aes = list(size = legend_size)),
           color = guide_legend(override.aes = list(size = legend_size))) +
    scale_fill_viridis() +
    theme(legend.title = element_text(size = text_size), 
          legend.text  = element_text(size = text_size),
          legend.key.size = unit(0.5, "lines"),legend.position = c(0.85, 0.8))
  
  assign(paste0("chr",chr,"_plot"),grid.arrange(p1, p2, p3, p4, nrow = 2))
  
  #rand2 <- sample(1:nrow(isec), 500000, replace = FALSE)
  list(chr_plots=get(paste0("chr",chr,"_plot")),SNPs_data=isec) #isec[rand2,]
}  

stopCluster(cl)
merge_chr_plots=foreach(fold.result=par_res, fold.num=icount()) %do%{
  fold.result$chr_plots
}
merge_SNPs_data=foreach(fold.result=par_res, fold.num=icount(), .combine = rbind) %do%{
  fold.result$SNPs_data
}

pdf(paste0("AF_afterImputation",pop,"_R2thr_",INFO_thr,"_v2.pdf"))
# grid.arrange(grobs=chr_plots,nrow=11)
for(chr in 1:length(merge_chr_plots)){
  print(as.ggplot(merge_chr_plots[[chr]]))
}
dev.off()
fwrite(merge_SNPs_data,paste0(unlist(strsplit(pop,"_"))[2],"_allINFO",INFO_thr,"_v2.txt"),row.names=F,quote=F,col.names =T)


text_size=8
legend_size=2
# rand1 <- sample(1:nrow(merge_SNPs_data), 500000, replace = FALSE)
# temp1 <- merge_SNPs_data[rand1,]
temp1 <- merge_SNPs_data
temp1$AF_GROUP <- factor(temp1$AF_GROUP, levels = c("1", "2","3"), labels = c("MAF > 5%", "MAF 0.5-5%", "MAF < 0.5%"))

# Plot INFO-value distributions
p1 = ggplot(temp1, aes(x=INFO, color=AF_GROUP, fill=AF_GROUP)) + 
  geom_histogram(aes(y=..density..), alpha=0.4,position="identity",bins = 100) + 
  geom_density(alpha=.2) +
  labs(title=paste0("Distribution of R2"),x="R2", y = "Density") +
  guides(shape = guide_legend(override.aes = list(size = legend_size)),
        color = guide_legend(override.aes = list(size = legend_size))) +
  scale_color_viridis(discrete=TRUE) +
  scale_fill_viridis(discrete=TRUE) +
  theme(legend.title = element_text(size = text_size), 
        legend.text  = element_text(size = text_size),
        legend.key.size = unit(0.3, "lines"),legend.position = c(0.5, 0.5))

# Plot AF of panel vs. imputed variants
p2 = ggplot(merge_SNPs_data,aes(x=AF_ref,y=AF_target)) + stat_binhex(bins = 500) +
  labs(title=paste0("Imputed AF vs. reference panel AF"),x="AF in ref data", y = "AF in target data") +
  guides(shape = guide_legend(override.aes = list(size = legend_size)),
        color = guide_legend(override.aes = list(size = legend_size))) +
  scale_fill_viridis() +
  theme(legend.title = element_text(size = text_size), 
        legend.text  = element_text(size = text_size),
        legend.key.size = unit(0.5, "lines"),legend.position = c(0.85, 0.3))

# Imputed data AF histogram for intersecting variants
p3 = ggplot(merge_SNPs_data, aes(x=AF_target)) + 
  geom_histogram(color="black", fill="white",bins=300)+
  labs(title=paste0("AF distribution of imputed variants"),x="Imputed AF", y = "Count")

merge_SNPs_data=merge_SNPs_data[order(merge_SNPs_data$CHR,merge_SNPs_data$POS),]
x_coords=sapply(13:22, function(chr){
  sub_chr=merge_SNPs_data[which(merge_SNPs_data$CHR==chr),]
  pos=sub_chr$POS[nrow(sub_chr)]
  return(which(merge_SNPs_data$POS==pos & merge_SNPs_data$CHR==chr))
})
df=data.frame(x=x_coords,y=rep(-0.55,length(x_coords)),xend=x_coords,yend=rep(0.55,length(x_coords)))
p4 = ggplot(merge_SNPs_data,aes(x=1:nrow(merge_SNPs_data),y=AF_target-AF_ref)) + stat_binhex(bins = 500) +
  labs(title=paste0("Absolute AF differences along the chromosome"),x="Chromosome position", y = "AF difference (imputed - panel)") +
  guides(shape = guide_legend(override.aes = list(size = legend_size)),
        color = guide_legend(override.aes = list(size = legend_size))) +
  scale_fill_viridis() +
  theme(legend.title = element_text(size = text_size), 
        legend.text  = element_text(size = text_size),
        legend.key.size = unit(0.5, "lines"),legend.position = c(0.85, 0.8)) + 
  geom_segment(aes(x = x, y = y, xend = xend , yend = yend),data=df) 


pdf(paste0("AF_afterImputation",pop,"_R2thr_",INFO_thr,"_allCHR.pdf"))
grid.arrange(p1, p2, p3, p4, nrow = 2)
dev.off()

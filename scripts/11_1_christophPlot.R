# Pathfinder
libs <- c("tidyverse", "DESeq2") # list libraries here
lapply(libs, require, character.only = T)
rm(libs)
#These are the outputs from DESEQ, keeping everything whether significant or not
raw.res <- readRDS("data/processed/models/unadjusted_de_interaction.RDS")
res1 <- results(raw.res, name = "genotype_KO_vs_WT") |> as.data.frame()
res2 <- results(raw.res, name = "treatment_TAC_vs_Sham") |> as.data.frame()

res1$X <- row.names(res1)
res2$X <- row.names(res2)
#match gene names
res2=res2[match(res1$X,res2$X),]
tokeep=!is.na(res2$X)
res1=res1[tokeep,]
res2=res2[tokeep,]

#get out fold changes
res1_fc=res1$log2FoldChange
res2_fc=res2$log2FoldChange

#how many bins do you want? (so 20 bins every 5% of the DE genes)
desired_bins=20

partition_size=floor(length(res1_fc)/(desired_bins))
seq=c(1:(desired_bins) *partition_size)

#This give the boundaries for the partitions
ranges1=c(min(res1_fc),sort(res1_fc)[seq],max(res1_fc))
ranges2=c(min(res2_fc),sort(res2_fc)[seq],max(res2_fc))

#quick function to assign a specific FC to a specific bin
get_bin <- function(x,ranges){
  return(min(which(x<ranges))-1)
}

#apply it, obviously
bin1=sapply(res1_fc,get_bin,ranges1)
bin2=sapply(res2_fc,get_bin,ranges2)

#I forget why I had to do this... something about the function I wrote
bin1[bin1=="Inf"]=desired_bins
bin2[bin2=="Inf"]=desired_bins

bin1[bin1==(desired_bins+1)]=desired_bins
bin2[bin2==(desired_bins+1)]=desired_bins


#make a dataframe.
positions=data.frame(symbol=res1$X, bin1=bin1,bin2=bin2)

#put everything into a NxN matrix, where N is the number of bins you want
big_matrix=matrix(0,nrow=desired_bins,ncol=desired_bins)
for(i in 1:nrow(positions)){
  big_matrix[positions$bin1[i],positions$bin2[i]]=big_matrix[positions$bin1[i],positions$bin2[i]]+1
}
positions2=melt(big_matrix) 

#This is figuring out exactly how many genes are in each bin so we can color it
temp=paste(bin1,bin2,sep="_")
temp2=paste(positions2$Var1,positions2$Var2,sep="_")
z=match(temp,temp2)
number=positions2[z,3]

positions$number=number
#This one is better
ggplot(positions,aes(x=bin1,y=bin2,color=number)) + geom_jitter() +scale_color_gradient(low="blue",high="red")+
  xlab("<--  Down    KO_vs_WT      Up  -->") + ylab("<--  Down      TAC_vs_Sham      Up -->")+theme_bw()+
  theme(axis.text = element_blank(),text=element_text(size=16,family="serif"))

#get an output.  The ones you care about are the ones with the top left and bottom right, so 
# Bin 1 = 1-5 & Bin 2 = 16-20  or Bin 1 = 16-20 & Bin 2 = 1-5
output = cbind(res2$x,res1_fc,res2_fc,bin1,bin2)  
colnames(output)=c("Gene","Fold Change Condition 1", "Fold Change Condition 2",
                   "Bin Condition 1","Bin Condition 2")



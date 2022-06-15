# load functions
# src_dir="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/Scripts/" # cluster
src_dir="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/"

source(paste (src_dir,"Microbial_Profiles_Functions.R",sep="")) 

sub_sample = function (df,sample_name,total_count) # df is an output of Tax_profile function, {"Value","Count","Percent","sample_name"}
{
  subset=df[as.character(df$sample_name)==as.character(sample_name),]
  subset_expanded=data.frame()
  for (i in 1:nrow(subset)){
    tmp=data.frame(Value=rep(subset$Value[i],subset$Count[i]),count=rep(1,subset$Count[i]))
    subset_expanded=rbind(subset_expanded,tmp)
  }
  subset_subsampled_tmp=subset_expanded[sample(1:nrow(subset_expanded),size = total_count,replace = F),]
  subset_subsampled=data.frame()
  sampled_tax=unique(subset_subsampled_tmp$Value)
  sampled_sum=sum(subset_subsampled_tmp$count)
  for (tax in sampled_tax) {
    tax_count=sum(subset_subsampled_tmp$count[subset_subsampled_tmp$Value==tax])
    tax_percent=tax_count/sampled_sum
    subset_subsampled=rbind(subset_subsampled,data.frame(Value=tax,Count=tax_count,Percent=tax_percent,sample_name=sample_name))
  }
  subset_subsampled=subset_subsampled[order(subset_subsampled$Count,decreasing = T),]
  return(subset_subsampled)
}
# Note: This script is more oganzied and complete script of  /Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Analyze_taxa_profile_only_4_probes_array.R

#### CLEAN START ANALYZE WITH NORMALIZED DATA and SUBSAMPLES ON GENUS AND REPEATS  25/10/2021 ####
#############################################################################
# repeat subsample 100 times (only on the genus level) and plot the (1) NMDS (2) UpSetR and record the values of the intersections and subsampled profiles
## FOR THE FULL SAMPLES OMNI4 and filtering rare species from each sample -> here!
set.seed(639245)
# load functions
# src_dir="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/Scripts/" # cluster
src_dir="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/"

source(paste (src_dir,"Microbial_Profiles_Functions.R",sep="")) 

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenMatrix <- function(mat) {
  ut <- upper.tri(mat)
  data.frame(
    row = rownames(mat)[row(mat)[ut]],
    column = rownames(mat)[col(mat)[ut]],
    val=(mat)[ut]
  )
}

sub_sample = function (df,sample_name,total_count) # df is an output of Tax_profile function, {"Value","Count","Percent","sample_name"}
{
  subset=df[as.character(df$sample_name)==as.character(sample_name),]
  subset_expanded=data.frame()
  for (i in 1:nrow(subset)){
    tmp=data.frame(Value=rep(subset$Value[i],subset$Count[i]),count=rep(1,subset$Count[i]))
    subset_expanded=rbind(subset_expanded,tmp)
  }
  subset_subsampled_tmp=subset_expanded[sample(1:nrow(subset_expanded),size = total_count,replace = F),]
  subset_subsampled=data.frame()
  sampled_tax=unique(subset_subsampled_tmp$Value)
  sampled_sum=sum(subset_subsampled_tmp$count)
  for (tax in sampled_tax) {
    tax_count=sum(subset_subsampled_tmp$count[subset_subsampled_tmp$Value==tax])
    tax_percent=tax_count/sampled_sum
    subset_subsampled=rbind(subset_subsampled,data.frame(Value=tax,Count=tax_count,Percent=tax_percent,sample_name=sample_name))
  }
  subset_subsampled=subset_subsampled[order(subset_subsampled$Count,decreasing = T),]
  return(subset_subsampled)
}
# Genus levels stackplots and heatmaps
# Another option: https://cran.r-project.org/web/packages/Polychrome/vignettes/testgg.html
library(Polychrome)
library(gplots)

# Global vars and data reading
path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/multimaped_filtered/MMSEQS2/UMI_FILTERED/"
library(openxlsx)
library(Polychrome)
library(gplots)
library(ggplot2)
DB_type="MMSEQS2"
num_of_subsamples_repeats=100
out_dir=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/",DB_type,"/UMI_FILTERED/","concise_4_probesOMNI","/subsample/Rp10M/repeats_",num_of_subsamples_repeats,"/",sep="")
dir.create(out_dir, showWarnings = FALSE,recursive = TRUE)

# The data itself
# The full sample
list_of_files_unsplitted=list.files(path = path,
                                    pattern = ".*_Ar_multimap_bacteria_R2.usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt")
list_of_files_unsplitted=paste(path,list_of_files_unsplitted,sep="")
sample_name_unsplitted=basename(list_of_files_unsplitted)

# The splitted by prob
list_of_files_splitted=list.files(path = path,
                                  pattern = ".*_Ar_multimap_bacteria_R2.PROB_.*.usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt")
list_of_files_splitted=paste(path,list_of_files_splitted,sep="")
sample_name_splitted=basename(list_of_files_splitted)

# Derek primer 515
Derek_AmpSeq_515_path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/AmpSeq/515_806/"
list_of_Derek_Amp_seq_515_files=list.files(path = Derek_AmpSeq_515_path,
                                           pattern = ".*uniq_lineage_count.txt")

list_of_Derek_Amp_seq_515_files=paste(Derek_AmpSeq_515_path,list_of_Derek_Amp_seq_515_files,sep="")
sample_name_Derek_Amp_seq_515=basename(list_of_Derek_Amp_seq_515_files)
sample_name_Derek_Amp_seq_515=paste(sample_name_Derek_Amp_seq_515,".515_806",sep="")
# Derek primer 799
Derek_AmpSeq_799_path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/AmpSeq/799_1192/"
list_of_Derek_Amp_seq_799_files=list.files(path = Derek_AmpSeq_799_path,
                                           pattern = ".*uniq_lineage_count.txt")

list_of_Derek_Amp_seq_799_files=paste(Derek_AmpSeq_799_path,list_of_Derek_Amp_seq_799_files,sep="")
sample_name_Derek_Amp_seq_799=basename(list_of_Derek_Amp_seq_799_files)
sample_name_Derek_Amp_seq_799=paste(sample_name_Derek_Amp_seq_799,".799_1192",sep="")

list_of_files=c(list_of_files_unsplitted,list_of_files_splitted,list_of_Derek_Amp_seq_799_files,list_of_Derek_Amp_seq_515_files)
sample_names=c(sample_name_unsplitted,sample_name_splitted,sample_name_Derek_Amp_seq_799,sample_name_Derek_Amp_seq_515)

sample_names=gsub(".usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt",replacement = "",fixed = T,x = sample_names)
sample_names=gsub(".usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.uniq_lineage_count.txt",replacement = "",fixed = T,x = sample_names)

sample_names=gsub("_Ar_multimap_bacteria_R2",replacement = "",fixed = T,x = sample_names)
sample_names=gsub("191029-069_",replacement = "",fixed = T,x = sample_names)

tmp_df=data.frame(file=list_of_files,sample_name=sample_names)
lineages_files_and_tag=tmp_df

# change the names of the AmpSeq samples
lineages_files_and_tag$sample_name=gsub("F5p515F2F4F6",replacement = "AmpSeq_A",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("C5p515F2F4F6",replacement = "AmpSeq_B",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("A2p515F2F4F6",replacement = "AmpSeq_C",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("E3p515F1F3F5",replacement = "AmpSeq_D",fixed = T,x = lineages_files_and_tag$sample_name)

# remove the bad splitted ones
lineages_files_and_tag=lineages_files_and_tag[!grepl("1.PROB_479",lineages_files_and_tag$sample_name),]
lineages_files_and_tag=lineages_files_and_tag[!grepl("1.PROB_1265",lineages_files_and_tag$sample_name),]

lineages_files_and_tag_x2=lineages_files_and_tag[!grepl(lineages_files_and_tag$sample_name,pattern = "([ABCD]1)",perl=T),]
# Should we exclude another one --> YES
# AmpSeq_C.799_1192 603 reads
# lineages_files_and_tag_x2=lineages_files_and_tag_x2[!grepl(lineages_files_and_tag_x2$sample_name,pattern = "AmpSeq_C.799_1192",fixed = T),]

# plot definitions
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
P40 <- sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
swatch(P40)
names(P40) <- NULL

col1<- colorRampPalette(c("white", "blue"))(50)
breaks1 = c(seq(0,0.5,length=51))

# out files
XLSX_matrix_file=paste(out_dir,"datasets_tables_Bacteria_AmpSeq_vs_Array_4probes_subsampled_equal_size_",DB_type,".UMI_FILTERED.xlsx",sep="")
XSLX_Obj=createWorkbook(XLSX_matrix_file)
pdf(file=paste(out_dir,"Array_4probes_Bacteria_all_samples_RpM_subsampled_equal_size_MMSEQS2_NT_stacks_labeled_NMDS_and_intersects.UMI_FILTERED.pdf",sep=""),height = 16.5,width = 23.4) # paper = "a4r", height = 8.3 , width = 11.7

Level="genus" # for simplicity we only analyze genus level, we can do a loop for all... for (Level in c ("genus")) # c ("phylum", "class","order","family","genus") #=c("Count", "superkingdom","phylum", "class","order","family","genus","species")

# read all data and decide the size to sample
all=Tax_profile(lineage_files_and_tag_df=lineages_files_and_tag_x2,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
                levels_to_account_for_exclude = c("genus"))

# we don't use the prob samples so lets remove them
all=all[!grepl(all$sample_name,pattern = "PROB"),]

# subset
# count the total number of reads in each sample
sample_count=data.frame()
for (sample_i in unique(all$sample_name)) {
  sample_total=sum(all$Count[all$sample_name==sample_i])
  sample_count=rbind(sample_count,data.frame(sample_name=sample_i,total_size=sample_total))
}

subset_size=5000 #min(sample_count$total_size)

# RpM
all_RpM=data.frame()
for (sample in sample_count$sample_name)
{
  sample_depth=sample_count$total_size[sample_count$sample_name==sample]
  rpm_sample=all[all$sample_name==sample,]
  rpm_sample$Rp10M=round(rpm_sample$Count/sample_depth*10^6,0)
  rpm_sample=data.frame(Value=rpm_sample$Value,Count=rpm_sample$Rp10M,Percent=rpm_sample$Rp10M/sum(rpm_sample$Rp10M),sample_name=rpm_sample$sample_name)
  all_RpM=rbind(all_RpM,rpm_sample)
}

sample_count_rpm=data.frame()
for (sample_i in unique(all_RpM$sample_name)) {
  sample_total=sum(all_RpM$Count[all_RpM$sample_name==sample_i])
  sample_count_rpm=rbind(sample_count_rpm,data.frame(sample_name=sample_i,total_size=sample_total))
}

intersects_per_subset=data.frame()
Bray_Curtis=data.frame()
# To do: aggregate results...
# add the subsampled profile to xls
# df for intersect size
# 
for (i_repeat in 1:num_of_subsamples_repeats) # the i_repeat should be accounted later in figures and aggregator
{
  message (paste("- Sstart subsample repeat",i_repeat))
  all_subsample=data.frame()
  for (sample_i in unique(all_RpM$sample_name)) {
    sub_sample_df=sub_sample(df = all_RpM,sample_name = sample_i, subset_size)
    message(paste ("\t-- subsample ",subset_size," reads out of ",sample_count_rpm$total_size[sample_count_rpm$sample_name==sample_i]," reads of sample ",sample_i))
    all_subsample=rbind(all_subsample,sub_sample_df)
  }
  # save subsample
  write.table(file = paste(out_dir,"subsample_",i_repeat,"_",Level,".UMI_FILTERED.profile.csv",sep=""),x =all_subsample, quote = F,row.names = F,sep = ";",col.names = T)
  
  ### HERE! --> Continue to update and replace all df with all_subsample
  all_with_others=data.frame()
  other_cutoff=0.025
  for (sample in unique (all_subsample$sample_name))
  {
    subset=all_subsample[all_subsample$sample_name==sample,]
    sabset_to_add=subset[subset$Percent>other_cutoff,]
    other_percent=sum(subset$Percent[subset$Percent<=other_cutoff])
    other_count=sum(subset$Count[subset$Percent<=other_cutoff])
    sabset_to_add=rbind(sabset_to_add,data.frame(Value=paste("Other (<=",other_cutoff*100,"%)",sep=""),Count=other_count,Percent=other_percent,sample_name=sample))
    all_with_others=rbind(all_with_others,sabset_to_add)
  }
  # no need for the PROB samples here --> before the change the data frame used was: all_with_others
  all_with_others_NO_PROB=all_with_others[!grepl(all_with_others$sample_name,pattern = "PROB"),]
  all_with_others_NO_PROB$sample_name <- factor(all_with_others_NO_PROB$sample_name,levels = c("A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                                               "B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                                               "C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                                               "D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  p=ggplot(all_with_others_NO_PROB, aes(fill=Value, y=Percent, x=sample_name,label=paste(signif(Percent*100,2),"% ",Value,sep=""))) + 
    geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
    theme(legend.text=element_text(size=15)) +
    scale_fill_manual(values = P40) + 
    ggtitle(paste("Subsample",i_repeat,"composition",Level)) +
    theme(plot.title = element_text(hjust = 0.5))
  print (p)
  
  # scale_fill_manual(values = as.vector(palette36))
  # scale_fill_brewer( palette = "Paired")
  
  # the heatmap - is it needed?!?
  # percent_matrix=make_table(flat_df = all_with_others_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  # count_matrix=make_table(flat_df = all_with_others_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  ## TODO: take into account that it may be that some fall under the other category...
  ## we can fix it with the same proceudre we did for missing data (we need to replace the NA with number and reduce it from the Other category)
  # percent_matrix_NA_0=percent_matrix
  # percent_matrix_NA_0[is.na(percent_matrix_NA_0)]=0
  # heatmap.2(x = as.matrix(percent_matrix_NA_0), col = col1, symm = FALSE,breaks=breaks1,trace = "none",margins=c(15,15),cexRow=0.9, main=paste(Level," - NA as 0, DB: ",DB_type, sep="")) # main = paste("Abundance of ",tax_level, " (out of the all reads)",sep="")
  
  # sheetName = paste(Level,"_wOt_",other_cutoff,"_percent_NA_0",sep="")
  # addWorksheet(XSLX_Obj, sheetName)
  # writeDataTable(x=percent_matrix_NA_0, wb = XSLX_Obj,sheet = sheetName,rowNames = TRUE)
  
  # Do NMDS between the full samples profiles - https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/ / https://jkzorz.github.io/2019/06/06/NMDS.html
  ## CONSIDER TO FILTER OUT THE RARE SPECIES.... TO AVOID DOUBLE ZERO AND OTHER PROBLEMS?
  
  # From each sample filter out species below a cutoff
  # table(all[all$Count<2,"Value"])
  # table(all$sample_name)
  
  ## length(subset$Percent[subset$Percent<=0.0005]) sample="AmpSeq_A.515_806"
  library(vegan)
  all_NO_PROB=all_subsample[!grepl(all_subsample$sample_name,pattern = "PROB"),]
  all_NO_PROB$sample_name <- factor(all_NO_PROB$sample_name,levels = c("A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                       "B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                       "C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                       "D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  for (sample in unique(all_NO_PROB$sample_name))
  {
    s_sum=sum(all_NO_PROB$Count[all_NO_PROB$sample_name==sample])
    message (paste ("QA: ",sample," total ",s_sum,sep=""))
  }
  
  percent_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  count_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  
  # percent_matrix_all_sub_samples=make_table(flat_df = all_NO_PROB_sub_samples,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  # count_matrix_all_sub_samples=make_table(flat_df = all_NO_PROB_sub_samples,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  
  community_matrix=t(as.matrix(percent_matrix_all_samples))
  #community_matrix=t(as.matrix(percent_matrix_all_sub_samples))
  community_matrix[is.na(community_matrix)]=0 # replace NA with 0
  example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                       k=2,trymax = 100,distance = "bray") # The number of reduced dimensions -> ,, autotransform = FALSE,distance ="euclidean"
  stressplot(example_NMDS)
  # plot (example_NMDS)
  
  #extract NMDS scores (x and y coordinates)
  data.scores = as.data.frame(scores(example_NMDS))
  data.scores$Sample=row.names(data.scores)
  data.scores$Sample=gsub(pattern=".799_1192",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".515_806",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="AmpSeq_",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="1|2",replacement = "",x=data.scores$Sample)
  
  data.scores$Type=row.names(data.scores) 
  data.scores$Type[grepl("^[ABCD]1$",perl=T,rownames(data.scores))]="Array2Probes"
  data.scores$Type[grepl("^[ABCD]2$",perl=T,rownames(data.scores))]="Array4Probes"
  data.scores$Type[grepl("515_806",perl=T,rownames(data.scores))]="AmpSeq_p515"
  data.scores$Type[grepl("799_1192",perl=T,rownames(data.scores))]="AmpSeq_p799"
  
  nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
    geom_point(size = 6, aes( shape = Sample, fill = Type,color=Type))+ 
    # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
    geom_text(aes(label=""),hjust=-0.1, vjust=-0.1) +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    # labs(x = "NMDS1", colour = "Plant", y = "NMDS2", shape = "Type")  + 
    # labs(x = "NMDS1", colour = "Method", y = "NMDS2", shape = "Plant")  + 
    labs(x = "NMDS1", fill = "Method", y = "NMDS2", shape = "Plant",colour = "Method")  + 
    scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    scale_shape_manual(values=c(21, 22, 23,24))+ 
    scale_fill_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    ggtitle(paste("Subsample ",i_repeat," - NMDS of the full profile\nsubsampled to equal size (",subset_size,") on the ",Level, " level (Stress=",signif(example_NMDS$stress,3),")",sep=""))
  print (nmds_plot,width = 4, height = 4)
  # dev.off()
  # colors for the methods and shape for plants ?
  
  # hirarchial clustering and distances heatmaps
  # calculate Bray-Curtis distance among samples
  comm.bc.dist <- vegdist(community_matrix, method = "bray")
  # cluster communities using average-linkage algorithm
  comm.bc.clust <- hclust(comm.bc.dist, method = "average")
  # plot cluster diagram
  plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")
  library(circlize)
  col_fun = colorRamp2(c(0, 1), c("red","white"))
  col_fun(seq(0, 1))
  
  library(ComplexHeatmap)
  dist_mat=as.matrix(comm.bc.dist)
  pHeatmapSim=Heatmap(dist_mat,col=col_fun,column_title = paste("Subsample ",i_repeat,"- Bray-Curtis dissimilarity",sep=""),name="Distance",
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%.2f", dist_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                      })
  print(pHeatmapSim) 
  dist_mat_flat=flattenMatrix(dist_mat)
  names(dist_mat_flat)=c("sample1","sample2",paste("BC_Diss","_sample_",i_repeat,sep=""))
  if (nrow(Bray_Curtis)==0)
  {
    Bray_Curtis=dist_mat_flat
  } else {
    Bray_Curtis=merge(Bray_Curtis,dist_mat_flat,by=c("sample1","sample2"),all=T)
  }
  # # clustering with k mean -> consider this: https://rstudio-pubs-static.s3.amazonaws.com/274963_9e95868f44d64a008e4bf197cdef079e.html
  # find_k <- function(data, plot_title = ""){
  #   max_clusters=dim(as.matrix(data))[1]
  #   #create data for within group sum of square
  #   wss = numeric(max_clusters)
  #   for (k in 1:(max_clusters-1)) wss[k] <- kmeans(data, k, nstart = 20)$tot.withinss
  #   
  #   #create data for average silhouette distance
  #   library(cluster)
  #   asw <- numeric(max_clusters)
  #   for (k in 2:(max_clusters-1)) asw[k] <- pam(data, k)$silinfo$avg.width
  #   
  #   res=data.frame(k=1:max_clusters,WSS=wss,asw=asw)
  #   #create s cree plot
  #   par(mar=c(5, 4, 4, 6))
  #   plot(1:max_clusters, wss, type = "b", main = plot_title, xlab = "Number of Clusters", ylab = "Within groups sum of squares")
  #   par(new = T)
  #   plot(1:max_clusters, asw, type = "l", lty = 2, col = "red", axes = F, xlab = NA, ylab = NA)
  #   axis(side = 4)
  #   mtext("Average Silhouette width", side = 4, line = 3)
  #   legend("topright", legend = c("WSS", "Si"), lty = c(1,2), col = c("black", "red"))
  #   return(res)
  # }
  # res=find_k(comm.bc.dist,"bray")
  # #conduct k-means to find clusters
  # km=pam(x = comm.bc.dist,k = 2)
  # plot(silhouette(km, comm.bc.dist), main = paste("Silhoulette plot using Bray-Curtis distance") , col = 1:2, border = NA, cex = 0.6)
  # 
  
  # Do the upSetR plots
  unique(all_NO_PROB$sample_name)
  tax_list=list(A2=all_subsample$Value[all_subsample$sample_name=="A2"],
                AmpSeq_A_799=all_subsample$Value[all_subsample$sample_name=="AmpSeq_A.799_1192"],
                AmpSeq_A_515=all_subsample$Value[all_subsample$sample_name=="AmpSeq_A.515_806"],
                
                B2=all_subsample$Value[all_subsample$sample_name=="B2"],
                AmpSeq_B_799=all_subsample$Value[all_subsample$sample_name=="AmpSeq_B.799_1192"],
                AmpSeq_B_515=all_subsample$Value[all_subsample$sample_name=="AmpSeq_B.515_806"],
                
                C2=all_subsample$Value[all_subsample$sample_name=="C2"],
                AmpSeq_C_799=all_subsample$Value[all_subsample$sample_name=="AmpSeq_C.799_1192"],
                AmpSeq_C_515=all_subsample$Value[all_subsample$sample_name=="AmpSeq_C.515_806"],
                
                D2=all_subsample$Value[all_subsample$sample_name=="D2"],
                AmpSeq_D_799=all_subsample$Value[all_subsample$sample_name=="AmpSeq_D.799_1192"],
                AmpSeq_D_515=all_subsample$Value[all_subsample$sample_name=="AmpSeq_D.515_806"]
  )
  
  m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
  cs = comb_size(m)
  ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),column_title = paste("A2 (Omni4Probs)",Level))
  ht = draw(ht)
  co = column_order(ht)
  nc = ncol(m)
  decorate_annotation("Intersection\nsize", {
    grid.text(cs[co], 
              x = 1:nc, 
              y = unit(cs[co], "native") + unit(1, "mm"), 
              gp = gpar(fontsize = 8), 
              just = "bottom",
              default.units = "native")
  })
  # separately for each plant 
  # pdf(file=paste(out_dir,"AmpSeq_Array4probes_Bacteria_all_samples_MMSEQS2_NT_upSetR.pdf",sep=""),height = 3,width = 6.5) 
  i_subset_intersects=data.frame()
  for (p in c("A","B","C","D"))
  {
    tax_list=list(Array=all_subsample$Value[all_subsample$sample_name==paste(p,"2",sep="")],
                  AmpSeq_p799=all_subsample$Value[all_subsample$sample_name==paste("AmpSeq_",p,".799_1192",sep="")],
                  AmpSeq_p515=all_subsample$Value[all_subsample$sample_name==paste("AmpSeq_",p,".515_806",sep="")])
    # for the sub samples
    # tax_list=list(Array=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste(p,"2",sep="")]),
    #               AmpSeq_799=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste("AmpSeq_",p,".799_1192",sep="")]),
    #               AmpSeq_515=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste("AmpSeq_",p,".515_806",sep="")]))
    m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
    cs = comb_size(m)
    row_size = set_size(m)
    ht = UpSet(m, 
               top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),
               right_annotation = upset_right_annotation(m, ylim = c(0, 1.1*max(set_size(m)))),
               column_title = paste("Subsample",i_repeat,"- Plant",p,Level,"subsampled to",subset_size,"reads"))
    ht = draw(ht)
    co = column_order(ht)
    row_od = row_order(ht)
    nc = ncol(m)
    decorate_annotation("Intersection\nsize", {
      grid.text(cs[co], 
                x = 1:nc, 
                y = unit(cs[co], "native") + unit(1, "mm"), 
                gp = gpar(fontsize = 10), 
                just = "bottom",
                default.units = "native")
    })
    decorate_annotation("Set size", {
      grid.text(row_size[row_od], 
                unit(row_size[row_od], "native") + unit(1, "mm"), 
                rev(seq_len(length(row_size))), 
                default.units = "native", just = "bottom", rot = -90,
                gp = gpar(fontsize = 10))
    })
    # library(eulerr)
    
    
    intersects=c("Array" = 0, "AmpSeq_p799" = 0, "AmpSeq_p515" = 0,
                 "Array&AmpSeq_p799&AmpSeq_p515"=0,  
                 "Array&AmpSeq_p799"=0,
                 "Array&AmpSeq_p515"=0,
                 "AmpSeq_p799&AmpSeq_p515"=0)
    for (i in names(cs)) { # assign
      if (i=="100") {intersects[["Array"]]=cs[[i]]}
      if (i=="010") {intersects[["AmpSeq_p799"]]=cs[[i]]}
      if (i=="001") {intersects[["AmpSeq_p515"]]=cs[[i]]}
      if (i=="111") {intersects[["Array&AmpSeq_p799&AmpSeq_p515"]]=cs[[i]]}
      if (i=="110") {intersects[["Array&AmpSeq_p799"]]=cs[[i]]}
      if (i=="101") {intersects[["Array&AmpSeq_p515"]]=cs[[i]]}
      if (i=="011") {intersects[["AmpSeq_p799&AmpSeq_p515"]]=cs[[i]]}
    }
    intersect_df=as.data.frame(intersects)
    row.names(intersect_df)=paste(p,row.names(intersect_df),sep="_")
    names(intersect_df)=paste("subset",i_repeat,"intersect",sep="_")
    
    i_subset_intersects=rbind(i_subset_intersects,intersect_df)
    names(i_subset_intersects)=paste("subset",i_repeat,"intersect",sep="_")
    # fit <- euler(intersects)
    # # Customize colors, remove borders, bump alpha, color labels white
    # plot(fit,
    #      fills = list(fill = c("#E69F00", "#56B4E9","#009E73","#CC79A7"), alpha = 0.5),
    #      labels = list(col = "white", font = 4))
  }
  if (nrow(intersects_per_subset)==0) {
    intersects_per_subset=i_subset_intersects
  } else {
    intersects_per_subset=merge(intersects_per_subset,i_subset_intersects,by="row.names",all=T)
    row.names(intersects_per_subset)=intersects_per_subset$Row.names
    intersects_per_subset=intersects_per_subset[,-1]
  }
  # dev.off()
  
  ## plot the subsamples profiles
  # all_subsamples_with_others=data.frame()
  # other_cutoff=0.02
  # for (sample in unique (all_NO_PROB_sub_samples$sample_name))
  # {
  #  subset=all_NO_PROB_sub_samples[all_NO_PROB_sub_samples$sample_name==sample,]
  #  names(subset)=c("Value","Count","Percent","sample_name")
  #  sabset_to_add=subset[subset$Percent>other_cutoff,]
  #  other_percent=sum(subset$Percent[subset$Percent<=other_cutoff])
  #  other_count=sum(subset$Count[subset$Percent<=other_cutoff])
  #  sabset_to_add=rbind(sabset_to_add,data.frame(Value=paste("Other (<=",other_cutoff*100,"%)",sep=""),Count=other_count,Percent=other_percent,sample_name=sample))
  #  all_subsamples_with_others=rbind(all_subsamples_with_others,sabset_to_add)
  # }
  ## no need for the PROB samples here --> before the change the data frame used was: all_with_others
  # all_subsamples_with_others_NO_PROB=all_subsamples_with_others[!grepl(all_with_others$sample_name,pattern = "PROB"),]
  # all_subsamples_with_others_NO_PROB$sample_name <- factor(all_subsamples_with_others_NO_PROB$sample_name,levels = c("A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
  #                                                                                              "B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
  #                                                                                              "C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
  #                                                                                              "D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  # p=ggplot(all_subsamples_with_others_NO_PROB, aes(fill=Value, y=Percent, x=sample_name,label=paste(signif(Percent*100,2),"% ",Value,sep=""))) + 
  #   geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
  #   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
  #   theme(legend.text=element_text(size=15)) +
  #   scale_fill_manual(values = P40) + 
  #   ggtitle(paste("Composition",Level)) +
  #   theme(plot.title = element_text(hjust = 0.5))
  # print (p)
  
  # TO SEE IF FURTHER UPDATE IS NEEDED....
  # # plot NMDS for each of the OMNI4 OMNI2 and its probes
  # for (sample in c("Ax","Bx","Cx","Dx"))
  # {
  #   samples_regex=NA
  #   if (sample=="Ax") {samples_regex="F5p515F2F4F6|A1|A2|AmpSeq_A"}
  #   if (sample=="Bx") {samples_regex="C5p515F2F4F6|B1|B2|AmpSeq_B"}
  #   if (sample=="Cx") {samples_regex="A2p515F2F4F6|C1|C2|AmpSeq_C"}
  #   if (sample=="Dx") {samples_regex="E3p515F1F3F5|D1|D2|AmpSeq_D"}
  #   sample_subset=lineages_files_and_tag[grepl(samples_regex,lineages_files_and_tag$sample_name),]
  #   if (sample=="Ax") 
  #   {
  #     sample_subset=sample_subset[!grepl("A2p515F2F4F6",sample_subset$sample_name),]
  #   }
  #   
  #   all_sample_subset=Tax_profile(lineage_files_and_tag_df=sample_subset,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
  #                                 levels_to_account_for_exclude = c("genus"))
  #   # subset
  #   # count the total number of reads in each sample
  #   sample_count=data.frame()
  #   for (sample_i in unique(all_sample_subset$sample_name)) {
  #     sample_total=sum(all_sample_subset$Count[all_sample_subset$sample_name==sample_i])
  #     sample_count=rbind(sample_count,data.frame(sample_name=sample_i,total_size=sample_total))
  #   }
  #   all_subsample=data.frame()
  #   subset_size=min(sample_count$total_size)
  #   for (sample_i in unique(all_sample_subset$sample_name)) {
  #     sub_sample_df=sub_sample(df = all_sample_subset,sample_name = sample_i, subset_size)
  #     message(paste ("-- subsample ",subset_size," reads out of ",sample_count$total_size[sample_count$sample_name==sample_i]," reads of sample ",sample_i))
  #     all_subsample=rbind(all_subsample,sub_sample_df)
  #   }
  #   ### HERE! --> Continue to update and replace all df with all_subsample
  #   percent_matrix_all_samples=make_table(flat_df = all_subsample,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  #   count_matrix_all_samples=make_table(flat_df = all_subsample,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  #   
  #   community_matrix=t(as.matrix(percent_matrix_all_samples))
  #   #community_matrix=t(as.matrix(percent_matrix_all_sub_samples))
  #   community_matrix[is.na(community_matrix)]=0 # replace NA with 0
  #   example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
  #                        k=2,trymax = 100,distance = "bray") # The number of reduced dimensions -> ,, autotransform = FALSE,distance ="euclidean"
  #   stressplot(example_NMDS)
  #   # plot (example_NMDS)
  #   
  #   #extract NMDS scores (x and y coordinates)
  #   data.scores = as.data.frame(scores(example_NMDS))
  #   data.scores$Sample=row.names(data.scores) # AmpSeq_p515; AmpSeq_p799; pr_1265; pr_unk; pr_479; pr_902; pr_799 -> these will have colors
  #   data.scores$Sample=gsub(pattern="PROB_799F",replacement = "pr_799",x=data.scores$Sample)
  #   data.scores$Sample=gsub(pattern="PROB_902R",replacement = "pr_902",x=data.scores$Sample)
  #   data.scores$Sample=gsub(pattern="PROB_479",replacement = "pr_479",x=data.scores$Sample)
  #   data.scores$Sample=gsub(pattern="PROB_1265",replacement = "pr_1265",x=data.scores$Sample)
  #   data.scores$Sample=gsub(pattern="PROB_UNKNOWN",replacement = "pr_unk",x=data.scores$Sample)
  #   data.scores$Sample=gsub(pattern="^[ABCD][12]$",replacement = "Full",x=data.scores$Sample,perl = T)
  #   data.scores$Sample=gsub(pattern="^[ABCD][12].",replacement = "",x=data.scores$Sample,perl = T)
  #   data.scores$Sample=gsub(pattern="515_806",replacement = "AmpSeq_p515",x=data.scores$Sample)
  #   data.scores$Sample=gsub(pattern="799_1192",replacement = "AmpSeq_p799",x=data.scores$Sample)
  #   data.scores$Sample=gsub(pattern="^AmpSeq_[ABCD].",replacement = "",x=data.scores$Sample,perl=T)
  #   
  #   # type -> will be symbol {AmpSeq | Omni2Probes | Omni4Probes}
  #   data.scores$Type=row.names(data.scores) 
  #   data.scores$Type[grepl("^[ABCD]1",perl=T,rownames(data.scores))]="Omni2Probes"
  #   data.scores$Type[grepl("^[ABCD]2",perl=T,rownames(data.scores))]="Omni4Probes"
  #   data.scores$Type[grepl("AmpSeq",perl=T,rownames(data.scores))]="AmpSeq"
  #   
  #   pdf(file=paste(out_dir,"AmpSeq_Array4probes_Bacteria_all_samples_MMSEQS2_NT_NMDS.",sample,".pdf",sep=""),height = 7,width = 8)
  #   nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  #     # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
  #     geom_point(size = 6, aes( shape = Type, fill = Sample,color=Sample))+ 
  #     # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
  #     geom_text(aes(label=Sample),hjust=-0.3, vjust=-0.3) +
  #     theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
  #           axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
  #           legend.text = element_text(size = 12, face ="bold", colour ="black"), 
  #           legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
  #           axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
  #           legend.title = element_text(size = 14, colour = "black", face = "bold"), 
  #           panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
  #           legend.key=element_blank()) + 
  #     # labs(x = "NMDS1", colour = "Plant", y = "NMDS2", shape = "Type")  + 
  #     # labs(x = "NMDS1", colour = "Method", y = "NMDS2", shape = "Plant")  + 
  #     labs(x = "NMDS1", fill = "Probe/Primer", y = "NMDS2", shape = "Method",colour = "Probe/Primer")  + 
  #     scale_colour_manual(values = c("#56B4E9","#0072B2","#000000","#E69F00","#009E73","#F0E442","#D55E00","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
  #     scale_shape_manual(values=c(21,24,22))+ 
  #     scale_fill_manual(values = c("#56B4E9","#0072B2","#000000","#E69F00","#009E73","#F0E442","#D55E00","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
  #     ggtitle(paste("NMDS of the full profile subset to equal size (",subset_size,") on the ",Level, " level (Stress=",signif(example_NMDS$stress,3),")",sep=""))
  #   print (nmds_plot)
  #   dev.off()
  # }
  # 
  # # drafts
  # km_result = transpose(km)
  # 
  # km = map(comm.bc.dist, pam, k = 3)
  # 
  # # extract cluster data
  # km_result = transpose(km)
  # km_cluster = as.data.frame(km_result$cluster)
  # colnames(km_cluster) <- dist_meas
  # 
  # plot(silhouette(km[[i]], comm.bc.dist), main = paste("Silhoulette plot using Bray-Curtis distance") , col = 1:3, border = NA, cex = 0.6)
  # 
  # 
  # 
  # # hirarchial clustering and distances heatmaps
  # # calculate Euclidean distance among samples
  # comm.ec.dist <- vegdist(community_matrix, method = "euclidean")
  # 
  # # cluster communities using average-linkage algorithm
  # comm.ec.clust <- hclust(comm.ec.dist, method = "average")
  # # plot cluster diagram
  # plot(comm.ec.clust, ylab = "Euclidean dissimilarity")
  # library(circlize)
  # col_fun = colorRamp2(c(0, 1), c("red","white"))
  # col_fun(seq(0, 1))
  # 
  # library(ComplexHeatmap)
  # Heatmap(as.matrix(comm.ec.dist),col=col_fun,column_title = "Euclidean dissimilarity",name="Distance")
  # 
  # 
  # library(qgraph)
  # qgraph(comm.bc.dist, layout='spring', vsize=3)
  # 
}
# continue here...
dev.off()
write.table(file=paste(out_dir,"intersects_all_subsamples.UMI_FILTERED.csv",sep=""),x=intersects_per_subset,quote = F,sep = ";",row.names = T)
write.table(file=paste(out_dir,"Bray_Curtis_dissimilarity_all_subsamples.UMI_FILTERED.csv",sep=""),x=Bray_Curtis,quote = F,sep = ";",row.names = F)

# Plot just the NMDS
library(vegan)
pdf(file=paste(out_dir,"Array_4probes_Bacteria_all_samples_RpM_subsampled_equal_size_MMSEQS2_NT_NMDS.UMI_FILTERED.pdf",sep=""),height = 7,width = 8) # paper = "a4r", height = 8.3 , width = 11.7

for (i_repeat in 1:num_of_subsamples_repeats)
{
  sample_file = paste(out_dir,"subsample_",i_repeat,"_",Level,".UMI_FILTERED.profile.csv",sep="")
  all_subsample = read.delim(file=sample_file,sep = ";",stringsAsFactors = F)
  all_NO_PROB=all_subsample[!grepl(all_subsample$sample_name,pattern = "PROB"),]
  all_NO_PROB$sample_name <- factor(all_NO_PROB$sample_name,levels = c("A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                       "B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                       "C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                       "D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  for (sample in unique(all_NO_PROB$sample_name))
  {
    s_sum=sum(all_NO_PROB$Count[all_NO_PROB$sample_name==sample])
    message (paste ("QA: ",sample," total ",s_sum,sep=""))
  }
  
  percent_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  count_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  
  # percent_matrix_all_sub_samples=make_table(flat_df = all_NO_PROB_sub_samples,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  # count_matrix_all_sub_samples=make_table(flat_df = all_NO_PROB_sub_samples,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  
  community_matrix=t(as.matrix(percent_matrix_all_samples))
  #community_matrix=t(as.matrix(percent_matrix_all_sub_samples))
  community_matrix[is.na(community_matrix)]=0 # replace NA with 0
  example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                       k=2,trymax = 100,distance = "bray") # The number of reduced dimensions -> ,, autotransform = FALSE,distance ="euclidean"
  # stressplot(example_NMDS)
  # plot (example_NMDS)
  
  #extract NMDS scores (x and y coordinates)
  data.scores = as.data.frame(scores(example_NMDS))
  data.scores$Sample=row.names(data.scores)
  data.scores$Sample=gsub(pattern=".799_1192",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".515_806",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="AmpSeq_",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="1|2",replacement = "",x=data.scores$Sample)
  
  data.scores$Type=row.names(data.scores) 
  data.scores$Type[grepl("^[ABCD]1$",perl=T,rownames(data.scores))]="Array2Probes"
  data.scores$Type[grepl("^[ABCD]2$",perl=T,rownames(data.scores))]="Array4Probes"
  data.scores$Type[grepl("515_806",perl=T,rownames(data.scores))]="AmpSeq_p515"
  data.scores$Type[grepl("799_1192",perl=T,rownames(data.scores))]="AmpSeq_p799"
  
  nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
    geom_point(size = 6, aes( shape = Sample, fill = Type,color=Type))+ 
    # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
    geom_text(aes(label=""),hjust=-0.1, vjust=-0.1) +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    # labs(x = "NMDS1", colour = "Plant", y = "NMDS2", shape = "Type")  + 
    # labs(x = "NMDS1", colour = "Method", y = "NMDS2", shape = "Plant")  + 
    labs(x = "NMDS1", fill = "Method", y = "NMDS2", shape = "Plant",colour = "Method")  + 
    scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    scale_shape_manual(values=c(21, 22, 23,24))+ 
    scale_fill_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    ggtitle(paste("Subsample ",i_repeat," - NMDS of the full profile\nsubsampled to equal size (",subset_size,") on the ",Level, " level (Stress=",signif(example_NMDS$stress,3),")",sep=""))
  print (nmds_plot,width = 4, height = 4)
}
dev.off()




# saveWorkbook(XSLX_Obj, file = XLSX_matrix_file, overwrite = TRUE)

#### END CLEAN SUBSAMPLES WITH NORMALIZED ####
##############################################



#### VERSION 2 ####
#### CLEAN START ANALYZE WITH SUBSAMPLES ON GENUS AND REPEATS (IGNORE AmpSeq C2.799 as it is too samll) 25/10/2021 ####
#############################################################################
# repeat subsample 100 times (only on the genus level) and plot the (1) NMDS (2) UpSetR and record the values of the intersections and subsampled profiles
## FOR THE FULL SAMPLES OMNI4 and filtering rare species from each sample -> here!
set.seed(639245)
# load functions
# src_dir="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/Scripts/" # cluster
src_dir="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/"

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenMatrix <- function(mat) {
  ut <- upper.tri(mat)
  data.frame(
    row = rownames(mat)[row(mat)[ut]],
    column = rownames(mat)[col(mat)[ut]],
    val=(mat)[ut]
  )
}


source(paste (src_dir,"Microbial_Profiles_Functions.R",sep="")) 

sub_sample = function (df,sample_name,total_count) # df is an output of Tax_profile function, {"Value","Count","Percent","sample_name"}
{
  subset=df[as.character(df$sample_name)==as.character(sample_name),]
  subset_expanded=data.frame()
  for (i in 1:nrow(subset)){
    tmp=data.frame(Value=rep(subset$Value[i],subset$Count[i]),count=rep(1,subset$Count[i]))
    subset_expanded=rbind(subset_expanded,tmp)
  }
  subset_subsampled_tmp=subset_expanded[sample(1:nrow(subset_expanded),size = total_count,replace = F),]
  subset_subsampled=data.frame()
  sampled_tax=unique(subset_subsampled_tmp$Value)
  sampled_sum=sum(subset_subsampled_tmp$count)
  for (tax in sampled_tax) {
    tax_count=sum(subset_subsampled_tmp$count[subset_subsampled_tmp$Value==tax])
    tax_percent=tax_count/sampled_sum
    subset_subsampled=rbind(subset_subsampled,data.frame(Value=tax,Count=tax_count,Percent=tax_percent,sample_name=sample_name))
  }
  subset_subsampled=subset_subsampled[order(subset_subsampled$Count,decreasing = T),]
  return(subset_subsampled)
}
# Genus levels stackplots and heatmaps
# Another option: https://cran.r-project.org/web/packages/Polychrome/vignettes/testgg.html
library(Polychrome)
library(gplots)

# Global vars and data reading
path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/multimaped_filtered/MMSEQS2/UMI_FILTERED/"
library(openxlsx)
library(Polychrome)
library(gplots)
library(ggplot2)
DB_type="MMSEQS2"
num_of_subsamples_repeats=100
out_dir=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/",DB_type,"/UMI_FILTERED/","concise_4_probesOMNI","/subsample/repeats_",num_of_subsamples_repeats,"/",sep="")
dir.create(out_dir, showWarnings = FALSE,recursive = TRUE)

# The data itself
# The full sample
list_of_files_unsplitted=list.files(path = path,
                                    pattern = ".*_Ar_multimap_bacteria_R2.usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt")
list_of_files_unsplitted=paste(path,list_of_files_unsplitted,sep="")
sample_name_unsplitted=basename(list_of_files_unsplitted)

# The splitted by prob
list_of_files_splitted=list.files(path = path,
                                  pattern = ".*_Ar_multimap_bacteria_R2.PROB_.*.usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt")
list_of_files_splitted=paste(path,list_of_files_splitted,sep="")
sample_name_splitted=basename(list_of_files_splitted)

# Derek primer 515
Derek_AmpSeq_515_path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/AmpSeq/515_806/"
list_of_Derek_Amp_seq_515_files=list.files(path = Derek_AmpSeq_515_path,
                                           pattern = ".*uniq_lineage_count.txt")

list_of_Derek_Amp_seq_515_files=paste(Derek_AmpSeq_515_path,list_of_Derek_Amp_seq_515_files,sep="")
sample_name_Derek_Amp_seq_515=basename(list_of_Derek_Amp_seq_515_files)
sample_name_Derek_Amp_seq_515=paste(sample_name_Derek_Amp_seq_515,".515_806",sep="")
# Derek primer 799
Derek_AmpSeq_799_path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/AmpSeq/799_1192/"
list_of_Derek_Amp_seq_799_files=list.files(path = Derek_AmpSeq_799_path,
                                           pattern = ".*uniq_lineage_count.txt")

list_of_Derek_Amp_seq_799_files=paste(Derek_AmpSeq_799_path,list_of_Derek_Amp_seq_799_files,sep="")
sample_name_Derek_Amp_seq_799=basename(list_of_Derek_Amp_seq_799_files)
sample_name_Derek_Amp_seq_799=paste(sample_name_Derek_Amp_seq_799,".799_1192",sep="")

list_of_files=c(list_of_files_unsplitted,list_of_files_splitted,list_of_Derek_Amp_seq_799_files,list_of_Derek_Amp_seq_515_files)
sample_names=c(sample_name_unsplitted,sample_name_splitted,sample_name_Derek_Amp_seq_799,sample_name_Derek_Amp_seq_515)

sample_names=gsub(".usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt",replacement = "",fixed = T,x = sample_names)
sample_names=gsub(".usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.uniq_lineage_count.txt",replacement = "",fixed = T,x = sample_names)

sample_names=gsub("_Ar_multimap_bacteria_R2",replacement = "",fixed = T,x = sample_names)
sample_names=gsub("191029-069_",replacement = "",fixed = T,x = sample_names)

tmp_df=data.frame(file=list_of_files,sample_name=sample_names)
lineages_files_and_tag=tmp_df

# change the names of the AmpSeq samples
lineages_files_and_tag$sample_name=gsub("F5p515F2F4F6",replacement = "AmpSeq_A",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("C5p515F2F4F6",replacement = "AmpSeq_B",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("A2p515F2F4F6",replacement = "AmpSeq_C",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("E3p515F1F3F5",replacement = "AmpSeq_D",fixed = T,x = lineages_files_and_tag$sample_name)

# remove the bad splitted ones
lineages_files_and_tag=lineages_files_and_tag[!grepl("1.PROB_479",lineages_files_and_tag$sample_name),]
lineages_files_and_tag=lineages_files_and_tag[!grepl("1.PROB_1265",lineages_files_and_tag$sample_name),]

lineages_files_and_tag_x2=lineages_files_and_tag[!grepl(lineages_files_and_tag$sample_name,pattern = "([ABCD]1)",perl=T),]
# Should we exclude another one --> YES
# AmpSeq_C.799_1192 603 reads
lineages_files_and_tag_x2=lineages_files_and_tag_x2[!grepl(lineages_files_and_tag_x2$sample_name,pattern = "AmpSeq_C.799_1192",fixed = T),]

# plot definitions
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
P40 <- sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
swatch(P40)
names(P40) <- NULL

col1<- colorRampPalette(c("white", "blue"))(50)
breaks1 = c(seq(0,0.5,length=51))

# out files
XLSX_matrix_file=paste(out_dir,"datasets_tables_Bacteria_AmpSeq_vs_Array_4probes_subsampled_equal_size_",DB_type,".UMI_FILTERED.xlsx",sep="")
XSLX_Obj=createWorkbook(XLSX_matrix_file)
pdf(file=paste(out_dir,"Array_4probes_Bacteria_all_samples_subsampled_equal_size_MMSEQS2_NT_stacks_labeled_NMDS_and_intersects.UMI_FILTERED.pdf",sep=""),height = 16.5,width = 23.4) # paper = "a4r", height = 8.3 , width = 11.7

Level="genus" # for simplicity we only analyze genus level, we can do a loop for all... for (Level in c ("genus")) # c ("phylum", "class","order","family","genus") #=c("Count", "superkingdom","phylum", "class","order","family","genus","species")

# read all data and decide the size to sample
all=Tax_profile(lineage_files_and_tag_df=lineages_files_and_tag_x2,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
                levels_to_account_for_exclude = c("genus"))

# we don't use the prob samples so lets remove them, if we want to use them we should sample each time the probs data of each sample and sum them up to give the full sample!
all=all[!grepl(all$sample_name,pattern = "PROB"),]

# subset
# count the total number of reads in each sample
sample_count=data.frame()
for (sample_i in unique(all$sample_name)) {
  sample_total=sum(all$Count[all$sample_name==sample_i])
  sample_count=rbind(sample_count,data.frame(sample_name=sample_i,total_size=sample_total))
}

subset_size=min(sample_count$total_size)

Bray_Curtis=data.frame()
intersects_per_subset=data.frame()
# To do: aggregate results...
# add the subsampled profile to xls
# df for intersect size
# 
for (i_repeat in 1:num_of_subsamples_repeats) # the i_repeat should be accounted later in figures and aggregator
{
  message (paste("- Sstart subsample repeat",i_repeat))
  all_subsample=data.frame()
  for (sample_i in unique(all$sample_name)) {
    sub_sample_df=sub_sample(df = all,sample_name = sample_i, subset_size)
    message(paste ("\t-- subsample ",subset_size," reads out of ",sample_count$total_size[sample_count$sample_name==sample_i]," reads of sample ",sample_i))
    all_subsample=rbind(all_subsample,sub_sample_df)
  }
  # save subsample
  write.table(file = paste(out_dir,"subsample_",i_repeat,"_",Level,".UMI_FILTERED.profile.csv",sep=""),x =all_subsample, quote = F,row.names = F,sep = ";",col.names = T)
  
  ### HERE! --> Continue to update and replace all df with all_subsample
  all_with_others=data.frame()
  other_cutoff=0.02
  for (sample in unique (all_subsample$sample_name))
  {
    subset=all_subsample[all_subsample$sample_name==sample,]
    sabset_to_add=subset[subset$Percent>other_cutoff,]
    other_percent=sum(subset$Percent[subset$Percent<=other_cutoff])
    other_count=sum(subset$Count[subset$Percent<=other_cutoff])
    sabset_to_add=rbind(sabset_to_add,data.frame(Value=paste("Other (<=",other_cutoff*100,"%)",sep=""),Count=other_count,Percent=other_percent,sample_name=sample))
    all_with_others=rbind(all_with_others,sabset_to_add)
  }
  # no need for the PROB samples here --> before the change the data frame used was: all_with_others
  all_with_others_NO_PROB=all_with_others[!grepl(all_with_others$sample_name,pattern = "PROB"),]
  all_with_others_NO_PROB$sample_name <- factor(all_with_others_NO_PROB$sample_name,levels = c("A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                                               "B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                                               "C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                                               "D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  p=ggplot(all_with_others_NO_PROB, aes(fill=Value, y=Percent, x=sample_name,label=paste(signif(Percent*100,2),"% ",Value,sep=""))) + 
    geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
    theme(legend.text=element_text(size=15)) +
    scale_fill_manual(values = P40) + 
    ggtitle(paste("Subsample",i_repeat,"composition",Level)) +
    theme(plot.title = element_text(hjust = 0.5))
  print (p)
    
  # scale_fill_manual(values = as.vector(palette36))
  # scale_fill_brewer( palette = "Paired")
  
  # the heatmap - is it needed?!?
  # percent_matrix=make_table(flat_df = all_with_others_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  # count_matrix=make_table(flat_df = all_with_others_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  ## TODO: take into account that it may be that some fall under the other category...
  ## we can fix it with the same proceudre we did for missing data (we need to replace the NA with number and reduce it from the Other category)
  # percent_matrix_NA_0=percent_matrix
  # percent_matrix_NA_0[is.na(percent_matrix_NA_0)]=0
  # heatmap.2(x = as.matrix(percent_matrix_NA_0), col = col1, symm = FALSE,breaks=breaks1,trace = "none",margins=c(15,15),cexRow=0.9, main=paste(Level," - NA as 0, DB: ",DB_type, sep="")) # main = paste("Abundance of ",tax_level, " (out of the all reads)",sep="")
  
  # sheetName = paste(Level,"_wOt_",other_cutoff,"_percent_NA_0",sep="")
  # addWorksheet(XSLX_Obj, sheetName)
  # writeDataTable(x=percent_matrix_NA_0, wb = XSLX_Obj,sheet = sheetName,rowNames = TRUE)
  
  # Do NMDS between the full samples profiles - https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/ / https://jkzorz.github.io/2019/06/06/NMDS.html
  ## CONSIDER TO FILTER OUT THE RARE SPECIES.... TO AVOID DOUBLE ZERO AND OTHER PROBLEMS?
  
  # From each sample filter out species below a cutoff
  # table(all[all$Count<2,"Value"])
  # table(all$sample_name)
  
  ## length(subset$Percent[subset$Percent<=0.0005]) sample="AmpSeq_A.515_806"
  library(vegan)
  all_NO_PROB=all_subsample[!grepl(all_subsample$sample_name,pattern = "PROB"),]
  all_NO_PROB$sample_name <- factor(all_NO_PROB$sample_name,levels = c("A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                       "B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                       "C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                       "D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  for (sample in unique(all_NO_PROB$sample_name))
  {
    s_sum=sum(all_NO_PROB$Count[all_NO_PROB$sample_name==sample])
    message (paste ("QA: ",sample," total ",s_sum,sep=""))
  }

  percent_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  count_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  
  # percent_matrix_all_sub_samples=make_table(flat_df = all_NO_PROB_sub_samples,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  # count_matrix_all_sub_samples=make_table(flat_df = all_NO_PROB_sub_samples,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  
  community_matrix=t(as.matrix(percent_matrix_all_samples))
  #community_matrix=t(as.matrix(percent_matrix_all_sub_samples))
  community_matrix[is.na(community_matrix)]=0 # replace NA with 0
  example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                       k=2,trymax = 100,distance = "bray") # The number of reduced dimensions -> ,, autotransform = FALSE,distance ="euclidean"
  stressplot(example_NMDS)
  # plot (example_NMDS)
  
  #extract NMDS scores (x and y coordinates)
  data.scores = as.data.frame(scores(example_NMDS))
  data.scores$Sample=row.names(data.scores)
  data.scores$Sample=gsub(pattern=".799_1192",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".515_806",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="AmpSeq_",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="1|2",replacement = "",x=data.scores$Sample)
  
  data.scores$Type=row.names(data.scores) 
  data.scores$Type[grepl("^[ABCD]1$",perl=T,rownames(data.scores))]="Array2Probes"
  data.scores$Type[grepl("^[ABCD]2$",perl=T,rownames(data.scores))]="Array4Probes"
  data.scores$Type[grepl("515_806",perl=T,rownames(data.scores))]="AmpSeq_p515"
  data.scores$Type[grepl("799_1192",perl=T,rownames(data.scores))]="AmpSeq_p799"
  
  nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
    geom_point(size = 6, aes( shape = Sample, fill = Type,color=Type))+ 
    # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
    geom_text(aes(label=""),hjust=-0.1, vjust=-0.1) +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    # labs(x = "NMDS1", colour = "Plant", y = "NMDS2", shape = "Type")  + 
    # labs(x = "NMDS1", colour = "Method", y = "NMDS2", shape = "Plant")  + 
    labs(x = "NMDS1", fill = "Method", y = "NMDS2", shape = "Plant",colour = "Method")  + 
    scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    scale_shape_manual(values=c(21, 22, 23,24))+ 
    scale_fill_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    ggtitle(paste("Subsample ",i_repeat," - NMDS of the full profile\nsubsampled to equal size (",subset_size,") on the ",Level, " level (Stress=",signif(example_NMDS$stress,3),")",sep=""))
  print (nmds_plot,width = 4, height = 4)
  # dev.off()
  # colors for the methods and shape for plants ?
  
  # hirarchial clustering and distances heatmaps
  # calculate Bray-Curtis distance among samples
  comm.bc.dist <- vegdist(community_matrix, method = "bray")
  # cluster communities using average-linkage algorithm
  comm.bc.clust <- hclust(comm.bc.dist, method = "average")
  # plot cluster diagram
  plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")
  library(circlize)
  col_fun = colorRamp2(c(0, 1), c("red","white"))
  col_fun(seq(0, 1))
  
  library(ComplexHeatmap)
  dist_mat=as.matrix(comm.bc.dist)
  pHeatmapSim=Heatmap(dist_mat,col=col_fun,column_title = paste("Subsample ",i_repeat,"- Bray-Curtis dissimilarity",sep=""),name="Distance",
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f", dist_mat[i, j]), x, y, gp = gpar(fontsize = 10))
          })
   print(pHeatmapSim) 
   dist_mat_flat=flattenMatrix(dist_mat)
   names(dist_mat_flat)=c("sample1","sample2",paste("BC_Diss","_sample_",i_repeat,sep=""))
   if (nrow(Bray_Curtis)==0)
   {
     Bray_Curtis=dist_mat_flat
   } else {
     Bray_Curtis=merge(Bray_Curtis,dist_mat_flat,by=c("sample1","sample2"),all=T)
   }  
    # # clustering with k mean -> consider this: https://rstudio-pubs-static.s3.amazonaws.com/274963_9e95868f44d64a008e4bf197cdef079e.html
    # find_k <- function(data, plot_title = ""){
    #   max_clusters=dim(as.matrix(data))[1]
    #   #create data for within group sum of square
    #   wss = numeric(max_clusters)
    #   for (k in 1:(max_clusters-1)) wss[k] <- kmeans(data, k, nstart = 20)$tot.withinss
    #   
    #   #create data for average silhouette distance
    #   library(cluster)
    #   asw <- numeric(max_clusters)
    #   for (k in 2:(max_clusters-1)) asw[k] <- pam(data, k)$silinfo$avg.width
    #   
    #   res=data.frame(k=1:max_clusters,WSS=wss,asw=asw)
    #   #create s cree plot
    #   par(mar=c(5, 4, 4, 6))
    #   plot(1:max_clusters, wss, type = "b", main = plot_title, xlab = "Number of Clusters", ylab = "Within groups sum of squares")
    #   par(new = T)
    #   plot(1:max_clusters, asw, type = "l", lty = 2, col = "red", axes = F, xlab = NA, ylab = NA)
    #   axis(side = 4)
    #   mtext("Average Silhouette width", side = 4, line = 3)
    #   legend("topright", legend = c("WSS", "Si"), lty = c(1,2), col = c("black", "red"))
    #   return(res)
    # }
    # res=find_k(comm.bc.dist,"bray")
    # #conduct k-means to find clusters
    # km=pam(x = comm.bc.dist,k = 2)
    # plot(silhouette(km, comm.bc.dist), main = paste("Silhoulette plot using Bray-Curtis distance") , col = 1:2, border = NA, cex = 0.6)
    # 
    
    # Do the upSetR plots
    unique(all_NO_PROB$sample_name)
    tax_list=list(A2=all_subsample$Value[all_subsample$sample_name=="A2"],
                  AmpSeq_A_799=all_subsample$Value[all_subsample$sample_name=="AmpSeq_A.799_1192"],
                  AmpSeq_A_515=all_subsample$Value[all_subsample$sample_name=="AmpSeq_A.515_806"],
                  
                  B2=all_subsample$Value[all_subsample$sample_name=="B2"],
                  AmpSeq_B_799=all_subsample$Value[all_subsample$sample_name=="AmpSeq_B.799_1192"],
                  AmpSeq_B_515=all_subsample$Value[all_subsample$sample_name=="AmpSeq_B.515_806"],
                  
                  C2=all_subsample$Value[all_subsample$sample_name=="C2"],
                  AmpSeq_C_799=all_subsample$Value[all_subsample$sample_name=="AmpSeq_C.799_1192"],
                  AmpSeq_C_515=all_subsample$Value[all_subsample$sample_name=="AmpSeq_C.515_806"],
                  
                  D2=all_subsample$Value[all_subsample$sample_name=="D2"],
                  AmpSeq_D_799=all_subsample$Value[all_subsample$sample_name=="AmpSeq_D.799_1192"],
                  AmpSeq_D_515=all_subsample$Value[all_subsample$sample_name=="AmpSeq_D.515_806"]
    )
 
    m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
    cs = comb_size(m)
    ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),column_title = paste("A2 (Omni4Probs)",Level))
    ht = draw(ht)
    co = column_order(ht)
    nc = ncol(m)
    decorate_annotation("Intersection\nsize", {
      grid.text(cs[co], 
                x = 1:nc, 
                y = unit(cs[co], "native") + unit(1, "mm"), 
                gp = gpar(fontsize = 8), 
                just = "bottom",
                default.units = "native")
    })
    # separately for each plant 
    # pdf(file=paste(out_dir,"AmpSeq_Array4probes_Bacteria_all_samples_MMSEQS2_NT_upSetR.pdf",sep=""),height = 3,width = 6.5) 
    i_subset_intersects=data.frame()
    for (p in c("A","B","C","D"))
    {
      tax_list=list(Array=all_subsample$Value[all_subsample$sample_name==paste(p,"2",sep="")],
                    AmpSeq_p799=all_subsample$Value[all_subsample$sample_name==paste("AmpSeq_",p,".799_1192",sep="")],
                    AmpSeq_p515=all_subsample$Value[all_subsample$sample_name==paste("AmpSeq_",p,".515_806",sep="")])
      # for the sub samples
      # tax_list=list(Array=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste(p,"2",sep="")]),
      #               AmpSeq_799=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste("AmpSeq_",p,".799_1192",sep="")]),
      #               AmpSeq_515=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste("AmpSeq_",p,".515_806",sep="")]))
      m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
      cs = comb_size(m)
      row_size = set_size(m)
      ht = UpSet(m, 
                 top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),
                 right_annotation = upset_right_annotation(m, ylim = c(0, 1.1*max(set_size(m)))),
                 column_title = paste("Subsample",i_repeat,"- Plant",p,Level,"subsampled to",subset_size,"reads"))
      ht = draw(ht)
      co = column_order(ht)
      row_od = row_order(ht)
      nc = ncol(m)
      decorate_annotation("Intersection\nsize", {
        grid.text(cs[co], 
                  x = 1:nc, 
                  y = unit(cs[co], "native") + unit(1, "mm"), 
                  gp = gpar(fontsize = 10), 
                  just = "bottom",
                  default.units = "native")
      })
      decorate_annotation("Set size", {
        grid.text(row_size[row_od], 
                  unit(row_size[row_od], "native") + unit(1, "mm"), 
                  rev(seq_len(length(row_size))), 
                  default.units = "native", just = "bottom", rot = -90,
                  gp = gpar(fontsize = 10))
      })
      # library(eulerr)
      
      
      intersects=c("Array" = 0, "AmpSeq_p799" = 0, "AmpSeq_p515" = 0,
                   "Array&AmpSeq_p799&AmpSeq_p515"=0,  
                   "Array&AmpSeq_p799"=0,
                   "Array&AmpSeq_p515"=0,
                   "AmpSeq_p799&AmpSeq_p515"=0)
      for (i in names(cs)) { # assign
        if (i=="100") {intersects[["Array"]]=cs[[i]]}
        if (i=="010") {intersects[["AmpSeq_p799"]]=cs[[i]]}
        if (i=="001") {intersects[["AmpSeq_p515"]]=cs[[i]]}
        if (i=="111") {intersects[["Array&AmpSeq_p799&AmpSeq_p515"]]=cs[[i]]}
        if (i=="110") {intersects[["Array&AmpSeq_p799"]]=cs[[i]]}
        if (i=="101") {intersects[["Array&AmpSeq_p515"]]=cs[[i]]}
        if (i=="011") {intersects[["AmpSeq_p799&AmpSeq_p515"]]=cs[[i]]}
      }
      intersect_df=as.data.frame(intersects)
      row.names(intersect_df)=paste(p,row.names(intersect_df),sep="_")
      names(intersect_df)=paste("subset",i_repeat,"intersect",sep="_")
      
      i_subset_intersects=rbind(i_subset_intersects,intersect_df)
      names(i_subset_intersects)=paste("subset",i_repeat,"intersect",sep="_")
      # fit <- euler(intersects)
      # # Customize colors, remove borders, bump alpha, color labels white
      # plot(fit,
      #      fills = list(fill = c("#E69F00", "#56B4E9","#009E73","#CC79A7"), alpha = 0.5),
      #      labels = list(col = "white", font = 4))
    }
    if (nrow(intersects_per_subset)==0) {
      intersects_per_subset=i_subset_intersects
    } else {
      intersects_per_subset=merge(intersects_per_subset,i_subset_intersects,by="row.names",all=T)
      row.names(intersects_per_subset)=intersects_per_subset$Row.names
      intersects_per_subset=intersects_per_subset[,-1]
    }
    # dev.off()
    
    ## plot the subsamples profiles
    # all_subsamples_with_others=data.frame()
    # other_cutoff=0.02
    # for (sample in unique (all_NO_PROB_sub_samples$sample_name))
    # {
    #  subset=all_NO_PROB_sub_samples[all_NO_PROB_sub_samples$sample_name==sample,]
    #  names(subset)=c("Value","Count","Percent","sample_name")
    #  sabset_to_add=subset[subset$Percent>other_cutoff,]
    #  other_percent=sum(subset$Percent[subset$Percent<=other_cutoff])
    #  other_count=sum(subset$Count[subset$Percent<=other_cutoff])
    #  sabset_to_add=rbind(sabset_to_add,data.frame(Value=paste("Other (<=",other_cutoff*100,"%)",sep=""),Count=other_count,Percent=other_percent,sample_name=sample))
    #  all_subsamples_with_others=rbind(all_subsamples_with_others,sabset_to_add)
    # }
    ## no need for the PROB samples here --> before the change the data frame used was: all_with_others
    # all_subsamples_with_others_NO_PROB=all_subsamples_with_others[!grepl(all_with_others$sample_name,pattern = "PROB"),]
    # all_subsamples_with_others_NO_PROB$sample_name <- factor(all_subsamples_with_others_NO_PROB$sample_name,levels = c("A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
    #                                                                                              "B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
    #                                                                                              "C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
    #                                                                                              "D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
    # p=ggplot(all_subsamples_with_others_NO_PROB, aes(fill=Value, y=Percent, x=sample_name,label=paste(signif(Percent*100,2),"% ",Value,sep=""))) + 
    #   geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
    #   theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
    #   theme(legend.text=element_text(size=15)) +
    #   scale_fill_manual(values = P40) + 
    #   ggtitle(paste("Composition",Level)) +
    #   theme(plot.title = element_text(hjust = 0.5))
    # print (p)
    
    # TO SEE IF FURTHER UPDATE IS NEEDED....
    # # plot NMDS for each of the OMNI4 OMNI2 and its probes
    # for (sample in c("Ax","Bx","Cx","Dx"))
    # {
    #   samples_regex=NA
    #   if (sample=="Ax") {samples_regex="F5p515F2F4F6|A1|A2|AmpSeq_A"}
    #   if (sample=="Bx") {samples_regex="C5p515F2F4F6|B1|B2|AmpSeq_B"}
    #   if (sample=="Cx") {samples_regex="A2p515F2F4F6|C1|C2|AmpSeq_C"}
    #   if (sample=="Dx") {samples_regex="E3p515F1F3F5|D1|D2|AmpSeq_D"}
    #   sample_subset=lineages_files_and_tag[grepl(samples_regex,lineages_files_and_tag$sample_name),]
    #   if (sample=="Ax") 
    #   {
    #     sample_subset=sample_subset[!grepl("A2p515F2F4F6",sample_subset$sample_name),]
    #   }
    #   
    #   all_sample_subset=Tax_profile(lineage_files_and_tag_df=sample_subset,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
    #                                 levels_to_account_for_exclude = c("genus"))
    #   # subset
    #   # count the total number of reads in each sample
    #   sample_count=data.frame()
    #   for (sample_i in unique(all_sample_subset$sample_name)) {
    #     sample_total=sum(all_sample_subset$Count[all_sample_subset$sample_name==sample_i])
    #     sample_count=rbind(sample_count,data.frame(sample_name=sample_i,total_size=sample_total))
    #   }
    #   all_subsample=data.frame()
    #   subset_size=min(sample_count$total_size)
    #   for (sample_i in unique(all_sample_subset$sample_name)) {
    #     sub_sample_df=sub_sample(df = all_sample_subset,sample_name = sample_i, subset_size)
    #     message(paste ("-- subsample ",subset_size," reads out of ",sample_count$total_size[sample_count$sample_name==sample_i]," reads of sample ",sample_i))
    #     all_subsample=rbind(all_subsample,sub_sample_df)
    #   }
    #   ### HERE! --> Continue to update and replace all df with all_subsample
    #   percent_matrix_all_samples=make_table(flat_df = all_subsample,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
    #   count_matrix_all_samples=make_table(flat_df = all_subsample,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
    #   
    #   community_matrix=t(as.matrix(percent_matrix_all_samples))
    #   #community_matrix=t(as.matrix(percent_matrix_all_sub_samples))
    #   community_matrix[is.na(community_matrix)]=0 # replace NA with 0
    #   example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
    #                        k=2,trymax = 100,distance = "bray") # The number of reduced dimensions -> ,, autotransform = FALSE,distance ="euclidean"
    #   stressplot(example_NMDS)
    #   # plot (example_NMDS)
    #   
    #   #extract NMDS scores (x and y coordinates)
    #   data.scores = as.data.frame(scores(example_NMDS))
    #   data.scores$Sample=row.names(data.scores) # AmpSeq_p515; AmpSeq_p799; pr_1265; pr_unk; pr_479; pr_902; pr_799 -> these will have colors
    #   data.scores$Sample=gsub(pattern="PROB_799F",replacement = "pr_799",x=data.scores$Sample)
    #   data.scores$Sample=gsub(pattern="PROB_902R",replacement = "pr_902",x=data.scores$Sample)
    #   data.scores$Sample=gsub(pattern="PROB_479",replacement = "pr_479",x=data.scores$Sample)
    #   data.scores$Sample=gsub(pattern="PROB_1265",replacement = "pr_1265",x=data.scores$Sample)
    #   data.scores$Sample=gsub(pattern="PROB_UNKNOWN",replacement = "pr_unk",x=data.scores$Sample)
    #   data.scores$Sample=gsub(pattern="^[ABCD][12]$",replacement = "Full",x=data.scores$Sample,perl = T)
    #   data.scores$Sample=gsub(pattern="^[ABCD][12].",replacement = "",x=data.scores$Sample,perl = T)
    #   data.scores$Sample=gsub(pattern="515_806",replacement = "AmpSeq_p515",x=data.scores$Sample)
    #   data.scores$Sample=gsub(pattern="799_1192",replacement = "AmpSeq_p799",x=data.scores$Sample)
    #   data.scores$Sample=gsub(pattern="^AmpSeq_[ABCD].",replacement = "",x=data.scores$Sample,perl=T)
    #   
    #   # type -> will be symbol {AmpSeq | Omni2Probes | Omni4Probes}
    #   data.scores$Type=row.names(data.scores) 
    #   data.scores$Type[grepl("^[ABCD]1",perl=T,rownames(data.scores))]="Omni2Probes"
    #   data.scores$Type[grepl("^[ABCD]2",perl=T,rownames(data.scores))]="Omni4Probes"
    #   data.scores$Type[grepl("AmpSeq",perl=T,rownames(data.scores))]="AmpSeq"
    #   
    #   pdf(file=paste(out_dir,"AmpSeq_Array4probes_Bacteria_all_samples_MMSEQS2_NT_NMDS.",sample,".pdf",sep=""),height = 7,width = 8)
    #   nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    #     # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
    #     geom_point(size = 6, aes( shape = Type, fill = Sample,color=Sample))+ 
    #     # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
    #     geom_text(aes(label=Sample),hjust=-0.3, vjust=-0.3) +
    #     theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
    #           axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
    #           legend.text = element_text(size = 12, face ="bold", colour ="black"), 
    #           legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
    #           axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
    #           legend.title = element_text(size = 14, colour = "black", face = "bold"), 
    #           panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
    #           legend.key=element_blank()) + 
    #     # labs(x = "NMDS1", colour = "Plant", y = "NMDS2", shape = "Type")  + 
    #     # labs(x = "NMDS1", colour = "Method", y = "NMDS2", shape = "Plant")  + 
    #     labs(x = "NMDS1", fill = "Probe/Primer", y = "NMDS2", shape = "Method",colour = "Probe/Primer")  + 
    #     scale_colour_manual(values = c("#56B4E9","#0072B2","#000000","#E69F00","#009E73","#F0E442","#D55E00","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    #     scale_shape_manual(values=c(21,24,22))+ 
    #     scale_fill_manual(values = c("#56B4E9","#0072B2","#000000","#E69F00","#009E73","#F0E442","#D55E00","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    #     ggtitle(paste("NMDS of the full profile subset to equal size (",subset_size,") on the ",Level, " level (Stress=",signif(example_NMDS$stress,3),")",sep=""))
    #   print (nmds_plot)
    #   dev.off()
    # }
    # 
    # # drafts
    # km_result = transpose(km)
    # 
    # km = map(comm.bc.dist, pam, k = 3)
    # 
    # # extract cluster data
    # km_result = transpose(km)
    # km_cluster = as.data.frame(km_result$cluster)
    # colnames(km_cluster) <- dist_meas
    # 
    # plot(silhouette(km[[i]], comm.bc.dist), main = paste("Silhoulette plot using Bray-Curtis distance") , col = 1:3, border = NA, cex = 0.6)
    # 
    # 
    # 
    # # hirarchial clustering and distances heatmaps
    # # calculate Euclidean distance among samples
    # comm.ec.dist <- vegdist(community_matrix, method = "euclidean")
    # 
    # # cluster communities using average-linkage algorithm
    # comm.ec.clust <- hclust(comm.ec.dist, method = "average")
    # # plot cluster diagram
    # plot(comm.ec.clust, ylab = "Euclidean dissimilarity")
    # library(circlize)
    # col_fun = colorRamp2(c(0, 1), c("red","white"))
    # col_fun(seq(0, 1))
    # 
    # library(ComplexHeatmap)
    # Heatmap(as.matrix(comm.ec.dist),col=col_fun,column_title = "Euclidean dissimilarity",name="Distance")
    # 
    # 
    # library(qgraph)
    # qgraph(comm.bc.dist, layout='spring', vsize=3)
    # 
}

dev.off()
write.table(file=paste(out_dir,"intersects_all_subsamples.UMI_FILTERED.csv",sep=""),x=intersects_per_subset,quote = F,sep = ";",row.names = T)
# saveWorkbook(XSLX_Obj, file = XLSX_matrix_file, overwrite = TRUE)
write.table(file=paste(out_dir,"Bray_Curtis_dissimilarity_all_subsamples.UMI_FILTERED.csv",sep=""),x=Bray_Curtis,quote = F,sep = ";",row.names = F)

# Plot just the NMDS
library(vegan)
pdf(file=paste(out_dir,"Array_4probes_Bacteria_all_samples_subsampled_equal_size_MMSEQS2_NT_NMDS.UMI_FILTERED.pdf",sep=""),height = 7,width = 8) # paper = "a4r", height = 8.3 , width = 11.7

for (i_repeat in 1:num_of_subsamples_repeats)
{
  sample_file = paste(out_dir,"subsample_",i_repeat,"_",Level,".UMI_FILTERED.profile.csv",sep="")
  all_subsample = read.delim(file=sample_file,sep = ";",stringsAsFactors = F)
  all_NO_PROB=all_subsample[!grepl(all_subsample$sample_name,pattern = "PROB"),]
  all_NO_PROB$sample_name <- factor(all_NO_PROB$sample_name,levels = c("A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                       "B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                       "C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                       "D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  for (sample in unique(all_NO_PROB$sample_name))
  {
    s_sum=sum(all_NO_PROB$Count[all_NO_PROB$sample_name==sample])
    message (paste ("QA: ",sample," total ",s_sum,sep=""))
  }
  
  percent_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  count_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  
  # percent_matrix_all_sub_samples=make_table(flat_df = all_NO_PROB_sub_samples,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  # count_matrix_all_sub_samples=make_table(flat_df = all_NO_PROB_sub_samples,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
  
  community_matrix=t(as.matrix(percent_matrix_all_samples))
  #community_matrix=t(as.matrix(percent_matrix_all_sub_samples))
  community_matrix[is.na(community_matrix)]=0 # replace NA with 0
  example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                       k=2,trymax = 100,distance = "bray") # The number of reduced dimensions -> ,, autotransform = FALSE,distance ="euclidean"
  # stressplot(example_NMDS)
  # plot (example_NMDS)
  
  #extract NMDS scores (x and y coordinates)
  data.scores = as.data.frame(scores(example_NMDS))
  data.scores$Sample=row.names(data.scores)
  data.scores$Sample=gsub(pattern=".799_1192",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".515_806",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="AmpSeq_",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="1|2",replacement = "",x=data.scores$Sample)
  
  data.scores$Type=row.names(data.scores) 
  data.scores$Type[grepl("^[ABCD]1$",perl=T,rownames(data.scores))]="Array2Probes"
  data.scores$Type[grepl("^[ABCD]2$",perl=T,rownames(data.scores))]="Array4Probes"
  data.scores$Type[grepl("515_806",perl=T,rownames(data.scores))]="AmpSeq_p515"
  data.scores$Type[grepl("799_1192",perl=T,rownames(data.scores))]="AmpSeq_p799"
  
  nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
    geom_point(size = 6, aes( shape = Sample, fill = Type,color=Type))+ 
    # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
    geom_text(aes(label=""),hjust=-0.1, vjust=-0.1) +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    # labs(x = "NMDS1", colour = "Plant", y = "NMDS2", shape = "Type")  + 
    # labs(x = "NMDS1", colour = "Method", y = "NMDS2", shape = "Plant")  + 
    labs(x = "NMDS1", fill = "Method", y = "NMDS2", shape = "Plant",colour = "Method")  + 
    scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    scale_shape_manual(values=c(21, 22, 23,24))+ 
    scale_fill_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    ggtitle(paste("Subsample ",i_repeat," - NMDS of the full profile\nsubsampled to equal size (",subset_size,") on the ",Level, " level (Stress=",signif(example_NMDS$stress,3),")",sep=""))
  print (nmds_plot,width = 4, height = 4)
}
dev.off()
#### END ANALYZE WITH SUBSAMPLES ON GENUS AND REPEATS  ####
############################################################

#### START clean For the full sample without sub-sampleming ####
####################################################
# load functions
library(ComplexHeatmap)
# src_dir="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/Scripts/" # cluster
src_dir="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/"

source(paste (src_dir,"Microbial_Profiles_Functions.R",sep="")) 
# Genus levels stackplots and heatmaps
# Another option: https://cran.r-project.org/web/packages/Polychrome/vignettes/testgg.html
library(Polychrome)
library(gplots)

# Global vars and data reading
path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/multimaped_filtered/MMSEQS2/UMI_FILTERED/"
library(openxlsx)
library(Polychrome)
library(gplots)
library(ggplot2)
DB_type="MMSEQS2"
out_dir=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/",DB_type,"/UMI_FILTERED/","concise_4_probesOMNI","/",sep="")
dir.create(out_dir, showWarnings = FALSE,recursive = TRUE)

# The data itself
# The full sample
list_of_files_unsplitted=list.files(path = path,
                                    pattern = ".*_Ar_multimap_bacteria_R2.usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt")
list_of_files_unsplitted=paste(path,list_of_files_unsplitted,sep="")
sample_name_unsplitted=basename(list_of_files_unsplitted)

# The splitted by prob
list_of_files_splitted=list.files(path = path,
                                  pattern = ".*_Ar_multimap_bacteria_R2.PROB_.*.usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt")
list_of_files_splitted=paste(path,list_of_files_splitted,sep="")
sample_name_splitted=basename(list_of_files_splitted)

# Derek primer 515
Derek_AmpSeq_515_path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/AmpSeq/515_806/"
list_of_Derek_Amp_seq_515_files=list.files(path = Derek_AmpSeq_515_path,
                                           pattern = ".*uniq_lineage_count.txt")

list_of_Derek_Amp_seq_515_files=paste(Derek_AmpSeq_515_path,list_of_Derek_Amp_seq_515_files,sep="")
sample_name_Derek_Amp_seq_515=basename(list_of_Derek_Amp_seq_515_files)
sample_name_Derek_Amp_seq_515=paste(sample_name_Derek_Amp_seq_515,".515_806",sep="")
# Derek primer 799
Derek_AmpSeq_799_path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/AmpSeq/799_1192/"
list_of_Derek_Amp_seq_799_files=list.files(path = Derek_AmpSeq_799_path,
                                           pattern = ".*uniq_lineage_count.txt")

list_of_Derek_Amp_seq_799_files=paste(Derek_AmpSeq_799_path,list_of_Derek_Amp_seq_799_files,sep="")
sample_name_Derek_Amp_seq_799=basename(list_of_Derek_Amp_seq_799_files)
sample_name_Derek_Amp_seq_799=paste(sample_name_Derek_Amp_seq_799,".799_1192",sep="")

list_of_files=c(list_of_files_unsplitted,list_of_files_splitted,list_of_Derek_Amp_seq_799_files,list_of_Derek_Amp_seq_515_files)
sample_names=c(sample_name_unsplitted,sample_name_splitted,sample_name_Derek_Amp_seq_799,sample_name_Derek_Amp_seq_515)

sample_names=gsub(".usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.UMI_FILTERED.uniq_lineage_count.txt",replacement = "",fixed = T,x = sample_names)
sample_names=gsub(".usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.uniq_lineage_count.txt",replacement = "",fixed = T,x = sample_names)

sample_names=gsub("_Ar_multimap_bacteria_R2",replacement = "",fixed = T,x = sample_names)
sample_names=gsub("191029-069_",replacement = "",fixed = T,x = sample_names)

tmp_df=data.frame(file=list_of_files,sample_name=sample_names)
lineages_files_and_tag=tmp_df

# change the names of the AmpSeq samples
lineages_files_and_tag$sample_name=gsub("F5p515F2F4F6",replacement = "AmpSeq_A",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("C5p515F2F4F6",replacement = "AmpSeq_B",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("A2p515F2F4F6",replacement = "AmpSeq_C",fixed = T,x = lineages_files_and_tag$sample_name)
lineages_files_and_tag$sample_name=gsub("E3p515F1F3F5",replacement = "AmpSeq_D",fixed = T,x = lineages_files_and_tag$sample_name)

# remove the bad splitted ones
lineages_files_and_tag=lineages_files_and_tag[!grepl("1.PROB_479",lineages_files_and_tag$sample_name),]
lineages_files_and_tag=lineages_files_and_tag[!grepl("1.PROB_1265",lineages_files_and_tag$sample_name),]

lineages_files_and_tag_x2=lineages_files_and_tag[!grepl(lineages_files_and_tag$sample_name,pattern = "([ABCD]1)",perl=T),]
# Should we exclude another one --> YES
# AmpSeq_C.799_1192 603 reads
lineages_files_and_tag_x2=lineages_files_and_tag_x2[!grepl(lineages_files_and_tag_x2$sample_name,pattern = "AmpSeq_C.799_1192",fixed = T),]

# plot definitions
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
P40 <- sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
swatch(P40)
names(P40) <- NULL

col1<- colorRampPalette(c("white", "blue"))(50)
breaks1 = c(seq(0,0.5,length=51))

# out files
XLSX_matrix_file=paste(out_dir,"datasets_tables_Bacteria_AmpSeq_vs_Array_4probes_",DB_type,".UMI_FILTERED.xlsx",sep="")
XSLX_Obj=createWorkbook(XLSX_matrix_file)
pdf(file=paste(out_dir,"AmpSeq_Array4probes_Bacteria_all_samples_MMSEQS2_NT_upSetR.UMI_FILTERED.pdf",sep=""),height = 3,width = 6.5) # paper = "a4r", height = 8.3 , width = 11.7

Level="genus" # for simplicity we only analyze genus level, we can do a loop for all... for (Level in c ("genus")) # c ("phylum", "class","order","family","genus") #=c("Count", "superkingdom","phylum", "class","order","family","genus","species")

# read all data and decide the size to sample
all=Tax_profile(lineage_files_and_tag_df=lineages_files_and_tag_x2,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
                levels_to_account_for_exclude = c("genus"))
# we don't use the prob samples so lets remove them, if we want to use them we should sample each time the probs data of each sample and sum them up to give the full sample!
all=all[!grepl(all$sample_name,pattern = "PROB"),]
plants_intersects=data.frame()
for (p in c("A","B","C","D"))
{
  tax_list=list(Array=all$Value[all$sample_name==paste(p,"2",sep="")],
                AmpSeq_p799=all$Value[all$sample_name==paste("AmpSeq_",p,".799_1192",sep="")],
                AmpSeq_p515=all$Value[all$sample_name==paste("AmpSeq_",p,".515_806",sep="")])
  # for the sub samples
  # tax_list=list(Array=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste(p,"2",sep="")]),
  #               AmpSeq_799=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste("AmpSeq_",p,".799_1192",sep="")]),
  #               AmpSeq_515=as.character(all_NO_PROB_sub_samples$tax[all_NO_PROB_sub_samples$sample_name==paste("AmpSeq_",p,".515_806",sep="")]))
  m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
  cs = comb_size(m)
  row_size = set_size(m)
  ht = UpSet(m, 
             top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),
             right_annotation = upset_right_annotation(m, ylim = c(0, 1.1*max(set_size(m)))),
             column_title = paste("Plant",p,Level,"UMI filtered"))
  ht = draw(ht)
  co = column_order(ht)
  row_od = row_order(ht)
  nc = ncol(m)
  decorate_annotation("Intersection\nsize", {
    grid.text(cs[co], 
              x = 1:nc, 
              y = unit(cs[co], "native") + unit(1, "mm"), 
              gp = gpar(fontsize = 10), 
              just = "bottom",
              default.units = "native")
  })
  decorate_annotation("Set size", {
    grid.text(row_size[row_od], 
              unit(row_size[row_od], "native") + unit(1, "mm"), 
              rev(seq_len(length(row_size))), 
              default.units = "native", just = "bottom", rot = -90,
              gp = gpar(fontsize = 10))
  })
  # library(eulerr)
  
  
  intersects=c("Array" = 0, "AmpSeq_p799" = 0, "AmpSeq_p515" = 0,
               "Array&AmpSeq_p799&AmpSeq_p515"=0,  
               "Array&AmpSeq_p799"=0,
               "Array&AmpSeq_p515"=0,
               "AmpSeq_p799&AmpSeq_p515"=0)
  for (i in names(cs)) { # assign
    if (i=="100") {intersects[["Array"]]=cs[[i]]}
    if (i=="010") {intersects[["AmpSeq_p799"]]=cs[[i]]}
    if (i=="001") {intersects[["AmpSeq_p515"]]=cs[[i]]}
    if (i=="111") {intersects[["Array&AmpSeq_p799&AmpSeq_p515"]]=cs[[i]]}
    if (i=="110") {intersects[["Array&AmpSeq_p799"]]=cs[[i]]}
    if (i=="101") {intersects[["Array&AmpSeq_p515"]]=cs[[i]]}
    if (i=="011") {intersects[["AmpSeq_p799&AmpSeq_p515"]]=cs[[i]]}
  }
  intersect_df=as.data.frame(intersects)
  row.names(intersect_df)=paste(p,row.names(intersect_df),sep="_")

  plants_intersects=rbind(plants_intersects,intersect_df)
}  
dev.off()

# Just as QA compare to the non-UMI filtered
path_WO_FILTER="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/multimaped_filtered/MMSEQS2/"
library(openxlsx)
library(Polychrome)
library(gplots)
library(ggplot2)
DB_type="MMSEQS2"

# The data itself
# The full sample
list_of_files_unsplitted_WO_FILTER=list.files(path = path_WO_FILTER,
                                    pattern = ".*_Ar_multimap_bacteria_R2.usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.uniq_lineage_count.txt")
list_of_files_unsplitted_WO_FILTER=paste(path_WO_FILTER,list_of_files_unsplitted_WO_FILTER,sep="")
sample_name_unsplitted_WO_FILTER=basename(list_of_files_unsplitted_WO_FILTER)

# The splitted by prob
list_of_files_splitted_WO_FILTER=list.files(path = path_WO_FILTER,
                                  pattern = ".*_Ar_multimap_bacteria_R2.PROB_.*.usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.uniq_lineage_count.txt")
list_of_files_splitted_WO_FILTER=paste(path_WO_FILTER,list_of_files_splitted_WO_FILTER,sep="")
sample_name_splitted_WO_FILTER=basename(list_of_files_splitted_WO_FILTER)

# Derek primer 515
Derek_AmpSeq_515_path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/AmpSeq/515_806/"
list_of_Derek_Amp_seq_515_files=list.files(path = Derek_AmpSeq_515_path,
                                           pattern = ".*uniq_lineage_count.txt")

list_of_Derek_Amp_seq_515_files=paste(Derek_AmpSeq_515_path,list_of_Derek_Amp_seq_515_files,sep="")
sample_name_Derek_Amp_seq_515=basename(list_of_Derek_Amp_seq_515_files)
sample_name_Derek_Amp_seq_515=paste(sample_name_Derek_Amp_seq_515,".515_806",sep="")
# Derek primer 799
Derek_AmpSeq_799_path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/New_Probs_Dec2020/data/AmpSeq/799_1192/"
list_of_Derek_Amp_seq_799_files=list.files(path = Derek_AmpSeq_799_path,
                                           pattern = ".*uniq_lineage_count.txt")

list_of_Derek_Amp_seq_799_files=paste(Derek_AmpSeq_799_path,list_of_Derek_Amp_seq_799_files,sep="")
sample_name_Derek_Amp_seq_799=basename(list_of_Derek_Amp_seq_799_files)
sample_name_Derek_Amp_seq_799=paste(sample_name_Derek_Amp_seq_799,".799_1192",sep="")

list_of_files_WO_FILTER=c(list_of_files_unsplitted_WO_FILTER,list_of_files_splitted_WO_FILTER,list_of_Derek_Amp_seq_799_files,list_of_Derek_Amp_seq_515_files)
sample_names_WO_FILTER=c(sample_name_unsplitted_WO_FILTER,sample_name_splitted_WO_FILTER,sample_name_Derek_Amp_seq_799,sample_name_Derek_Amp_seq_515)

sample_names_WO_FILTER=gsub(".usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.uniq_lineage_count.txt",replacement = "",fixed = T,x = sample_names_WO_FILTER)
sample_names_WO_FILTER=gsub(".usearch_unique_vs_NT_Jan2021.mmseq2.m8.best_hit_multiple_and_tax_With_tax.LCA.uniq_lineage_count.txt",replacement = "",fixed = T,x = sample_names_WO_FILTER)

sample_names_WO_FILTER=gsub("_Ar_multimap_bacteria_R2",replacement = "",fixed = T,x = sample_names_WO_FILTER)
sample_names_WO_FILTER=gsub("191029-069_",replacement = "",fixed = T,x = sample_names_WO_FILTER)

tmp_df_WO_FILTER=data.frame(file=list_of_files_WO_FILTER,sample_name=sample_names_WO_FILTER)
lineages_files_and_tag_WO_FILTER=tmp_df_WO_FILTER

# change the names of the AmpSeq samples
lineages_files_and_tag_WO_FILTER$sample_name=gsub("F5p515F2F4F6",replacement = "AmpSeq_A",fixed = T,x = lineages_files_and_tag_WO_FILTER$sample_name)
lineages_files_and_tag_WO_FILTER$sample_name=gsub("C5p515F2F4F6",replacement = "AmpSeq_B",fixed = T,x = lineages_files_and_tag_WO_FILTER$sample_name)
lineages_files_and_tag_WO_FILTER$sample_name=gsub("A2p515F2F4F6",replacement = "AmpSeq_C",fixed = T,x = lineages_files_and_tag_WO_FILTER$sample_name)
lineages_files_and_tag_WO_FILTER$sample_name=gsub("E3p515F1F3F5",replacement = "AmpSeq_D",fixed = T,x = lineages_files_and_tag_WO_FILTER$sample_name)

# remove the bad splitted ones
lineages_files_and_tag_WO_FILTER=lineages_files_and_tag_WO_FILTER[!grepl("1.PROB_479",lineages_files_and_tag_WO_FILTER$sample_name),]
lineages_files_and_tag_WO_FILTER=lineages_files_and_tag_WO_FILTER[!grepl("1.PROB_1265",lineages_files_and_tag_WO_FILTER$sample_name),]

lineages_files_and_tag_x2_WO_FILTER=lineages_files_and_tag_WO_FILTER[!grepl(lineages_files_and_tag_WO_FILTER$sample_name,pattern = "([ABCD]1)",perl=T),]
# Should we exclude another one --> YES
# AmpSeq_C.799_1192 603 reads
lineages_files_and_tag_x2_WO_FILTER=lineages_files_and_tag_x2_WO_FILTER[!grepl(lineages_files_and_tag_x2_WO_FILTER$sample_name,pattern = "AmpSeq_C.799_1192",fixed = T),]

# plot definitions
P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
P40 <- sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
swatch(P40)
names(P40) <- NULL

col1<- colorRampPalette(c("white", "blue"))(50)
breaks1 = c(seq(0,0.5,length=51))

Level="genus" # for simplicity we only analyze genus level, we can do a loop for all... for (Level in c ("genus")) # c ("phylum", "class","order","family","genus") #=c("Count", "superkingdom","phylum", "class","order","family","genus","species")

# read all data and decide the size to sample
all_WO_Filter=Tax_profile(lineage_files_and_tag_df=lineages_files_and_tag_x2_WO_FILTER,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
                levels_to_account_for_exclude = c("genus"))
# we don't use the prob samples so lets remove them, if we want to use them we should sample each time the probs data of each sample and sum them up to give the full sample!
all_WO_Filter=all_WO_Filter[!grepl(all_WO_Filter$sample_name,pattern = "PROB"),]

names(all_WO_Filter)[2:3]=paste("WO_UMI_FILTER",names(all_WO_Filter)[2:3],sep=".")

qa_joind=merge(all_WO_Filter,all,by=c("Value","sample_name"),keep.all=T)
View(qa_joind[is.na(qa_joind$Count),])

# count the total number of reads in each sample
sample_count_WO_FILTER=data.frame()
for (sample_i in unique(all_WO_Filter$sample_name)) {
  sample_total=sum(all_WO_Filter$WO_UMI_FILTER.Count[all_WO_Filter$sample_name==sample_i])
  sample_count_WO_FILTER=rbind(sample_count_WO_FILTER,data.frame(sample_name=sample_i,total_size=sample_total))
}

subset_size=min(sample_count$total_size)


#### END


#### TO UPDATE FROM HERE!




### full samples also with the Omni2Probs
XLSX_matrix_file=paste(out_dir,"datasets_tables_Bacteria_AmpSeq_vs_Array_2OR4probes",DB_type,".xlsx",sep="")
XSLX_Obj=createWorkbook(XLSX_matrix_file)


pdf(file=paste(out_dir,"Array_2OR4probes_Bacteria_all_samples_MMSEQS2_NT_stacks_labeled_NMDS_and_heatmaps.pdf",sep=""),height = 16.5,width = 23.4) # paper = "a4r", height = 8.3 , width = 11.7
col1<- colorRampPalette(c("white", "blue"))(50)
breaks1 = c(seq(0,0.5,length=51))
#lineages_files_and_tag_x2=lineages_files_and_tag[!grepl(lineages_files_and_tag$sample_name,pattern = "([ABCD]1)",perl=T),]
for (Level in c ("phylum", "class","order","family","genus")) # =c("Count", "superkingdom","phylum", "class","order","family","genus","species")
{
  all=Tax_profile(lineage_files_and_tag_df=lineages_files_and_tag,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
                  levels_to_account_for_exclude = c("genus"))
  all_with_others=data.frame()
  other_cutoff=0.02
  for (sample in unique (all$sample_name))
  {
    subset=all[all$sample_name==sample,]
    sabset_to_add=subset[subset$Percent>other_cutoff,]
    other_percent=sum(subset$Percent[subset$Percent<=other_cutoff])
    other_count=sum(subset$Count[subset$Percent<=other_cutoff])
    sabset_to_add=rbind(sabset_to_add,data.frame(Value=paste("Other (<=",other_cutoff*100,"%)",sep=""),Count=other_count,Percent=other_percent,sample_name=sample))
    all_with_others=rbind(all_with_others,sabset_to_add)
  }
  # no need for the PROB samples here --> before the change the data frame used was: all_with_others
  all_with_others_NO_PROB=all_with_others[!grepl(all_with_others$sample_name,pattern = "PROB"),]
  all_with_others_NO_PROB$sample_name <- factor(all_with_others_NO_PROB$sample_name,levels = c("A1","A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                                               "B1","B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                                               "C1","C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                                               "D1","D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  p=ggplot(all_with_others_NO_PROB, aes(fill=Value, y=Percent, x=sample_name,label=paste(signif(Percent*100,2),"% ",Value,sep=""))) + 
    geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
    theme(legend.text=element_text(size=15)) +
    scale_fill_manual(values = P40) + 
    ggtitle(paste("Composition",Level)) +
    theme(plot.title = element_text(hjust = 0.5))
  print (p)
  
  # scale_fill_manual(values = as.vector(palette36))
  # scale_fill_brewer( palette = "Paired")
  
  # the heatmap
  percent_matrix=make_table(flat_df = all_with_others_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  # TODO: take into account that it may be that some fall under the other category...
  # we can fix it with the same proceudre we did for missing data (we need to replace the NA with number and reduce it from the Other category)
  percent_matrix_NA_0=percent_matrix
  percent_matrix_NA_0[is.na(percent_matrix_NA_0)]=0
  heatmap.2(x = as.matrix(percent_matrix_NA_0), col = col1, symm = FALSE,breaks=breaks1,trace = "none",margins=c(15,15),cexRow=0.9, main=paste(Level," - NA as 0, DB: ",DB_type, sep="")) # main = paste("Abundance of ",tax_level, " (out of the all reads)",sep="")
  
  sheetName = paste(Level,"_wOt_",other_cutoff,"_percent_NA_0",sep="")
  addWorksheet(XSLX_Obj, sheetName)
  writeDataTable(x=percent_matrix_NA_0, wb = XSLX_Obj,sheet = sheetName,rowNames = TRUE)
  
  # Do NMDS between the full samples profiles - https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/ / https://jkzorz.github.io/2019/06/06/NMDS.html
  library(vegan)
  all_NO_PROB=all[!grepl(all$sample_name,pattern = "PROB"),]
  all_NO_PROB$sample_name <- factor(all_NO_PROB$sample_name,levels = c("A1","A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                       "B1","B2","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                       "C1","C2","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                       "D1","D2","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  percent_matrix_all_samples=make_table(flat_df = all_NO_PROB,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  
  community_matrix=t(as.matrix(percent_matrix_all_samples))
  community_matrix[is.na(community_matrix)]=0 # replace NA with 0
  example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                       k=2,trymax = 100) # The number of reduced dimensions -> ,distance = "bray"
  stressplot(example_NMDS)
  # plot (example_NMDS)
  
  #extract NMDS scores (x and y coordinates)
  data.scores = as.data.frame(scores(example_NMDS))
  data.scores$Sample=row.names(data.scores)
  data.scores$Sample=gsub(pattern=".799_1192",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".515_806",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="AmpSeq_",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="1|2",replacement = "",x=data.scores$Sample)
  
  data.scores$Type=row.names(data.scores) 
  data.scores$Type[grepl("^[ABCD]1$",perl=T,rownames(data.scores))]="Omni2Probes"
  data.scores$Type[grepl("^[ABCD]2$",perl=T,rownames(data.scores))]="Omni4Probes"
  data.scores$Type[grepl("AmpSeq",perl=T,rownames(data.scores))]="AmpSeq"
  
  
  nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
    geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    labs(x = "NMDS1", colour = "", y = "NMDS2", shape = "Type")  + 
    scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    ggtitle(paste("NMDS of full profile in level ",Level," (Stress=",signif(example_NMDS$stress,3),")",sep=""))
  print (nmds_plot)
  
}
dev.off()
saveWorkbook(XSLX_Obj, file = XLSX_matrix_file, overwrite = TRUE)

##### Comparing the probs and AmpSeq
XLSX_matrix_file=paste(out_dir,"datasets_tables_Bacteria_AmpSeq_vs_2OR4probes",DB_type,".xlsx",sep="")
XSLX_Obj=createWorkbook(XLSX_matrix_file)


pdf(file=paste(out_dir,"2OR4probes_Bacteria_all_samples_MMSEQS2_NT_stacks_labeled_NMDS_and_heatmaps.pdf",sep=""),height = 16.5,width = 23.4) # paper = "a4r", height = 8.3 , width = 11.7
col1<- colorRampPalette(c("white", "blue"))(50)
breaks1 = c(seq(0,0.5,length=51))
#lineages_files_and_tag_x2=lineages_files_and_tag[!grepl(lineages_files_and_tag$sample_name,pattern = "([ABCD]1)",perl=T),]
for (Level in c ("phylum", "class","order","family","genus")) # =c("Count", "superkingdom","phylum", "class","order","family","genus","species")
{
  all=Tax_profile(lineage_files_and_tag_df=lineages_files_and_tag,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
                  levels_to_account_for_exclude = c("genus"))
  all_with_others=data.frame()
  other_cutoff=0.035
  for (sample in unique (all$sample_name))
  {
    subset=all[all$sample_name==sample,]
    sabset_to_add=subset[subset$Percent>other_cutoff,]
    other_percent=sum(subset$Percent[subset$Percent<=other_cutoff])
    other_count=sum(subset$Count[subset$Percent<=other_cutoff])
    sabset_to_add=rbind(sabset_to_add,data.frame(Value=paste("Other (<=",other_cutoff*100,"%)",sep=""),Count=other_count,Percent=other_percent,sample_name=sample))
    all_with_others=rbind(all_with_others,sabset_to_add)
  }
  # no need for the PROB samples here --> before the change the data frame used was: all_with_others
  all_with_others_NO_ARRRAY=all_with_others[!grepl(all_with_others$sample_name,pattern = "^[ABCD][12]$"),]
  all_with_others_NO_ARRRAY$sample_name <- factor(all_with_others_NO_ARRRAY$sample_name,levels = c("A1.PROB_799F","A1.PROB_902R","A1.PROB_UNKNOWN","A2.PROB_1265","A2.PROB_479","A2.PROB_799F","A2.PROB_902R","A2.PROB_UNKNOWN","AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                                                   "B1.PROB_799F","B1.PROB_902R","B1.PROB_UNKNOWN","B2.PROB_1265","B2.PROB_479","B2.PROB_799F","B2.PROB_902R","B2.PROB_UNKNOWN","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                                                   "C1.PROB_799F","C1.PROB_902R","C1.PROB_UNKNOWN","C2.PROB_1265","C2.PROB_479","C2.PROB_799F","C2.PROB_902R","C2.PROB_UNKNOWN","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                                                   "D1.PROB_799F","D1.PROB_902R","D1.PROB_UNKNOWN","D2.PROB_1265","D2.PROB_479","D2.PROB_799F","D2.PROB_902R","D2.PROB_UNKNOWN","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  p=ggplot(all_with_others_NO_ARRRAY, aes(fill=Value, y=Percent, x=sample_name,label=paste(signif(Percent*100,2),"% ",Value,sep=""))) + 
    geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
    theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
    theme(legend.text=element_text(size=15)) +
    scale_fill_manual(values = P40) + 
    ggtitle(paste("Composition",Level)) +
    theme(plot.title = element_text(hjust = 0.5))
  print (p)
  
  # scale_fill_manual(values = as.vector(palette36))
  # scale_fill_brewer( palette = "Paired")
  
  # the heatmap
  percent_matrix=make_table(flat_df = all_with_others_NO_ARRRAY,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  # TODO: take into account that it may be that some fall under the other category...
  # we can fix it with the same proceudre we did for missing data (we need to replace the NA with number and reduce it from the Other category)
  percent_matrix_NA_0=percent_matrix
  percent_matrix_NA_0[is.na(percent_matrix_NA_0)]=0
  heatmap.2(x = as.matrix(percent_matrix_NA_0), col = col1, symm = FALSE,breaks=breaks1,trace = "none",margins=c(15,15),cexRow=0.9, main=paste(Level," - NA as 0, DB: ",DB_type, sep="")) # main = paste("Abundance of ",tax_level, " (out of the all reads)",sep="")
  
  sheetName = paste(Level,"_wOt_",other_cutoff,"_percent_NA_0",sep="")
  addWorksheet(XSLX_Obj, sheetName)
  writeDataTable(x=percent_matrix_NA_0, wb = XSLX_Obj,sheet = sheetName,rowNames = TRUE)
  
  # Do NMDS between the full samples profiles - https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/ / https://jkzorz.github.io/2019/06/06/NMDS.html
  library(vegan)
  all_NO_ARRRAY=all[!grepl(all$sample_name,pattern = "^[ABCD][12]$"),]
  all_NO_ARRRAY$sample_name <- factor(all_NO_ARRRAY$sample_name,levels = c("A1.PROB_799F","A1.PROB_902R","A1.PROB_UNKNOWN","A2.PROB_1265","A2.PROB_479","A2.PROB_799F","A2.PROB_902R","A2.PROB_UNKNOWN","AmpSeq_A.515_806", "AmpSeq_A.799_1192",
                                                                           "B1.PROB_799F","B1.PROB_902R","B1.PROB_UNKNOWN","B2.PROB_1265","B2.PROB_479","B2.PROB_799F","B2.PROB_902R","B2.PROB_UNKNOWN","AmpSeq_B.515_806","AmpSeq_B.799_1192",
                                                                           "C1.PROB_799F","C1.PROB_902R","C1.PROB_UNKNOWN","C2.PROB_1265","C2.PROB_479","C2.PROB_799F","C2.PROB_902R","C2.PROB_UNKNOWN","AmpSeq_C.515_806","AmpSeq_C.799_1192",
                                                                           "D1.PROB_799F","D1.PROB_902R","D1.PROB_UNKNOWN","D2.PROB_1265","D2.PROB_479","D2.PROB_799F","D2.PROB_902R","D2.PROB_UNKNOWN","AmpSeq_D.515_806","AmpSeq_D.799_1192"))
  percent_matrix_all_samples=make_table(flat_df = all_NO_ARRRAY,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
  
  community_matrix=t(as.matrix(percent_matrix_all_samples))
  community_matrix[is.na(community_matrix)]=0 # replace NA with 0
  example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                       k=2,trymax = 100) # The number of reduced dimensions -> ,distance = "bray"
  stressplot(example_NMDS)
  # plot (example_NMDS)
  
  #extract NMDS scores (x and y coordinates)
  data.scores = as.data.frame(scores(example_NMDS))
  data.scores$Sample=row.names(data.scores)
  data.scores$Sample=gsub(pattern=".799_1192",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".515_806",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="AmpSeq_",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".PROB_1265",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".PROB_799F",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".PROB_902R",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".PROB_479",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern=".PROB_UNKNOWN",replacement = "",x=data.scores$Sample)
  data.scores$Sample=gsub(pattern="1|2",replacement = "",x=data.scores$Sample)

  data.scores$Type=row.names(data.scores) 
  data.scores$Type[grepl("^[ABCD]1.",perl=T,rownames(data.scores))]="Omni2Probes"
  data.scores$Type[grepl("^[ABCD]2.",perl=T,rownames(data.scores))]="Omni4Probes"
  data.scores$Type[grepl("AmpSeq",perl=T,rownames(data.scores))]="AmpSeq"
  
  
  nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
    geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    labs(x = "NMDS1", colour = "Plant", y = "NMDS2", shape = "Type")  + 
    scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
    ggtitle(paste("NMDS of full profile in level ",Level," (Stress=",signif(example_NMDS$stress,3),")",sep=""))
  print (nmds_plot)
  # ordiplot(example_NMDS)
}
dev.off()
saveWorkbook(XSLX_Obj, file = XLSX_matrix_file, overwrite = TRUE)

# https://mb3is.megx.net/gustame/dissimilarity-based-methods/nmds
# Stress values equal to or below 0.1 are considered fair, while values equal to or below 0.05 indicate good fit.

# do upSetR
library(ComplexHeatmap)
# pdf(file = paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12/data/multimaped_filtered/MMSEQS2/spatial_info/OMNI12.All_Bacterial_ITS_Probs.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.Top50.","corr.upsetR.pdf",sep=""),
#    height=5, width=20)


# upset(fromList(nodes_list),empty.intersections = "on")#, order.by = "freq"
# movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
# upset(as.data.frame(test),empty.intersections = "on")
# with complexHeatMap: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
# make a list with entteties to compare
lineages_files_and_tag_x2_NO_PROBS=lineages_files_and_tag_x2[!grepl(lineages_files_and_tag_x2$sample_name,pattern = "PROB"),]
Level="genus"  
all=Tax_profile(lineage_files_and_tag_df=lineages_files_and_tag_x2_NO_PROBS,tax_level=Level,topX=NA,Filter_of_Interest="Bacteria",ToExclude=c("unclassified","Chloroplast","Mitochondria","uncultured"),
                levels_to_account_for_exclude = c("genus"))

pdf(file = paste(out_dir,"A2.all.upSetR.pdf",sep=""),height=3, width=9)
tax_list=list(Array=all$Value[all$sample_name=="A2"],
              AmpSeq_799=all$Value[all$sample_name=="AmpSeq_A.799_1192"],
              AmpSeq_515=all$Value[all$sample_name=="AmpSeq_A.515_806"])

m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
cs = comb_size(m)
ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),column_title = paste("A2 (Omni4Probs)",Level))
ht = draw(ht)
co = column_order(ht)
nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = 18), 
            just = "bottom",
            default.units = "native")
})
dev.off()


pdf(file = paste(out_dir,"B2.all.upSetR.pdf",sep=""),height=3, width=9)
tax_list=list(Array=all$Value[all$sample_name=="B2"],
              AmpSeq_799=all$Value[all$sample_name=="AmpSeq_B.799_1192"],
              AmpSeq_515=all$Value[all$sample_name=="AmpSeq_B.515_806"])

m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
cs = comb_size(m)
ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),column_title = paste("B2 (Omni4Probs)",Level))
ht = draw(ht)
co = column_order(ht)
nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = 18), 
            just = "bottom",
            default.units = "native")
})
dev.off()

pdf(file = paste(out_dir,"C2.all.upSetR.pdf",sep=""),height=3, width=9)
tax_list=list(Array=all$Value[all$sample_name=="C2"],
              AmpSeq_799=all$Value[all$sample_name=="AmpSeq_C.799_1192"],
              AmpSeq_515=all$Value[all$sample_name=="AmpSeq_C.515_806"])

m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
cs = comb_size(m)
ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),column_title = paste("C2 (Omni4Probs)",Level))
ht = draw(ht)
co = column_order(ht)
nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = 18), 
            just = "bottom",
            default.units = "native")
})
dev.off()

pdf(file = paste(out_dir,"D2.all.upSetR.pdf",sep=""),height=3, width=9)
tax_list=list(Array=all$Value[all$sample_name=="D2"],
              AmpSeq_799=all$Value[all$sample_name=="AmpSeq_D.799_1192"],
              AmpSeq_515=all$Value[all$sample_name=="AmpSeq_D.515_806"])

m = make_comb_mat(tax_list) # https://support.bioconductor.org/p/118557/
cs = comb_size(m)
ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),column_title = paste("D2 (Omni4Probs)",Level))
ht = draw(ht)
co = column_order(ht)
nc = ncol(m)
decorate_annotation("Intersection\nsize", {
  grid.text(cs[co], 
            x = 1:nc, 
            y = unit(cs[co], "native") + unit(1, "mm"), 
            gp = gpar(fontsize = 18), 
            just = "bottom",
            default.units = "native")
})
dev.off()


# Compare the Omni4 with the AmpSeq


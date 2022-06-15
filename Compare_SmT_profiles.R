# load functions
# src_dir="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/Scripts/" # cluster
src_dir="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/"

source(paste (src_dir,"Microbial_Profiles_Functions.R",sep="")) 

all_samples_flat=data.frame()
all_samples_matrix=data.frame()
experiments=c("OMNI12","OMNI13")
# agregate all samples
for (exp in experiments)
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (exp=="OMNI13") {samples=c(samples,"D1")}
  for (s in samples)
  {
    message(paste("==",exp,s))
    data_matrix_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/Bacterial_and_UNKNOWN/",exp,"_",s,".All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos.csv",sep="")
    data=read.delim(file = data_matrix_file,sep = ";",stringsAsFactors = F)
    
    # extract the 'under the tissue'
    # Apr 2022
    # local version - Apr2022
    expression_data_file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/RNA_data/NewestData/220405_semiwild_dataset_rawcounts_filtered.",exp,".",s,".tsv",sep="") # the latest
    # cluster version - Apr2022
    # expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/Apr2022/220405_semiwild_dataset_rawcounts_filtered.",Exp,".",smpl,".tsv",sep="") # the latest
    message(paste("[INFO] Expression data: '",expression_data_file,"'",sep=""))
    expression_data=read.delim(expression_data_file,sep="\t",stringsAsFactors = F)
    expression_data=expression_data[,colSums(expression_data)>0]
    # old
    # expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",exp,".",s,".tsv",sep="")
    # expression_data=read.delim(file = expression_data_file,sep = "\t",stringsAsFactors = F)
    under_the_tissue_spots=names(expression_data)
    message(paste("[INFO] Total 'under the tissue' pixels: ",length(under_the_tissue_spots),sep=""))
    
    data_under_the_tissue=data[,names(data) %in% under_the_tissue_spots]
    row.names(data_under_the_tissue)=data$tax
    message(paste("[INFO] Total 'under the tissue' pixels with Bacterial data: ",NCOL(data_under_the_tissue),sep=""))
    # View(data_under_the_tissue[1:10,1:10])
    sample_profile=as.data.frame(rowSums(data_under_the_tissue))
    sample_name=paste(exp,s,sep=".")
    names(sample_profile)=sample_name
    if (nrow(all_samples_flat)==0) {
      all_samples_flat=data.frame(tax=row.names(sample_profile),count=sample_profile[[1]],percent=sample_profile[[1]]/sum(sample_profile[[1]]),sample=sample_name)
      all_samples_matrix=sample_profile
    } else {
      all_samples_flat=rbind(all_samples_flat,data.frame(tax=row.names(sample_profile),count=sample_profile[[1]],percent=sample_profile[[1]]/sum(sample_profile[[1]]),sample=sample_name))
      all_samples_matrix=merge(all_samples_matrix,sample_profile,by="row.names",all=T)
      row.names(all_samples_matrix)=all_samples_matrix$Row.names
      all_samples_matrix=all_samples_matrix[,-1]
    }
  }
}

# add the samples labels with sections annotations etc
samples_label=data.frame(sample= c("OMNI12.A1","OMNI12.A2","OMNI12.B1","OMNI12.B2","OMNI12.C1","OMNI12.C2",
                                   "OMNI13.A1","OMNI13.A2","OMNI13.B1","OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1"))

P1L1=c("OMNI12.A1","OMNI12.A2","OMNI12.B2")
P1L2=c("OMNI12.B1","OMNI12.C1","OMNI12.C2")
P2L1=c("OMNI13.A1","OMNI13.A2","OMNI13.B1")
P2L2=c("OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")

samples_label$PlantLeaf[samples_label$sample %in% P1L1]="P1.L1"
samples_label$PlantLeaf[samples_label$sample %in% P1L2]="P1.L2"
samples_label$PlantLeaf[samples_label$sample %in% P2L1]="P2.L1"
samples_label$PlantLeaf[samples_label$sample %in% P2L2]="P2.L2"

section1=c("OMNI12.A1","OMNI12.B1","OMNI13.A1","OMNI13.D1")
section2=c("OMNI12.A2","OMNI12.C1","OMNI13.C1")
section3=c("OMNI13.A2","OMNI13.B2")
section4=c("OMNI12.B2","OMNI12.C2","OMNI13.B1","OMNI13.C2")

samples_label$Section[samples_label$sample %in% section1]="1"
samples_label$Section[samples_label$sample %in% section2]="2"
samples_label$Section[samples_label$sample %in% section3]="3"
samples_label$Section[samples_label$sample %in% section4]="4"

samples_label$lable=paste(samples_label$PlantLeaf,samples_label$Section,sep=".")

samples_label_line=as.data.frame(t(samples_label)[c(1,4),])
names(samples_label_line)=samples_label_line[1,]
samples_label_line=samples_label_line[-1,]
# save all the samples matrix
# OMNI12_OMNI13_Bacterial_profiles_file="~/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacterial_and_UNKNOWN.UMI_filtered.Under_the_tissue.genus_profiles.counts.tsv"
OMNI12_OMNI13_Bacterial_profiles_file="~/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacterial_and_UNKNOWN.UMI_filtered.Under_the_tissue_Apr2022_220405.genus_profiles.counts.tsv"

all_samples_matrix_for_file=data.frame(genus=row.names(all_samples_matrix),all_samples_matrix)
all_samples_matrix_for_file=rbind(all_samples_matrix_for_file,data.frame(genus="sample_name",samples_label_line)) # add the section name we used in the paper

write.table(file=OMNI12_OMNI13_Bacterial_profiles_file,x = all_samples_matrix_for_file,quote = F,sep = "\t",row.names = F)
rm(all_samples_matrix_for_file)

# Percent matrix
samples_sum=colSums(all_samples_matrix,na.rm = T)
all_samples_matrix_p=all_samples_matrix
for (i_sample in 1:ncol(all_samples_matrix_p))
{
  all_samples_matrix_p[[i_sample]]=all_samples_matrix_p[[i_sample]]/samples_sum[i_sample]
}
all_samples_matrix_p_for_file=data.frame(genus=row.names(all_samples_matrix_p),all_samples_matrix_p)
all_samples_matrix_p_for_file=rbind(all_samples_matrix_p_for_file,data.frame(genus="sample_name",samples_label_line)) # add the section name we used in the paper

# OMNI12_OMNI13_Bacterial_Pprofiles_file="~/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacterial_and_UNKNOWN.UMI_filtered.Under_the_tissue.genus_profiles.percents.tsv"
OMNI12_OMNI13_Bacterial_Pprofiles_file="~/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacterial_and_UNKNOWN.UMI_filtered.Under_the_tissue_Apr2022_220405.genus_profiles.percents.tsv"


write.table(file=OMNI12_OMNI13_Bacterial_Pprofiles_file,x = all_samples_matrix_p_for_file,quote = F,sep = "\t",row.names = F)
rm(all_samples_matrix_p_for_file)



# do the comparisions between samples
all_with_others=data.frame()
other_cutoff=0.01
for (sample in unique (all_samples_flat$sample))
{
  subset=all_samples_flat[all_samples_flat$sample==sample,]
  sabset_to_add=subset[subset$percent>other_cutoff,]
  other_percent=sum(subset$percent[subset$percent<=other_cutoff])
  other_count=sum(subset$count[subset$percent<=other_cutoff])
  sabset_to_add=rbind(sabset_to_add,data.frame(tax=paste("Other (<=",other_cutoff*100,"%)",sep=""),count=other_count,percent=other_percent,sample=sample))
  all_with_others=rbind(all_with_others,sabset_to_add)
}

# improve all with others to not exclude taxa that pass the cutoff from other samples
all_with_others_relaxed=data.frame()
all_tax_considered=unique(all_with_others$tax)
all_tax_considered_without_other_bin=all_tax_considered[!grepl(x=all_tax_considered,pattern = "Other")]
Other_bin_name=all_tax_considered[grepl(x=all_tax_considered,pattern = "Other")]

for (sample in unique (all_samples_flat$sample))
{
  subset=all_with_others[all_with_others$sample==sample&all_with_others$tax!=Other_bin_name,]
  Other_bin_line=all_with_others[all_with_others$sample==sample&all_with_others$tax==Other_bin_name,]
  for (tax in all_tax_considered_without_other_bin) # see if any tax included in the final set need to be added and correct the Other bin
  {
      if (sum(subset$tax==tax)==0) # was below the other cutoff and not included
      {
        if (nrow(all_samples_flat[all_samples_flat$sample==sample&all_samples_flat$tax==tax,])>0) { # exists
          low_abundance_line=all_samples_flat[all_samples_flat$sample==sample&all_samples_flat$tax==tax,]
          
          # substract from the Other bin
          Other_bin_line$count=Other_bin_line$count-low_abundance_line$count
          Other_bin_line$percent=Other_bin_line$percent-low_abundance_line$percent
          
          # add the low abundance line
          all_with_others_relaxed=rbind(all_with_others_relaxed,low_abundance_line)
          rm (low_abundance_line)
        } 
      }
  }
  # add all the other tax and the updated Other bin
  all_with_others_relaxed=rbind(all_with_others_relaxed,subset,Other_bin_line)
  rm (subset,Other_bin_line,tax)
}  

all_with_others$name <- factor(all_with_others$sample,levels = c("OMNI12.A1","OMNI12.A2","OMNI12.B2",
                                                                 "OMNI12.B1","OMNI12.C1","OMNI12.C2",
                                                                 "OMNI13.A1","OMNI13.A2","OMNI13.B1",
                                                                 "OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")
)
all_with_others_relaxed$name <- factor(all_with_others_relaxed$sample,levels = c("OMNI12.A1","OMNI12.A2","OMNI12.B2",
                                                                                 "OMNI12.B1","OMNI12.C1","OMNI12.C2",
                                                                                 "OMNI13.A1","OMNI13.A2","OMNI13.B1",
                                                                                 "OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")
)

# stack plots
# plot definitions
# Genus levels stackplots and heatmaps
# Another option: https://cran.r-project.org/web/packages/Polychrome/vignettes/testgg.html
library(Polychrome)
library(ggplot2)

P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
P40 <- sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
swatch(P40)
names(P40) <- NULL

# The strict version
all_with_others=merge(x=all_with_others,y=samples_label,by="sample")

## pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Bacteria_Genus.profiles.pdf",width = 28,height = 14)
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.pdf",width = 28,height = 14) # did it because was not sure the above is the UMI filterd...
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.pdf",width = 28,height = 14) # To make sure identical with previous and make sure all in one place with the data-files

# Apr2022 under the tissue version
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.pdf",width = 28,height = 14) # To make sure identical with previous and make sure all in one place with the data-files

p=ggplot(all_with_others, aes(fill=tax, y=percent, x=lable,label=paste(signif(percent*100,2),"% ",tax,sep=""))) + 
  geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
  theme(legend.text=element_text(size=15)) +
  scale_fill_manual(values = P40) + 
  # ggtitle(paste("Subsample",i_repeat,"composition",Level)) +
  theme(plot.title = element_text(hjust = 0.5))
print (p)                                      
dev.off()

# save the flat data
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.data.tsv",
#             x=all_with_others,quote = F,sep = "\t",row.names = F)

# Apr2022
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.data.tsv",
            x=all_with_others,quote = F,sep = "\t",row.names = F)



# The relaxed version
all_with_others_relaxed=merge(x=all_with_others_relaxed,y=samples_label,by="sample")

# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.relaxed.pdf",width = 28,height = 14) # To make sure identical with previous and make sure all in one place with the data-files
# Apr2022 under the tissue
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.pdf",width = 28,height = 14) # To make sure identical with previous and make sure all in one place with the data-files
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.25x14.pdf",width = 25,height = 14) # To make sure identical with previous and make sure all in one place with the data-files

p=ggplot(all_with_others_relaxed, aes(fill=tax, y=percent, x=lable,label=paste(signif(percent*100,2),"% ",tax,sep=""))) + 
  geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
  theme(legend.text=element_text(size=15)) +
  scale_fill_manual(values = P40) + 
  # ggtitle(paste("Subsample",i_repeat,"composition",Level)) +
  theme(plot.title = element_text(hjust = 0.5))
print (p)                                      
dev.off()

# save the flat data
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.relaxed.data.tsv",
#             x=all_with_others_relaxed,quote = F,sep = "\t",row.names = F)

# Apr2022 under the tissue
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.data.tsv",
               x=all_with_others_relaxed,quote = F,sep = "\t",row.names = F)


# Just as QA a file with both relaxed an non relaxed
tmp_both=merge(all_with_others_relaxed,all_with_others,by=c("sample","tax"),all=T)
names(tmp_both)[3:14]=gsub(x=names(tmp_both)[3:14],pattern = ".x$",replacement = ".relaxed",perl = T)
names(tmp_both)[3:14]=gsub(x=names(tmp_both)[3:14],pattern = ".y$",replacement = ".strict",perl = T)
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.relaxed_and_strict.data.tsv",
#             x=tmp_both,quote = F,sep = "\t",row.names = F)

# Apr2022 under the tissue
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed_and_strict.data.tsv",
               x=tmp_both,quote = F,sep = "\t",row.names = F)

rm(tmp_both)
# some summary of average profile

average_profile_df=data.frame()
for (g in unique(all_with_others_relaxed$tax))
{
  r=data.frame(tax=g,
               mean=signif(x = mean(all_with_others_relaxed$percent[all_with_others_relaxed$tax==g]),3),
               sd=signif(x = sd(all_with_others_relaxed$percent[all_with_others_relaxed$tax==g]),3),
               var=signif(x = var(all_with_others_relaxed$percent[all_with_others_relaxed$tax==g]),3),
               median=signif(x=median(all_with_others_relaxed$percent[all_with_others_relaxed$tax==g]),3),
               n=sum(all_with_others_relaxed$tax==g)
               )
  average_profile_df=rbind(average_profile_df,r)
}
average_profile_df[order(average_profile_df$median,decreasing = T),]

# save the average table
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.relaxed_average_profile.tsv",
#             x=average_profile_df,quote = F,sep = "\t",row.names = F)
# Apr2022 under the tissue
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed_average_profile.tsv",
               x=average_profile_df,quote = F,sep = "\t",row.names = F)

# Back to table across samples
all_with_others_relaxed_table=make_table(flat_df = all_with_others_relaxed[,c("tax","count","percent","lable")],column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.relaxed.samples_profile.tsv",
#            x=all_with_others_relaxed_table,quote = F,sep = "\t",row.names = T)

# Apr2022 under the tissue
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.samples_profile.tsv",
            x=all_with_others_relaxed_table,quote = F,sep = "\t",row.names = T)

# NMDS
library(vegan)

# percent_matrix_all_samples=make_table(flat_df = all_with_others,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)

percent_matrix_all_samples=make_table(flat_df = all_samples_flat,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
count_matrix_all_samples=make_table(flat_df = all_samples_flat,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)

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

P1L1=c("OMNI12.A1","OMNI12.A2","OMNI12.B2")
P1L2=c("OMNI12.B1","OMNI12.C1","OMNI12.C2")
P2L1=c("OMNI13.A1","OMNI13.A2","OMNI13.B1")
P2L2=c("OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")

data.scores$PlantLeaf=row.names(data.scores) 
data.scores$PlantLeaf[row.names(data.scores) %in% P1L1]="P1.L1"
data.scores$PlantLeaf[row.names(data.scores) %in% P1L2]="P1.L2"
data.scores$PlantLeaf[row.names(data.scores) %in% P2L1]="P2.L1"
data.scores$PlantLeaf[row.names(data.scores) %in% P2L2]="P2.L2"

section1=c("OMNI12.A1","OMNI12.B1","OMNI13.A1","OMNI13.D1")
section2=c("OMNI12.A2","OMNI12.C1","OMNI13.C1")
section3=c("OMNI13.A2","OMNI13.B2")
section4=c("OMNI12.B2","OMNI12.C2","OMNI13.B1","OMNI13.C2")

data.scores$Section=row.names(data.scores) 
data.scores$Section[row.names(data.scores) %in% section1]="1"
data.scores$Section[row.names(data.scores) %in% section2]="2"
data.scores$Section[row.names(data.scores) %in% section3]="3"
data.scores$Section[row.names(data.scores) %in% section4]="4"

nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
  geom_point(size = 9, aes( shape = PlantLeaf, fill = Section,color=Section))+ 
  # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
  geom_text(aes(label=paste(PlantLeaf,Section,sep=".")),hjust=-0.1, vjust=-0.1) +
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
  labs(x = "NMDS1", shape = "Plant.Leaf",fill = "Section", y = "NMDS2",colour = "Section")  + #  shape = "Sample",
  scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
  scale_shape_manual(values=c(21, 22, 23,24))+ 
  scale_fill_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7")) +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
  ggtitle(paste("NMDS of the full Bacterial profiles on the Genus level (Stress=",signif(example_NMDS$stress,3),")",sep=""))
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13.Bacterial_genus.NMDS.pdf",width = 7,height = 7)
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.NMDS.pdf",width = 7,height = 7) # Just to make more organized
# Apr2022 under the tissue
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.under_the_tissue_Apr2022_220405.NMDS.pdf",width = 7,height = 7) # Just to make more organized

print (nmds_plot,width = 4, height = 4)
dev.off()

# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.BaryCurtis_Dist.pdf",width = 7,height = 7) # Just to make more organized
# Apr2022 under the tissue
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.under_the_tissue_Apr2022_220405.BaryCurtis_Dist.pdf",width = 7,height = 7) # Just to make more organized

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
pHeatmapSim=Heatmap(dist_mat,col=col_fun,name="Distance", #column_title = paste("Subsample ",i_repeat,"- Bray-Curtis dissimilarity",sep="")
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.2f", dist_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                    })
print(pHeatmapSim) 
dev.off()

# save the tables
###
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.NMDS_community_matrix.tsv",
#             x=community_matrix,
#             quote = F,sep = "\t",row.names = T
#             )
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.NMDS_scores.tsv",
#             x=data.scores,
#             quote = F,sep = "\t",row.names = T
# )
# saveRDS(example_NMDS, file = "/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.NMDS_obj.rds")

# Apr2022 under the tissue version
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.under_the_tissue_Apr2022_220405.NMDS_community_matrix.tsv",
            x=community_matrix,
            quote = F,sep = "\t",row.names = T
)
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.under_the_tissue_Apr2022_220405.NMDS_scores.tsv",
            x=data.scores,
            quote = F,sep = "\t",row.names = T
)
saveRDS(example_NMDS, file = "/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Bacterial_genus.under_the_tissue_Apr2022_220405.NMDS_obj.rds")


rm("community_matrix" , "data", "data_matrix_file","data_under_the_tissue","data.scores","example_NMDS","exp",
  "experiments","expression_data","expression_data_file","nmds_plot","other_count","other_cutoff","other_percent","p",                        
  "s","sabset_to_add","sample","sample_name","sample_profile","samples","under_the_tissue_spots","subset")

# Only the Top20, including those missing without 'Others' category
TopX=20
all_samples_flat_TopX=data.frame()
for (sample in unique (all_samples_flat$sample))
{
  subset=all_samples_flat[all_samples_flat$sample==sample,]
  subset=subset[order(subset$count,decreasing = T),]
  sabset_to_add=subset[1:TopX,]
  all_samples_flat_TopX=rbind(all_samples_flat_TopX,sabset_to_add)
}
count_matrix_all_samples_TopX=make_table(flat_df = all_samples_flat_TopX,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)
message(paste("[INFO] Before filling, still have total of",sum(is.na(count_matrix_all_samples_TopX)),"missing values"))

# fill in the missing
for (i in 1:nrow(count_matrix_all_samples_TopX)) {
  for (j in 1:ncol(count_matrix_all_samples_TopX)) {
    if (is.na(count_matrix_all_samples_TopX[i,j])) # missing - i.e. this species is not on the TopX of that sample
    {
      sample_name=names(count_matrix_all_samples_TopX)[j]
      tax=row.names(count_matrix_all_samples_TopX)[i]
      # check if it is found out of the TopX
      if (length(all_samples_flat$count[all_samples_flat$sample==sample_name&
                                  all_samples_flat$tax==tax])>0)
      {
        if (length(all_samples_flat$count[all_samples_flat$sample==sample_name&
                                          all_samples_flat$tax==tax])==1)
        {
          count_matrix_all_samples_TopX[i,j]=all_samples_flat$count[all_samples_flat$sample==sample_name&all_samples_flat$tax==tax]
          message(paste("[INFO] found a filling for tax",tax,"in sample",sample_name,":",all_samples_flat$count[all_samples_flat$sample==sample_name&all_samples_flat$tax==tax]))
        } else {
          message (paste ("[WARNING] more than one value found for",sample_name, "and tax",tax,"- Don't take the fill!!"))
        }
      }
    }
  }
}

message(paste("[INFO] After filling, still have total of",sum(is.na(count_matrix_all_samples_TopX)),"missing values"))

P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
P40 <- sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
swatch(P40)
names(P40) <- NULL

# back to flat db
all_samples_TopX_count=data.frame()
all_samples_TopX=data.frame()

for (i in 1:nrow(count_matrix_all_samples_TopX)) {
  for (j in 1:ncol(count_matrix_all_samples_TopX)) {
    if (!is.na(count_matrix_all_samples_TopX[i,j])) {
      sample_name=names(count_matrix_all_samples_TopX)[j]
      tax=row.names(count_matrix_all_samples_TopX)[i]
      record=data.frame(tax=tax,
                        count=count_matrix_all_samples_TopX[i,j],
                        sample=sample_name)
        
      all_samples_TopX_count=rbind(all_samples_TopX_count,record) 
    }
  }
}

for (sample in unique (all_samples_TopX_count$sample))
{
  subset=all_samples_TopX_count[all_samples_TopX_count$sample==sample,]
  subset$percent=subset$count/sum(subset$count)
  all_samples_TopX=rbind(all_samples_TopX,subset)
}

#pdf(file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Top",TopX,"_Bacterial_Genus.profiles.pdf",sep=""),width = 28,height = 14)
# Apr2022 under the tissue version
pdf(file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Top",TopX,"_Bacterial_Genus.under_the_tissue_Apr2022_220405.profiles.pdf",sep=""),width = 28,height = 14)
p=ggplot(all_samples_TopX, aes(fill=tax, y=percent, x=sample,label=paste(signif(percent*100,2),"% ",tax,sep=""))) + 
  geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
  theme(legend.text=element_text(size=15)) +
  scale_fill_manual(values = P40) + 
  # ggtitle(paste("Subsample",i_repeat,"composition",Level)) +
  theme(plot.title = element_text(hjust = 0.5))
print (p)                                      
dev.off()

# Convert it to percentage matrix
percent_matrix_all_samples_TopX=data.frame(row.names = row.names(count_matrix_all_samples_TopX))
# data.frame(matrix(NA, nrow = nrow(count_matrix_all_samples_TopX), ncol = ncol(count_matrix_all_samples_TopX)))
for (i in 1:ncol(count_matrix_all_samples_TopX))
{
  sample_name=names(count_matrix_all_samples_TopX)[i]
  percent_matrix_all_samples_TopX[[sample_name]]=count_matrix_all_samples_TopX[[sample_name]]/sum(count_matrix_all_samples_TopX[[sample_name]])
}


# Do The NMDS
library(vegan)
community_matrix=t(as.matrix(percent_matrix_all_samples_TopX))
#community_matrix=t(as.matrix(percent_matrix_all_sub_samples))
community_matrix[is.na(community_matrix)]=0 # replace NA with 0
example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                     k=2,trymax = 100,distance = "bray") # The number of reduced dimensions -> ,, autotransform = FALSE,distance ="euclidean"
stressplot(example_NMDS)
# plot (example_NMDS)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(example_NMDS))
data.scores$Sample=row.names(data.scores)

P1L1=c("OMNI12.A1","OMNI12.A2","OMNI12.B2")
P1L2=c("OMNI12.B1","OMNI12.C1","OMNI12.C2")
P2L1=c("OMNI13.A1","OMNI13.A2","OMNI13.B1")
P2L2=c("OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")

data.scores$PlantLeaf=row.names(data.scores) 
data.scores$PlantLeaf[row.names(data.scores) %in% P1L1]="P1.L1"
data.scores$PlantLeaf[row.names(data.scores) %in% P1L2]="P1.L2"
data.scores$PlantLeaf[row.names(data.scores) %in% P2L1]="P2.L1"
data.scores$PlantLeaf[row.names(data.scores) %in% P2L2]="P2.L2"

section1=c("OMNI12.A1","OMNI12.B1","OMNI13.A1","OMNI13.D1")
section2=c("OMNI12.A2","OMNI12.C1","OMNI13.C1")
section3=c("OMNI13.A2","OMNI13.B2")
section4=c("OMNI12.B2","OMNI12.C2","OMNI13.B1","OMNI13.C2")

data.scores$Section=row.names(data.scores) 
data.scores$Section[row.names(data.scores) %in% section1]="1"
data.scores$Section[row.names(data.scores) %in% section2]="2"
data.scores$Section[row.names(data.scores) %in% section3]="3"
data.scores$Section[row.names(data.scores) %in% section4]="4"

nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
  geom_point(size = 9, aes( shape = PlantLeaf, fill = Section,color=Section))+ 
  # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
  geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
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
  labs(x = "NMDS1", shape =  "Plant.Leaf" ,fill = "Section", y = "NMDS2",colour = "Section")  + #  shape = "Sample",
  scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
  scale_shape_manual(values=c(21, 22, 23,24))+ 
  scale_fill_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7")) +   # https://jrnold.github.io/ggthemes/reference/colorblind.html
  ggtitle(paste("NMDS of the Top ",TopX," Bacterial genera (Stress=",signif(example_NMDS$stress,3),")",sep=""))

# pdf(file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Top",TopX,"_Bacterial_genus.NMDS.pdf",sep=""),width = 7,height = 7)

# Apr2022 uner the tissue version
pdf(file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Top",TopX,"_Bacterial_genus.under_the_tissue_Apr2022_220405.NMDS.pdf",sep=""),width = 7,height = 7)
print (nmds_plot,width = 4, height = 4)
dev.off()


# hirarchial clustering and distances heatmaps
# calculate Bray-Curtis distance among samples
comm.bc.dist <- vegdist(community_matrix, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
# plot cluster diagram
library(dendextend)

# order.dendrogram(as.dendrogram(comm.bc.clust))
# comm.bc.clust$labels
leafs_colors=c(1,1,2,1,2,2,3,3,3,4,4,4,4)
dend=as.dendrogram(comm.bc.clust)
## [1] 1 2 3 1 2 3
# But sort them based on their order in dend:
colors_to_use <- leafs_colors[order.dendrogram(dend)]
colors_to_use
## [1] 1 1 2 2 3 3
# Now we can use them
labels_colors(dend) <- colors_to_use
# Now each state has a color
labels_colors(dend) 
pdf(file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Top",TopX,"_Bacterial_genus.under_the_tissue_Apr2022_220405.Bray_Curtis_sim.pdf",sep=""),width = 7,height = 7)
plot(dend, main = "Bray-Curtis dissimilarity")

# plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")
labels_colors(comm.bc.clust) <- 1:5
# Now each state has a color
labels_colors(dend) 

# profiles with the dendograms using ggtree
# https://yulab-smu.top/treedata-book/chapter9.html?q=dend#dendrogram
# https://yulab-smu.top/treedata-book/chapter7.html


library(circlize)
col_fun = colorRamp2(c(0, 1), c("red","white"))
col_fun(seq(0, 1))

library(ComplexHeatmap)
dist_mat=as.matrix(comm.bc.dist)
pHeatmapSim=Heatmap(dist_mat,col=col_fun,name="Distance", #column_title = paste("Subsample ",i_repeat,"- Bray-Curtis dissimilarity",sep="")
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.2f", dist_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                    })
print(pHeatmapSim) 

dev.off()

### FUNGI
# Start all clean -> rm()

# load functions
# src_dir="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/Scripts/" # cluster
src_dir="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/"

source(paste (src_dir,"Microbial_Profiles_Functions.R",sep="")) 

all_samples_flat_Fungi=data.frame()
all_samples_matrix_Fungi=data.frame()
experiments=c("OMNI12","OMNI13")
# agregate all samples
for (exp in experiments)
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (exp=="OMNI13") {samples=c(samples,"D1")}
  for (s in samples)
  {
    message(paste("--",exp,s))
    data_matrix_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/Fungi_and_UNKNOWN/",exp,"_",s,".All_ITS_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos.csv",sep="")
    data=read.delim(file = data_matrix_file,sep = ";",stringsAsFactors = F)
    
    
    # extract the 'under the tissue'
    # Apr 2022
    # local version - Apr2022
    expression_data_file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/RNA_data/NewestData/220405_semiwild_dataset_rawcounts_filtered.",exp,".",s,".tsv",sep="") # the latest
    # cluster version - Apr2022
    # expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/Apr2022/220405_semiwild_dataset_rawcounts_filtered.",Exp,".",smpl,".tsv",sep="") # the latest
    message(paste("[INFO] Expression data: '",expression_data_file,"'",sep=""))
    expression_data=read.delim(expression_data_file,sep="\t",stringsAsFactors = F)
    expression_data=expression_data[,colSums(expression_data)>0]
    # old
    # expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",exp,".",s,".tsv",sep="")
    # expression_data=read.delim(file = expression_data_file,sep = "\t",stringsAsFactors = F)
    under_the_tissue_spots=names(expression_data)
    message(paste("[INFO] Total 'under the tissue' pixels: ",length(under_the_tissue_spots),sep=""))
    
    data_under_the_tissue=data[,names(data) %in% under_the_tissue_spots]
    row.names(data_under_the_tissue)=data$tax
    message(paste("[INFO] Total 'under the tissue' pixels with Fungi data: ",NCOL(data_under_the_tissue),sep=""))
    
    ## old
    # # extract the 'under the tissue'
    # expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",exp,".",s,".tsv",sep="")
    # expression_data=read.delim(file = expression_data_file,sep = "\t",stringsAsFactors = F)
    # under_the_tissue_spots=names(expression_data)
    
    # data_under_the_tissue=data[,names(data) %in% under_the_tissue_spots]
    # row.names(data_under_the_tissue)=data$tax
    # View(data_under_the_tissue[1:10,1:10])
    
    sample_profile=as.data.frame(rowSums(data_under_the_tissue))
    sample_name=paste(exp,s,sep=".")
    names(sample_profile)=sample_name
    # if (nrow(all_samples_flat)==0) { ### IS IT A BUG? Should it be: all_samples_flat_Fungi
    if (nrow(all_samples_flat_Fungi)==0) {
      all_samples_flat_Fungi=data.frame(tax=row.names(sample_profile),count=sample_profile[[1]],percent=sample_profile[[1]]/sum(sample_profile[[1]]),sample=sample_name)
      all_samples_matrix_Fungi=sample_profile
    } else {
      all_samples_flat_Fungi=rbind(all_samples_flat_Fungi,data.frame(tax=row.names(sample_profile),count=sample_profile[[1]],percent=sample_profile[[1]]/sum(sample_profile[[1]]),sample=sample_name))
      all_samples_matrix_Fungi=merge(all_samples_matrix_Fungi,sample_profile,by="row.names",all=T)
      row.names(all_samples_matrix_Fungi)=all_samples_matrix_Fungi$Row.names
      all_samples_matrix_Fungi=all_samples_matrix_Fungi[,-1]
    }
  }
}

# add the samples labels with sections annotations etc
samples_label=data.frame(sample= c("OMNI12.A1","OMNI12.A2","OMNI12.B1","OMNI12.B2","OMNI12.C1","OMNI12.C2",
                                   "OMNI13.A1","OMNI13.A2","OMNI13.B1","OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1"))

P1L1=c("OMNI12.A1","OMNI12.A2","OMNI12.B2")
P1L2=c("OMNI12.B1","OMNI12.C1","OMNI12.C2")
P2L1=c("OMNI13.A1","OMNI13.A2","OMNI13.B1")
P2L2=c("OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")

samples_label$PlantLeaf[samples_label$sample %in% P1L1]="P1.L1"
samples_label$PlantLeaf[samples_label$sample %in% P1L2]="P1.L2"
samples_label$PlantLeaf[samples_label$sample %in% P2L1]="P2.L1"
samples_label$PlantLeaf[samples_label$sample %in% P2L2]="P2.L2"

section1=c("OMNI12.A1","OMNI12.B1","OMNI13.A1","OMNI13.D1")
section2=c("OMNI12.A2","OMNI12.C1","OMNI13.C1")
section3=c("OMNI13.A2","OMNI13.B2")
section4=c("OMNI12.B2","OMNI12.C2","OMNI13.B1","OMNI13.C2")

samples_label$Section[samples_label$sample %in% section1]="1"
samples_label$Section[samples_label$sample %in% section2]="2"
samples_label$Section[samples_label$sample %in% section3]="3"
samples_label$Section[samples_label$sample %in% section4]="4"

samples_label$lable=paste(samples_label$PlantLeaf,samples_label$Section,sep=".")

samples_label_line=as.data.frame(t(samples_label)[c(1,4),])
names(samples_label_line)=samples_label_line[1,]
samples_label_line=samples_label_line[-1,]
# save all the samples matrix
# OMNI12_OMNI13_Fungi_profiles_file="~/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_and_UNKNOWN.UMI_filtered.Under_the_tissue.genus_profiles.counts.tsv"

# Apr2022 under the tissue
OMNI12_OMNI13_Fungi_profiles_file="~/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_and_UNKNOWN.UMI_filtered.under_the_tissue_Apr2022_220405.genus_profiles.counts.tsv"

all_samples_matrix_for_file=data.frame(genus=row.names(all_samples_matrix_Fungi),all_samples_matrix_Fungi)
all_samples_matrix_for_file=rbind(all_samples_matrix_for_file,data.frame(genus="sample_name",samples_label_line)) # add the section name we used in the paper

write.table(file=OMNI12_OMNI13_Fungi_profiles_file,x = all_samples_matrix_for_file,quote = F,sep = "\t",row.names = F)
rm(all_samples_matrix_for_file)

# Percent matrix
samples_sum=colSums(all_samples_matrix_Fungi,na.rm = T)
all_samples_matrix_Fungi_p=all_samples_matrix_Fungi
for (i_sample in 1:ncol(all_samples_matrix_Fungi_p))
{
  all_samples_matrix_Fungi_p[[i_sample]]=all_samples_matrix_Fungi_p[[i_sample]]/samples_sum[i_sample]
}
all_samples_matrix_p_for_file=data.frame(genus=row.names(all_samples_matrix_Fungi_p),all_samples_matrix_Fungi_p)
all_samples_matrix_p_for_file=rbind(all_samples_matrix_p_for_file,data.frame(genus="sample_name",samples_label_line)) # add the section name we used in the paper

# old
# OMNI12_OMNI13_Fungi_Pprofiles_file="~/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_and_UNKNOWN.UMI_filtered.Under_the_tissue.genus_profiles.percents.tsv"
# Apr2022 under the tissue
OMNI12_OMNI13_Fungi_Pprofiles_file="~/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_and_UNKNOWN.UMI_filtered.under_the_tissue_Apr2022_220405.genus_profiles.percents.tsv"


write.table(file=OMNI12_OMNI13_Fungi_Pprofiles_file,x = all_samples_matrix_p_for_file,quote = F,sep = "\t",row.names = F)
rm(all_samples_matrix_p_for_file)


# do the comparisions between samples
all_with_others_Fungi=data.frame()
other_cutoff=0.01
for (sample in unique (all_samples_flat_Fungi$sample))
{
  subset=all_samples_flat_Fungi[all_samples_flat_Fungi$sample==sample,]
  sabset_to_add=subset[subset$percent>other_cutoff,]
  other_percent=sum(subset$percent[subset$percent<=other_cutoff])
  other_count=sum(subset$count[subset$percent<=other_cutoff])
  sabset_to_add=rbind(sabset_to_add,data.frame(tax=paste("Other (<=",other_cutoff*100,"%)",sep=""),count=other_count,percent=other_percent,sample=sample))
  all_with_others_Fungi=rbind(all_with_others_Fungi,sabset_to_add)
}

# improve all with others to not exclude taxa that pass the cutoff from other samples
all_with_others_Fungi_relaxed=data.frame()
all_tax_considered=unique(all_with_others_Fungi$tax)
all_tax_considered_without_other_bin=all_tax_considered[!grepl(x=all_tax_considered,pattern = "Other")]
Other_bin_name=all_tax_considered[grepl(x=all_tax_considered,pattern = "Other")]

for (sample in unique (all_samples_flat_Fungi$sample))
{
  subset=all_with_others_Fungi[all_with_others_Fungi$sample==sample&all_with_others_Fungi$tax!=Other_bin_name,]
  Other_bin_line=all_with_others_Fungi[all_with_others_Fungi$sample==sample&all_with_others_Fungi$tax==Other_bin_name,]
  for (tax in all_tax_considered_without_other_bin) # see if any tax included in the final set need to be added and correct the Other bin
  {
    if (sum(subset$tax==tax)==0) # was below the other cutoff and not included
    {
      if (nrow(all_samples_flat_Fungi[all_samples_flat_Fungi$sample==sample&all_samples_flat_Fungi$tax==tax,])>0) { # exists
        low_abundance_line=all_samples_flat_Fungi[all_samples_flat_Fungi$sample==sample&all_samples_flat_Fungi$tax==tax,]
        
        # substract from the Other bin
        Other_bin_line$count=Other_bin_line$count-low_abundance_line$count
        Other_bin_line$percent=Other_bin_line$percent-low_abundance_line$percent
        
        # add the low abundance line
        all_with_others_Fungi_relaxed=rbind(all_with_others_Fungi_relaxed,low_abundance_line)
        rm (low_abundance_line)
      } 
    }
  }
  # add all the other tax and the updated Other bin
  all_with_others_Fungi_relaxed=rbind(all_with_others_Fungi_relaxed,subset,Other_bin_line)
  rm (subset,Other_bin_line,tax)
}  

all_with_others_Fungi$name <- factor(all_with_others_Fungi$sample,levels = c("OMNI12.A1","OMNI12.A2","OMNI12.B2",
                                                                             "OMNI12.B1","OMNI12.C1","OMNI12.C2",
                                                                             "OMNI13.A1","OMNI13.A2","OMNI13.B1",
                                                                             "OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")
)
all_with_others_Fungi_relaxed$name <- factor(all_with_others_Fungi_relaxed$sample,levels = c("OMNI12.A1","OMNI12.A2","OMNI12.B2",
                                                                                             "OMNI12.B1","OMNI12.C1","OMNI12.C2",
                                                                                             "OMNI13.A1","OMNI13.A2","OMNI13.B1",
                                                                                             "OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")
)

# stack plots
# plot definitions
# Genus levels stackplots and heatmaps
# Another option: https://cran.r-project.org/web/packages/Polychrome/vignettes/testgg.html
library(Polychrome)
library(ggplot2)

P40 <- createPalette(40, c("#FF0000", "#00FF00", "#0000FF"), range = c(30, 80))
swatch(P40)
P40 <- sortByHue(P40)
P40 <- as.vector(t(matrix(P40, ncol=4)))
swatch(P40)
names(P40) <- NULL

all_with_others_Fungi=merge(all_with_others_Fungi,samples_label,by="sample")
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_ITS_Genus.profiles.pdf",width = 28,height = 14)

# old
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.UMI_filtered.pdf",width = 28,height = 14) # To make it more organized
# Apr2022 under the tissue
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.pdf",width = 28,height = 14) # To make it more organized

p=ggplot(all_with_others_Fungi, aes(fill=tax, y=percent, x=lable,label=paste(signif(percent*100,2),"% ",tax,sep=""))) + 
  geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
  theme(legend.text=element_text(size=15)) +
  scale_fill_manual(values = P40) + 
  # ggtitle(paste("Subsample",i_repeat,"composition",Level)) +
  theme(plot.title = element_text(hjust = 0.5))
print (p)                                      
dev.off()

# save the flat data
# old
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.UMI_filtered.data.tsv",
#             x=all_with_others_Fungi,quote = F,sep = "\t",row.names = F)
# Apr2022 under the tissue
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.data.tsv",
            x=all_with_others_Fungi,quote = F,sep = "\t",row.names = F)



# The relaxed version
all_with_others_Fungi_relaxed=merge(x=all_with_others_Fungi_relaxed,y=samples_label,by="sample")
# old undeer the tissue
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.UMI_filtered.relaxed.pdf",width = 28,height = 14) # To make sure identical with previous and make sure all in one place with the data-files

# Apr2022 under the tissue
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.pdf",width = 28,height = 14) # To make sure identical with previous and make sure all in one place with the data-files
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.25x14.pdf",width = 25,height = 14) # To make sure identical with previous and make sure all in one place with the data-files

p=ggplot(all_with_others_Fungi_relaxed, aes(fill=tax, y=percent, x=lable,label=paste(signif(percent*100,2),"% ",tax,sep=""))) + 
  geom_bar(position="fill", stat="identity") + geom_text(size = 3, position = position_stack(vjust = 0.5)) + theme_bw() +
  theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5,size = 15)) + 
  theme(legend.text=element_text(size=15)) +
  scale_fill_manual(values = P40) + 
  # ggtitle(paste("Subsample",i_repeat,"composition",Level)) +
  theme(plot.title = element_text(hjust = 0.5))
print (p)                                      
dev.off()

# save the flat data
# old under the tissue
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.UMI_filtered.relaxed.data.tsv",
#             x=all_with_others_Fungi_relaxed,quote = F,sep = "\t",row.names = F)

# Apr2022 under the tissue
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.data.tsv",
            x=all_with_others_Fungi_relaxed,quote = F,sep = "\t",row.names = F)

# Just as QA a file with both relaxed an non relaxed
tmp_both=merge(all_with_others_Fungi_relaxed,all_with_others_Fungi,by=c("sample","tax"),all=T)
names(tmp_both)[3:14]=gsub(x=names(tmp_both)[3:14],pattern = ".x$",replacement = ".relaxed",perl = T)
names(tmp_both)[3:14]=gsub(x=names(tmp_both)[3:14],pattern = ".y$",replacement = ".strict",perl = T)
# old under the tissue
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.UMI_filtered.relaxed_and_strict.data.tsv",
#             x=tmp_both,quote = F,sep = "\t",row.names = F)
# Apr2022 under the tissue

write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed_and_strict.data.tsv",
            x=tmp_both,quote = F,sep = "\t",row.names = F)
rm(tmp_both)

# some summary of average profile
average_Fungi_profile_df=data.frame()
for (g in unique(all_with_others_Fungi_relaxed$tax))
{
  r=data.frame(tax=g,
               mean=signif(x = mean(all_with_others_Fungi_relaxed$percent[all_with_others_Fungi_relaxed$tax==g]),3),
               sd=signif(x = sd(all_with_others_Fungi_relaxed$percent[all_with_others_Fungi_relaxed$tax==g]),3),
               var=signif(x = var(all_with_others_Fungi_relaxed$percent[all_with_others_Fungi_relaxed$tax==g]),3),
               median=signif(x=median(all_with_others_Fungi_relaxed$percent[all_with_others_Fungi_relaxed$tax==g]),3),
               n=sum(all_with_others_Fungi_relaxed$tax==g)
  )
  average_Fungi_profile_df=rbind(average_Fungi_profile_df,r)
}
average_Fungi_profile_df[order(average_Fungi_profile_df$median,decreasing = T),]

# save the average table
# old under the tissue
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.UMI_filtered.relaxed_average_profile.tsv",
#             x=average_Fungi_profile_df,quote = F,sep = "\t",row.names = F)

# Apr2022 under the tissue
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed_average_profile.tsv",
            x=average_Fungi_profile_df,quote = F,sep = "\t",row.names = F)

# make a table
all_with_others_Fungi_relaxed_table=make_table(flat_df = all_with_others_Fungi_relaxed[,c("tax","count","percent","lable")],column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
# old under the tissue
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.UMI_filtered.relaxed.samples_profile.tsv",
#             x=all_with_others_Fungi_relaxed_table,quote = F,sep = "\t",row.names = T)

# Apr2022 under the tissue 
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Fungi_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.samples_profile.tsv",
            x=all_with_others_Fungi_relaxed_table,quote = F,sep = "\t",row.names = T)

# NMDS -- here
library(vegan)

# percent_matrix_all_samples=make_table(flat_df = all_with_others,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
percent_matrix_all_samples_Fungi=make_table(flat_df = all_samples_flat_Fungi,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 3)
count_matrix_all_samples_Fungi=make_table(flat_df = all_samples_flat_Fungi,column_number_sample_name = 4,column_number_of_value = 1,column_number_of_count = 2)

community_matrix_Fungi=t(as.matrix(percent_matrix_all_samples_Fungi))
#community_matrix=t(as.matrix(percent_matrix_all_sub_samples))
community_matrix_Fungi[is.na(community_matrix_Fungi)]=0 # replace NA with 0
example_NMDS=metaMDS(community_matrix_Fungi, # Our community-by-species matrix
                     k=2,trymax = 100,distance = "bray") # The number of reduced dimensions -> ,, autotransform = FALSE,distance ="euclidean"
stressplot(example_NMDS)
# plot (example_NMDS)

#extract NMDS scores (x and y coordinates)
data.scores = as.data.frame(scores(example_NMDS))
data.scores$Sample=row.names(data.scores)

P1L1=c("OMNI12.A1","OMNI12.A2","OMNI12.B2")
P1L2=c("OMNI12.B1","OMNI12.C1","OMNI12.C2")
P2L1=c("OMNI13.A1","OMNI13.A2","OMNI13.B1")
P2L2=c("OMNI13.B2","OMNI13.C1","OMNI13.C2","OMNI13.D1")

data.scores$PlantLeaf=row.names(data.scores) 
data.scores$PlantLeaf[row.names(data.scores) %in% P1L1]="P1.L1"
data.scores$PlantLeaf[row.names(data.scores) %in% P1L2]="P1.L2"
data.scores$PlantLeaf[row.names(data.scores) %in% P2L1]="P2.L1"
data.scores$PlantLeaf[row.names(data.scores) %in% P2L2]="P2.L2"

section1=c("OMNI12.A1","OMNI12.B1","OMNI13.A1","OMNI13.D1")
section2=c("OMNI12.A2","OMNI12.C1","OMNI13.C1")
section3=c("OMNI13.A2","OMNI13.B2")
section4=c("OMNI12.B2","OMNI12.C2","OMNI13.B1","OMNI13.C2")

data.scores$Section=row.names(data.scores) 
data.scores$Section[row.names(data.scores) %in% section1]="1"
data.scores$Section[row.names(data.scores) %in% section2]="2"
data.scores$Section[row.names(data.scores) %in% section3]="3"
data.scores$Section[row.names(data.scores) %in% section4]="4"

nmds_plot = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  # geom_point(size = 4, aes( shape = Type, colour = Sample))+ 
  geom_point(size = 9, aes( shape = PlantLeaf, fill = Section,color=Section))+ 
  # geom_text(aes(label=row.names(data.scores)),hjust=-0.1, vjust=-0.1) +
  geom_text(aes(label=paste(PlantLeaf,Section,sep=".")),hjust=-0.1, vjust=-0.1) +
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
  labs(x = "NMDS1", shape = "Plant.Leaf",fill = "Section", y = "NMDS2",colour = "Section")  + #  shape = "Sample",
  scale_colour_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7"))  +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
  scale_shape_manual(values=c(21, 22, 23,24))+ 
  scale_fill_manual(values = c("#E69F00", "#56B4E9","#009E73","#CC79A7")) +  # https://jrnold.github.io/ggthemes/reference/colorblind.html
  ggtitle(paste("NMDS of the full ITS profiles on the Genus level (Stress=",signif(example_NMDS$stress,3),")",sep=""))
#pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13.ITS_genus.NMDS.pdf",width = 7,height = 7)

# old under the tissue
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.NMDS.pdf",width = 7,height = 7) # Just to make more organized

# Apr2022 under the tissue
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.under_the_tissue_Apr2022_220405.NMDS.pdf",width = 7,height = 7) # Just to make more organized

print (nmds_plot,width = 4, height = 4)
dev.off()

# old under the tissue
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.BaryCurtis_Dist.pdf",width = 7,height = 7) # Just to make more organized

# Apr2022 under the tissue
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.under_the_tissue_Apr2022_220405.BaryCurtis_Dist.pdf",width = 7,height = 7) # Just to make more organized

# hirarchial clustering and distances heatmaps
# calculate Bray-Curtis distance among samples
comm.bc.dist <- vegdist(community_matrix_Fungi, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust <- hclust(comm.bc.dist, method = "average")
# plot cluster diagram
plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")
library(circlize)
col_fun = colorRamp2(c(0, 1), c("red","white"))
col_fun(seq(0, 1))

library(ComplexHeatmap)
dist_mat=as.matrix(comm.bc.dist)
pHeatmapSim=Heatmap(dist_mat,col=col_fun,name="Distance", #column_title = paste("Subsample ",i_repeat,"- Bray-Curtis dissimilarity",sep="")
                    cell_fun = function(j, i, x, y, width, height, fill) {
                      grid.text(sprintf("%.2f", dist_mat[i, j]), x, y, gp = gpar(fontsize = 10))
                    })
print(pHeatmapSim) 
dev.off()

# save the tables
# old under the tissue
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.NMDS_community_matrix.tsv",
#             x=community_matrix_Fungi,
#             quote = F,sep = "\t",row.names = T
# )
# write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.NMDS_scores.tsv",
#             x=data.scores,
#             quote = F,sep = "\t",row.names = T
# )
# saveRDS(example_NMDS, file = "/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.NMDS_obj.rds")

# Apr2022 under the tissue
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.under_the_tissue_Apr2022_220405.NMDS_community_matrix.tsv",
            x=community_matrix_Fungi,
            quote = F,sep = "\t",row.names = T
)
write.table(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.Fungi_genus.under_the_tissue_Apr2022_220405.NMDS_scores.tsv",
            x=data.scores,
            quote = F,sep = "\t",row.names = T
)
saveRDS(example_NMDS, file = "/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13.under_the_tissue_Apr2022_220405.Fungi_genus.NMDS_obj.rds")



### Compare and validate if anything change as result of Apr2022 under the tissue

# profile data
# Bacteria
Bacteria_old_profile_all_samples=read.delim(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.UMI_filtered.relaxed.data.tsv",sep="\t",stringsAsFactors = F)
Bacteria_new_profile_all_samples=read.delim(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12_OMNI13_Profiles/OMNI12_OMNI13_Bacteria_Genus.profiles.under_the_tissue_Apr2022_220405.UMI_filtered.relaxed.data.tsv",sep="\t",stringsAsFactors = F)
names(Bacteria_old_profile_all_samples)=paste(names(Bacteria_old_profile_all_samples),"old",sep=".")
names(Bacteria_new_profile_all_samples)=paste(names(Bacteria_new_profile_all_samples),"new",sep=".")
joined=merge(Bacteria_old_profile_all_samples,Bacteria_new_profile_all_samples,by.x=c("sample.old","tax.old"),by.y=c("sample.new","tax.new"),all=T)
joined$percent_diff=joined$percent.new-joined$percent.old
joined$percent_diff=signif(joined$percent_diff*100,3)


#### LEFTOVERS IGNORE ALL BELOW THIS LINE ####
##############################################
dist_mat_flat=flattenMatrix(dist_mat)
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
  
#### LEFTOVERS  
  


"A2", "AmpSeq_A.515_806", "AmpSeq_A.799_1192",
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





# Global vars and data reading





path="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12/data/multimaped_filtered/


New_Probs_Dec2020/data/multimaped_filtered/MMSEQS2/UMI_FILTERED/"
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

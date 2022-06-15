args = commandArgs(trailingOnly=TRUE)

in_RNA_spatial=args[1]   # in_RNA_spatial="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.OMNI13.A1.tsv"
                         # in_RNA_spatial="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/Apr2022/220405_semiwild_dataset_rawcounts_filtered.OMNI12.A1.tsv"
                         # in_RNA_spatial="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/Apr2022/hotspots/RNASumCutoff100/220405_semiwild_dataset_rawcounts_filtered.PerGene_localG_2x2.Genes_with_total_sum_grt100.OMNI13.C2.tsv"
in_Taxon_spatial=args[2] # in_Taxon_spatial="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/Fungi_and_UNKNOWN/OMNI13_A1.All_ITS_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top_50_genus.spatial_pos.csv"
                         # in_Taxon_spatial="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/Bacterial_and_UNKNOWN/Under_tissue/Apr22_expression_220405/OMNI12_A1.All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top_50_genus.under_tissue.based_on_joined_Bacteria_Fungi_Expression_220405.spatial_pos.csv"
                         # in_Taxon_spatial="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/Bacterial_and_UNKNOWN/Under_tissue/Apr22_expression_220405/hotspots/OMNI13_C2.All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top_50_genus.2x2_Getis_Ord_localG_hotspots.under_tissue_based_on_joined_Bacteria_Fungi_Expression_220405.spatial_pos.csv"
out_prefix=args[3]       # out_prefix="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Apr2022/220405_semiwild_dataset_rawcounts_filtered.RNA.with_Fungi.A1"
                         # out_prefix="test.prefix."
top_tax_corr=args[4]     # top_tax_corr=1
top_tax_RF=args[5]       # top_tax_RF=1
num.threads=args[6]      # num.threads=1
Bacterial_cutoff=args[7] # Bacterial_cutoff=-99999 for hotspots
RNA_cutoff=args[8]       # RNA_cutoff=-99999 for hotspots
Boruta_RNAsum_cutoff=args[9]  # Boruta_RNAsum_cutoff=-99999 for hotspots
Pixels_subset_file=args[10] # optional, if provided would take anly these pixels from the datasets
                            # Pixels_subset_file="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/OMNI12_A1.Bacteria_Fungi_expression_220405.under_the_tissue_pixels.txt"
                            # Pixels_subset_file=NA

# defaults for Boruta
if (is.na(num.threads) | num.threads=="NA"){num.threads=1} else {
  num.threads=as.numeric(num.threads)
}

# defaults For the correlation
if (is.na(Bacterial_cutoff)) {
  Bacterial_cutoff=5
} else {
  Bacterial_cutoff=as.numeric(Bacterial_cutoff)
}
if (is.na(RNA_cutoff)) {
  RNA_cutoff=10
} else {
  RNA_cutoff=as.numeric(RNA_cutoff)
}

if (is.na(Boruta_RNAsum_cutoff)) {
  Boruta_RNAsum_cutoff=0
} else {
  Boruta_RNAsum_cutoff=as.numeric(Boruta_RNAsum_cutoff)
}

Sprot_annotation_file="/ebio/abt6_projects8/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/DB/uniprot_reviewed_ATh_3702.tab"
# Sprot_annotation_file="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/DB/uniprot_reviewed_ATh_3702.tab"
# top_tax_corr=1
# top_tax_RF=1
n_trees=1000
message(paste("in_RNA:",in_RNA_spatial))
message(paste("in_Taxon_spatial:",in_Taxon_spatial))
message(paste("out_prefix:",out_prefix))

message(paste("Sprot_annotation_file:",Sprot_annotation_file))
message(paste("top_tax_corr:",top_tax_corr))
message(paste("top_tax_RF:",top_tax_RF))
message(paste("Bacterial_cutoff:",Bacterial_cutoff))
message(paste("RNA_cutoff:",RNA_cutoff))
message(paste("Boruta_RNAsum_cutoff:",Boruta_RNAsum_cutoff))
message(paste("Boruta_threads:",num.threads))

message(paste("Pixels_subset_file:",Pixels_subset_file))

Pixels_subset_list=data.frame()
if (!is.na(Pixels_subset_file))
{
  Pixels_subset_list=read.delim(file=Pixels_subset_file,stringsAsFactors = F,sep = "\t")
}  
RNA_data=read.delim(file=in_RNA_spatial,sep = "\t",stringsAsFactors = F)
message("== Sample data: RNA_data[1:10,1:10]")
print(RNA_data[1:min(10,NROW(RNA_data)),1:min(10,NCOL(RNA_data))])
# dim(RNA_data)
if (!is.na(Pixels_subset_file)) # exptraxt only the requested pixels
{
  RNA_data=RNA_data[,paste("X",Pixels_subset_list$xy,sep="")]
}
# dim(RNA_data)

# remove pixels with more than 75% of the genes are NA
gene_with_NA_per_pixel=colSums(is.na(RNA_data))
pixels_with_more_than75perNA=names(gene_with_NA_per_pixel)[gene_with_NA_per_pixel/NROW(RNA_data)>0.75]
if (length(pixels_with_more_than75perNA)>0)
{
  message(paste("[WARNING] Filtering ",length(pixels_with_more_than75perNA)," pixels with more than 75% missing data: ",paste(pixels_with_more_than75perNA,collapse = ", ")," {before filtering: ",NCOL(RNA_data)," pixels}",sep=""))
  RNA_data=RNA_data[,names(gene_with_NA_per_pixel)[gene_with_NA_per_pixel/NROW(RNA_data)<=0.75]]
  message(paste("[INFO] Pixels left after filtering: ",NCOL(RNA_data),sep=""))
}
# remove genes with NA in one of the pixels
num_genes_with_NA=sum(is.na(rowSums(RNA_data)))
if (num_genes_with_NA>0)
{
  names_genes_with_NA=row.names(RNA_data)[is.na(rowSums(RNA_data))]
  message(paste("[WARNING] Filtering ",num_genes_with_NA," genes with some missing data: ",paste(names_genes_with_NA,collapse = ", ")," {before filtering: ",NROW(RNA_data),"}",sep=""))
  RNA_data=RNA_data[!is.na(rowSums(RNA_data)),]
  message(paste("[INFO] Left with ",NROW(RNA_data)," genes after filtering",sep=""))
}

Taxon_Spatial=read.delim(file = in_Taxon_spatial,sep = ";",header = T,stringsAsFactors = F)
row.names(Taxon_Spatial)=Taxon_Spatial[,1]
Taxon_Spatial=Taxon_Spatial[,-1]
# dim (Taxon_Spatial)
if (!is.na(Pixels_subset_file)) # exptraxt only the requested pixels
{
  Taxon_Spatial=Taxon_Spatial[,paste("X",Pixels_subset_list$xy,sep="")]
}
# dim (Taxon_Spatial)
# For the correlation
RNA_microbial_corr=data.frame(gene=c(),genus=c(),r=c(),p_val=c(),genus_sum=c(),gene_sum=c(),i_genus=c(),i_gene=c())

# correlate all top microbial profiles with the RNA seq 
for (i_genus in 1:top_tax_corr)
{
  i_tax_profile=Taxon_Spatial[i_genus,]
  # dim(i_tax_profile)
  common_pos=intersect(names(i_tax_profile),names(RNA_data))  #names(i_tax_profile)[names(i_tax_profile) %in% names(RNA_data)] # pos where we have both RNA and Microbial data --> under the tissue
  RNA_for_tax_profile=RNA_data[,common_pos]
  i_tax_profile_corresponding=i_tax_profile[,common_pos]
  # dim(i_tax_profile_corresponding)
  genus_sum=sum(unlist(i_tax_profile_corresponding))
  genus_name=row.names(i_tax_profile_corresponding)
  if (genus_sum>=Bacterial_cutoff) {
    for (i_gene in 1:nrow(RNA_for_tax_profile)) {
      gene_sum=sum(unlist(RNA_for_tax_profile[i_gene,]))
      if (gene_sum>RNA_cutoff)
      {
        res=cor.test(unlist(i_tax_profile_corresponding),unlist(RNA_for_tax_profile[i_gene,]))
        rec=data.frame(gene=row.names(RNA_for_tax_profile)[i_gene],genus=genus_name,r=res$estimate,p_val=res$p.value,genus_sum=genus_sum,gene_sum=gene_sum,i_genus=i_genus,i_gene=i_gene)
        RNA_microbial_corr=rbind(RNA_microbial_corr,rec)
      } else {message (paste("\t-- skipped correlation with",row.names(RNA_for_tax_profile)[i_gene],"total expression:",gene_sum))}
    } 
  }  else {message (paste("-- skipped correlation with tax",row.names(Taxon_Spatial)[i_genus]))}
}
RNA_microbial_corr$FDR.BH=p.adjust(RNA_microbial_corr$p_val,"BH")
RNA_microbial_corr$gene.orig=RNA_microbial_corr$gene # to remove extra suffix with "_.*$"
RNA_microbial_corr$gene=gsub(x=RNA_microbial_corr$gene,pattern = "_.*$",replacement = "",perl=T)
corr_file=paste(out_prefix,".corr.M",Bacterial_cutoff,"_R",RNA_cutoff,".tsv",sep="")
write.table(x=RNA_microbial_corr,file = corr_file,sep = "\t",row.names=F)

# add annotations
# annotation=read.delim(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12/data/slide12_raw_counts.RNA.gene_names.UNIPROT_map.tab",sep="\t")
# annotation$Cross.reference..Araport.=gsub(x=annotation$Cross.reference..Araport.,pattern = ";$",replacement = "",perl = T)
# corr_file_with_annotation=merge(RNA_microbial_corr,annotation,by.x="gene",by.y="Cross.reference..Araport.",all.x=T)
# corr_file_annotated=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12/data/slide12_raw_counts.RNA.C2.corr.top50_genus_M5_R10.annotated.tsv",sep="")
# write.table(x=corr_file_with_annotation,file = corr_file_annotated,sep = "\t",row.names=F)

Sprot_annotation=read.delim(file=Sprot_annotation_file,sep="\t")
Sprot_annotation$Cross.reference..Araport.=gsub(x=Sprot_annotation$Cross.reference..Araport.,pattern = ";$",replacement = "",perl = T)
corr_file_with_annotation=merge(RNA_microbial_corr,Sprot_annotation,by.x="gene",by.y="Cross.reference..Araport.",all.x=T)
corr_file_annotated=paste(out_prefix,".corr.M",Bacterial_cutoff,"_R",RNA_cutoff,".Sprot_annotated.tsv",sep="")
write.table(x=corr_file_with_annotation,file = corr_file_annotated,sep = "\t",row.names=F)


# names of genes with abs(r)>=0.2
# unique(RNA_microbial_corr$gene[abs(RNA_microbial_corr$r)>=0.05])
message("== Pairs with r>=0.05")
print(paste(unique(RNA_microbial_corr$gene[RNA_microbial_corr$r>=0.05]),sep="\n"))


## Prepare the dataset for RandomForest
# each genus is a dataset for which we try to identify the genes that are correlated with
# each spot on the array is a sample with ~14,000 features (genes) trying to predict the abundance of the genus on that spot
# each spot should be normalized
# all gene names

library(randomForest)
library(Boruta)
expression_profile_per_spot=as.data.frame(t(RNA_data)) # each row is a spot and each col is a gene
# dim(expression_profile_per_spot)
message("== Sample data: expression_profile_per_spot[1:10,1:10]")
print(expression_profile_per_spot[1:min(10,NROW(expression_profile_per_spot)),1:min(10,NCOL(expression_profile_per_spot))])
expression_profile_per_spot.non_ambiguous=expression_profile_per_spot[,!grepl (pattern = "ambiguous",names (expression_profile_per_spot))]
# dim (expression_profile_per_spot.non_ambiguous)
tax_profile_per_spot=as.data.frame(t(Taxon_Spatial))
# dim(tax_profile_per_spot)
message("== Sample data: tax_profile_per_spot[1:10,1:10]")
print(tax_profile_per_spot[1:min(10,NROW(tax_profile_per_spot)),1:min(10,top_tax_RF,NCOL(tax_profile_per_spot))])

overlaping_spots_RNA_microbial=intersect(row.names(tax_profile_per_spot),row.names(expression_profile_per_spot))
# length(overlaping_spots_RNA_microbial)
expression_profile_per_spot_overlap_microbial=expression_profile_per_spot.non_ambiguous[overlaping_spots_RNA_microbial,]
# dim (expression_profile_per_spot_overlap_microbial)
# out_prefix=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI6/data/multimaped_filtered/210630_lableaf062_raw_counts_mtrbcpfilt.noCp.A1.",sep="")
# look on this parralel option: https://www.listendata.com/2015/05/speeding-up-random-forest.html
for (i_tax in 1:top_tax_RF){
  tax_name=names(tax_profile_per_spot)[i_tax]
  message(paste ("=== new Boruta starts: ",Sys.time(),sep=""))
  message(paste("Boruta tax name:",tax_name))
  i_tax_profile_per_spot=data.frame(row.names=overlaping_spots_RNA_microbial,tax=tax_profile_per_spot[overlaping_spots_RNA_microbial,i_tax])
  subset.df=merge(expression_profile_per_spot_overlap_microbial,i_tax_profile_per_spot,by="row.names")
  row.names(subset.df)=subset.df$Row.names
  # names(subset.df)=="Row.names"
  subset.df=subset.df[,!names(subset.df)=="Row.names"]
  # filter genes WO expression
  subset.df_above_min_expression=subset.df[,colSums(subset.df)>Boruta_RNAsum_cutoff]
  
  # if by mistake the Boruta_RNAsum_cutoff remove the tax col add it manually...
  if (!("tax" %in% names(subset.df_above_min_expression)))
  {
    subset.df_above_min_expression$tax=subset.df$tax
  }
  
  # report how many were filtered by Boruta_RNAsum_cutoff
  after_filter=NCOL(subset.df_above_min_expression)
  before_filter=NCOL(subset.df)
  message(paste ("[INFO] when applyeing the filter of colSums>",Boruta_RNAsum_cutoff," on the dataset, ",before_filter-after_filter," features were removed [",after_filter,":",before_filter,"]",sep=""))
  
  # remove keep.forest = FALSE, -> issue a warning in ranger:
  # In ranger::ranger(data = x, dependent.variable.name = "shadow.Boruta.decision",  ... :
  #                    Unused arguments: keep.forest
  message("Boruta(tax~.,data=subset.df_above_min_expression,ntree=n_trees,num.threads=num.threads)") # keep.forest = FALSE,
  Boruta_results=Boruta(tax~.,data=subset.df_above_min_expression,ntree=n_trees,num.threads=num.threads) # keep.forest = FALSE,
  
  # save files
  gene_importance_table_file=paste(out_prefix,".",i_tax,".",tax_name,".ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".per_spot_RF.Boruta_varImp.tsv",sep="")
  Boruta_obj_file=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".Boruta_Obj.R",sep="")
  gene_importance_plot=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".Boruta_varImp.pdf",sep="")
  # learning_curve_plot=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".learning_curve.pdf",sep="")
  
  save(Boruta_results,file=Boruta_obj_file)
  
  pdf (gene_importance_plot,width=60,height = 15)
  plot(Boruta_results)
  # varImpPlot(RF) # ,cex = 0.8)
  dev.off()
  
  print(Boruta_results)      
  
  # join annotation and var importance
  # 
  Boruta_importance=attStats(Boruta_results)
  Boruta_importance$orig_name=row.names(Boruta_importance)
  row.names(Boruta_importance)=gsub(x=row.names(Boruta_importance),pattern = "_.*$",replacement = "",perl=T)
  Boruta_importance_and_annotation=merge(Boruta_importance,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
  write.table(file = gene_importance_table_file,x=Boruta_importance_and_annotation,sep = "\t",quote = FALSE, row.names = FALSE)
  
  # Do the TentativeRoughFix
  Boruta_TentativeRoughFix_obj_file=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".Boruta_TentativeRoughFix_Obj.R",sep="")
  Boruta_importance_with_RoughFix_decision_gene_importance_table_file=paste(out_prefix,".",i_tax,".",tax_name,".ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="")
    
  Boruta_results_TentativeRoughFix=TentativeRoughFix(Boruta_results)
  save(Boruta_results_TentativeRoughFix, file = Boruta_TentativeRoughFix_obj_file)
  print(Boruta_results_TentativeRoughFix)
  
  # Boruta_importance=attStats(Boruta_results)
  Boruta_importance_TentativeRoughFix=attStats(Boruta_results_TentativeRoughFix)
  
  
  decision_TentativeRoughFix=data.frame(row.names=row.names(Boruta_importance_TentativeRoughFix),decision_RoughFix=Boruta_importance_TentativeRoughFix$decision)
  
  Boruta_importance_with_RoughFix_decision=merge(Boruta_importance,decision_TentativeRoughFix,by.x="orig_name",by.y="row.names")
  row.names(Boruta_importance_with_RoughFix_decision)=gsub(x=Boruta_importance_with_RoughFix_decision$orig_name,pattern = "_.*$",replacement = "",perl=T)
  # Boruta_importance_with_RoughFix_decision=merge(Boruta_importance,decision_TentativeRoughFix,by="row.names")
  # row.names(Boruta_importance_with_RoughFix_decision)=Boruta_importance_with_RoughFix_decision$Row.names
  # Boruta_importance_with_RoughFix_decision=Boruta_importance_with_RoughFix_decision[,-1]
  
  Boruta_importance_with_RoughFix_decision$final_decision=as.character(Boruta_importance_with_RoughFix_decision$decision)
  
  Boruta_importance_with_RoughFix_decision$final_decision[Boruta_importance_with_RoughFix_decision$decision!=Boruta_importance_with_RoughFix_decision$decision_RoughFix]=paste(Boruta_importance_with_RoughFix_decision$decision_RoughFix[Boruta_importance_with_RoughFix_decision$decision!=Boruta_importance_with_RoughFix_decision$decision_RoughFix],".RoughFix",sep="")
  
  # Boruta_importance_and_annotation=merge(Boruta_results,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
  # Boruta_importance_and_annotation=merge(Boruta_importance,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
  Boruta_importance_with_RoughFix_decision_and_annotation=merge(Boruta_importance_with_RoughFix_decision,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
  names(Boruta_importance_with_RoughFix_decision_and_annotation)[1]=c("gene")
  write.table(file = Boruta_importance_with_RoughFix_decision_gene_importance_table_file,x= Boruta_importance_with_RoughFix_decision_and_annotation,sep="\t",quote = FALSE,row.names = FALSE)
  
  
  message(paste ("=== Boruta ended: ",Sys.time(),sep=""))
  
  # Do the RF step to get the accuracy of the model  
  message(paste ("=== RF model with all features started: ",Sys.time(),sep=""))
  message("randomForest(tax ~ ., data=subset.df_above_min_expression,ntree=n_trees,importance=TRUE,keep.forest = FALSE)")
  RF=randomForest(tax ~ ., data=subset.df_above_min_expression,ntree=n_trees,importance=TRUE,keep.forest = FALSE)
  
  # save files
  gene_importance_table_file=paste(out_prefix,".",i_tax,".",tax_name,".ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".per_spot_RF.RF_varImp.tsv",sep="")
  RF_obj_file=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".RF_Obj.R",sep="")
  gene_importance_plot=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".RF_varImp.pdf",sep="")
  learning_curve_plot=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".learning_curve.pdf",sep="")
  
  save(RF,file=RF_obj_file)
  
  # plot the RF var importance
  pdf (gene_importance_plot,width=60,height = 15)
  varImpPlot(RF) # ,cex = 0.8)
  dev.off()
  
  # text version of the model
  print(RF)      
  
  # plot the learning curve
  pdf (learning_curve_plot,width=20,height = 15)
  plot (RF)
  dev.off()
  
  # join annotation and var importance
  RF_importance=as.data.frame(RF$importance)
  RF_importance$orig_name=row.names(RF_importance)
  row.names(RF_importance)=gsub(x=row.names(RF_importance),pattern = "_.*$",replacement = "",perl=T)
  RF_importance_and_annotation=merge(RF_importance,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
  write.table(file = gene_importance_table_file,x= RF_importance_and_annotation,sep="\t",quote = FALSE,row.names = FALSE)

  message(paste ("=== RF with all features ended: ",Sys.time(),sep=""))

  message(paste ("=== RF model with BORUTA selected features started: ",Sys.time(),sep=""))
  Boruta_final_selected_genes=Boruta_importance_with_RoughFix_decision_and_annotation$orig_name[grepl(x = Boruta_importance_with_RoughFix_decision_and_annotation$final_decision,pattern="Confirmed")]
  message(paste("= Total selected features by BORUTA:",length(Boruta_final_selected_genes)))
  message("randomForest(tax ~ ., data=subset.df_above_min_expression[,c(Boruta_final_selected_genes,\"tax\")],ntree=n_trees,importance=TRUE,keep.forest = FALSE)")
  RF_BORUTA_Selected=randomForest(tax ~ ., data=subset.df_above_min_expression[,c(Boruta_final_selected_genes,"tax")],ntree=n_trees,importance=TRUE,keep.forest = FALSE)
  
  # save files
  gene_importance_table_file=paste(out_prefix,".",i_tax,".",tax_name,".ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".BORUTA_selected_features.per_spot_RF.RF_varImp.tsv",sep="")
  RF_obj_file=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".BORUTA_selected_features.RF_Obj.R",sep="")
  gene_importance_plot=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".BORUTA_selected_features.RF_varImp.pdf",sep="")
  learning_curve_plot=paste(out_prefix,".",i_tax,".",tax_name,".per_spot_RF.ntree",n_trees,".RNAsum",Boruta_RNAsum_cutoff,".BORUTA_selected_features.learning_curve.pdf",sep="")
  
  save(RF_BORUTA_Selected,file=RF_obj_file)
  
  # plot the RF var importance
  pdf (gene_importance_plot,width=60,height = 15)
  varImpPlot(RF_BORUTA_Selected) # ,cex = 0.8)
  dev.off()
  
  # text version of the model
  print(RF_BORUTA_Selected)      
  
  # plot the learning curve
  pdf (learning_curve_plot,width=20,height = 15)
  plot (RF_BORUTA_Selected)
  dev.off()
  
  # join annotation and var importance
  RF_importance=as.data.frame(RF_BORUTA_Selected$importance)
  RF_importance$orig_name=row.names(RF_importance)
  row.names(RF_importance)=gsub(x=row.names(RF_importance),pattern = "_.*$",replacement = "",perl=T)
  RF_importance_and_annotation=merge(RF_importance,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
  write.table(file = gene_importance_table_file,x= RF_importance_and_annotation,sep="\t",quote = FALSE,row.names = FALSE)
  
  message(paste ("=== RF with BORUTA selected features ended: ",Sys.time(),sep=""))
  
}

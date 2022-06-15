library('openxlsx')
library (ComplexHeatmap)

plot_upSetR=function (elements_list,title=NA) {
  if (is.na(title)) {title=""}
  m = make_comb_mat(elements_list) # https://support.bioconductor.org/p/118557/
  cs = comb_size(m)
  row_size = set_size(m)
  
  ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))),column_title = title)
  ht = draw(ht)
  co = column_order(ht)
  row_od = row_order(ht)
  
  nc = ncol(m)
  
  decorate_annotation("Intersection\nsize", {
    grid.text(cs[co], 
              x = 1:nc, 
              y = unit(cs[co], "native") + unit(1, "mm"), 
              gp = gpar(fontsize = 9), 
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
}


# create xls file
# wb <- createWorkbook()

# check which are missing
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (S in samples)
  {
    # reads based
    # Boruta_reads_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.",S,".1.Top50_genus.per_spot_RFntree1000.per_spot_RF.Boruta_varImp.tsv",sep="");
    # Boruta_reads_Bacteria_File=paste("/Volumes//spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.",S,".1.Top50_genus.per_spot_RFntree1000.per_spot_RF.Boruta_varImp.tsv",sep="");
    # Fungi_reads_cor_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.",S,".corr.M5_R10.tsv",sep="")
    # Bacteria_reads_cor_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.",S,".corr.M5_R10.tsv",sep="");
    
    # hotspots based
    # Boruta_localG_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/hotspots/Boruta/joint_RNA_16S_ITS_Files/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.2x2.Getis_Ord.joint_spots.",S,".1.Getis_Ord_2x2_G.per_spot_RFntree1000.per_spot_RF.Boruta_varImp.tsv",sep="");
    # Boruta_localG_Bacteria_File=paste("/Volumes//spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/hotspots/Boruta/joint_RNA_16S_ITS_Files/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.2x2.Getis_Ord.joint_spots.",S,".1.Getis_Ord_2x2_G.per_spot_RFntree1000.per_spot_RF.Boruta_varImp.tsv",sep="");
    # Fungi_localG_cor_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/hotspots/Boruta/joint_RNA_16S_ITS_Files/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.2x2.Getis_Ord.joint_spots.",S,".corr.M-99999_R-99999.Sprot_annotated.tsv",sep="")
    # Bacteria_localG_cor_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/hotspots/Boruta/joint_RNA_16S_ITS_Files/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.2x2.Getis_Ord.joint_spots.",S,".corr.M-99999_R-99999.Sprot_annotated.tsv",sep="")
    
    ## Apr2022 expression data
    # reads based
    Boruta_reads_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum10/","220405_semiwild_dataset_rawcounts_filtered.with_ITS_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.",Exp,".",S,".1.Top50_ITS_genus.ntree1000.RNAsum10.per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="");
    Boruta_reads_Bacteria_File=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum10/","220405_semiwild_dataset_rawcounts_filtered.with_Bacteria_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.",Exp,".",S,".1.Top50_Bacteria_genus.ntree1000.RNAsum10.per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="");
    Fungi_reads_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum10/","220405_semiwild_dataset_rawcounts_filtered.with_ITS_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.",Exp,".",S,".corr.M5_R10.tsv",sep="")
    Bacteria_reads_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum10/","220405_semiwild_dataset_rawcounts_filtered.with_Bacteria_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.",Exp,".",S,".corr.M5_R10.tsv",sep="");
    
    # hotspots based
    Boruta_localG_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum100/hotspots/220405_semiwild_dataset_rawcounts_filtered.Genes_with_total_sum_grt100.with_ITS_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.localG_2x2.",Exp,".",S,".1.all_Fungi_LocalG.ntree1000.RNAsum-99999.per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="")
    Boruta_localG_Bacteria_File=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum100/hotspots/220405_semiwild_dataset_rawcounts_filtered.Genes_with_total_sum_grt100.with_Bacteria_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.localG_2x2.",Exp,".",S,".1.all_Bacteria_LocalG.ntree1000.RNAsum-99999.per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="")
    Fungi_localG_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum100/hotspots/220405_semiwild_dataset_rawcounts_filtered.Genes_with_total_sum_grt100.with_ITS_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.localG_2x2.",Exp,".",S,".corr.M-99999_R-99999.Sprot_annotated.tsv",sep="")
    Bacteria_localG_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum100/hotspots/220405_semiwild_dataset_rawcounts_filtered.Genes_with_total_sum_grt100.with_Bacteria_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.localG_2x2.",Exp,".",S,".corr.M-99999_R-99999.Sprot_annotated.tsv",sep="")
    
    for (file in c(Boruta_reads_Fungi_file,Boruta_reads_Bacteria_File,Fungi_reads_cor_file,Bacteria_reads_cor_file,
                  Boruta_localG_Fungi_file,Boruta_localG_Bacteria_File,Fungi_localG_cor_file,Bacteria_localG_cor_file))
    {
      if (!file.exists(file)) {
        message (paste (Exp,S,file," is missing!",sep="\t"))
      }
    }
  }
}
rm ("Bacteria_localG_cor_file","Bacteria_reads_cor_file","Boruta_localG_Bacteria_File","Boruta_localG_Fungi_file","Boruta_reads_Bacteria_File", 
   "Boruta_reads_Fungi_file","Exp","file","Fungi_localG_cor_file","Fungi_reads_cor_file"       
   ,"S")    
      

# analyze

Sprot_annotation_file="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/DB/uniprot_reviewed_ATh_3702.tab"
Sprot_annotation=read.delim(file=Sprot_annotation_file,sep="\t")
Sprot_annotation$Cross.reference..Araport.=gsub(x=Sprot_annotation$Cross.reference..Araport.,pattern = ";$",replacement = "",perl = T)


selected_genes_from_all_sections=data.frame()
selected_genes_from_all_sections_with_all_details=data.frame()
selected_global_stats=data.frame() # how many were selected by each type out of the possible candidates

for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (S in samples)
  {
    # reads based
    # Boruta_reads_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.",S,".1.Top50_genus.per_spot_RFntree1000.per_spot_RF.Boruta_varImp.tsv",sep="");
    # Boruta_reads_Bacteria_File=paste("/Volumes//spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.",S,".1.Top50_genus.per_spot_RFntree1000.per_spot_RF.Boruta_varImp.tsv",sep="");
    # Fungi_reads_cor_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.",S,".corr.M5_R10.tsv",sep="")
    # Bacteria_reads_cor_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.",S,".corr.M5_R10.tsv",sep="");
    
    # hotspots based
    # Boruta_localG_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/hotspots/Boruta/joint_RNA_16S_ITS_Files/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.2x2.Getis_Ord.joint_spots.",S,".1.Getis_Ord_2x2_G.per_spot_RFntree1000.per_spot_RF.Boruta_varImp.tsv",sep="");
    # Boruta_localG_Bacteria_File=paste("/Volumes//spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/hotspots/Boruta/joint_RNA_16S_ITS_Files/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.2x2.Getis_Ord.joint_spots.",S,".1.Getis_Ord_2x2_G.per_spot_RFntree1000.per_spot_RF.Boruta_varImp.tsv",sep="");
    # Fungi_localG_cor_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/hotspots/Boruta/joint_RNA_16S_ITS_Files/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.2x2.Getis_Ord.joint_spots.",S,".corr.M-99999_R-99999.Sprot_annotated.tsv",sep="")
    # Bacteria_localG_cor_file=paste("/Volumes/spatial_array_metatranscriptomics//data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/hotspots/Boruta/joint_RNA_16S_ITS_Files/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.2x2.Getis_Ord.joint_spots.",S,".corr.M-99999_R-99999.Sprot_annotated.tsv",sep="")
    
    ## Apr2022 expression data
    # reads based
    Boruta_reads_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum10/","220405_semiwild_dataset_rawcounts_filtered.with_ITS_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.",Exp,".",S,".1.Top50_ITS_genus.ntree1000.RNAsum10.per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="");
    Boruta_reads_Bacteria_File=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum10/","220405_semiwild_dataset_rawcounts_filtered.with_Bacteria_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.",Exp,".",S,".1.Top50_Bacteria_genus.ntree1000.RNAsum10.per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="");
    Fungi_reads_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum10/","220405_semiwild_dataset_rawcounts_filtered.with_ITS_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.",Exp,".",S,".corr.M5_R10.tsv",sep="")
    Bacteria_reads_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum10/","220405_semiwild_dataset_rawcounts_filtered.with_Bacteria_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.",Exp,".",S,".corr.M5_R10.tsv",sep="");
    
    # hotspots based
    Boruta_localG_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum100/hotspots/220405_semiwild_dataset_rawcounts_filtered.Genes_with_total_sum_grt100.with_ITS_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.localG_2x2.",Exp,".",S,".1.all_Fungi_LocalG.ntree1000.RNAsum-99999.per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="")
    Boruta_localG_Bacteria_File=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum100/hotspots/220405_semiwild_dataset_rawcounts_filtered.Genes_with_total_sum_grt100.with_Bacteria_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.localG_2x2.",Exp,".",S,".1.all_Bacteria_LocalG.ntree1000.RNAsum-99999.per_spot_RF.Boruta_with_TentativeRoughFix.varImp.tsv",sep="")
    Fungi_localG_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum100/hotspots/220405_semiwild_dataset_rawcounts_filtered.Genes_with_total_sum_grt100.with_ITS_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.localG_2x2.",Exp,".",S,".corr.M-99999_R-99999.Sprot_annotated.tsv",sep="")
    Bacteria_localG_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta/RNASum100/hotspots/220405_semiwild_dataset_rawcounts_filtered.Genes_with_total_sum_grt100.with_Bacteria_and_UNKNOWN.UMI_filtered.SUM_Top_50_genus.localG_2x2.",Exp,".",S,".corr.M-99999_R-99999.Sprot_annotated.tsv",sep="")
    
    # Continue only if all output files exists...
    if (file.exists(Boruta_reads_Fungi_file) & file.exists(Boruta_reads_Fungi_file) & file.exists(Fungi_reads_cor_file) & file.exists(Bacteria_reads_cor_file) &
        file.exists(Boruta_localG_Fungi_file) & file.exists(Boruta_localG_Fungi_file) & file.exists(Fungi_localG_cor_file) & file.exists(Boruta_localG_Fungi_file)) 
    {
      # reads based
      Boruta_reads_Fungi_data=read.delim(file = Boruta_reads_Fungi_file,sep="\t",stringsAsFactors = F)
      Boruta_reads_Bacterial_data=read.delim(file = Boruta_reads_Bacteria_File,sep="\t",stringsAsFactors = F)
      Fungi_reads_cor_data=read.delim(file = Fungi_reads_cor_file,sep="\t",stringsAsFactors = F)
      Bacterial_reads_cor_data=read.delim(file = Bacteria_reads_cor_file,sep="\t",stringsAsFactors = F)

      # localG based
      Boruta_localG_Fungi_data=read.delim(file = Boruta_localG_Fungi_file,sep="\t",stringsAsFactors = F)
      Boruta_localG_Bacterial_data=read.delim(file = Boruta_localG_Bacteria_File,sep="\t",stringsAsFactors = F)
      Fungi_localG_cor_data=read.delim(file = Fungi_localG_cor_file,sep="\t",stringsAsFactors = F)
      Bacterial_localG_cor_data=read.delim(file = Bacteria_localG_cor_file,sep="\t",stringsAsFactors = F)
      
      # Boruta
      Boruta_reads_Bacterial_selected=Boruta_reads_Bacterial_data[grepl(x = Boruta_reads_Bacterial_data$final_decision,pattern = "Confirmed",fixed = T),]
      Boruta_reads_Fungi_selected=Boruta_reads_Fungi_data[grepl(x = Boruta_reads_Fungi_data$final_decision ,pattern = "Confirmed",fixed = T),]
      
      Boruta_localG_Bacterial_selected=Boruta_localG_Bacterial_data[grepl(x = Boruta_localG_Bacterial_data$final_decision,pattern = "Confirmed",fixed = T),]
      Boruta_localG_Fungi_selected=Boruta_localG_Fungi_data[grepl(x = Boruta_localG_Fungi_data$final_decision ,pattern = "Confirmed",fixed = T),]
      
      # cor based
      signif_tresh=0.01
      Bacterial_reads_cor_selected=Bacterial_reads_cor_data[Bacterial_reads_cor_data$FDR.BH<signif_tresh,]
      Fungi_reads_cor_selected=Fungi_reads_cor_data[Fungi_reads_cor_data$FDR.BH<signif_tresh,]
      
      Bacterial_localG_cor_selected=Bacterial_localG_cor_data[Bacterial_localG_cor_data$FDR.BH<signif_tresh,]
      Fungi_localG_cor_selected=Fungi_localG_cor_data[Fungi_localG_cor_data$FDR.BH<signif_tresh,]
      
      ### for the BORUTA selected genes filter those not significant by cor -> NEW FOR THIS ANALYSES
      
      Boruta_reads_Bacterial_selected_and_significant_cor_gene_id=Boruta_reads_Bacterial_selected$gene[(Boruta_reads_Bacterial_selected$gene %in% Bacterial_reads_cor_selected$gene)]
      Boruta_reads_Fungi_selected_and_significant_cor_gene_id=Boruta_reads_Fungi_selected$gene[(Boruta_reads_Fungi_selected$gene %in% Fungi_reads_cor_selected$gene)]
      
      Boruta_localG_Bacterial_selected_and_significant_cor_gene_id=Boruta_localG_Bacterial_selected$gene[(Boruta_localG_Bacterial_selected$gene %in% Bacterial_localG_cor_selected$gene)]
      Boruta_localG_Fungi_selected_and_significant_cor_gene_id=Boruta_localG_Fungi_selected$gene[(Boruta_localG_Fungi_selected$gene %in% Fungi_localG_cor_selected$gene)]
      
      # keep the global stats
      #QA
      if (NROW(Boruta_reads_Bacterial_data)!=NROW(Boruta_reads_Fungi_data)) {message (paste("[WARN] for",Exp,S,"the number of genes considered for Bacteria is not equal to that considered for Fungi - Boruta reads [",NROW(Boruta_reads_Bacterial_data),"!=",NROW(Boruta_reads_Fungi_data),"]",sep=" "))}
      if (NROW(Boruta_localG_Bacterial_data)!=NROW(Boruta_localG_Fungi_data)) {message (paste("[WARN] for",Exp,S,"the number of genes considered for Bacteria is not equal to that considered for Fungi - Boruta localG [",NROW(Boruta_localG_Bacterial_data),"!=",NROW(Boruta_localG_Fungi_data),"]",sep=" "))}
      
      if (NROW(Bacterial_reads_cor_data)!=NROW(Fungi_reads_cor_data)) {message (paste("[WARN] for",Exp,S,"the number of genes considered for Bacteria is not equal to that considered for Fungi - Spearman reads [",NROW(Bacterial_reads_cor_data),"!=",NROW(Fungi_reads_cor_data),"]",sep=" "))}
      if (NROW(Bacterial_localG_cor_data)!=NROW(Fungi_localG_cor_data)) {message (paste("[WARN] for",Exp,S,"the number of genes considered for Bacteria is not equal to that considered for Fungi - Spearman localG [",NROW(Bacterial_localG_cor_data),"!=",NROW(Fungi_localG_cor_data),"]",sep=" "))}
      
      selected_stat_record=data.frame(experiment=Exp,
                                      sample=S,
                                      tested_genes_read=NROW(Boruta_reads_Bacterial_data),
                                      tested_genes_localG=NROW(Boruta_localG_Bacterial_data),
                                      BORUTA_localG_selected_Bacteria=NROW(Boruta_localG_Bacterial_selected),
                                      BORUTA_localG_selected_Fungi=NROW(Boruta_localG_Fungi_selected),
                                      SPEARMAN_0.01_localG_selected_Bacteria=NROW(Bacterial_localG_cor_selected),
                                      SPEARMAN_0.01_localG_selected_Fungi=NROW(Fungi_localG_cor_selected),
                                      BORUTA_and_SPEARMAN_localG_selected_Bacteria=length(Boruta_localG_Bacterial_selected_and_significant_cor_gene_id),
                                      BORUTA_and_SPEARMAN_locolG_selected_Fungi=length(Boruta_localG_Fungi_selected_and_significant_cor_gene_id),
                                      
                                      BORUTA_reads_selected_Bacteria=NROW(Boruta_reads_Bacterial_selected),
                                      BORUTA_reads_selected_Fungi=NROW(Boruta_reads_Fungi_selected),
                                      SPEARMAN_0.01_reads_selected_Bacteria=NROW(Bacterial_reads_cor_selected),
                                      SPEARMAN_0.01_reads_selected_Fungi=NROW(Fungi_reads_cor_selected),
                                      BORUTA_and_SPEARMAN_reads_selected_Bacteria=length(Boruta_reads_Bacterial_selected_and_significant_cor_gene_id),
                                      BORUTA_and_SPEARMAN_reads_selected_Fungi=length(Boruta_reads_Fungi_selected_and_significant_cor_gene_id),
                                      
                                      # count of the old ones were selected
                                      BORUTA_reads_or_localG_selected_Bacteria=length(unique(c(Boruta_reads_Bacterial_selected$gene,Boruta_localG_Bacterial_selected$gene))),
                                      BORUTA_reads_or_localG_selected_Fungi=length(unique(c(Boruta_reads_Fungi_selected$gene,Boruta_localG_Fungi_selected$gene))),
                                      
                                      # count the new one selected
                                      BORUTA_and_SPEARMAN_reads_or_localG_selected_Bacteria=length(unique(c(Boruta_reads_Bacterial_selected_and_significant_cor_gene_id,Boruta_localG_Bacterial_selected_and_significant_cor_gene_id))),
                                      BORUTA_and_SPEARMAN_reads_or_localG_selected_Fungi=length(unique(c(Boruta_reads_Fungi_selected_and_significant_cor_gene_id,Boruta_localG_Fungi_selected_and_significant_cor_gene_id)))
                                      )
      selected_global_stats=rbind(selected_global_stats,selected_stat_record)
      ### END NEW FOR THIS ANALYSES
      
      # intersects the eight
      all_Bacterial_genes=list(Cor_reads_Bacterial=Bacterial_reads_cor_selected$gene,
                               Cor_localG_Bacterial=Bacterial_localG_cor_selected$gene,
                               Boruta_reads_Bacterial=Boruta_reads_Bacterial_selected$gene,
                               Boruta_localG_Bacterial=Boruta_localG_Bacterial_selected$gene
      )
      all_Bacterial_genes_PA=as.data.frame(list_to_matrix(all_Bacterial_genes))                
      # all_Bacterial_genes_PA$sum=rowSums(all_Bacterial_genes_PA)
      all_Bacterial_genes_PA$Bacteria_CorBoruta_agreement=rowSums(all_Bacterial_genes_PA)
      all_Bacterial_genes_PA_annotated=merge(all_Bacterial_genes_PA,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
      # old
      # intersects_upsetR_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/Intersects_Cor_and_Boruta_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Bacteria_Fungi_sum.",S,".pdf",sep="");
      
      # Apr2022
      intersects_upsetR_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/Intersects_Cor_and_Boruta_reads_or_2x2_localG.Apr2022_220405_semiwild_dataset_rawcounts_filtered.RNA.with_Top50genus_Bacteria_Fungi_sum.",S,".pdf",sep="");
      pdf (intersects_upsetR_file,height = 7,width = 14)
      plot_upSetR(all_Bacterial_genes)
      
      # cp /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/Intersects_Cor_and_Boruta_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Bacteria_Fungi_sum.* ~/Dropbox/PostDoc/Projects/16S_array/Host_response/Intersect_BORUTA_Cor/OMNI13
      # cp /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/Intersects_Cor_and_Boruta_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Bacteria_Fungi_sum.* ~/Dropbox/PostDoc/Projects/16S_array/Host_response/Intersect_BORUTA_Cor/OMNI12        
      
      all_Fungi_genes=list(Cor_reads_Fungi=Fungi_reads_cor_selected$gene,
                           Cor_localG_Fungi=Fungi_localG_cor_selected$gene,
                           Boruta_reads_Fungi=Boruta_reads_Fungi_selected$gene,
                           Boruta_localG_Fungi=Boruta_localG_Fungi_selected$gene
                          )
      all_Fungi_genes_PA=as.data.frame(list_to_matrix(all_Fungi_genes))                
      # all_Fungi_genes_PA$sum=rowSums(all_Fungi_genes_PA)
      all_Fungi_genes_PA$Fungi_CorBoruta_agreement=rowSums(all_Fungi_genes_PA)
      all_Fungi_genes_PA_annotated=merge(all_Fungi_genes_PA,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
      
      plot_upSetR(all_Fungi_genes)
      
      
      all_genes=c(all_Bacterial_genes,all_Fungi_genes)
      all_genes_PA=as.data.frame(list_to_matrix(all_genes))                
      # all_genes_PA$sum=rowSums(all_genes_PA)
      all_genes_PA$Bacteria_Fungi_CorBoruta_agreement=rowSums(all_genes_PA)
      all_genes_PA_annotated=merge(all_genes_PA,Sprot_annotation,by.x="row.names",by.y="Cross.reference..Araport.",all.x=T)
      
      plot_upSetR(all_genes)
      
      dev.off()
      
      # The final genes for the section -> Boruta_reads or Boruta_localG ONLY IF THEIR COR VALUE SIGNIFICANT
      # Fungi_genes_names=unique(c(Boruta_reads_Fungi_selected$gene,Boruta_localG_Fungi_selected$gene)) # The old not considering the cor value
      Fungi_genes_names=unique(c(Boruta_reads_Fungi_selected_and_significant_cor_gene_id,Boruta_localG_Fungi_selected_and_significant_cor_gene_id))
      #Fungi_genes_final_selected=all_genes_PA_annotated[all_genes_PA_annotated$Row.names %in% Fungi_genes_names,]
      Fungi_genes_final_selected=all_Fungi_genes_PA_annotated[all_Fungi_genes_PA_annotated$Row.names %in% Fungi_genes_names,]
      
      names(Fungi_genes_final_selected)[1]=c("gene")
      
      # Bacterial_genes_names=unique(c(Boruta_reads_Bacterial_selected$gene,Boruta_localG_Bacterial_selected$gene)) # The old not considering the cor value
      Bacterial_genes_names=unique(c(Boruta_reads_Bacterial_selected_and_significant_cor_gene_id,Boruta_localG_Bacterial_selected_and_significant_cor_gene_id))
      # Bacterial_genes_final_selected=all_genes_PA_annotated[all_genes_PA_annotated$Row.names %in% Bacterial_genes_names,]
      Bacterial_genes_final_selected=all_Bacterial_genes_PA_annotated[all_Bacterial_genes_PA_annotated$Row.names %in% Bacterial_genes_names,]
      names(Bacterial_genes_final_selected)[1]=c("gene")
      
      
      # Add the correlations values
      names(Bacterial_reads_cor_data)[3:9]=paste(names(Bacterial_reads_cor_data)[3:9],"_Reads",sep="")
      Bacterial_genes_final_selected=merge(Bacterial_genes_final_selected,Bacterial_reads_cor_data[,c("gene","gene_sum_Reads","i_gene_Reads","r_Reads","FDR.BH_Reads")],by="gene",all.x = T)
      
      names(Bacterial_localG_cor_data)[3:9]=paste(names(Bacterial_localG_cor_data)[3:9],"_LocalG",sep="")
      Bacterial_genes_final_selected=merge(Bacterial_genes_final_selected,Bacterial_localG_cor_data[,c("gene","gene_sum_LocalG","i_gene_LocalG","r_LocalG","FDR.BH_LocalG")],by="gene",all.x = T)
      
      names(Fungi_reads_cor_data)[3:9]=paste(names(Fungi_reads_cor_data)[3:9],"_Reads",sep="")
      Fungi_genes_final_selected=merge(Fungi_genes_final_selected,Fungi_reads_cor_data[,c("gene","gene_sum_Reads","i_gene_Reads","r_Reads","FDR.BH_Reads")],by="gene",all.x = T)
      
      names(Fungi_localG_cor_data)[3:9]=paste(names(Fungi_localG_cor_data)[3:9],"_LocalG",sep="")
      Fungi_genes_final_selected=merge(Fungi_genes_final_selected,Fungi_localG_cor_data[,c("gene","gene_sum_LocalG","i_gene_LocalG","r_LocalG","FDR.BH_LocalG")],by="gene",all.x = T)
      
      # save the final files
      # old
      # Boruta_based_selected_Bacteria_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Bacteria_sum.",S,".tdl",sep="");
      # Boruta_based_selected_Bacteria_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Bacteria_sum.",S,".tdl",sep="");
      # Apr2022
      Boruta_based_selected_Bacteria_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.Apr2022_220405_semiwild_dataset_rawcounts_filtered.RNA.with_Top50genus_Bacteria_sum.",S,".tdl",sep="");
      write.table(x=Bacterial_genes_final_selected,file = Boruta_based_selected_Bacteria_file,quote = F,sep = "\t",row.names = F)
      
      # old
      # Boruta_based_selected_Fungi_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Fungi_sum.",S,".tdl",sep="");
      # Boruta_based_selected_Fungi_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Fungi_sum.",S,".tdl",sep="");
      # Apr2022
      Boruta_based_selected_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.Apr2022_220405_semiwild_dataset_rawcounts_filtered.RNA.with_Top50genus_Fungi_sum.",S,".tdl",sep="");
      write.table(x=Fungi_genes_final_selected,file = Boruta_based_selected_Fungi_file,quote = F,sep = "\t",row.names = F)
      
      
      # Agregate compact
      unique_sample_name=paste(Exp,S,sep="_")
      
      Bacteria_df_to_add=data.frame(row.names=Bacterial_genes_final_selected$gene,
                                    V1=rowSums(Bacterial_genes_final_selected[,grepl(x=names(Bacterial_genes_final_selected),pattern = "Bacterial",fixed = T)]))
      names(Bacteria_df_to_add)=paste(unique_sample_name,"Bacteria_CorBoruta_agreement",sep=".")
      
      Fungi_df_to_add=data.frame(row.names=Fungi_genes_final_selected$gene,
                                 V1=rowSums(Fungi_genes_final_selected[,grepl(x=names(Fungi_genes_final_selected),pattern = "Fungi",fixed = T)]))
      names(Fungi_df_to_add)=paste(unique_sample_name,"Fungi_CorBoruta_agreement",sep=".")
      
      if (nrow(selected_genes_from_all_sections)==0) 
      {
        selected_genes_from_all_sections=Bacteria_df_to_add
        selected_genes_from_all_sections=merge(selected_genes_from_all_sections,Fungi_df_to_add,by="row.names",all=T)
      } else {
        selected_genes_from_all_sections=merge(selected_genes_from_all_sections,Bacteria_df_to_add,by.x="Row.names",by.y="row.names",all=T)
        selected_genes_from_all_sections=merge(selected_genes_from_all_sections,Fungi_df_to_add,by.x="Row.names",by.y="row.names",all=T)
      }
      rm (Bacteria_df_to_add,Fungi_df_to_add)
      # agregate_detailed
      #Fungi_corr_df_to_add=Fungi_genes_final_selected[,c("gene","sum","gene_sum_Reads","i_gene_Reads","r_Reads","FDR.BH_Reads",
      #                                                   "gene_sum_LocalG","i_gene_LocalG","r_LocalG","FDR.BH_LocalG")]
      Fungi_corr_df_to_add=Fungi_genes_final_selected[,c("gene","Fungi_CorBoruta_agreement","gene_sum_Reads","i_gene_Reads","r_Reads","FDR.BH_Reads",
                                                         "gene_sum_LocalG","i_gene_LocalG","r_LocalG","FDR.BH_LocalG")]
      names(Fungi_corr_df_to_add)[2:NCOL(Fungi_corr_df_to_add)]=paste(unique_sample_name,names(Fungi_corr_df_to_add)[2:NCOL(Fungi_corr_df_to_add)],sep=".")
      names(Fungi_corr_df_to_add)[3:NCOL(Fungi_corr_df_to_add)]=paste(names(Fungi_corr_df_to_add)[3:NCOL(Fungi_corr_df_to_add)],"Fungi",sep=".")
      
      Fungi_corr_df_to_add[[paste(unique_sample_name,"selected_Fungi",sep=".")]]=1
 
      # Bacteria_corr_df_to_add=Bacterial_genes_final_selected[,c("gene","sum","gene_sum_Reads","i_gene_Reads","r_Reads","FDR.BH_Reads",
      #                                                          "gene_sum_LocalG","i_gene_LocalG","r_LocalG","FDR.BH_LocalG")]
      Bacteria_corr_df_to_add=Bacterial_genes_final_selected[,c("gene","Bacteria_CorBoruta_agreement","gene_sum_Reads","i_gene_Reads","r_Reads","FDR.BH_Reads",
                                                                "gene_sum_LocalG","i_gene_LocalG","r_LocalG","FDR.BH_LocalG")]
      names(Bacteria_corr_df_to_add)[2:NCOL(Bacteria_corr_df_to_add)]=paste(unique_sample_name,names(Bacteria_corr_df_to_add)[2:NCOL(Bacteria_corr_df_to_add)],sep=".")
      names(Bacteria_corr_df_to_add)[3:NCOL(Bacteria_corr_df_to_add)]=paste(names(Bacteria_corr_df_to_add)[3:NCOL(Bacteria_corr_df_to_add)],"Bacteria",sep=".")
      
      Bacteria_corr_df_to_add[[paste(unique_sample_name,"selected_Bacteria",sep=".")]]=1
                                 
      # add tp aggregate
      if (nrow(selected_genes_from_all_sections_with_all_details)==0) 
      {
        selected_genes_from_all_sections_with_all_details=Bacteria_corr_df_to_add
        selected_genes_from_all_sections_with_all_details=merge(selected_genes_from_all_sections_with_all_details,Fungi_corr_df_to_add,by="gene",all=T)
      } else {
        selected_genes_from_all_sections_with_all_details=merge(selected_genes_from_all_sections_with_all_details,Bacteria_corr_df_to_add,by="gene",all=T)
        selected_genes_from_all_sections_with_all_details=merge(selected_genes_from_all_sections_with_all_details,Fungi_corr_df_to_add,by="gene",all=T)
      }
      rm (Bacteria_corr_df_to_add,Fungi_corr_df_to_add)
      
      rm (Boruta_reads_Bacterial_selected_and_significant_cor_gene_id,Boruta_reads_Fungi_selected_and_significant_cor_gene_id,
          Boruta_localG_Bacterial_selected_and_significant_cor_gene_id, Boruta_localG_Fungi_selected_and_significant_cor_gene_id,
          selected_stat_record) # created in the new analyses
      
      rm ("all_Bacterial_genes","all_Bacterial_genes_PA","all_Bacterial_genes_PA_annotated","all_Fungi_genes","all_Fungi_genes_PA",
          "all_Fungi_genes_PA_annotated","all_genes","all_genes_PA","all_genes_PA_annotated",
          "Bacteria_localG_cor_file","Bacteria_reads_cor_file","Bacterial_genes_final_selected","Bacterial_genes_names","Bacterial_localG_cor_data",
          "Bacterial_localG_cor_selected","Bacterial_reads_cor_data","Bacterial_reads_cor_selected","Boruta_based_selected_Bacteria_file","Boruta_based_selected_Fungi_file",  
          "Boruta_localG_Bacteria_File","Boruta_localG_Bacterial_data","Boruta_localG_Bacterial_selected","Boruta_localG_Fungi_data","Boruta_localG_Fungi_file",         
          "Boruta_localG_Fungi_selected","Boruta_reads_Bacteria_File","Boruta_reads_Bacterial_data","Boruta_reads_Bacterial_selected","Boruta_reads_Fungi_data",           
          "Boruta_reads_Fungi_file","Boruta_reads_Fungi_selected","Fungi_genes_final_selected","Fungi_genes_names","Fungi_localG_cor_data","Fungi_localG_cor_file",
          "Fungi_localG_cor_selected","Fungi_reads_cor_data","Fungi_reads_cor_file","Fungi_reads_cor_selected","intersects_upsetR_file","signif_tresh","unique_sample_name")
    } # if exists closure
  } # Sample closure
} # Exp closure

# write global stats of selected numbers
# old
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_stats_per_section_Bacteria_Fungi_WO_HS_separation.txt"
# Apr2022

# Add the leaf section labels
samples_label=data.frame(sample= c("OMNI12_A1","OMNI12_A2","OMNI12_B1","OMNI12_B2","OMNI12_C1","OMNI12_C2",
                                   "OMNI13_A1","OMNI13_A2","OMNI13_B1","OMNI13_B2","OMNI13_C1","OMNI13_C2","OMNI13_D1"))

P1L1=c("OMNI12_A1","OMNI12_A2","OMNI12_B2")
P1L2=c("OMNI12_B1","OMNI12_C1","OMNI12_C2")
P2L1=c("OMNI13_A1","OMNI13_A2","OMNI13_B1")
P2L2=c("OMNI13_B2","OMNI13_C1","OMNI13_C2","OMNI13_D1")

samples_label$PlantLeaf[samples_label$sample %in% P1L1]="P1.L1"
samples_label$PlantLeaf[samples_label$sample %in% P1L2]="P1.L2"
samples_label$PlantLeaf[samples_label$sample %in% P2L1]="P2.L1"
samples_label$PlantLeaf[samples_label$sample %in% P2L2]="P2.L2"

section1=c("OMNI12_A1","OMNI12_B1","OMNI13_A1","OMNI13_D1")
section2=c("OMNI12_A2","OMNI12_C1","OMNI13_C1")
section3=c("OMNI13_A2","OMNI13_B2")
section4=c("OMNI12_B2","OMNI12_C2","OMNI13_B1","OMNI13_C2")

samples_label$Section[samples_label$sample %in% section1]="1"
samples_label$Section[samples_label$sample %in% section2]="2"
samples_label$Section[samples_label$sample %in% section3]="3"
samples_label$Section[samples_label$sample %in% section4]="4"

samples_label$lable=paste(samples_label$PlantLeaf,samples_label$Section,sep=".")
names(selected_global_stats)[2]="SectionSample"
selected_global_stats$sample=paste(selected_global_stats$experiment,selected_global_stats$SectionSample,sep="_")
selected_global_stats=merge(x = selected_global_stats,y = samples_label,by="sample")
file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_selected_stats_per_section_Bacteria_Fungi_Apr2022_220405_semiwild_dataset_rawcounts_filtered.WO_HS_separation.txt"
write.table(x=selected_global_stats,file = file,quote = F,sep = "\t",row.names = F)

# Sum data per plant leaf and do some intersetcs ploting
# pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Whole_array_gene_selected_by_BORUTA_per_leaf_upSetR.pdf")

# old
# pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Whole_array_gene_selected_by_BORUTA_AND_significant_0.01_spearman_per_leaf_upSetR.pdf")
# Apr2022

pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/Whole_array_Apr2022_220405_semiwild_dataset_rawcounts_filtered_gene_selected_by_BORUTA_AND_significant_0.01_spearman_per_leaf_upSetR.pdf")

# do the per plant section counting
for (type in c("Bacteria","Fungi")) {
  Pleaf=c()
  # t=names(selected_genes_from_all_sections_with_all_details)[grepl(pattern = "selected",x=names(selected_genes_from_all_sections_with_all_details))]
  if (type=="Bacteria") {Pleaf=c("OMNI12_A1.selected_Bacteria|OMNI12_A2.selected_Bacteria|OMNI12_B2.selected_Bacteria",
                                 "OMNI12_B1.selected_Bacteria|OMNI12_C1.selected_Bacteria|OMNI12_C2.selected_Bacteria",
                                 "OMNI13_A1.selected_Bacteria|OMNI13_A2.selected_Bacteria|OMNI13_B1.selected_Bacteria",
                                 "OMNI13_B2.selected_Bacteria|OMNI13_C1.selected_Bacteria|OMNI13_C2.selected_Bacteria|OMNI13_D1.selected_Bacteria")}
  if (type=="Fungi") {Pleaf=c("OMNI12_A1.selected_Fungi|OMNI12_A2.selected_Fungi|OMNI12_B2.selected_Fungi",
                              "OMNI12_B1.selected_Fungi|OMNI12_C1.selected_Fungi|OMNI12_C2.selected_Fungi",
                              "OMNI13_A1.selected_Fungi|OMNI13_A2.selected_Fungi|OMNI13_B1.selected_Fungi",
                              "OMNI13_B2.selected_Fungi|OMNI13_C1.selected_Fungi|OMNI13_C2.selected_Fungi|OMNI13_D1.selected_Fungi")}
  for (i in 1:4)
  {
    # label=paste("PL",i,type,sep="_")
    label="NA"
    if (i==1) {label=paste("P1.L1")}
    if (i==2) {label=paste("P1.L2")}
    if (i==3) {label=paste("P2.L1")}
    if (i==4) {label=paste("P2.L2")}
    label=paste(label,type,sep=".")
    PL_elected_data=selected_genes_from_all_sections_with_all_details[,grepl(x=names(selected_genes_from_all_sections_with_all_details),pattern = Pleaf[i])]
    selected_genes_from_all_sections_with_all_details[[label]]=rowSums(PL_elected_data,na.rm = T)
    PL_elected_data=cbind(selected_genes_from_all_sections_with_all_details[,1],PL_elected_data)
    names(PL_elected_data)[1]="gene"
    # filter all NA lines
    PL_elected_data_nonNA=PL_elected_data[rowSums(PL_elected_data[,2:NCOL(PL_elected_data)],na.rm = T)>0,]
    PL_elected_data_nonNA[is.na(PL_elected_data_nonNA)]=0
    names(PL_elected_data_nonNA)=gsub(x=names(PL_elected_data_nonNA),pattern = paste(".selected",type,sep="_"),replacement = "")
    
    # replace the name with the correct label
    for (j in 2:length(names(PL_elected_data_nonNA)))
    {
      sample_name=names(PL_elected_data_nonNA)[j]
      label1=samples_label$lable[samples_label$sample==sample_name]
      message (paste ("==",sample_name,label))
      names(PL_elected_data_nonNA)[j]=label1
    }
    # plot_upSetR(PL_elected_data_nonNA,paste ("Selected",type,"leaf",i))
    
    # save each plot separetley to make it easier to join in inkScape
    # pdf(paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/UpSetR_sep/Whole_array_Apr2022_220405_semiwild_dataset_rawcounts_filtered_gene_selected_by_BORUTA_AND_significant_0.01_spearman_per_leaf_upSetR",".",label,".pdf",sep=""),width=4,height = 3)
    # plot_upSetR(PL_elected_data_nonNA,"")
    # dev.off()
    
    plot_upSetR(PL_elected_data_nonNA,paste ("Selected",label))
    rm(PL_elected_data_nonNA,PL_elected_data)
  }
}

# intersect per leaf [1-4] selected Bacteria/Fungi
all_selected=data.frame()
figures_list=list()
for (i in (1:4))
{
  label="NA"
  if (i==1) {label=paste("P1.L1")}
  if (i==2) {label=paste("P1.L2")}
  if (i==3) {label=paste("P2.L1")}
  if (i==4) {label=paste("P2.L2")}
  
  selected_per_leaf=data.frame()
  for (type in c("Bacteria","Fungi"))
  {
    label1=paste(label,type,sep=".")
    # tmp=selected_genes_from_all_sections_with_all_details[,c("gene",paste("PL",i,type,sep="_"))]
    tmp=selected_genes_from_all_sections_with_all_details[,c("gene",label1)]
    if (NROW(selected_per_leaf)==0)
    {
      selected_per_leaf=tmp[tmp[[2]]>=2,]
    } else {
      selected_per_leaf=merge(selected_per_leaf,tmp[tmp[[2]]>=2,],by="gene",all=T)
    }
  }
  selected_per_leaf[is.na(selected_per_leaf)]=0
  selected_per_leaf[[2]][selected_per_leaf[[2]]>0]=1
  selected_per_leaf[[3]][selected_per_leaf[[3]]>0]=1
  
  # plot_upSetR(selected_per_leaf,paste("Selected leaf",i))
  plot_upSetR(elements_list = selected_per_leaf,title = label)
  if (NROW(all_selected)==0)
  {
    all_selected=selected_per_leaf
  } else {
    all_selected=merge(all_selected,selected_per_leaf,by="gene",all=T)
  }
  rm (selected_per_leaf)
}
dev.off()
# old
# pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Whole_array_gene_selected_by_BORUTA_per_leaf_upSetR_intersect_all_selected.pdf",width = 18,height = 10)
# pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Whole_array_gene_selected_by_BORUTA_AND_significant_0.01_spearman_per_leaf_upSetR_intersect_all_selected.pdf",width = 18,height = 10)

# Apr2022
# pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/Whole_array_Apr2022_220405_semiwild_dataset_rawcounts_filtered_gene_selected_by_BORUTA_AND_significant_0.01_spearman_per_leaf_upSetR_intersect_all_selected.pdf",width = 18,height = 10)
pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/Whole_array_Apr2022_220405_semiwild_dataset_rawcounts_filtered_gene_selected_by_BORUTA_AND_significant_0.01_spearman_per_leaf_upSetR_intersect_all_selected.A4R.pdf",width=11,height = 6)
all_selected[is.na(all_selected)]=0
plot_upSetR(all_selected)
dev.off()

# IMPORTANT NEED TO KEEP THE NAMES AND LABELS COORDINATED
# intersect l1-l4 - selected Fungi; selected Bacteria - (1) per leaf (2) across all
selected_genes_from_all_sections_with_all_details$sum_sections_Bacteria=rowSums(selected_genes_from_all_sections_with_all_details[,grepl(x=names(selected_genes_from_all_sections_with_all_details),pattern = "PL_\\d_Bacteria",perl = T)],na.rm = T)
selected_genes_from_all_sections_with_all_details$sum_sections_Fungi=rowSums(selected_genes_from_all_sections_with_all_details[,grepl(x=names(selected_genes_from_all_sections_with_all_details),pattern = "PL_\\d_Fungi",perl = T)],na.rm = T)
selected_genes_from_all_sections_with_all_details$Bacterial_response=rowSums(selected_genes_from_all_sections_with_all_details[,grepl(x=names(selected_genes_from_all_sections_with_all_details),
                                                                                                                                  pattern = "PL_1_Bacteria|PL_2_Bacteria|PL_3_Bacteria|PL_4_Bacteria")]>=2,na.rm = T)
selected_genes_from_all_sections_with_all_details$Fungi_response=rowSums(selected_genes_from_all_sections_with_all_details[,grepl(x=names(selected_genes_from_all_sections_with_all_details),
                                                                                                                                      pattern = "PL_1_Fungi|PL_2_Fungi|PL_3_Fungi|PL_4_Fungi")]>=2,na.rm = T)
selected_genes_from_all_sections_with_all_details$response_type=NA
selected_genes_from_all_sections_with_all_details$response_type[selected_genes_from_all_sections_with_all_details$Fungi_response>0&
                                                                selected_genes_from_all_sections_with_all_details$Bacterial_response==0]="Fungi_Unique?"
selected_genes_from_all_sections_with_all_details$response_type[selected_genes_from_all_sections_with_all_details$Bacterial_response>0&
                                                                selected_genes_from_all_sections_with_all_details$Fungi_response==0]="Bacteria_Unique?"
selected_genes_from_all_sections_with_all_details$response_type[selected_genes_from_all_sections_with_all_details$Fungi_response>0&
                                                                  selected_genes_from_all_sections_with_all_details$sum_sections_Bacteria==0]="Fungi_Unique?"
selected_genes_from_all_sections_with_all_details$response_type[selected_genes_from_all_sections_with_all_details$Bacterial_response>0&
                                                                  selected_genes_from_all_sections_with_all_details$sum_sections_Fungi==0]="Bacteria_Unique"

selected_genes_from_all_sections_with_all_details$response_type[selected_genes_from_all_sections_with_all_details$Bacterial_response>0&
                                                                  selected_genes_from_all_sections_with_all_details$Fungi_response>0]="Shared"
# add annotation to the agregate
selected_genes_from_all_sections_annotated=merge(selected_genes_from_all_sections,Sprot_annotation,by.x="Row.names",by.y="Cross.reference..Araport.",all.x=T)
selected_genes_from_all_sections_with_all_details_annotated=merge(selected_genes_from_all_sections_with_all_details,Sprot_annotation,by.x="gene",by.y="Cross.reference..Araport.",all.x=T)



# save

# tmp_file="~/all_selected_by_Boruta_all_sections_ell_exp_Bacteria_Fungi.txt"
#file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.txt"
# old
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.txt"
# Apr2022
file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_Apr2022_220405_semiwild_dataset_rawcounts_filtered.WO_HS_separation.txt"
write.table(x=selected_genes_from_all_sections_annotated,file = file,quote = F,sep = "\t",row.names = F)

# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_by_plant_Bacteria_Fungi_WO_HS_separation.txt"
# old
file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_AND_significant_0.01_spearman_by_plant_Bacteria_Fungi_WO_HS_separation.txt"
# Apr2022
file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_selected_by_Boruta_AND_significant_0.01_spearman_by_plant_Bacteria_Fungi_Apr2022_220405_semiwild_dataset_rawcounts_filtered.WO_HS_separation.txt"
write.table(x=all_selected,file = file,quote = F,sep = "\t",row.names = F)

selected_names_and_type=selected_genes_from_all_sections_with_all_details_annotated[!is.na(selected_genes_from_all_sections_with_all_details_annotated$response_type),c("gene","response_type")]
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_Bacteria_Fungi_and_type_WO_HS_separation.txt"
# old 
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_AND_significant_0.01_spearman_Bacteria_Fungi_and_type_WO_HS_separation.txt"
# Apr2022
file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_selected_by_Boruta_AND_significant_0.01_spearman_Bacteria_Fungi_and_type.Apr2022_220405_semiwild_dataset_rawcounts_filtered.WO_HS_separation.txt"
write.table(x=selected_names_and_type,file = file,quote = F,sep = "\t",row.names = F)

# add median r per plant per type
library(matrixStats)
for (r_type in c("r_Reads.Fungi","r_LocalG.Fungi","r_Reads.Bacteria","r_LocalG.Bacteria"))
{
  r_data=selected_genes_from_all_sections_with_all_details_annotated[,grepl(x=names(selected_genes_from_all_sections_with_all_details_annotated),pattern = r_type)]
  r_matrix=as.matrix(r_data)
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".all_Median",sep="")]]=rowMedians(x = r_matrix,na.rm = T)
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".all_Mean",sep="")]]=rowMeans2(r_matrix,na.rm = T)
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".all_std",sep="")]]=rowSds(r_matrix,na.rm = T)
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".total_pos_r",sep="")]]=rowSums(as.matrix(r_data>0),na.rm = T)
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".total_neg_r",sep="")]]=rowSums(as.matrix(r_data<0),na.rm = T)
  
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".max_value",sep="")]]=do.call(pmax, c(r_data, na.rm=TRUE))
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".min_value",sep="")]]=do.call(pmin, c(r_data, na.rm=TRUE))
  
  r_matrix[is.na(r_matrix)]=-Inf
  max_col=max.col(r_matrix)
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".max_sample",sep="")]]=gsub(x = colnames(r_matrix)[max_col],pattern = paste(".",r_type,sep=""),replacement = "")
  
  r_matrix[is.infinite(r_matrix)]=Inf  
  min_col=max.col(r_matrix*-1)
  selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".min_sample",sep="")]]=gsub(x = colnames(r_matrix)[min_col],pattern = paste(".",r_type,sep=""),replacement = "")

  rm (r_matrix,r_data,max_col,min_col)
}

# QA
# for (r_type in c("r_Reads.Fungi","r_LocalG.Fungi","r_Reads.Bacteria","r_LocalG.Bacteria"))
# {
#   r_data.tmp=selected_genes_from_all_sections_with_all_details_annotated[,grepl(x=names(selected_genes_from_all_sections_with_all_details_annotated),pattern = r_type)]
#   r_data.tmp1=cbind(selected_genes_from_all_sections_with_all_details_annotated$gene,r_data.tmp)
#   r_data.tmp2=r_data.tmp[,1:13]
#   r_matrix=as.matrix(r_data.tmp2)
#   selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".all_Median",sep="")]]=rowMedians(x = r_matrix,na.rm = T)
#   selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".all_Mean",sep="")]]=rowMeans2(r_matrix,na.rm = T)
#   selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".all_std",sep="")]]=rowSds(r_matrix,na.rm = T)
#   selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".total_pos_r",sep="")]]=rowSums(as.matrix(r_data>0),na.rm = T)
#   selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".total_neg_r",sep="")]]=rowSums(as.matrix(r_data<0),na.rm = T)
#   
#   tmp_max=do.call(pmax, c(r_data.tmp2, na.rm=TRUE))
#   tmp_max=cbind(selected_genes_from_all_sections_with_all_details_annotated$gene,tmp_max)
#   selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".max_value",sep="")]]=do.call(pmax, c(r_data, na.rm=TRUE))
#   selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".min_value",sep="")]]=do.call(pmin, c(r_data, na.rm=TRUE))
#   
#   r_matrix[is.na(r_matrix)]=-Inf
#   error_example=r_matrix[c(1030,1031),] # AT4G20260 = 1030
#   max_col=max.col(error_example)
#   colnames(error_example)[max_col]
#   # selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".max_sample",sep="")]]=gsub(x = colnames(r_matrix)[max_col],pattern = paste(".",r_type,sep=""),replacement = "")
#   # min_col=max.col(r_matrix*-1)
#   
#   error_example[is.infinite(error_example)]=Inf  
#   min_col=max.col(error_example*-1)
#   colnames(error_example)[min_col]
#   
#   # selected_genes_from_all_sections_with_all_details_annotated[[paste(r_type,".min_sample",sep="")]]=gsub(x = colnames(r_matrix)[min_col],pattern = paste(".",r_type,sep=""),replacement = "")
#   
#   rm (r_matrix,r_data,max_col,min_col)
# }


# add the agreement info to the selected details
# Seems not to add much... Validate for Apr2022 version
# selected_genes_from_all_sections_with_all_details_annotated_extended=merge(selected_genes_from_all_sections_with_all_details_annotated,selected_genes_from_all_sections,by.y="Row.names",by.x="gene",all=T)

# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.detailed.with_agreement_info.txt"

# old
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.detailed.with_agreement_info.txt"
# Apr2022 
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_Apr2022_220405_semiwild_dataset_rawcounts_filtered.WO_HS_separation.detailed.with_agreement_info.txt"
# write.table(x=selected_genes_from_all_sections_with_all_details_annotated_extended,file = file,quote = F,sep = "\t",row.names = F)

# old
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.detailed.txt"
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.detailed.txt"
# Apr2022
file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_Apr2022_220405_semiwild_dataset_rawcounts_filtered.WO_HS_separation.detailed.txt"

write.table(x=selected_genes_from_all_sections_with_all_details_annotated,file = file,quote = F,sep = "\t",row.names = F)

# make another df which we handle for ploting (i.e ignoring those selected only in few sections but not in 2/3 of a leaf)
final_selected_genes_from_all_sections_with_all_details_annotated=selected_genes_from_all_sections_with_all_details_annotated[!is.na(selected_genes_from_all_sections_with_all_details_annotated$response_type),]

# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/FINAL_selected_by_Boruta_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.detailed.txt"
# file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/FINAL_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.detailed.txt"
# Apr2022
file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/FINAL_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_Apr2022_220405_semiwild_dataset_rawcounts_filtered.WO_HS_separation.detailed.txt"
write.table(x=final_selected_genes_from_all_sections_with_all_details_annotated,file = file,quote = F,sep = "\t",row.names = F)


# Compare with old expression results
old_results=read.delim(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/FINAL_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.detailed.txt",stringsAsFactors = F,sep="\t")
compare_old_new_results=merge(data.frame(gene=old_results$gene,old=1),data.frame(gene=final_selected_genes_from_all_sections_with_all_details_annotated$gene,new=1,final_selected_genes_from_all_sections_with_all_details_annotated[,288:323]),by="gene",all=T)
compare_old_new_results[is.na(compare_old_new_results)]=0
compare_old_new_results=merge(compare_old_new_results,Sprot_annotation,by.x="gene",by.y="Cross.reference..Araport.",all.x=T)

message(paste("== newly selected genes: ",sum(compare_old_new_results$old==0&compare_old_new_results$new==1),sep=""))
message(paste("== not selected genes: ",sum(compare_old_new_results$old==1&compare_old_new_results$new==0),sep=""))
not_selected_now=compare_old_new_results[compare_old_new_results$old==1&compare_old_new_results$new==0,]
newly_selected_now=compare_old_new_results[compare_old_new_results$old==0&compare_old_new_results$new==1,]


### DO THE GO ENRICHMENT Apr2022
#### Plot GO enrichment - Shared genes
GO_enrichment_level0_Shared_file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/shared_response_selected_genes_Apr2022_GO_PANTHER.level0.txt"
GO_enrichment_level0_Shared=read.delim(file = GO_enrichment_level0_Shared_file,sep = "\t",header = T,stringsAsFactors = F)
GO_enrichment_level0_Shared_sorted=GO_enrichment_level0_Shared[with(GO_enrichment_level0_Shared, order(number_in_list, fold_enrichment,decreasing = T)),]

library(ggplot2)
pdf("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/shared_response_selected_genes_Apr2022_GO_PANTHER.level0.Top25_byN.pdf",width = 10)
ggplot(GO_enrichment_level0_Shared_sorted[1:25,], aes(x=reorder(title,fold_enrichment), y=fold_enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",number_in_list,sep="")), nudge_y=4, color="black",size=3) + 
  theme(axis.text = element_text(size = 13)) +
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

pdf("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/shared_response_selected_genes_Apr2022_GO_PANTHER.level0.all_terms.pdf",width = 10,height = 10)
ggplot(GO_enrichment_level0_Shared_sorted, aes(x=reorder(title,fold_enrichment), y=fold_enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",number_in_list,sep="")), nudge_y=4, color="black",size=3) + 
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

# Plot GO enrichment - all genes
GO_enrichment_level0_all_file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/all_selected_genes_Apr2022_GO_PANTHER.level0.txt"
GO_enrichment_level0_all=read.delim(file = GO_enrichment_level0_all_file,sep = "\t",header = T,stringsAsFactors = F)
GO_enrichment_level0_all_sorted=GO_enrichment_level0_all[with(GO_enrichment_level0_all, order(number_in_list, fold_enrichment,decreasing = T)),]

library(ggplot2)
pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/all_selected_genes_Apr2022_GO_PANTHER.level0.Top29_byN.pdf",width = 10)
# ggplot(GO_enrichment_level0_all_sorted[1:25,], aes(x=reorder(title,number_in_list), y=fold_enrichment,fill=FDR)) + 
sum(GO_enrichment_level0_all_sorted$number_in_list>=10)
# 29
ggplot(GO_enrichment_level0_all_sorted[1:29,], aes(x=reorder(title,fold_enrichment), y=fold_enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",number_in_list,sep="")), nudge_y=3, color="black",size=3) + 
  theme(axis.text = element_text(size = 13)) +
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/all_selected_genes_Apr2022_GO_PANTHER.level0.all_terms.pdf",width = 10,height = 15)
ggplot(GO_enrichment_level0_all_sorted, aes(x=reorder(title,fold_enrichment), y=fold_enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",number_in_list,sep="")), nudge_y=4, color="black",size=3) + 
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()


### GO enrichment based on DAVID
GO_enrichment_DAVID_file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/DAVID/all_selected_genes_Apr2022_David_Functional_Annotation_Chart.GO_BP.txt"
GO_enrichment_DAVID=read.delim(file = GO_enrichment_DAVID_file,sep = "\t",header = T,stringsAsFactors = F)
GO_enrichment_DAVID_sorted=GO_enrichment_DAVID[with(GO_enrichment_DAVID, order(Count, Fold.Enrichment,decreasing = T)),]

library(ggplot2)
pdf("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/DAVID/all_selected_genes_Apr2022_David_Functional_Annotation_Chart.GO_BP.N_grteq20.pdf",width = 10)
# ggplot(GO_enrichment_level0_all_sorted[1:25,], aes(x=reorder(title,number_in_list), y=fold_enrichment,fill=FDR)) + 
sum(GO_enrichment_DAVID_sorted$Count>=20&GO_enrichment_DAVID_sorted$FDR<=0.05)
# 15
ggplot(GO_enrichment_DAVID_sorted[GO_enrichment_DAVID_sorted$Count>=20&GO_enrichment_DAVID_sorted$FDR<=0.05,], aes(x=reorder(Term,Fold.Enrichment), y=Fold.Enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",Count,sep="")), nudge_y=3, color="black",size=3) + 
  theme(axis.text = element_text(size = 13)) +
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

# All significant
pdf("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/GO/DAVID/all_selected_genes_Apr2022_David_Functional_Annotation_Chart.GO_BP.pdf",width = 10,height = 10)
ggplot(GO_enrichment_DAVID_sorted[GO_enrichment_DAVID_sorted$FDR<=0.05,], aes(x=reorder(Term,Fold.Enrichment), y=Fold.Enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",Count,sep="")), nudge_y=4, color="black",size=3) + 
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

##### Do some ploting...
# new test based on: https://michaelminn.net/tutorials/r-point-analysis/index.html
library(sp)
all_samples_spatial_data=readRDS("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.and_Top50_Bacteria_Fungi_Sum.spatial_objects_list.rds")
# START OLD SECTION
# # read all data files to a list of spatial object indexed by their sample names
# all_samples_spatial_data=list()
# message ("## Read all data into spatial objects")
# for (Exp in c("OMNI12","OMNI13"))
# {
#   samples=c("A1","A2","B1","B2","C1","C2")
#   if (Exp=="OMNI13") {samples=c(samples,"D1")}
#   for (smpl in samples)
#   {
#     sample_name=paste(Exp,smpl,sep="_")
#     message(paste("--",sample_name))
#     hotspost_df=data.frame()
#     data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
#     data=read.delim(data_file,sep=";",stringsAsFactors = F)
#     data$all_Bacteria=rowSums(data[,1:50])
#     data$all_Fungi=rowSums(data[,51:100])
#     
#     xy_list=strsplit(row.names(data), "x")
#     xy=as.data.frame(t(as.data.frame(xy_list)))
#     xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
#     names(xy)=c("x","y")
#     data_with_xy=cbind(xy,data)
#     rownames(data_with_xy)=rownames(data)
#     data_with_xy$x=as.numeric(data_with_xy$x)
#     data_with_xy$y=as.numeric(data_with_xy$y)
#     
#     # Now with the expression
#     expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",Exp,".",smpl,".tsv",sep="")
#     expression_data=read.delim(expression_data_file,sep="\t",stringsAsFactors = F)
#     expression_data=t(expression_data) # row = position; col = gene
#     # Only genes with more than overall 10 raeds
#     gene_sum=data.frame(count=colSums(expression_data))
#     gene_sum$percent_of_data=gene_sum$count/(sum(gene_sum$count))
#     # sum(gene_sum$percent[gene_sum$count>=10])
#     # row.names(gene_sum)
#     gene_sum_ordered=gene_sum[order(gene_sum$count,decreasing = T),]
#     # create the spatial object
#     xy_list=strsplit(row.names(expression_data), "x")
#     xy=as.data.frame(t(as.data.frame(xy_list)))
#     xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
#     names(xy)=c("x","y")
#     expression_data_with_xy=cbind(xy,expression_data)
#     expression_data_with_xy$x=as.numeric(expression_data_with_xy$x)
#     expression_data_with_xy$y=as.numeric(expression_data_with_xy$y)
#     
#     rownames(expression_data_with_xy)=rownames(expression_data)
#     # subset genes
#     RNA_cutoff=100 # Only genes with more than 100 reads are accounted
#     expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
#     
#     # Join the expression and all_Fungi all_Bacteria
#     joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
#     
#     # create the spatial object
#     spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
#     
#     #store it in the a list
#     all_samples_spatial_data[[sample_name]]=spatial_data
#     
#     # clean up
#     rm (sample_name,hotspost_df,data_file,data,xy_list,xy,data_with_xy,
#         expression_data_file,expression_data,gene_sum,gene_sum_ordered,expression_data_with_xy,RNA_cutoff,expression_data_with_xy_subset,
#         joint_Expression_Bacteria_Fungi_data,spatial_data)
#   }
# }
# rm (Exp,smpl) 
# END OLD SECTION

### RNA Cutoff=10 --> For Moran stats
# all_samples_spatial_data_RNACutoff10=list()
# message ("## Read all data into spatial objects")
# for (Exp in c("OMNI12","OMNI13"))
# {
#   samples=c("A1","A2","B1","B2","C1","C2")
#   if (Exp=="OMNI13") {samples=c(samples,"D1")}
#   for (smpl in samples)
#   {
#     sample_name=paste(Exp,smpl,sep="_")
#     message(paste("--",sample_name))
#     hotspost_df=data.frame()
#     data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
#     data=read.delim(data_file,sep=";",stringsAsFactors = F)
#     data$all_Bacteria=rowSums(data[,1:50])
#     data$all_Fungi=rowSums(data[,51:100])
#     
#     xy_list=strsplit(row.names(data), "x")
#     xy=as.data.frame(t(as.data.frame(xy_list)))
#     xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
#     names(xy)=c("x","y")
#     data_with_xy=cbind(xy,data)
#     rownames(data_with_xy)=rownames(data)
#     data_with_xy$x=as.numeric(data_with_xy$x)
#     data_with_xy$y=as.numeric(data_with_xy$y)
#     
#     # Now with the expression
#     expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",Exp,".",smpl,".tsv",sep="")
#     expression_data=read.delim(expression_data_file,sep="\t",stringsAsFactors = F)
#     expression_data=t(expression_data) # row = position; col = gene
#     # Only genes with more than overall 10 raeds
#     gene_sum=data.frame(count=colSums(expression_data))
#     gene_sum$percent_of_data=gene_sum$count/(sum(gene_sum$count))
#     # sum(gene_sum$percent[gene_sum$count>=10])
#     # row.names(gene_sum)
#     gene_sum_ordered=gene_sum[order(gene_sum$count,decreasing = T),]
#     # create the spatial object
#     xy_list=strsplit(row.names(expression_data), "x")
#     xy=as.data.frame(t(as.data.frame(xy_list)))
#     xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
#     names(xy)=c("x","y")
#     expression_data_with_xy=cbind(xy,expression_data)
#     expression_data_with_xy$x=as.numeric(expression_data_with_xy$x)
#     expression_data_with_xy$y=as.numeric(expression_data_with_xy$y)
#     
#     rownames(expression_data_with_xy)=rownames(expression_data)
#     # subset genes
#     RNA_cutoff=10 # Only genes with more than 100 reads are accounted
#     expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
#     
#     # Join the expression and all_Fungi all_Bacteria
#     joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
#     
#     # create the spatial object
#     spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
#     
#     #store it in the a list
#     all_samples_spatial_data_RNACutoff10[[sample_name]]=spatial_data
#     
#     # clean up
#     rm (sample_name,hotspost_df,data_file,data,xy_list,xy,data_with_xy,
#         expression_data_file,expression_data,gene_sum,gene_sum_ordered,expression_data_with_xy,RNA_cutoff,expression_data_with_xy_subset,
#         joint_Expression_Bacteria_Fungi_data,spatial_data)
#   }
# }
# rm (Exp,smpl) 


# calculate the localG per layer
library(raster)
library(spdep)  # poly2nb
library(classInt) # needed?
Getis_Ord_GI_per_grid = function (grid_size,spatial_data,layer_name)
{
  pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
  box = round(extent(spatial_data) / pixelsize) * pixelsize
  template = raster(box, crs = spatial_data,
                    nrows = (box@ymax - box@ymin) / pixelsize, 
                    ncols = (box@xmax - box@xmin) / pixelsize)
  getisraster = rasterize(spatial_data, template, field = layer_name, fun = sum)
  getisgrid = rasterToPolygons(getisraster)
  # plot(getisgrid)
  # Create the list of neighbors
  neighbors = poly2nb(getisgrid)
  weighted_neighbors = nb2listw(neighbors, zero.policy=T)
  
  # plot(getisgrid, border = 'lightgrey')
  # plot(neighbors, coordinates(getisgrid), add=TRUE, col = 'red')
  
  # Perform the local G analysis (Getis-Ord GI*)
  # local_g=localG(getisgrid$layer, weighted_neighbors)
  # local_g1=cbind(getisgrid, as.matrix(local_g))
  # names(local_g1)[2]="gstat"
  
  getisgrid$HOTSPOT = as.vector(localG(getisgrid$layer, weighted_neighbors))
  
  # tm_shape(local_g1) + tm_fill("gstat", palette = "RdBu", style = "pretty") +
  #  tm_borders(alpha=.4)
  
  # calculate the multiple testing adjusted p-value based on the number of nighbours+1
  getisgrid$HOTSPOT.p=2*pnorm(-abs((getisgrid$HOTSPOT)))
  getisgrid$HOTSPOT.p.SP_FDR=p.adjustSP(getisgrid$HOTSPOT.p, neighbors, "BH")
  
  # globalMoran <- moran.test(getisgrid$layer, weighted_neighbors)
  return (getisgrid)
}

Getis_Ord_GI_per_grid_and_Moran = function (grid_size,spatial_data,layer_name)
{
  pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
  box = round(extent(spatial_data) / pixelsize) * pixelsize
  template = raster(box, crs = spatial_data,
                    nrows = (box@ymax - box@ymin) / pixelsize, 
                    ncols = (box@xmax - box@xmin) / pixelsize)
  getisraster = rasterize(spatial_data, template, field = layer_name, fun = sum)
  getisgrid = rasterToPolygons(getisraster)
  # plot(getisgrid)
  # Create the list of neighbors
  neighbors = poly2nb(getisgrid)
  weighted_neighbors = nb2listw(neighbors, zero.policy=T)
  
  # plot(getisgrid, border = 'lightgrey')
  # plot(neighbors, coordinates(getisgrid), add=TRUE, col = 'red')
  
  # Perform the local G analysis (Getis-Ord GI*)
  # local_g=localG(getisgrid$layer, weighted_neighbors)
  # local_g1=cbind(getisgrid, as.matrix(local_g))
  # names(local_g1)[2]="gstat"
  
  getisgrid$HOTSPOT = as.vector(localG(getisgrid$layer, weighted_neighbors))
  
  # tm_shape(local_g1) + tm_fill("gstat", palette = "RdBu", style = "pretty") +
  #  tm_borders(alpha=.4)
  
  # calculate the multiple testing adjusted p-value based on the number of nighbours+1
  getisgrid$HOTSPOT.p=2*pnorm(-abs((getisgrid$HOTSPOT)))
  getisgrid$HOTSPOT.p.SP_FDR=p.adjustSP(getisgrid$HOTSPOT.p, neighbors, "BH")
  
  globalMoran <- moran.test(getisgrid$layer, weighted_neighbors)
  results=list()
  results[[1]]=getisgrid
  results[[2]]=globalMoran
  return (results)
}

# ploting functions
library(RColorBrewer)
library(tmap)
get_hotspots_maps=function (spatial_data,layer,grid_size=2,title=NA)
{
  if (is.na(title)) {title=layer}
  # pdf
  # pdf (file = out_pdf,height = 7,width = 7)
  
  # calculate the localG stas
  if (length(spatial_data[[layer]])>0)
  {
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = layer)
    pmap=tm_shape(getisgrid) + 
      tm_fill("HOTSPOT", 
              palette = "-RdBu",
              style = "pretty", title="G stat",n=5,legend.reverse=T, midpoint=0) +
      tm_borders(alpha=.4) + tm_legend(legend.text.size=0.5,legend.title.size=0.6) +
      tm_layout(title = title,title.size=0.8)
  } else {
    pmap=NULL
  }
  
  # plot
  #  brewer.pal(n = 3, name = "RdBu")
  #  breaks = c(-20, -1.96, -1, 1, 1.96, 20)
  #  palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
  #  col = palette[cut(getisgrid$HOTSPOT, breaks)]
  #  plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
  #  legend("bottom", inset=.02, title="Z score (local spatial G[i] statistic)",
  #         legend =c("<-1.96","-1","1","1.96",">1.96"), 
  #         fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
  #         xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')

 
  
  return (pmap)
}

get_hotspots_maps_fixed=function (spatial_data,layer,grid_size=2,title=NA,only_significant=F)
{
  if (is.na(title)) {title=layer}
  # pdf
  # pdf (file = out_pdf,height = 7,width = 7)
  
  # calculate the localG stas
  if (length(spatial_data[[layer]])>0)
  {
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = layer)
    if (only_significant) { # make the non significant hotspots value as NA
      getisgrid$HOTSPOT[getisgrid$HOTSPOT.p.SP_FDR>0.05]=NA
    }
    # https://geocompr.github.io/post/2019/tmap-color-scales/
    if (min(getisgrid$HOTSPOT,na.rm = T) < -12 | max(getisgrid$HOTSPOT,na.rm=T)>12) {
      message(paste("== [WARNING]",title,"min:",min(getisgrid$HOTSPOT,na.rm = T),"max:",max(getisgrid$HOTSPOT,na.rm=T)))
    }
    pmap=tm_shape(getisgrid) + 
      tm_fill("HOTSPOT", 
              palette = "-RdBu",
              style = "fixed",title="G stat",legend.reverse=T, midpoint=0,breaks = c(-12,-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10,12),textNA = "NS") +
      # style = "cont", breaks = c(-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10) ,title="G stat",n=5,legend.reverse=T, midpoint=0) +
      tm_borders(alpha=.4) + tm_legend(legend.text.size=0.5,legend.title.size=0.6) +
      tm_layout(title = title,title.size=0.8)
  } else {
    pmap=NULL
  }
  
  # plot
  #  brewer.pal(n = 3, name = "RdBu")
  #  breaks = c(-20, -1.96, -1, 1, 1.96, 20)
  #  palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
  #  col = palette[cut(getisgrid$HOTSPOT, breaks)]
  #  plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
  #  legend("bottom", inset=.02, title="Z score (local spatial G[i] statistic)",
  #         legend =c("<-1.96","-1","1","1.96",">1.96"), 
  #         fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
  #         xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
  return (pmap)
}

library(sf)
get_dot_map=function(spatial_data,layer,title=NA,palette="OrRd") # take spatialDataFrame and return tm_map for the specified layer
                                                                 # "GnBu","PuBuGn"
{
  if (is.na(title)) {title=layer}
  p=NULL
  if (length(spatial_data[[layer]])>0) # validate layer exists
  {
    tmp.sp=spatial_data[layer]
    # tmp.sp[[layer]][tmp.sp[[layer]]<10]=NA
    tmp.sf=st_as_sf(tmp.sp)
    
    p <- tm_shape(tmp.sf) +
      tm_dots(layer, shape = 19, alpha = 0.5, size = 0.1, 
              palette = palette) + # for legend title title=title) +
      tm_layout(title = title,title.size=0.8) 
  }
   return (p)
}

# dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated
# dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Bacteria+dataset_to_plot$sum_sections_Fungi,
#                                       pmax(dataset_to_plot$r_LocalG.Bacteria.max_value,dataset_to_plot$r_LocalG.Fungi.max_value),
#                                       decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_gene.pdf"))
# old
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_gene_BORUTA_AND_significant_0.01_spearman.pdf"))

# Apr2022
# OnlySignificant=TRUE
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_selected_gene_Apr2022_220405_Expression_BORUTA_AND_significant_0.01_spearman.OnlySignificant_HS.pdf"))
# pdf(pdf_file,width = 21,height = 14)

# defense related genes
dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[grepl(x=final_selected_genes_from_all_sections_with_all_details_annotated$Gene.ontology..biological.process.,pattern = "defense"),]
dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Bacteria+dataset_to_plot$sum_sections_Fungi,
                                      pmax(dataset_to_plot$r_LocalG.Bacteria.max_value,dataset_to_plot$r_LocalG.Fungi.max_value),
                                      decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/defense_response_selected_gene.pdf"))
# pdf(pdf_file,width = 21,height = 14) 
# Apr 2022
 OnlySignificant=TRUE
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/defense_response_selected_gene_Apr2022_220405_Expression_BORUTA_AND_significant_0.01_spearman.OnlySignificant_HS.pdf.pdf"))
 pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/defense_response_selected_gene_Apr2022_220405_Expression_BORUTA_AND_significant_0.01_spearman.OMNI13_B1.OnlySignificant_HS.pdf.pdf"))
 # OnlySignificant=FALSE
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/defense_response_selected_gene_Apr2022_220405_Expression_BORUTA_AND_significant_0.01_spearman.All_HS.pdf.pdf"))

pdf(pdf_file,width = 21,height = 14) 

# Fungi response
# dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[grepl(x=final_selected_genes_from_all_sections_with_all_details_annotated$response_type,pattern = "Fungi"),]
# dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Fungi,
#                                       dataset_to_plot$r_LocalG.Fungi.max_value,
#                                       decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Fungi_response_selected_gene.pdf"))
# pdf(pdf_file,width = 21,height = 14) 

# Bacteria response
# dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[grepl(x=final_selected_genes_from_all_sections_with_all_details_annotated$response_type,pattern = "Bacteria"),]
# dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Bacteria,
#                                       dataset_to_plot$r_LocalG.Bacteria.max_value,
#                                       decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Bacteria_response_selected_gene.pdf"))
# pdf(pdf_file,width = 21,height = 14) 

# do the ploting
for (i in 1:NROW(dataset_to_plot))
{
  # final_selected_genes_from_all_sections_with_all_details_annotated[1,]
  # nPosFungi=NA
  # nPosBac=NA
  nFungi=dataset_to_plot$sum_sections_Fungi[i]
  nBac=dataset_to_plot$sum_sections_Bacteria[i]
  gene_name=dataset_to_plot$Entry.name[i]
  response_type=dataset_to_plot$response_type[i]
  
  gene_id=dataset_to_plot[i,"gene"]
  plots_list=list()
  for (cor_type in c("r_LocalG","r_Reads"))
  {
    for (tax in c("Bacteria","Fungi")) {
      plot_name=paste(tax,cor_type,sep=".")
      
      cor_sign="max" # for positive correlation
      
      if (dataset_to_plot[[paste(cor_type,tax,"total_pos_r",sep=".")]][i]<
          dataset_to_plot[[paste(cor_type,tax,"total_neg_r",sep=".")]][i]) {
        cor_sign="min" # mostly negative correlations
      }
      # take the sample name to plot
      # the best
      # sample_to_plot=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
      #                                        pattern = paste(cor_type,".",tax,".",cor_sign,"_sample",sep=""))]
      # fixed
      sample_to_plot="OMNI13_B1"
      cor_value=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                        pattern = paste(cor_type,".",tax,".",cor_sign,"_value",sep=""))]
      cor_value=signif(cor_value,3)
      
      
      layer_to_plot=NA
      if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
      if (tax=="Fungi"){layer_to_plot="all_Fungi"}
      
      if (cor_type == "r_LocalG")
      {
        # plots_list[[plot_name]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        # plots_list[[paste(plot_name,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,sep="")) # the expression
      
        # with fixed scale and significant
        plots_list[[plot_name]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""),only_significant = OnlySignificant)
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,sep=""),only_significant = OnlySignificant) # the expression
      }
      if (cor_type== "r_Reads")
      {
        plots_list[[plot_name]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,sep="")) # the expression
      }
    } # end {Baccteria|Fungi}
  } # end {LocalG | Reads}
  # plots_list
  
  # https://github.com/r-tmap/tmap/issues/511
  
  # current.mode <- tmap_mode("plot")
  ## tmap_arrange(plots_list, widths = c(.75, .75))
  
  plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
  plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
  plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
  plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
  
  # for (i in 3:4) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
  print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 4, nrow = 2))
  if (i%%10==0) {message(paste("Finished ploting",i,"out of",NROW(dataset_to_plot),sep=" "))}
  flush.console()
  # tmap_mode(current.mode)
}

dev.off()

# for each section plot all the expression, Fungi and Bacteria

##### Final set of figures
gene_to_plot=c("AT3G01500","AT4G14400","AT2G14560") # BCA1_ARATH, ACD6, LURP1, 
sample_to_plot="OMNI13_B1"
# sample_to_plot="OMNI12_A1"

dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[final_selected_genes_from_all_sections_with_all_details_annotated$gene %in% gene_to_plot,]
dataset_to_plot=dataset_to_plot[order(dataset_to_plot[[paste(sample_to_plot,"r_LocalG","Bacteria",sep=".")]],decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_gene_for_paper.OMNI13_B1.pdf"))
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_gene_for_paper.",sample_to_plot,".BORUTA_AND_SPEARMAN.pdf",sep=""))
# Apr2022
# OnlySignificant=TRUE
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/selected_gene_for_paper.",sample_to_plot,".BORUTA_AND_SPEARMAN_OnlySignificant_HS.pdf",sep=""))
# pdf(pdf_file,width = 21,height = 16) 

# Apr2022
OnlySignificant=FALSE
pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/selected_gene_for_paper.",sample_to_plot,".BORUTA_AND_SPEARMAN_ALL_HS.pdf",sep=""))
pdf(pdf_file,width = 21,height = 16) 

plots_list=list()
for (cor_type in c("r_LocalG","r_Reads")) 
{
  for (tax in c("Bacteria","Fungi")) 
  {
    layer_to_plot=NA
    # the marker (16S/ITS)
    if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
    if (tax=="Fungi"){layer_to_plot="all_Fungi"}
    if (cor_type == "r_LocalG")
    {
      # plots_list[[paste(sample_to_plot,tax,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(sample_to_plot," - ",tax,sep=""))
      plots_list[[paste(sample_to_plot,tax,sep=".")]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(sample_to_plot," - ",tax,sep=""),only_significant = OnlySignificant)
    }
    if (cor_type== "r_Reads")
    {
      plots_list[[paste(sample_to_plot,tax,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(sample_to_plot," - ",tax,sep=""))
    }
  }
  # the expression
  for (i in 1:NROW(dataset_to_plot))
  {
    nFungi=dataset_to_plot$sum_sections_Fungi[i]
    nBac=dataset_to_plot$sum_sections_Bacteria[i]
    gene_name=dataset_to_plot$Entry.name[i]
    response_type=dataset_to_plot$response_type[i]
    
    gene_id=dataset_to_plot[i,"gene"]
    gene_sum_reads=sum(all_samples_spatial_data[[sample_to_plot]][[gene_id]])
    Bcor_value=dataset_to_plot[i,paste(sample_to_plot,cor_type,"Bacteria",sep=".")]
    Bcor_value=signif(Bcor_value,3)
    Fcor_value=dataset_to_plot[i,paste(sample_to_plot,cor_type,"Fungi",sep=".")]
    Fcor_value=signif(Fcor_value,3)
    if (cor_type == "r_LocalG")
    {
          # plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ntotal reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,sep="")) # the exression
          plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ntotal reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,sep=""),only_significant = OnlySignificant) # the exression
    }
    if (cor_type== "r_Reads")
    {
        plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,"\ntotal reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,sep="")) # the expression
    }
  }
  # plots_list
  # https://github.com/r-tmap/tmap/issues/511
  
  # current.mode <- tmap_mode("plot")
  ## tmap_arrange(plots_list, widths = c(.75, .75))
#  if (length(samples_with_expression)==length(all_samples))
#  {
#    plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
#    plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
#    plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
#    plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
#    for (i in 5:length(all_samples)) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
#  }
  print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 3, nrow = 2))
  # tmap_mode(current.mode)
#  if (i%%10==0) {message(paste("Finished ploting",i,"out of",NROW(dataset_to_plot),sep=" "))}
  flush.console()
}

dev.off()

# save the plotted df
# write.table(dataset_to_plot,file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_gene_for_paper.OMNI13_B1.txt",sep="\t",row.names = F,quote = F)
# write.table(dataset_to_plot,file=paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_gene_for_paper.",sample_to_plot,".BORUTA_AND_SPEARMAN.txt",sep=""),sep="\t",row.names = F,quote = F)
write.table(dataset_to_plot,file=paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/selected_gene_for_paper.Apr2022.",sample_to_plot,".BORUTA_AND_SPEARMAN.txt",sep=""),sep="\t",row.names = F,quote = F)

### plot the selected genes in all sections
gene_to_plot=c("AT3G01500","AT4G14400","AT2G14560") # BCA1_ARATH, ACD6, LURP1, 
# sample_to_plot="OMNI13_B1"
# sample_to_plot="OMNI12_A1"

dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[final_selected_genes_from_all_sections_with_all_details_annotated$gene %in% gene_to_plot,]


cor_type="r_LocalG"
for (OnlySignificant in c(TRUE,FALSE))
{
  plots_list=list()
  for (sample_to_plot in names(all_samples_spatial_data))
  {
    # Bacteria, Fungi
    for (tax in c("Bacteria","Fungi")) 
    {
      layer_to_plot=NA
      if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
      if (tax=="Fungi"){layer_to_plot="all_Fungi"}
      plots_list[[paste(sample_to_plot,tax,sep=".")]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(sample_to_plot," - ",tax,sep=""),only_significant = OnlySignificant)
    } 
    
    for (i in 1:NROW(dataset_to_plot))
    {
      nFungi=dataset_to_plot$sum_sections_Fungi[i]
      nBac=dataset_to_plot$sum_sections_Bacteria[i]
      gene_name=dataset_to_plot$Entry.name[i]
      response_type=dataset_to_plot$response_type[i]
      
      gene_id=dataset_to_plot[i,"gene"]
      gene_sum_reads=sum(all_samples_spatial_data[[sample_to_plot]][[gene_id]])
      Bcor_value=dataset_to_plot[i,paste(sample_to_plot,cor_type,"Bacteria",sep=".")]
      Bcor_value=signif(Bcor_value,3)
      Fcor_value=dataset_to_plot[i,paste(sample_to_plot,cor_type,"Fungi",sep=".")]
      Fcor_value=signif(Fcor_value,3)
      if (cor_type == "r_LocalG")
      {
        # plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ntotal reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,sep="")) # the exression
        plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ntotal reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,sep=""),only_significant = OnlySignificant) # the exression
      }
    }
  }
  
  if (OnlySignificant)
  {
    pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/selected_gene_for_paper.all_samples.OnlySignificant_HS.pdf",sep=""))
    pdf(pdf_file,width = 35,height = 30) 
    print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 5, nrow = 13))
    dev.off()
  }
  if (!OnlySignificant)
  {
    pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/selected_gene_for_paper.all_samples.All_HS.pdf",sep=""))
    pdf(pdf_file,width = 35,height = 16) 
    print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 5, nrow = 13))
    dev.off()
  }
}

# plot reads
cor_type="r_Reads"
plots_list=list()
for (sample_to_plot in names(all_samples_spatial_data))
{
  # Bacteria, Fungi
  for (tax in c("Bacteria","Fungi")) 
  {
    layer_to_plot=NA
    if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
    if (tax=="Fungi"){layer_to_plot="all_Fungi"}
    plots_list[[paste(sample_to_plot,tax,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(sample_to_plot," - ",tax,sep=""))
  } 
  
  for (i in 1:NROW(dataset_to_plot))
  {
    nFungi=dataset_to_plot$sum_sections_Fungi[i]
    nBac=dataset_to_plot$sum_sections_Bacteria[i]
    gene_name=dataset_to_plot$Entry.name[i]
    response_type=dataset_to_plot$response_type[i]
    
    gene_id=dataset_to_plot[i,"gene"]
    gene_sum_reads=sum(all_samples_spatial_data[[sample_to_plot]][[gene_id]])
    Bcor_value=dataset_to_plot[i,paste(sample_to_plot,cor_type,"Bacteria",sep=".")]
    Bcor_value=signif(Bcor_value,3)
    Fcor_value=dataset_to_plot[i,paste(sample_to_plot,cor_type,"Fungi",sep=".")]
    Fcor_value=signif(Fcor_value,3)
    if (cor_type == "r_Reads")
    {
      # plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ntotal reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,sep="")) # the exression
      plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,"\ntotal reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,sep="")) # the expression
    }
  }
}
pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/selected_gene_for_paper.all_samples.reads.pdf",sep=""))
pdf(pdf_file,width = 35,height = 30) 
print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 5, nrow = 13))
dev.off()
  
  

# -- CONTINUE UPDATE HERE!
#### As a control, plot the localG of all the genes on that sections
OnlySignificant=TRUE
RNASum_Cutoff=100
for (Exp_to_plot in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp_to_plot=="OMNI13") {samples=c(samples,"D1")}
  for (S_to_plot in samples)
  {
    # Exp_to_plot="OMNI13"
    # S_to_plot="B1"
    sample_to_plot=paste(Exp_to_plot,S_to_plot,sep="_")
    
    sample_data=all_samples_spatial_data[[sample_to_plot]]@data
    sample_sums=(colSums(sample_data))
    gene_ids=names(sample_sums)[order(sample_sums,decreasing = T)]
    
    # load the info about the genes of these section
    # Boruta_based_selected_Bacteria_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp_to_plot,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Bacteria_sum.",S_to_plot,".tdl",sep="");
    # Apr2022
    Boruta_based_selected_Bacteria_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp_to_plot,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/","Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.Apr2022_220405_semiwild_dataset_rawcounts_filtered.RNA.with_Top50genus_Bacteria_sum.",S_to_plot,".tdl",sep="")
    Boruta_based_selected_Bacteria=read.delim(file = Boruta_based_selected_Bacteria_file,stringsAsFactors = F)
    # Boruta_based_selected_Fungi_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Fungi_sum.",S,".tdl",sep="");
    # Boruta_based_selected_Fungi_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp_to_plot,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Fungi_sum.",S_to_plot,".tdl",sep="");
    # Apr2022
    Boruta_based_selected_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp_to_plot,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/","Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.Apr2022_220405_semiwild_dataset_rawcounts_filtered.RNA.with_Top50genus_Fungi_sum.",S_to_plot,".tdl",sep="")
    Boruta_based_selected_Fungi=read.delim(file = Boruta_based_selected_Fungi_file,stringsAsFactors = F)
    
    ## TO DO: BoxPlot of the expression sum
    
    # old
    # pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_genes",sample_to_plot,"pdf",sep="."))
    # Apr2022
    pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_genes",sample_to_plot,"OnlySignificant_HS","pdf",sep="."))
    
    pdf(pdf_file,width = 35,height = 16) 
    
    ncol=5
    plots_list=rep(NULL, ncol*2)
    for (i in 1:min(length(gene_ids),1300))
    {
      if (i%%ncol==1&i>1) # flush the aggregated plots
      {
        print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = ncol, nrow = 2))
        plots_list=rep(NULL, ncol*2)
      }
      for (cor_type in c("r_LocalG","r_Reads")) 
      {
        gene_id=gene_ids[i]
        gene_sum_reads=sum(all_samples_spatial_data[[sample_to_plot]][[gene_id]])
        
        # TO DO add if in final set of selected, pvalue for r
        Bcor_value=Boruta_based_selected_Bacteria[Boruta_based_selected_Bacteria$gene==gene_id,cor_type]
        Bcor_value=signif(Bcor_value,3)
        Fcor_value=Boruta_based_selected_Fungi[Boruta_based_selected_Fungi$gene==gene_id,cor_type]
        Fcor_value=signif(Fcor_value,3)
        
        # FDR.BH_LocalG FDR.BH_Reads
        Final_selected=gene_id%in%final_selected_genes_from_all_sections_with_all_details_annotated$gene
        plot_index=(i-1)%%ncol+1
        if (cor_type == "r_LocalG")
        {
          if (gene_sum_reads>=RNASum_Cutoff)
          {
            # plots_list[[plot_index]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ni=",i," total reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,"\nFinal_selection:",Final_selected,sep="")) # the exression
            plots_list[[plot_index]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ni=",i," total reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,"\nFinal_selection:",Final_selected,sep=""),only_significant = OnlySignificant) # the exression
          } else {
            plots_list[[plot_index]]=NULL
          }
        }
        if (cor_type== "r_Reads")
        {
          plots_list[[plot_index+ncol]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,"\ni=",i," total reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,"\nFinal_selection:",Final_selected,sep="")) # the expression
        }
      }
      if (i%%(ncol*2)==0) {message(paste("Finished ploting",i,"out of",length(gene_ids),sep=" "))}
      flush.console()
    }
    # plot the reminder
    print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = ncol, nrow = 2))
    dev.off()
    sum(gene_ids%in%final_selected_genes_from_all_sections_with_all_details_annotated$gene)
    
    rm (sample_data,gene_ids,sample_to_plot,
        Boruta_based_selected_Bacteria_file,Boruta_based_selected_Bacteria,
        Boruta_based_selected_Fungi_file,Boruta_based_selected_Fungi,
        plot_index,col)
  }
}
rm(Exp_to_plot,S_to_plot)

### For single
# Exp_to_plot="OMNI13"
# S_to_plot="B1"
sample_to_plot=paste(Exp_to_plot,S_to_plot,sep="_")

sample_data=all_samples_spatial_data[[sample_to_plot]]@data
sample_sums=(colSums(sample_data))
gene_ids=names(sample_sums)[order(sample_sums,decreasing = T)]

# load the info about the genes of these section
# Old
# Boruta_based_selected_Bacteria_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp_to_plot,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Bacteria_sum.",S_to_plot,".tdl",sep="");
# Apr2022
Boruta_based_selected_Bacteria_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp_to_plot,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Expression_220405/","Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.Apr2022_220405_semiwild_dataset_rawcounts_filtered.RNA.with_Top50genus_Bacteria_sum.",S_to_plot,".tdl",sep="")

Boruta_based_selected_Bacteria=read.delim(file = Boruta_based_selected_Bacteria_file,stringsAsFactors = F)
# Boruta_based_selected_Fungi_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Fungi_sum.",S,".tdl",sep="");
# Apr20222
Boruta_based_selected_Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp_to_plot,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Expression_220405/","Boruta_selected_reads_or_2x2_localG_AND_significant_0.01_spearman.Apr2022_220405_semiwild_dataset_rawcounts_filtered.RNA.with_Top50genus_Fungi_sum.",S_to_plot,".tdl",sep="")
Boruta_based_selected_Fungi=read.delim(file = Boruta_based_selected_Fungi_file,stringsAsFactors = F)

## TO DO: BoxPlot of the expression sum

# old 
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_genes",sample_to_plot,"pdf",sep="."))
# Ap2022
OnlySignificant=TRUE
pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_genes",sample_to_plot,"OnlySignificant_HS","pdf",sep="."))
pdf(pdf_file,width = 35,height = 16) 

ncol=5
plots_list=rep(NULL, ncol*2)
for (i in 1:length(gene_ids))
{
  if (i%%ncol==1&i>1) # flush the aggregated plots
  {
    print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = ncol, nrow = 2))
    plots_list=rep(NULL, ncol*2)
  }
  for (cor_type in c("r_LocalG","r_Reads")) 
  {
    gene_id=gene_ids[i]
    gene_sum_reads=sum(all_samples_spatial_data[[sample_to_plot]][[gene_id]])
    
    # TO DO add if in final set of selected, pvalue for r
    Bcor_value=Boruta_based_selected_Bacteria[Boruta_based_selected_Bacteria$gene==gene_id,cor_type]
    Bcor_value=signif(Bcor_value,3)
    Fcor_value=Boruta_based_selected_Fungi[Boruta_based_selected_Fungi$gene==gene_id,cor_type]
    Fcor_value=signif(Fcor_value,3)
    
    # FDR.BH_LocalG FDR.BH_Reads
    Final_selected=gene_id%in%final_selected_genes_from_all_sections_with_all_details_annotated$gene
    
    plot_index=(i-1)%%ncol+1
    if (cor_type == "r_LocalG")
    {
      # plots_list[[plot_index]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ni=",i," total reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,"\nFinal_selection:",Final_selected,sep="")) # the exression
      plots_list[[plot_index]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,"\ni=",i," total reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,"\nFinal_selection:",Final_selected,sep=""),only_significant = OnlySignificant) # the exression
      
    }
    if (cor_type== "r_Reads")
    {
      plots_list[[plot_index+ncol]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,"\ni=",i," total reads:",gene_sum_reads,"\nBac_r=",Bcor_value," Fun_r=",Fcor_value,"\nFinal_selection:",Final_selected,sep="")) # the expression
    }
  }
  if (i%%(ncol*2)==0) {message(paste("Finished ploting",i,"out of",length(gene_ids),sep=" "))}
  flush.console()
}
# plot the reminder
print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = ncol, nrow = 2))
dev.off()
sum(gene_ids%in%final_selected_genes_from_all_sections_with_all_details_annotated$gene)

rm (sample_data,gene_ids,Exp_to_plot,S_to_plot,sample_to_plot,
    Boruta_based_selected_Bacteria_file,Boruta_based_selected_Bacteria,
    Boruta_based_selected_Fungi_file,Boruta_based_selected_Fungi,
    plot_index,col)

#### collect the Moran stats for all smaples
rm (tmp,tax,type,template,spatial_data,S_to_plot,r_type,results,S,box,sample_sums,sample_to_plot,localG_Obj,Moran_test_Obj)
all_genes_all_samples_expression_and_Moran=data.frame()

for (Exp_to_plot in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp_to_plot=="OMNI13") {samples=c(samples,"D1")}
  for (S_to_plot in samples)
  {
    # Exp_to_plot="OMNI13"
    # S_to_plot="B1"
    sample_name=paste(Exp_to_plot,S_to_plot,sep="_")
    message (paste ("==",sample_name))
    sample_data=all_samples_spatial_data[[sample_name]]@data
    sample_sums=(colSums(sample_data))
    sample_stats=data.frame(gene_id=names(sample_sums))
    sample_stats[[paste(sample_name,"readsCount",sep=".")]]=sample_sums
    sample_stats[[paste(sample_name,"MoranStat",sep=".")]]=NA
    sample_stats[[paste(sample_name,"MoranPval",sep=".")]]=NA
    sample_stats=sample_stats[order(sample_stats[[paste(sample_name,"readsCount",sep=".")]],decreasing = T),]
    for (i in 1:NROW(sample_stats))
    {
        gene_id=sample_stats$gene_id[i]
        if (sum(gene_id %in% names(all_samples_spatial_data[[sample_name]]))>0)
        {
          results=Getis_Ord_GI_per_grid_and_Moran(grid_size = 2,spatial_data = all_samples_spatial_data[[sample_name]],layer_name = gene_id)
          localG_Obj=results[[1]]
          Moran_test_Obj=results[[2]]
          sample_stats[[paste(sample_name,"MoranStat",sep=".")]][i]=Moran_test_Obj$statistic
          sample_stats[[paste(sample_name,"MoranPval",sep=".")]][i]=Moran_test_Obj$p.value
        } else {
          message (paste("[WARNNING] Could not find spatial data for ",gene_id," in ",sample_name," skip...",sep=""))
        }
          
        if (i%%100==0) {message(paste("Finished ",i,"out of",NROW(sample_stats),sep=" "))}
        flush.console()
    }
    if (NROW(all_genes_all_samples_expression_and_Moran)>0)
    {
      all_genes_all_samples_expression_and_Moran=merge(all_genes_all_samples_expression_and_Moran,sample_stats,by="gene_id",all=T)
    } else {
      all_genes_all_samples_expression_and_Moran=sample_stats
    }
    rm (sample_name,sample_data,sample_stats)
  } # finish sample
}
all_genes_all_samples_expression_and_Moran_annotated=merge(all_genes_all_samples_expression_and_Moran,Sprot_annotation,by.x="gene_id",by.y="Cross.reference..Araport.",all.x=T)
# save!    
#write.table(file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_samples_all_genes_total_reads_and_Moran.txt",
#            x = all_genes_all_samples_expression_and_Moran_annotated,quote = F,sep = "\t",row.names = F)
write.table(file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/Apr2022/all_samples_all_genes_Apr2022_220405_total_reads_and_Moran.txt",
               x = all_genes_all_samples_expression_and_Moran_annotated,quote = F,sep = "\t",row.names = F)
            
# -- CONTINUE UPDATE HERE!


# find examples of genes with high expression but low Moran
all_genes_all_samples_expression_and_Moran_annotated=read.delim("~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_samples_all_genes_total_reads_and_Moran.txt",sep = "\t",stringsAsFactors = F)
View(all_genes_all_samples_expression_and_Moran_annotated)

house_keeping_gene_to_plot=c("AT3G18780",
                             "AT1G75780",
                             "AT1G49240",
                             "AT5G62690",
                             "AT2G28390")


# plot the gene expression dist
# library(tidyverse)

 library(tidyr)
# library(ggplot2)
# library(dplyr)
all_genes_all_samples_expression_and_Moran_annotated[!grepl(all_genes_all_samples_expression_and_Moran_annotated$gene_id,pattern = "all",fixed = T),
                                                     grepl(x=names(all_genes_all_samples_expression_and_Moran_annotated),pattern = "readsCount")] %>% 
  gather(key="MesureType", value="Val") %>%
  ggplot( aes(x=MesureType, y=log2(Val), fill=MesureType)) +
  geom_violin()+ theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

all_genes_all_samples_expression_and_Moran_annotated[!grepl(all_genes_all_samples_expression_and_Moran_annotated$gene_id,pattern = "all",fixed = T),
                                                     grepl(x=names(all_genes_all_samples_expression_and_Moran_annotated),pattern = "MoranStat")] %>% 
  gather(key="MesureType", value="Val") %>%
  ggplot( aes(x=MesureType, y=Val, fill=MesureType)) +
  # geom_violin()
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("sample") + ylab("Moran I stat")

# plot the dist of pValues [SHULD IT BE CORRECTED AGAIN FOR Multiple testing?]
significnat_Moran_genes=colSums(all_genes_all_samples_expression_and_Moran_annotated[!grepl(all_genes_all_samples_expression_and_Moran_annotated$gene_id,pattern = "all",fixed = T),
                                                                                     grepl(x=names(all_genes_all_samples_expression_and_Moran_annotated),pattern = "MoranPval")]<0.05,na.rm = T)
total_genes=colSums(!is.na(all_genes_all_samples_expression_and_Moran_annotated[!grepl(all_genes_all_samples_expression_and_Moran_annotated$gene_id,pattern = "all",fixed = T),
                                                                                grepl(x=names(all_genes_all_samples_expression_and_Moran_annotated),pattern = "MoranPval")]))
all_genes_all_samples_expression_and_Moran_annotated[!grepl(all_genes_all_samples_expression_and_Moran_annotated$gene_id,pattern = "all",fixed = T),
                                                     grepl(x=names(all_genes_all_samples_expression_and_Moran_annotated),pattern = "MoranPval")] %>% 
  # gather(key="MesureType", value="Val") %>%
  pivot_longer(cols = ends_with("MoranPval"),names_to = "sample", values_to = "Val") %>%
  ggplot( aes(x=sample, y=-log10(Val), fill=sample)) +
  # geom_violin()
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("sample") + ylab("Moran I -log10(Pval)") +
  geom_hline(yintercept = -log10(0.05),color="red") + 
  
  # https://community.rstudio.com/t/how-to-add-text-annotation-over-each-boxplot-of-grouped-data-in-facet-grids-in-ggplot2/93869/2
  # geom_text(data = pans2, aes(y = yloc, label = label,color=trt), 
  #                          position = position_dodge(width = .75))
signif(((significnat_Moran_genes/total_genes)*100),3)



#### MORAN Stats for RNA Cutoff=10
#### collect the Moran stats for all smaples
rm (tmp,tax,type,template,spatial_data,S_to_plot,r_type,results,S,box,sample_sums,sample_to_plot,localG_Obj,Moran_test_Obj)
all_genes_all_samples_expression_and_Moran_RNA_Cutoff10=data.frame()

for (Exp_to_plot in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp_to_plot=="OMNI13") {samples=c(samples,"D1")}
  for (S_to_plot in samples)
  {
    # Exp_to_plot="OMNI13"
    # S_to_plot="B1"
    sample_name=paste(Exp_to_plot,S_to_plot,sep="_")
    message (paste ("==",sample_name))
    sample_data=all_samples_spatial_data_RNACutoff10[[sample_name]]@data
    sample_sums=(colSums(sample_data))
    sample_stats=data.frame(gene_id=names(sample_sums))
    sample_stats[[paste(sample_name,"readsCount",sep=".")]]=sample_sums
    sample_stats[[paste(sample_name,"MoranStat",sep=".")]]=NA
    sample_stats[[paste(sample_name,"MoranPval",sep=".")]]=NA
    sample_stats=sample_stats[order(sample_stats[[paste(sample_name,"readsCount",sep=".")]],decreasing = T),]
    for (i in 1:NROW(sample_stats))
    {
      gene_id=sample_stats$gene_id[i]
      if (sum(gene_id %in% names(all_samples_spatial_data_RNACutoff10[[sample_name]]))>0)
      {
        results=Getis_Ord_GI_per_grid_and_Moran(grid_size = 2,spatial_data = all_samples_spatial_data_RNACutoff10[[sample_name]],layer_name = gene_id)
        localG_Obj=results[[1]]
        Moran_test_Obj=results[[2]]
        sample_stats[[paste(sample_name,"MoranStat",sep=".")]][i]=Moran_test_Obj$statistic
        sample_stats[[paste(sample_name,"MoranPval",sep=".")]][i]=Moran_test_Obj$p.value
      } else {
        message (paste("[WARNNING] Could not fing spatial data for ",gene_id," in ",sample_name," skip...",sep=""))
      }
      
      if (i%%100==0) {message(paste("Finished ",i,"out of",NROW(sample_stats),sep=" "))}
      flush.console()
    }
    if (NROW(all_genes_all_samples_expression_and_Moran_RNA_Cutoff10)>0)
    {
      all_genes_all_samples_expression_and_Moran_RNA_Cutoff10=merge(all_genes_all_samples_expression_and_Moran_RNA_Cutoff10,sample_stats,by="gene_id",all=T)
    } else {
      all_genes_all_samples_expression_and_Moran_RNA_Cutoff10=sample_stats
    }
    rm (sample_name,sample_data,sample_stats)
  } # finish sample
}
all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated=merge(all_genes_all_samples_expression_and_Moran_RNA_Cutoff10,Sprot_annotation,by.x="gene_id",by.y="Cross.reference..Araport.",all.x=T)
# save!    
write.table(file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_samples_all_genes_total_reads_and_Moran_RNA_Cutoff10.txt",
            x = all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated,quote = F,sep = "\t",row.names = F)

# plot the gene expression dist
# library(tidyverse)

library(tidyr)
# library(ggplot2)
# library(dplyr)
all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated[!grepl(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated$gene_id,pattern = "all",fixed = T),
                                                     grepl(x=names(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated),pattern = "readsCount")] %>% 
  gather(key="MesureType", value="Val") %>%
  ggplot( aes(x=MesureType, y=log2(Val), fill=MesureType)) +
  geom_violin()+ theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1))

all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated[!grepl(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated$gene_id,pattern = "all",fixed = T),
                                                     grepl(x=names(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated),pattern = "MoranStat")] %>% 
  gather(key="MesureType", value="Val") %>%
  ggplot( aes(x=MesureType, y=Val, fill=MesureType)) +
  # geom_violin()
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("sample") + ylab("Moran I stat")

# plot the dist of pValues [SHULD IT BE CORRECTED AGAIN FOR Multiple testing?]
significnat_Moran_genes_RNACutoff10=colSums(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated[!grepl(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated$gene_id,pattern = "all",fixed = T),
                                                                                     grepl(x=names(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated),pattern = "MoranPval")]<0.05,na.rm = T)
total_genes_RNACutoff10=colSums(!is.na(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated[!grepl(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated$gene_id,pattern = "all",fixed = T),
                                                                                grepl(x=names(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated),pattern = "MoranPval")]))
all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated[!grepl(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated$gene_id,pattern = "all",fixed = T),
                                                     grepl(x=names(all_genes_all_samples_expression_and_Moran_RNACutoff10_annotated),pattern = "MoranPval")] %>% 
  # gather(key="MesureType", value="Val") %>%
  pivot_longer(cols = ends_with("MoranPval"),names_to = "sample", values_to = "Val") %>%
  ggplot( aes(x=sample, y=-log10(Val), fill=sample)) +
  # geom_violin()
  geom_boxplot() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  xlab("sample") + ylab("Moran I -log10(Pval)") +
  geom_hline(yintercept = -log10(0.05),color="red") + 
  
  # https://community.rstudio.com/t/how-to-add-text-annotation-over-each-boxplot-of-grouped-data-in-facet-grids-in-ggplot2/93869/2
  # geom_text(data = pans2, aes(y = yloc, label = label,color=trt), 
  #                          position = position_dodge(width = .75))
  signif(((significnat_Moran_genes_RNACutoff10/total_genes_RNACutoff10)*100),3)
###




### HERE!
#### Plot GO enrichment - Shared genes
GO_enrichment_level0_Shared_file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_genes_BORUTA_and_SPEARMAN.shared.PANTHER.GO.level0.txt"
GO_enrichment_level0_Shared=read.delim(file = GO_enrichment_level0_Shared_file,sep = "\t",header = T,stringsAsFactors = F)
GO_enrichment_level0_Shared_sorted=GO_enrichment_level0_Shared[with(GO_enrichment_level0_Shared, order(number_in_list, fold_enrichment,decreasing = T)),]

library(ggplot2)
pdf("~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_genes_BORUTA_and_SPEARMAN.shared.PANTHER.GO.level0.Top25_byN.pdf",width = 10)
ggplot(GO_enrichment_level0_Shared_sorted[1:25,], aes(x=reorder(title,fold_enrichment), y=fold_enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",number_in_list,sep="")), nudge_y=4, color="black",size=3) + 
  theme(axis.text = element_text(size = 13)) +
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

pdf("~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_genes_BORUTA_and_SPEARMAN.shared.PANTHER.GO.level0.all_terms.pdf",width = 10,height = 10)
ggplot(GO_enrichment_level0_Shared_sorted, aes(x=reorder(title,fold_enrichment), y=fold_enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",number_in_list,sep="")), nudge_y=4, color="black",size=3) + 
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

# Plot GO enrichment - all genes
GO_enrichment_level0_all_file="~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_genes_BORUTA_and_SPEARMAN.all.PANTHER.GO.level0.txt"
GO_enrichment_level0_all=read.delim(file = GO_enrichment_level0_all_file,sep = "\t",header = T,stringsAsFactors = F)
GO_enrichment_level0_all_sorted=GO_enrichment_level0_all[with(GO_enrichment_level0_all, order(number_in_list, fold_enrichment,decreasing = T)),]

library(ggplot2)
pdf("~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_genes_BORUTA_and_SPEARMAN.all.PANTHER.GO.level0.Top25_byN.pdf",width = 10)
# ggplot(GO_enrichment_level0_all_sorted[1:25,], aes(x=reorder(title,number_in_list), y=fold_enrichment,fill=FDR)) + 
ggplot(GO_enrichment_level0_all_sorted[1:25,], aes(x=reorder(title,fold_enrichment), y=fold_enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",number_in_list,sep="")), nudge_y=3, color="black",size=3) + 
  theme(axis.text = element_text(size = 13)) +
 # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

pdf("~/Dropbox/PostDoc/Projects/16S_array/Host_response/selected_genes_BORUTA_and_SPEARMAN.all.PANTHER.GO.level0.all_terms.pdf",width = 10,height = 15)
ggplot(GO_enrichment_level0_all_sorted, aes(x=reorder(title,fold_enrichment), y=fold_enrichment,fill=FDR)) + 
  geom_bar(stat = "identity") +
  coord_flip() + theme_bw() +
  geom_text(aes(label = paste("n=",number_in_list,sep="")), nudge_y=4, color="black",size=3) + 
  # geom_text(aes(label = signif(fold_enrichment,3)), nudge_y=-3, color="white",size=2) + 
  ylab("Fold enrichment") + xlab("GO term")
dev.off()

####### plot all sections for expression localG
## all
dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated
dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Bacteria+dataset_to_plot$sum_sections_Fungi,
                                      pmax(dataset_to_plot$r_LocalG.Bacteria.max_value,dataset_to_plot$r_LocalG.Fungi.max_value),
                                      decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_gene.all_sections.pdf"))
pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_gene_BORUTA_AND_significant_0.01_spearman.all_sections.pdf"))
pdf(pdf_file,width = 63,height = 21) 


## defense response
# dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[grepl(x=final_selected_genes_from_all_sections_with_all_details_annotated$Gene.ontology..biological.process.,pattern = "defense"),]
# dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Bacteria+dataset_to_plot$sum_sections_Fungi,
#                                      pmax(dataset_to_plot$r_LocalG.Bacteria.max_value,dataset_to_plot$r_LocalG.Fungi.max_value),
#                                      decreasing = T),]
## pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/defense_response_selected_gene.all_sections.pdf"))
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/defense_response_selected_gene_BORUTA_AND_significant_0.01_spearman.all_sections.pdf"))
# pdf(pdf_file,width = 63,height = 21)

cor_type="r_LocalG"
all_samples=c("OMNI12_A1","OMNI12_A2","OMNI12_B2","OMNI12_B1","OMNI12_C1","OMNI12_C2",
              "OMNI13_A1","OMNI13_A2","OMNI13_B1","OMNI13_B2","OMNI13_C1","OMNI13_C2","OMNI13_D1")

# we need to plot the Bacteria/Fungi data only once, the only thing that change is the r with the given gene
# so we store it in a different list
Bacteria_Fungi_plots=list()
for (tax in c("Bacteria","Fungi")) {
  for (sample_to_plot in all_samples)
  {
    layer_to_plot=NA
    if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
    if (tax=="Fungi"){layer_to_plot="all_Fungi"}
    if (cor_type == "r_LocalG")
    {
      Bacteria_Fungi_plots[[paste(sample_to_plot,tax,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(sample_to_plot," - ",tax,sep=""))
    }
    if (cor_type== "r_Reads")
    {
      Bacteria_Fungi_plots[[paste(sample_to_plot,tax,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(sample_to_plot," - ",tax,sep=""))
    }
  } # end samples 
}   # end {Baccteria|Fungi} 


for (i in 1:NROW(dataset_to_plot))
{
  # final_selected_genes_from_all_sections_with_all_details_annotated[1,]
  # nPosFungi=NA
  # nPosBac=NA
  nFungi=dataset_to_plot$sum_sections_Fungi[i]
  nBac=dataset_to_plot$sum_sections_Bacteria[i]
  gene_name=dataset_to_plot$Entry.name[i]
  response_type=dataset_to_plot$response_type[i]
  
  gene_id=dataset_to_plot[i,"gene"]
  # dataset_to_plot[1,c("OMNI12_A1.sum.Bacteria","OMNI12_A1.sum.Fungi","OMNI12_A1.gene_sum_Reads.Bacteria","OMNI12_A1.gene_sum_Reads.Fungi")] 
  plots_list=list()
  
  # first all expression plots
  # ONLY THE SECTIONS HAVING EXPRESSION FOR THE GENE
  samples_with_expression=c()
  for (sample_to_plot in all_samples)
  {
    gene_sum_reads=sum(all_samples_spatial_data[[sample_to_plot]][[gene_id]])
    if (gene_sum_reads>0)
    {
      if (cor_type == "r_LocalG")
      {
        plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot," total reads:",gene_sum_reads,sep="")) # the exression
      }
      if (cor_type== "r_Reads")
      {
        plots_list[[paste(sample_to_plot,gene_id,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot," total reads:",gene_sum_reads,sep="")) # the expression
      }
      samples_with_expression=c(samples_with_expression,sample_to_plot)
    }
  }
# ONLY THE SECTIONS HAVING EXPRESSION FOR THE GENE
  # Now pull the data for Bacteria / Fungi
  
  for (tax in c("Bacteria","Fungi")) {
    # for (sample_to_plot in all_samples)
    for (sample_to_plot in samples_with_expression)
    {
      cor_value=dataset_to_plot[i,paste(sample_to_plot,cor_type,tax,sep=".")]
      boruta_selected="YES"
      if (is.na(cor_value)) 
      {
        boruta_selected="NO"
      } else {
        cor_value=signif(cor_value,3)
      }
      plots_list[[paste(sample_to_plot,tax,sep=".")]]=Bacteria_Fungi_plots[[paste(sample_to_plot,tax,sep=".")]] + tm_layout(main.title = paste ("r=",cor_value," BORUTA: ",boruta_selected))
    } # end samples 
  }   # end {Baccteria|Fungi} 
  # plots_list
  # https://github.com/r-tmap/tmap/issues/511
  
  # current.mode <- tmap_mode("plot")
  ## tmap_arrange(plots_list, widths = c(.75, .75))
  if (length(samples_with_expression)==length(all_samples))
  {
    plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
    plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
    plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
    plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
    for (i in 5:length(all_samples)) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
  }
  print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = length(samples_with_expression), nrow = 3))
  # tmap_mode(current.mode)
  if (i%%10==0) {message(paste("Finished ploting",i,"out of",NROW(dataset_to_plot),sep=" "))}
  flush.console()
}

dev.off()
# here

### plot the same section in all layers - determined by the hotspot Fungi based -> tmp - later choose the best
# defense related genes
dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[grepl(x=final_selected_genes_from_all_sections_with_all_details_annotated$Gene.ontology..biological.process.,pattern = "defense"),]
dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Bacteria+dataset_to_plot$sum_sections_Fungi,
                                      pmax(dataset_to_plot$r_LocalG.Bacteria.max_value,dataset_to_plot$r_LocalG.Fungi.max_value),
                                      decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/defense_response_selected_gene.max_Fungi_localG_in_all.TESTS.pdf"))
pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/defense_response_selected_gene_BORUTA_AND_significant_0.01_spearman.max_Fungi_localG_in_all.TESTS.pdf"))
pdf(pdf_file,width = 21,height = 14) 

for (i in 1:NROW(dataset_to_plot))
{
  # final_selected_genes_from_all_sections_with_all_details_annotated[1,]
  # nPosFungi=NA
  # nPosBac=NA
  nFungi=dataset_to_plot$sum_sections_Fungi[i]
  nBac=dataset_to_plot$sum_sections_Bacteria[i]
  gene_name=dataset_to_plot$Entry.name[i]
  response_type=dataset_to_plot$response_type[i]
  
  gene_id=dataset_to_plot[i,"gene"]
  plots_list=list()
  
  cor_sign="max" # for positive correlation
  
  if (dataset_to_plot[[paste("r_LocalG.Fungi","total_pos_r",sep=".")]][i]<
      dataset_to_plot[[paste("r_LocalG.Fungi","total_neg_r",sep=".")]][i]) {
    cor_sign="min" # mostly negative correlations
  }
  # take the sample name to plot
  sample_to_plot=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                         pattern = paste("r_LocalG.Fungi",".",cor_sign,"_sample",sep=""))]
  
  for (cor_type in c("r_LocalG","r_Reads"))
  {
    for (tax in c("Bacteria","Fungi")) {
      plot_name=paste(tax,cor_type,sep=".")
      
      
      cor_value=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                        pattern = paste(cor_type,".",tax,".",cor_sign,"_value",sep=""))]
      cor_value=signif(cor_value,3)
      
      
      layer_to_plot=NA
      if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
      if (tax=="Fungi"){layer_to_plot="all_Fungi"}
      
      if (cor_type == "r_LocalG")
      {
        plots_list[[plot_name]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,sep="")) # the exression
      }
      if (cor_type== "r_Reads")
      {
        plots_list[[plot_name]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,sep="")) # the expression
      }
    } # end {Baccteria|Fungi}
  } # end {LocalG | Reads}
  # plots_list
  
  # https://github.com/r-tmap/tmap/issues/511
  
  # current.mode <- tmap_mode("plot")
  ## tmap_arrange(plots_list, widths = c(.75, .75))
  
  plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
  plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
  plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
  plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
  
  # for (i in 3:4) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
  print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 4, nrow = 2))
  # tmap_mode(current.mode)
}

dev.off()

# plot the same section in all layers - determined by the hotspot Bacteria based -> tmp - later choose the best
# defense related genes
dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[grepl(x=final_selected_genes_from_all_sections_with_all_details_annotated$Gene.ontology..biological.process.,pattern = "defense"),]
dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Bacteria+dataset_to_plot$sum_sections_Fungi,
                                      pmax(dataset_to_plot$r_LocalG.Bacteria.max_value,dataset_to_plot$r_LocalG.Fungi.max_value),
                                      decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/defense_response_selected_gene.max_Bacteria_localG_in_all.TESTS.pdf"))
pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/defense_response_selected_gene_BORUTA_AND_significant_0.01_spearman.max_Bacteria_localG_in_all.TESTS.pdf"))


pdf(pdf_file,width = 21,height = 14) 

for (i in 1:NROW(dataset_to_plot))
{
  # final_selected_genes_from_all_sections_with_all_details_annotated[1,]
  # nPosFungi=NA
  # nPosBac=NA
  nFungi=dataset_to_plot$sum_sections_Fungi[i]
  nBac=dataset_to_plot$sum_sections_Bacteria[i]
  gene_name=dataset_to_plot$Entry.name[i]
  response_type=dataset_to_plot$response_type[i]
  
  gene_id=dataset_to_plot[i,"gene"]
  plots_list=list()
  
  cor_sign="max" # for positive correlation
  
  if (dataset_to_plot[[paste("r_LocalG.Bacteria","total_pos_r",sep=".")]][i]<
      dataset_to_plot[[paste("r_LocalG.Bacteria","total_neg_r",sep=".")]][i]) {
    cor_sign="min" # mostly negative correlations
  }
  # take the sample name to plot
  sample_to_plot=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                         pattern = paste("r_LocalG.Bacteria",".",cor_sign,"_sample",sep=""))]
  
  for (cor_type in c("r_LocalG","r_Reads"))
  {
    for (tax in c("Bacteria","Fungi")) {
      plot_name=paste(tax,cor_type,sep=".")
      
      
      cor_value=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                        pattern = paste(cor_type,".",tax,".",cor_sign,"_value",sep=""))]
      cor_value=signif(cor_value,3)
      
      
      layer_to_plot=NA
      if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
      if (tax=="Fungi"){layer_to_plot="all_Fungi"}
      
      if (cor_type == "r_LocalG")
      {
        plots_list[[plot_name]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,sep="")) # the exression
      }
      if (cor_type== "r_Reads")
      {
        plots_list[[plot_name]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,sep="")) # the expression
      }
    } # end {Baccteria|Fungi}
  } # end {LocalG | Reads}
  # plots_list
  
  # https://github.com/r-tmap/tmap/issues/511
  
  # current.mode <- tmap_mode("plot")
  ## tmap_arrange(plots_list, widths = c(.75, .75))
  
  plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
  plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
  plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
  plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
  
  # for (i in 3:4) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
  print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 4, nrow = 2))
  # tmap_mode(current.mode)
}

dev.off()


### NOT TESTED.... NEED TO UPDATE FROM HERE also _AND_significant_0.01_spearman

# select the section with the highest value among options and plot the same section for Bacteria/Fungi/Expression
for (i in 1:NROW(dataset_to_plot))
{
  # final_selected_genes_from_all_sections_with_all_details_annotated[1,]
  # nPosFungi=NA
  # nPosBac=NA
  nFungi=dataset_to_plot$sum_sections_Fungi[i]
  nBac=dataset_to_plot$sum_sections_Bacteria[i]
  gene_name=dataset_to_plot$Entry.name[i]
  response_type=dataset_to_plot$response_type[i]
  gene_id=dataset_to_plot[i,"gene"]
  
  # select the max section among all options {localG/reads in Bacteria/Fungi}
  max_cor_type=NA
  max_tax=NA
  cor_sign=NA
  # compare the value among all types
  # are there more positive or negative correlation? chooose the max
  all_cor_signs=dataset_to_plot[i,grepl(x=names(dataset_to_plot),pattern = paste("total_pos_r|total_neg_r",sep=""))] 
  max_sign=which.max(all_cor_signs)
  if (grepl(x=names(max_sign),"pos")){cor_sign="max"}
  else {cor_sign="min"}
      
  # compare the r values
  all_cor_values=dataset_to_plot[i,grepl(x=names(dataset_to_plot),pattern = paste(cor_sign,"_value",sep=""))] 
  max_cor_value=which.max(abs(all_cor_values))
  tax_and_type=gsub(x=names(max_cor_value),pattern = paste("_value",sep=""),replacement = "") # will equalent to: cor_type,".",tax,".",cor_sign,
  sample_to_plot=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                         pattern = paste(tax_and_type,"_sample",sep=""))]
  
  
  # TO CONTINUE....
  
  
  
  plots_list=list()
  for (cor_type in c("r_LocalG","r_Reads"))
  {
    for (tax in c("Bacteria","Fungi")) {
      plot_name=paste(tax,cor_type,sep=".")
      # 
      # cor_sign="max" # for positive correlation
      # 
      # if (dataset_to_plot[[paste(cor_type,tax,"total_pos_r",sep=".")]][i]<
      #     dataset_to_plot[[paste(cor_type,tax,"total_neg_r",sep=".")]][i]) {
      #   cor_sign="min" # mostly negative correlations
      # }
      # # take the sample name to plot
      # sample_to_plot=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
      #                                        pattern = paste(cor_type,".",tax,".",cor_sign,"_sample",sep=""))]
      # cor_value=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
      #                                   pattern = paste(cor_type,".",tax,".",cor_sign,"_value",sep=""))]
      # cor_value=signif(cor_value,3)
      
      
      layer_to_plot=NA
      if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
      if (tax=="Fungi"){layer_to_plot="all_Fungi"}
      
      if (cor_type == "r_LocalG")
      {
        plots_list[[plot_name]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,sep="")) # the exression
      }
      if (cor_type== "r_Reads")
      {
        plots_list[[plot_name]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,sep="")) # the expression
      }
    } # end {Baccteria|Fungi}
  } # end {LocalG | Reads}
  # plots_list
  
  # https://github.com/r-tmap/tmap/issues/511
  
  # current.mode <- tmap_mode("plot")
  ## tmap_arrange(plots_list, widths = c(.75, .75))
  
  plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
  plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
  plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
  plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
  
  # for (i in 3:4) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
  print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 4, nrow = 2))
  # tmap_mode(current.mode)
}
  
  
  
  
  tax_and_type=strsplit(tax_and_type,".",fixed = T)[[1]]
  sample_to_plot=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                         pattern = paste(cor_type,".",tax,".",cor_sign,"_sample",sep=""))]
  
  
  
                    names(dataset_to_plot)[grepl(x=names(dataset_to_plot),pattern = paste("_value",sep=""))]  
      # take the sample name to plot
      sample_to_plot=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                             pattern = paste(cor_type,".",tax,".",cor_sign,"_sample",sep=""))]
      cor_value=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                        pattern = paste(cor_type,".",tax,".",cor_sign,"_value",sep=""))]
      cor_value=signif(cor_value,3)
    }
  }
  
  plots_list=list()
  for (cor_type in c("r_LocalG","r_Reads"))
  {
    for (tax in c("Bacteria","Fungi")) {
      plot_name=paste(tax,cor_type,sep=".")
      
      cor_sign="max" # for positive correlation
      
      if (dataset_to_plot[[paste(cor_type,tax,"total_pos_r",sep=".")]][i]<
          dataset_to_plot[[paste(cor_type,tax,"total_neg_r",sep=".")]][i]) {
        cor_sign="min" # mostly negative correlations
      }
      # take the sample name to plot
      sample_to_plot=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                             pattern = paste(cor_type,".",tax,".",cor_sign,"_sample",sep=""))]
      cor_value=dataset_to_plot[i,grepl(x=names(dataset_to_plot),
                                        pattern = paste(cor_type,".",tax,".",cor_sign,"_value",sep=""))]
      cor_value=signif(cor_value,3)
      
      
      layer_to_plot=NA
      if (tax=="Bacteria"){layer_to_plot="all_Bacteria"}
      if (tax=="Fungi"){layer_to_plot="all_Fungi"}
      
      if (cor_type == "r_LocalG")
      {
        plots_list[[plot_name]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,grid_size = 2,title = paste(gene_id," - ",sample_to_plot,sep="")) # the exression
      }
      if (cor_type== "r_Reads")
      {
        plots_list[[plot_name]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(plot_name,"\n",sample_to_plot," r=",cor_value,sep=""))
        plots_list[[paste(plot_name,gene_id,sep=".")]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = gene_id,title = paste(gene_id," - ",sample_to_plot,sep="")) # the expression
      }
    } # end {Baccteria|Fungi}
  } # end {LocalG | Reads}
  # plots_list
  
  # https://github.com/r-tmap/tmap/issues/511
  
  # current.mode <- tmap_mode("plot")
  ## tmap_arrange(plots_list, widths = c(.75, .75))
  
  plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
  plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
  plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
  plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
  
  # for (i in 3:4) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
  print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 4, nrow = 2))
  # tmap_mode(current.mode)
}

dev.off()



# add some random set of proteins / defense response not selected


# DRAFTS...
  
  # plot the expression
  
  # https://bookdown.org/lexcomber/brunsdoncomber2e/Ch3.html
  coords.tmp <- cbind(quakes$long, quakes$lat)
  # create the SpatialPointsDataFrame
  quakes.sp <- SpatialPointsDataFrame(coords.tmp, 
                                      data = data.frame(quakes), 
                                      proj4string = CRS("+proj=longlat "))
  # convert to sf
  quakes_sf <- st_as_sf(quakes.sp)
  
  
  
  
  tm_shape(all_samples_spatial_data[[sample_to_plot]]) + 
    tm_fill(layer_to_plot, 
            palette = "-RdBu",
            style = "pretty", title="G stat",n=5,legend.reverse=T) +
    tm_borders(alpha=.4) + tm_legend(legend.text.size=0.5,legend.title.size=0.6) +
    tm_layout(title = title)
  
  ...
} 

# current.mode <- tmap_mode("plot")
# tmap_arrange(w1, w2, w3, w4, widths = c(.25, .75))
# tmap_mode(current.mode)
# is it positive or negative correlation
{
  
}




plot_hotspots_maps=function (spatial_data,out_pdf)
{
  # pdf
  pdf (file = out_pdf,height = 7,width = 7)
  
  # calculate the localG stas
  for (i in 1:length(names(spatial_data)))
  {
    layer=names(spatial_data)[i]
    if (i%%10==0) {
      message(paste("\t--",i,"out of",length(names(spatial_data))))
    }

    # plot
    grid_size=2
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
    #  brewer.pal(n = 3, name = "RdBu")
    #  breaks = c(-20, -1.96, -1, 1, 1.96, 20)
    #  palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
    #  col = palette[cut(getisgrid$HOTSPOT, breaks)]
    #  plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
    #  legend("bottom", inset=.02, title="Z score (local spatial G[i] statistic)",
    #         legend =c("<-1.96","-1","1","1.96",">1.96"), 
    #         fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
    #         xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    
    pmap=tm_shape(getisgrid) + 
      tm_fill("HOTSPOT", 
              palette = "-RdBu",
              style = "pretty", title="G stat",n=5,legend.reverse=T) +
      tm_borders(alpha=.4) + tm_legend(legend.text.size=0.5,legend.title.size=0.6) +
      tm_layout(title = layer)
    print (pmap)
  }
  dev.off()
}
exp="OMNI12"
smpl="A1"
genes_to_plot=c("AT1G67090","AT5G54770","AT3G01500","AT5G38410","AT5G38430","AT2G42540","AT1G60950","AT1G42970","AT2G21660",
                "AT5G38420","AT2G39730","AT2G06520","AT1G06680","AT5G01530","AT1G29930","AT1G20340","AT3G63160","AT2G26500",
                "AT2G42530","AT3G60750",'AT4G03280',"AT3G21055","AT1G55670","AT1G61520","AT1G32060","AT1G31330","AT5G66570",
                "AT2G37220","AT4G38970","AT1G12900","AT3G14210","AT2G30570","AT4G10340","AT2G25510","AT2G21330","AT1G09340",
                "AT3G14420","AT4G37930","AT3G54890","AT4G01150","AT1G79040","AT1G44575","AT3G12780","AT3G26650")


hotspost_df=data.frame()
data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
data=read.delim(data_file,sep=";",stringsAsFactors = F)
data$all_Bacteria=rowSums(data[,1:50])
data$all_Fungi=rowSums(data[,51:100])

xy_list=strsplit(row.names(data), "x")
xy=as.data.frame(t(as.data.frame(xy_list)))
xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
names(xy)=c("x","y")
data_with_xy=cbind(xy,data)
rownames(data_with_xy)=rownames(data)
data_with_xy$x=as.numeric(data_with_xy$x)
data_with_xy$y=as.numeric(data_with_xy$y)

# Now with the expression
expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,".tsv",sep="")
expression_data=read.delim(expression_data_file,sep="\t",stringsAsFactors = F)
expression_data=t(expression_data) # row = position; col = gene
# Only genes with more than overall 10 raeds
gene_sum=data.frame(count=colSums(expression_data))
gene_sum$percent_of_data=gene_sum$count/(sum(gene_sum$count))
# sum(gene_sum$percent[gene_sum$count>=10])
# row.names(gene_sum)
gene_sum_ordered=gene_sum[order(gene_sum$count,decreasing = T),]
# create the spatial object
xy_list=strsplit(row.names(expression_data), "x")
xy=as.data.frame(t(as.data.frame(xy_list)))
xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
names(xy)=c("x","y")
expression_data_with_xy=cbind(xy,expression_data)
expression_data_with_xy$x=as.numeric(expression_data_with_xy$x)
expression_data_with_xy$y=as.numeric(expression_data_with_xy$y)

rownames(expression_data_with_xy)=rownames(expression_data)

Boruta_based_selected_Bacteria_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Bacteria_sum.",smpl,".tdl",sep="");
Boruta_based_selected_Bacteria=read.delim(file=Boruta_based_selected_Bacteria_file,sep = "\t",stringsAsFactors = F)

Boruta_based_selected_Fungi_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/Boruta_selected_reads_or_2x2_localG.210708_raw_counts_mtrbcpfiltered.RNA.with_Top50genus_Fungi_sum.",smpl,".tdl",sep="")
Boruta_based_selected_Fungi=read.delim(file=Boruta_based_selected_Fungi_file,sep = "\t",stringsAsFactors = F)

####### all Bacteria
genes_to_plot=c("AT1G06680","AT1G08380","AT1G09340","AT1G12900","AT1G16880","AT1G20020","AT1G20340","AT1G29910","AT1G29930","AT1G31580",
                "AT1G42970","AT1G52230","AT1G61520","AT1G67090","AT1G68010","AT1G72610","AT2G01020","AT2G06520","AT2G10940","AT2G16600",
                "AT2G21170","AT2G21660","AT2G26500","AT2G29630","AT2G34420","AT2G35370","AT2G37220","AT2G42530","AT2G42540","AT2G45180",
                "AT2G46820","AT2G47400","AT3G01500","AT3G09390","AT3G11630","AT3G12780","AT3G13750","AT3G14210","AT3G14415","AT3G14420",
                "AT3G21055","AT3G22120","AT3G22840","AT3G26070","AT3G26650","AT3G27850","AT3G41768","AT3G41979","AT3G44310","AT3G45140",
                "AT3G47470","AT3G53460","AT3G54890","AT3G55800","AT3G60750","AT3G63160","AT4G01150","AT4G03280","AT4G03520","AT4G04640",
                "AT4G10340","AT4G12800","AT4G14400","AT4G21280","AT4G24770","AT4G25050","AT4G26530","AT4G27700","AT4G28750","AT4G32260",
                "AT4G33010","AT4G37930","AT5G01530","AT5G12860","AT5G13630","AT5G13650","AT5G14740","AT5G20630","AT5G26742","AT5G30510",
                "AT5G38410","AT5G38420","AT5G38430","AT5G42530","AT5G54770","AT5G64040","AT5G66190","AT5G66570")
# subset genes
# RNA_cutoff=100 # Only genes with more than 100 reads are accounted
# expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",genes_to_plot)]
# Join the expression and all_Fungi all_Bacteria
joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
# create the spatial object
spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
out_pdf=paste(exp,".",smpl,".","BACTERIA_in_ALL_SECTIONS.hotspots.pdf",sep="")
plot_hotspots_maps(out_pdf = out_pdf,spatial_data = spatial_data)

# extract the subset of genes data
subset_r=Boruta_based_selected_Bacteria[Boruta_based_selected_Bacteria$gene %in% genes_to_plot ,]
out_r_bacteria=paste(exp,".",smpl,".Bacteria.","BACTERIA_in_ALL_SECTIONS.hotspots.r.tdl",sep="")
write.table(x=subset_r,file = out_r_bacteria,quote = F,sep = "\t",row.names = F)

subset_r=Boruta_based_selected_Fungi[Boruta_based_selected_Fungi$gene %in% genes_to_plot ,]
out_r_Fungi=paste(exp,".",smpl,"Fungi.","BACTERIA_in_ALL_SECTIONS.hotspots.r.tdl",sep="")
write.table(x=subset_r,file = out_r_Fungi,quote = F,sep = "\t",row.names = F)
######## all Fungi
genes_to_plot=c("AT1G06680","AT1G09340","AT1G12900","AT1G15405","AT1G16880","AT1G20020","AT1G20340","AT1G23310","AT1G29910",
                "AT1G29930","AT1G31580","AT1G42970","AT1G44575","AT1G55490","AT1G55670","AT1G67090","AT2G01010","AT2G01020",
                "AT2G06520","AT2G10940","AT2G21660","AT2G26500","AT2G29630","AT2G34420","AT2G35370","AT2G39730","AT2G42530",
                "AT2G42540","AT2G45180","AT2G47400","AT3G01500","AT3G11630","AT3G12780","AT3G13750","AT3G14415","AT3G21055",
                "AT3G22120","AT3G22840","AT3G26650","AT3G41768","AT3G41979","AT3G45140","AT3G47470","AT3G54890","AT3G55800",
                "AT3G60750","AT3G63160","AT4G01150","AT4G04640","AT4G05180","AT4G10340","AT4G12800","AT4G21280","AT4G26530",
                "AT4G28750","AT4G32260","AT4G33010","AT4G37930","AT5G01530","AT5G13630","AT5G20630","AT5G38410","AT5G38430",
                "AT5G42530","AT5G64040","AT5G66190")
# subset genes
# RNA_cutoff=100 # Only genes with more than 100 reads are accounted
# expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",genes_to_plot)]
# Join the expression and all_Fungi all_Bacteria
joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
# create the spatial object
spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
out_pdf=paste(exp,".",smpl,".","FUNGI_in_ALL_SECTIONS.hotspots.pdf",sep="")
plot_hotspots_maps(out_pdf = out_pdf,spatial_data = spatial_data)

# extract the subset of genes data
subset_r=Boruta_based_selected_Bacteria[Boruta_based_selected_Bacteria$gene %in% genes_to_plot ,]
out_r_bacteria=paste(exp,".",smpl,".Bacteria.","FUNGI_in_ALL_SECTIONS.hotspots.r.tdl",sep="")
write.table(x=subset_r,file = out_r_bacteria,quote = F,sep = "\t",row.names = F)

subset_r=Boruta_based_selected_Fungi[Boruta_based_selected_Fungi$gene %in% genes_to_plot ,]
out_r_Fungi=paste(exp,".",smpl,"Fungi.","FUNGI_in_ALL_SECTIONS.hotspots.r.tdl",sep="")
write.table(x=subset_r,file = out_r_Fungi,quote = F,sep = "\t",row.names = F)


# all_Bacteria+Fungi
genes_to_plot=c("AT1G06680","AT1G09340","AT1G12900","AT1G16880","AT1G20020","AT1G20340","AT1G29910",
                "AT1G29930","AT1G31580","AT1G42970","AT1G67090","AT2G01020","AT2G06520","AT2G10940",
                "AT2G21660","AT2G26500","AT2G29630","AT2G34420","AT2G35370","AT2G42530","AT2G42540",
                "AT2G45180","AT2G47400","AT3G01500","AT3G11630","AT3G12780","AT3G13750","AT3G14415",
                "AT3G21055","AT3G22120","AT3G22840","AT3G26650","AT3G41768","AT3G41979","AT3G45140",
                "AT3G47470","AT3G54890","AT3G55800","AT3G60750","AT3G63160","AT4G01150","AT4G04640",
                "AT4G10340","AT4G12800","AT4G21280","AT4G26530","AT4G28750","AT4G32260","AT4G33010",
                "AT4G37930","AT5G01530","AT5G13630","AT5G20630","AT5G38410","AT5G38430","AT5G42530",
                "AT5G64040","AT5G66190")

# subset genes
# RNA_cutoff=100 # Only genes with more than 100 reads are accounted
# expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",genes_to_plot)]
# Join the expression and all_Fungi all_Bacteria
joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
# create the spatial object
spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
out_pdf=paste(exp,".",smpl,".","BACTERIA_and_FUNGI_in_ALL_SECTIONS.hotspots.pdf",sep="")
plot_hotspots_maps(out_pdf = out_pdf,spatial_data = spatial_data)

# extract the subset of genes data
subset_r=Boruta_based_selected_Bacteria[Boruta_based_selected_Bacteria$gene %in% genes_to_plot ,]
out_r_bacteria=paste(exp,".",smpl,".Bacteria.","BACTERIA_and_FUNGI_in_ALL_SECTIONS.hotspots.r.tdl",sep="")
write.table(x=subset_r,file = out_r_bacteria,quote = F,sep = "\t",row.names = F)

subset_r=Boruta_based_selected_Fungi[Boruta_based_selected_Fungi$gene %in% genes_to_plot ,]
out_r_Fungi=paste(exp,".",smpl,"Fungi.","BACTERIA_and_FUNGI_in_ALL_SECTIONS.hotspots.r.tdl",sep="")
write.table(x=subset_r,file = out_r_Fungi,quote = F,sep = "\t",row.names = F)




# test
genes_to_plot=c("AT1G67090","AT5G54770","AT3G01500","AT5G38410","AT5G38430","AT2G42540","AT1G60950","AT1G42970","AT2G21660",
                "AT5G38420","AT2G39730","AT2G06520","AT1G06680","AT5G01530","AT1G29930","AT1G20340","AT3G63160","AT2G26500",
                "AT2G42530","AT3G60750",'AT4G03280',"AT3G21055","AT1G55670","AT1G61520","AT1G32060","AT1G31330","AT5G66570",
                "AT2G37220","AT4G38970","AT1G12900","AT3G14210","AT2G30570","AT4G10340","AT2G25510","AT2G21330","AT1G09340",
                "AT3G14420","AT4G37930","AT3G54890","AT4G01150","AT1G79040","AT1G44575","AT3G12780","AT3G26650")



  
  
  if (nrow(hotspost_df)==0)
  {
    hotspost_df=getisgrid.df
  } else {
    hotspost_df=merge(hotspost_df,getisgrid.df[,1:4],by="row.names",all=T)
    row.names(hotspost_df)=hotspost_df$Row.names
    hotspost_df=hotspost_df[,-1]
  }
}








      
      
      
      
      all_hits_PA$sum=rowSums(all_hits_PA)
      missed=row.names(all_hits_PA)[all_hits_PA$sum==1&all_hits_PA$target==1]
      missed_genus=q_data[q_data$V1 %in% missed,]
      missed_genus$sample=s
      
      
      all_genes=list(Bacterial_reads_cor=Bacterial_reads_cor_selected$gene,
                     Fungi_reads_cor=Fungi_reads_cor_selected$gene,
                     Bacterial_localG_cor=Bacterial_localG_cor_selected$gene,
                     Fungi_localG_cor=Fungi_localG_cor_selected$gene,
                     
                     Boruta_reads_Bacterial=Boruta_reads_Bacterial_selected$gene,
                     Boruta_reads_Fungi=Boruta_reads_Fungi_selected$gene,
                     
                     Boruta_localG_Bacterial=Boruta_localG_Bacterial_selected$gene,
                     Boruta_localG_Fungi=Boruta_localG_Fungi_selected$gene
                     )

      intersect_file=paste("/Volumes/hashkenazy/Stefania_Eduardo/evaluate_probes/Eduardo_samples_reads/MMSeqs2_vs_SILVA_138.1/ID_0.99.hits.R1_or_R2/probes_grep/",s,".0.99.hits.R1_or_R2.uniq_count.t_seq.probs_intersects.pdf",sep="")
      pdf(intersect_file,width = 14,height = 7)
      
      dev.off()
    }

    
    
    
    
    
    Bacterial_data$IncMSE_z=(Bacterial_data$X.IncMSE - mean(Bacterial_data$X.IncMSE))/sd(Bacterial_data$X.IncMSE)
    Fungi_data$IncMSE_z=(Fungi_data$X.IncMSE - mean(Fungi_data$X.IncMSE))/sd(Fungi_data$X.IncMSE)
    
    Bacterial_signif_Z=Bacterial_data[Bacterial_data$IncMSE_z>1.96,]
    names(Bacterial_signif_Z)[c(2,17)]=paste(S,"Bacteria",names(Bacterial_signif_Z)[c(2,17)],sep="_")
    
    Fungi_signif_Z=Fungi_data[Fungi_data$IncMSE_z>1.96,]
    names(Fungi_signif_Z)[c(2,17)]=paste(S,"Fungi",names(Fungi_signif_Z)[c(2,17)],sep="_")
    
    Fungi_cor_001=Fungi_cor_data[Fungi_cor_data$FDR.BH<0.01,]
    names(Fungi_cor_001)[c(3,6,8,9)]=paste(S,"Fungi",names(Fungi_cor_001)[c(3,6,8,9)],sep="_")
    
    Bacteria_cor_001=Bacterial_cor_data[Bacterial_cor_data$FDR.BH<0.01,]
    names(Bacteria_cor_001)[c(3,6,8,9)]=paste(S,"Bacteria",names(Bacteria_cor_001)[c(3,6,8,9)],sep="_")
    
    
    if (nrow(all_SignifVarImp)==0) {
      all_SignifVarImp=Bacterial_signif_Z[,c(1,2,17)]
      all_SignifVarImp=merge(all_SignifVarImp,Fungi_signif_Z[,c(1,2,17)],all=T,by="Row.names")
      all_SignifVarImp=merge(all_SignifVarImp,Bacteria_cor_001[,c(1,3,6,8,9)],all=T,by.x="Row.names",by.y="gene")
      all_SignifVarImp=merge(all_SignifVarImp,Fungi_cor_001[,c(1,3,6,8,9)],all=T,by.x="Row.names",by.y="gene")
    }
    else
    {
      all_SignifVarImp=merge(all_SignifVarImp,Bacterial_signif_Z[,c(1,2,17)],all=T,by="Row.names")
      all_SignifVarImp=merge(all_SignifVarImp,Fungi_signif_Z[,c(1,2,17)],all=T,by="Row.names")
      all_SignifVarImp=merge(all_SignifVarImp,Bacteria_cor_001[,c(1,3,6,8,9)],all=T,by.x="Row.names",by.y="gene")
      all_SignifVarImp=merge(all_SignifVarImp,Fungi_cor_001[,c(1,3,6,8,9)],all=T,by.x="Row.names",by.y="gene")
    }
    
    annotation=Fungi_data[,c(1,4:16)]
    annotation=rbind(annotation,Bacterial_data[,c(1,4:16)])
    annotation=unique(annotation)
    
  }
  
  all_SignifVarImp=merge(all_SignifVarImp,annotation,by="Row.names")
  all_SignifVarImp$Bacterial_RF_agreement=rowSums(!is.na(all_SignifVarImp[,grepl(x = names(all_SignifVarImp),pattern = "Bacteria_IncMSE_z")]))
  all_SignifVarImp$Fungi_RF_agreement=rowSums(!is.na(all_SignifVarImp[,grepl(x = names(all_SignifVarImp),pattern = "Fungi_IncMSE_z")]))
  all_SignifVarImp$Bacterial_Fungi_RF_agreement=all_SignifVarImp$Fungi_RF_agreement+all_SignifVarImp$Bacterial_RF_agreement
  
  all_SignifVarImp$Bacterial_cor_agreement=rowSums(!is.na(all_SignifVarImp[,grepl(x = names(all_SignifVarImp),pattern = "Bacteria_r")]))
  all_SignifVarImp$Fungi_cor_agreement=rowSums(!is.na(all_SignifVarImp[,grepl(x = names(all_SignifVarImp),pattern = "Fungi_r")]))
  all_SignifVarImp$Bacterial_Fungi_cor_agreement=all_SignifVarImp$Fungi_cor_agreement+all_SignifVarImp$Bacterial_cor_agreement
  
  # save
  out_sum_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria_Fungi.WithUnknown_probs.SUM_Top_50_genus.ntree1000.per_spot_RF.varImp_intersect_all_sections.csv",sep="")
  write.table(file = out_sum_file,x=all_SignifVarImp,quote = F,sep = ";",row.names = F)
  
  # add_worsheet
  name=paste(Exp,"reads",sep="_")
  addWorksheet(wb, name)
  writeData(wb, name, all_SignifVarImp)
  
  
  ###### With the Hotspots data
  all_SignifVarImp_HotsPots=data.frame()
  
  for (S in samples)
  {
    Fungi_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/hotspots/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.2x2.Getis_Ord.",S,".1.Getis_Ord_2x2_G.ntree1000.per_spot_RF.varImp.tsv",sep="");
    Bacteria_File=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/hotspots/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.2x2.Getis_Ord.",S,".1.Getis_Ord_2x2_G.ntree1000.per_spot_RF.varImp.tsv",sep="");
    
    Fungi_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Fungi_corr/WithUnknown_probs/hotspots/210708_raw_counts_mtrbcpfiltered.RNA.with_Fungi.2x2.Getis_Ord.",S,".corr.M5_R10.Sprot_annotated.tsv",sep="");
    Bacteria_cor_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/RNA_Bacteria_corr/WithUnknown_probs/hotspots/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria.2x2.Getis_Ord.",S,".corr.M5_R10.Sprot_annotated.tsv",sep="");
    
    if (file.exists(Fungi_file)) {
      Fungi_data=read.delim(file = Fungi_file,sep="\t",stringsAsFactors = F)
    }
    if (file.exists(Bacteria_File)) {
      Bacterial_data=read.delim(file = Bacteria_File,sep="\t",stringsAsFactors = F)
    }
    if (file.exists(Fungi_cor_file)) {
      Fungi_cor_data=read.delim(file = Fungi_cor_file,sep="\t",stringsAsFactors = F)
    }
    if (file.exists(Bacteria_cor_file)) {
      Bacterial_cor_data=read.delim(file = Bacteria_cor_file,sep="\t",stringsAsFactors = F)
    }
    
    if (nrow (Bacterial_data)>0) {
      Bacterial_data$IncMSE_z=(Bacterial_data$X.IncMSE - mean(Bacterial_data$X.IncMSE))/sd(Bacterial_data$X.IncMSE)
      Bacterial_signif_Z=Bacterial_data[Bacterial_data$IncMSE_z>1.96,]
      names(Bacterial_signif_Z)[c(2,17)]=paste(S,"Bacteria",names(Bacterial_signif_Z)[c(2,17)],sep="_")
      annotation=rbind(annotation,Bacterial_data[,c(1,4:16)])
    }
    
    if (nrow (Fungi_data)>0) {
      Fungi_data$IncMSE_z=(Fungi_data$X.IncMSE - mean(Fungi_data$X.IncMSE))/sd(Fungi_data$X.IncMSE)
      Fungi_signif_Z=Fungi_data[Fungi_data$IncMSE_z>1.96,]
      names(Fungi_signif_Z)[c(2,17)]=paste(S,"Fungi",names(Fungi_signif_Z)[c(2,17)],sep="_")
      annotation=rbind(annotation,Fungi_data[,c(1,4:16)])
    }
    
    
    if (nrow (Fungi_cor_data)>0)
    {
      Fungi_cor_001=Fungi_cor_data[Fungi_cor_data$FDR.BH<0.01,]
      names(Fungi_cor_001)[c(3,6,8,9)]=paste(S,"Fungi",names(Fungi_cor_001)[c(3,6,8,9)],sep="_")
      names(Fungi_cor_data)[1]="Row.names"
      annotation=rbind(annotation,Fungi_cor_data[,c(1,10:22)])
    }
    if (nrow (Bacterial_cor_data)>0)
    {
      Bacteria_cor_001=Bacterial_cor_data[Bacterial_cor_data$FDR.BH<0.01,]
      names(Bacteria_cor_001)[c(3,6,8,9)]=paste(S,"Bacteria",names(Bacteria_cor_001)[c(3,6,8,9)],sep="_")
      names(Bacterial_cor_data)[1]="Row.names"
      annotation=rbind(annotation,Bacterial_cor_data[,c(1,10:22)])
    }
    
    if (nrow(all_SignifVarImp_HotsPots)==0) {
      all_SignifVarImp_HotsPots=Bacterial_signif_Z[,c(1,2,17)]
      all_SignifVarImp_HotsPots=merge(all_SignifVarImp_HotsPots,Fungi_signif_Z[,c(1,2,17)],all=T,by="Row.names")
      all_SignifVarImp_HotsPots=merge(all_SignifVarImp_HotsPots,Bacteria_cor_001[,c(1,3,6,8,9)],all=T,by.x="Row.names",by.y="gene")
      all_SignifVarImp_HotsPots=merge(all_SignifVarImp_HotsPots,Fungi_cor_001[,c(1,3,6,8,9)],all=T,by.x="Row.names",by.y="gene")
    } else {
      all_SignifVarImp_HotsPots=merge(all_SignifVarImp_HotsPots,Bacterial_signif_Z[,c(1,2,17)],all=T,by="Row.names")
      all_SignifVarImp_HotsPots=merge(all_SignifVarImp_HotsPots,Fungi_signif_Z[,c(1,2,17)],all=T,by="Row.names")
      all_SignifVarImp_HotsPots=merge(all_SignifVarImp_HotsPots,Bacteria_cor_001[,c(1,3,6,8,9)],all=T,by.x="Row.names",by.y="gene")
      all_SignifVarImp_HotsPots=merge(all_SignifVarImp_HotsPots,Fungi_cor_001[,c(1,3,6,8,9)],all=T,by.x="Row.names",by.y="gene")
    }
    
    annotation=unique(annotation)
    
    rm("Bacteria_cor_file","Bacteria_File","Bacterial_cor_data","Bacterial_data","Fungi_cor_data","Fungi_cor_file","Fungi_data","Fungi_file","Bacteria_cor_001","Bacterial_signif_Z","Fungi_cor_001","Fungi_signif_Z")
    
  }
}


all_SignifVarImp_HotsPots=merge(all_SignifVarImp_HotsPots,annotation,by="Row.names")
all_SignifVarImp_HotsPots$Bacterial_RF_agreement=rowSums(!is.na(all_SignifVarImp_HotsPots[,grepl(x = names(all_SignifVarImp_HotsPots),pattern = "Bacteria_IncMSE_z")]))
all_SignifVarImp_HotsPots$Fungi_RF_agreement=rowSums(!is.na(all_SignifVarImp_HotsPots[,grepl(x = names(all_SignifVarImp_HotsPots),pattern = "Fungi_IncMSE_z")]))
all_SignifVarImp_HotsPots$Bacterial_Fungi_RF_agreement=all_SignifVarImp_HotsPots$Fungi_RF_agreement+all_SignifVarImp_HotsPots$Bacterial_RF_agreement

all_SignifVarImp_HotsPots$Bacterial_cor_agreement=rowSums(!is.na(all_SignifVarImp_HotsPots[,grepl(x = names(all_SignifVarImp_HotsPots),pattern = "Bacteria_r")]))
all_SignifVarImp_HotsPots$Fungi_cor_agreement=rowSums(!is.na(all_SignifVarImp_HotsPots[,grepl(x = names(all_SignifVarImp_HotsPots),pattern = "Fungi_r")]))
all_SignifVarImp_HotsPots$Bacterial_Fungi_cor_agreement=all_SignifVarImp_HotsPots$Fungi_cor_agreement+all_SignifVarImp_HotsPots$Bacterial_cor_agreement

# save
hotspots_out_sum_file=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria_Fungi.WithUnknown_probs.SUM_Top_50_genus.2x2.Getis_Ord.ntree1000.per_spot_RF.varImp_intersect_all_sections.csv",sep="")
write.table(file = hotspots_out_sum_file,x=all_SignifVarImp_HotsPots,quote = F,sep = ";",row.names = F)

# add_worsheet
name=paste(Exp,"G.O_2x2",sep="_")
addWorksheet(wb, name)
writeData(wb, name, all_SignifVarImp_HotsPots)

# Join both
names(all_SignifVarImp_HotsPots)[c(2:73,87:92)]=paste("HS",names(all_SignifVarImp_HotsPots)[c(2:73,87:92)],sep="_")
names(all_SignifVarImp)[c(2:73,87:92)]=paste("Reads",names(all_SignifVarImp)[c(2:73,87:92)],sep="_")

all_all=merge(all_SignifVarImp[,c(1:73,87:92)],all_SignifVarImp_HotsPots[,c(1:73,87:92)],by="Row.names",all=T)
all_all=merge(all_all,annotation,by="Row.names",all=T)

name=paste(Exp,"G.O_2x2_and_reads",sep="_")
addWorksheet(wb, name)
writeData(wb, name, all_all)

# After running the 2 Exps save the excel file
out_xlsx=paste("/Volumes//spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/210708_raw_counts_mtrbcpfiltered.RNA.with_Bacteria_Fungi.WithUnknown_probs.SUM_Top_50_genus.2x2.Getis_Ord_and_reads.ntree1000.per_spot_RF.varImp_intersect.xlsx",sep="")
saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)

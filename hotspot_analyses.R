# new test based on: https://michaelminn.net/tutorials/r-point-analysis/index.html
library(sp)
library(raster)
library(spdep)
library(openxlsx)
library(RColorBrewer)
library(classInt)

Getis_Ord_GI_per_spot = function (spatial_data,layer) # layer=genus_name
{
  # based on grid
  # pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
  # box = round(extent(spatial_data) / pixelsize) * pixelsize
  # template = raster(box, crs = spatial_data,
  #                   nrows = (box@ymax - box@ymin) / pixelsize, 
  #                   ncols = (box@xmax - box@xmin) / pixelsize)
  # getisraster = rasterize(spatial_data, template, field = genus_name, fun = sum)
  # getisgrid = rasterToPolygons(getisraster)
  ## plot(getisgrid)
  ## Create the list of neighbors
  # neighbors = poly2nb(getisgrid)
  # weighted_neighbors = nb2listw(neighbors, zero.policy=T)
  
  ## Perform the local G analysis (Getis-Ord GI*)
  # getisgrid$HOTSPOT = as.vector(localG(getisgrid$layer, weighted_neighbors))
  # globalMoran <- moran.test(getisgrid$layer, weighted_neighbors)
  
  ### for each spot separately
  # Create the list of neighbors
  # based on https://gis.stackexchange.com/questions/262887/converting-spatial-points-to-neighbours-list-using-r
  xy=coordinates(spatial_data)
  neighbors=knn2nb(knearneigh(xy))
  weighted_neighbors = nb2listw(neighbors, zero.policy=T)
  
  plot(OA.Census, border = 'lightgrey')
  plot(nb, coordinates(OA.Census), add=TRUE, col = 'red')
  ## Perform the local G analysis (Getis-Ord GI*)
  HOTSPOT=as.vector(localG(spatial_data[[layer]], weighted_neighbors))
  
  spatial_hotspots=SpatialPointsDataFrame(xy, data=data.frame(HOTSPOT=HOTSPOT))
  return (spatial_hotspots)
}

Getis_Ord_GI_per_grid = function (grid_size,spatial_data,genus_name)
{
  pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
  box = round(extent(spatial_data) / pixelsize) * pixelsize
  template = raster(box, crs = spatial_data,
                    nrows = (box@ymax - box@ymin) / pixelsize, 
                    ncols = (box@xmax - box@xmin) / pixelsize)
  getisraster = rasterize(spatial_data, template, field = genus_name, fun = sum)
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


globalMoran_per_spot = function (grid_size,spatial_data,genus_name)
{
  # https://rpubs.com/quarcs-lab/spatial-autocorrelation
  # pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
  # box = round(extent(spatial_data) / pixelsize) * pixelsize
  # template = raster(box, crs = spatial_data,
  #                   nrows = (box@ymax - box@ymin) / pixelsize, 
  #                   ncols = (box@xmax - box@xmin) / pixelsize)
  # getisraster = rasterize(spatial_data, template, field = genus_name, fun = sum)
  # getisgrid = rasterToPolygons(getisraster)
  # plot(getisgrid)
  ## Create the list of neighbors
  # neighbors = poly2nb(getisgrid)
  # weighted_neighbors = nb2listw(neighbors, zero.policy=T)
  
  # based on https://gis.stackexchange.com/questions/262887/converting-spatial-points-to-neighbours-list-using-r
  xy=coordinates(spatial_data)
  neighbors=knn2nb(knearneigh(xy))
  weighted_neighbors = nb2listw(neighbors, zero.policy=T)
  
  ## Perform the globalMoran
  # globalMoran <- moran.test(getisgrid$layer, weighted_neighbors)
  # moran <- moran.plot(getisgrid$layer, listw = weighted_neighbors)
  
  # Perform the globalMoran
  globalMoran <- moran.test(spatial_data[[genus_name]], weighted_neighbors)
  return (globalMoran)
}
globalMoran_per_grid = function (grid_size,spatial_data,genus_name)
{
  # https://rpubs.com/quarcs-lab/spatial-autocorrelation
  pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
  box = round(extent(spatial_data) / pixelsize) * pixelsize
  template = raster(box, crs = spatial_data,
                    nrows = (box@ymax - box@ymin) / pixelsize, 
                    ncols = (box@xmax - box@xmin) / pixelsize)
  getisraster = rasterize(spatial_data, template, field = genus_name, fun = sum)
  getisgrid = rasterToPolygons(getisraster)
  # plot(getisgrid)
  # Create the list of neighbors
  neighbors = poly2nb(getisgrid)
  weighted_neighbors = nb2listw(neighbors, zero.policy=T)
  
  # Perform the globalMoran
  globalMoran <- moran.test(getisgrid$layer, weighted_neighbors)
  # moran <- moran.plot(getisgrid$layer, listw = weighted_neighbors)
  return (globalMoran)
}
# library("gridExtra") 
for (smpl in c("A1","A2","B1","B2","C1","C2"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.pdf",sep="")
  pdf(out_pdfs,width = 7,height = 7)
  out_xls=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.xlsx",sep="")
  XSLX_Obj=createWorkbook(out_xls)
  
  out_summary_overlaps=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.Bacteria_Fungi_overlap.csv",sep="")

  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  for (grid_size in c(2:15,24))
  {
    i=0
    hotspost_df=data.frame()
    global_moran_all_Bacteria=NA
    global_moran_all_Fungi=NA
    for (genus_name in names(spatial_data)) # loop all genus and specific grid size
    {
      i=i+1
      type="Bacteria"
      if (i>50&i<100) {type="Fungi"}
      if (i>100) {type=""}
      total=sum(spatial_data[[genus_name]])
      if(total>5)
      {
        if (grid_size==3) # plot only for the first time
        {
          # http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS3_MakingMaps_part1_mappingVectorData.html
          # spatial_data_for_plot=spatial_data
          # spatial_data_for_plot[[genus_name]][spatial_data_for_plot[[genus_name]]==0]=NA
          # breaks.qt <- classIntervals(spatial_data_for_plot[[genus_name]], n = 5, style = "equal", intervalClosure = "right")
          # my.palette <- brewer.pal(n = 8, name = "PuRd")
          # p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name," - ",type," ",i),cex=1),cuts=9,at=c(-1,1,breaks.qt$brks[2:length(breaks.qt$brks)]))
          p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name," - ",type," ",i),cex=1))
          plot(p)
          rm(p)
        }
  
        # Create a regular grid of 1/8-mile-square crime areas via a raster
        # library(raster)
        message(paste("\t- genus ",genus_name," grid size=",grid_size,sep=""))
        pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
        box = round(extent(spatial_data) / pixelsize) * pixelsize
        template = raster(box, crs = spatial_data,
                          nrows = (box@ymax - box@ymin) / pixelsize, 
                          ncols = (box@xmax - box@xmin) / pixelsize)
        getisraster = rasterize(spatial_data, template, field = genus_name, fun = sum)
        getisgrid = rasterToPolygons(getisraster)
        # plot(getisgrid)
        # Create the list of neighbors
        neighbors = poly2nb(getisgrid)
        weighted_neighbors = nb2listw(neighbors, zero.policy=T)
        
        # Perform the local G analysis (Getis-Ord GI*)
        getisgrid$HOTSPOT = as.vector(localG(getisgrid$layer, weighted_neighbors))
        # Color the grid cells based on the z-score
        # join with the others
        getisgrid.df=as.data.frame(getisgrid)
        names(getisgrid.df)=paste(names(getisgrid.df),genus_name,sep="_")
        getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
        names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
        getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
        getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
        row.names(getisgrid.df)=getisgrid.df$xy
        getisgrid.df=getisgrid.df[,-1]        
        if (genus_name=="all_Bacteria") {
          global_moran_all_Bacteria=moran.test(getisgrid$layer, weighted_neighbors)
        }
        if (genus_name=="all_Fungi") {
          global_moran_all_Fungi=moran.test(getisgrid$layer, weighted_neighbors)
        }
        if (nrow(hotspost_df)==0)
        {
          hotspost_df=getisgrid.df
        } else {
          hotspost_df=merge(hotspost_df,getisgrid.df[,1:2],by="row.names",all=T)
          row.names(hotspost_df)=hotspost_df$Row.names
          hotspost_df=hotspost_df[,-1]
        }
        breaks = c(-20, -1.96, -1, 1, 1.96, 20)
        palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
        col = palette[cut(getisgrid$HOTSPOT, breaks)]
        
        # Plot
        # spplot(getisgrid)
        # par(mfrow=c(2,1))
        # print (p)
        # par(mfrow=c(2,2))
        plot(getisgrid, col=col,main=paste("local G analysis",genus_name," - ",type," ",i,"\n",grid_size,"x",grid_size,"total read=",total))
        legend("bottomleft", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
               c("<-1.96","-1","1","1.96",">1.96"), fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"), horiz=TRUE, cex=0.8)
        rm(getisgrid)
        
      } # end condition on number of reads
      #sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
      #sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    } # end looping all genus with specific grid size
    # write the hotspot data for the spcific grid size
    sheetName = paste("gridSize=",grid_size,sep="")
    addWorksheet(XSLX_Obj, sheetName)
    writeDataTable(x=hotspost_df, wb = XSLX_Obj,sheet = sheetName,rowNames = TRUE)
    
    # summary stats per grid size
    total_Bacteria_hotspots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96)
    total_Fungi_hotspots=sum(hotspost_df$HOTSPOT_all_Fungi>=1.96)
    total_Bacteria_coldspots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96)
    total_Fungi_coldspots=sum(hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_Bacteria_Fungi_HotSpots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
    overlap_Bacteria_Fungi_ColdSpots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_HotSpot_Bacteria_ColdSpot_Fungi=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_ColdSpot_Bacteria_HotSpot_Fungi=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
    global_moran_all_Bacteria_stat=global_moran_all_Bacteria[["estimate"]][["Moran I statistic"]]
    global_moran_all_Bacteria_p=global_moran_all_Bacteria[["p.value"]]
    global_moran_all_Fungi_stat=global_moran_all_Fungi[["estimate"]][["Moran I statistic"]]
    global_moran_all_Fungi_p=global_moran_all_Fungi[["p.value"]]
    record=data.frame(grid_size=grid_size,
                      total_Bacteria_hotspots=total_Bacteria_hotspots,
                      total_Fungi_hotspots=total_Fungi_hotspots,
                      total_Bacteria_coldspots=total_Bacteria_coldspots,
                      total_Fungi_coldspots=total_Fungi_coldspots,
                      overlap_Bacteria_Fungi_HotSpots=overlap_Bacteria_Fungi_HotSpots,
                      overlap_Bacteria_Fungi_ColdSpots=overlap_Bacteria_Fungi_ColdSpots,
                      overlap_HotSpot_Bacteria_ColdSpot_Fungi=overlap_HotSpot_Bacteria_ColdSpot_Fungi,
                      overlap_ColdSpot_Bacteria_HotSpot_Fungi=overlap_ColdSpot_Bacteria_HotSpot_Fungi,
                      total_pixls=nrow(hotspost_df),
                      global_moran_all_Bacteria_stat=global_moran_all_Bacteria_stat,
                      global_moran_all_Bacteria_p=global_moran_all_Bacteria_p,
                      global_moran_all_Fungi_stat=global_moran_all_Fungi_stat,
                      global_moran_all_Fungi_p=global_moran_all_Fungi_p
                      )
    hotspots_overlaps_per_grid_size=rbind(hotspots_overlaps_per_grid_size,record)
    rm (hotspost_df,"getisgrid.df","getisraster",
        "overlap_Bacteria_Fungi_ColdSpots","overlap_Bacteria_Fungi_HotSpots",
        "record","sheetName","neighbors","template","weighted_neighbors")
    
    rm ("box","breaks","col")
    
  } # end loop of different grid sizes
  # write the hotspots_overlaps_per_grid_size data
  write.table(file = out_summary_overlaps,x = hotspots_overlaps_per_grid_size,quote = F,sep = ";",row.names = F)
  saveWorkbook(XSLX_Obj, file = out_xls, overwrite = TRUE)
  
  
  rm ("data" ,"data_file","data_with_xy","genus_name","grid_size","hotspots_overlaps_per_grid_size","i",
      "out_pdfs","out_summary_overlaps","out_xls","palette","pixelsize","spatial_data",
      "XSLX_Obj","xy","xy_list","total")

  dev.off()
}

# for easy comparision plot again now for each grid size
for (smpl in c("A1","A2","B1","B2","C1","C2"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots_compact.pdf",sep="")
  pdf(out_pdfs,width = 14,height = 14)
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  
  i=0
  for (genus_name in names(spatial_data)) # loop all genus of a sample
  {
    i=i+1
    type="Bacteria"
    if (i>50&i<100) {type="Fungi"}
    if (i>100) {type=""}
    total=sum(spatial_data[[genus_name]])
    if(total>5)
    {
      par(mfrow = c(1, 1))
      p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name," - ",type," ",i,"-","total reads",total),cex=1))
      plot(p)
      getisgrid_2=Getis_Ord_GI_per_grid(grid_size = 2,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_3=Getis_Ord_GI_per_grid(grid_size = 3,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_6=Getis_Ord_GI_per_grid(grid_size = 6,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_12=Getis_Ord_GI_per_grid(grid_size = 12,spatial_data = spatial_data,genus_name = genus_name)
      breaks = c(-20, -1.96, -1, 1, 1.96, 20)
      palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(4, 4, 1, 1))
      grid_size=2
      col = palette[cut(getisgrid_2$HOTSPOT, breaks)]
      plot(getisgrid_2, col=col,main=paste(genus_name,"-",type," ",i,"-",grid_size,"x",grid_size))
      grid_size=3
      col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
      plot(getisgrid_3, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      grid_size=6
      col = palette[cut(getisgrid_6$HOTSPOT, breaks)]
      plot(getisgrid_6, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      grid_size=12
      col = palette[cut(getisgrid_12$HOTSPOT, breaks)]
      plot(getisgrid_12, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
             legend =c("<-1.96","-1","1","1.96",">1.96"), 
             fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
             xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    }
  }
  dev.off()
}

# for easy comparision plot again now for each grid size only all Bacteria all Fungi
for (smpl in c("A1","A2","B1","B2","C1","C2"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.AllBacFun.UnderTissue.hotspots.pdf",sep="")
  pdf(out_pdfs,width = 14,height = 14)
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  
  i=0
  for (genus_name in names(spatial_data)[c(101,102)]) # loop all genus of a sample
  {
    i=i+1
    type="Bacteria"
    if (i>50&i<100) {type="Fungi"}
    if (i>100) {type=""}
    total=sum(spatial_data[[genus_name]])
    if(total>5)
    {
      par(mfrow = c(1, 1))
      p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name,"total reads",total),cex=1))
      plot(p)
      getisgrid_2=Getis_Ord_GI_per_grid(grid_size = 2,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_3=Getis_Ord_GI_per_grid(grid_size = 3,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_6=Getis_Ord_GI_per_grid(grid_size = 6,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_12=Getis_Ord_GI_per_grid(grid_size = 12,spatial_data = spatial_data,genus_name = genus_name)
      breaks = c(-20, -1.96, -1, 1, 1.96, 20)
      palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(4, 4, 1, 1))
      grid_size=2
      col = palette[cut(getisgrid_2$HOTSPOT, breaks)]
      plot(getisgrid_2, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_2),"regions"))
      grid_size=3
      col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
      plot(getisgrid_3, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_3),"regions"))
      grid_size=6
      col = palette[cut(getisgrid_6$HOTSPOT, breaks)]
      plot(getisgrid_6, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_6),"regions"))
      grid_size=12
      col = palette[cut(getisgrid_12$HOTSPOT, breaks)]
      plot(getisgrid_12, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_12),"regions"))
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
             legend =c("<-1.96","-1","1","1.96",">1.96"), 
             fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
             xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    }
  }
  dev.off()
}

###### ALL FOR OMNI13
for (smpl in c("A1","A2","B1","B2","C1","C2","D1"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.pdf",sep="")
  pdf(out_pdfs,width = 7,height = 7)
  out_xls=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.xlsx",sep="")
  XSLX_Obj=createWorkbook(out_xls)
  
  out_summary_overlaps=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.Bac_Fun_over.csv",sep="")
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  for (grid_size in c(2:15,24))
  {
    i=0
    hotspost_df=data.frame()
    global_moran_all_Bacteria=NA
    global_moran_all_Fungi=NA
    for (genus_name in names(spatial_data)) # loop all genus and specific grid size
    {
      i=i+1
      type="Bacteria"
      if (i>50&i<100) {type="Fungi"}
      if (i>100) {type=""}
      total=sum(spatial_data[[genus_name]])
      if(total>5)
      {
        if (grid_size==3) # plot only for the first time
        {
          p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name," - ",type," ",i),cex=1))
          plot(p)
          rm(p)
        }
        
        # Create a regular grid of 1/8-mile-square crime areas via a raster
        # library(raster)
        message(paste("\t- genus ",genus_name," grid size=",grid_size,sep=""))
        pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
        box = round(extent(spatial_data) / pixelsize) * pixelsize
        template = raster(box, crs = spatial_data,
                          nrows = (box@ymax - box@ymin) / pixelsize, 
                          ncols = (box@xmax - box@xmin) / pixelsize)
        getisraster = rasterize(spatial_data, template, field = genus_name, fun = sum)
        getisgrid = rasterToPolygons(getisraster)
        # plot(getisgrid)
        # Create the list of neighbors
        neighbors = poly2nb(getisgrid)
        weighted_neighbors = nb2listw(neighbors, zero.policy=T)
        
        # Perform the local G analysis (Getis-Ord GI*)
        getisgrid$HOTSPOT = as.vector(localG(getisgrid$layer, weighted_neighbors))
        # Color the grid cells based on the z-score
        # join with the others
        getisgrid.df=as.data.frame(getisgrid)
        names(getisgrid.df)=paste(names(getisgrid.df),genus_name,sep="_")
        getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
        names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
        getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
        getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
        row.names(getisgrid.df)=getisgrid.df$xy
        getisgrid.df=getisgrid.df[,-1]        
        if (genus_name=="all_Bacteria") {
          global_moran_all_Bacteria=moran.test(getisgrid$layer, weighted_neighbors)
        }
        if (genus_name=="all_Fungi") {
          global_moran_all_Fungi=moran.test(getisgrid$layer, weighted_neighbors)
        }
        if (nrow(hotspost_df)==0)
        {
          hotspost_df=getisgrid.df
        } else {
          hotspost_df=merge(hotspost_df,getisgrid.df[,1:2],by="row.names",all=T)
          row.names(hotspost_df)=hotspost_df$Row.names
          hotspost_df=hotspost_df[,-1]
        }
        breaks = c(-20, -1.96, -1, 1, 1.96, 20)
        palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
        col = palette[cut(getisgrid$HOTSPOT, breaks)]
        
        # Plot
        # spplot(getisgrid)
        # par(mfrow=c(2,1))
        # print (p)
        # par(mfrow=c(2,2))
        plot(getisgrid, col=col,main=paste("local G analysis",genus_name," - ",type," ",i,"\n",grid_size,"x",grid_size,"total read=",total))
        legend("bottomleft", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
               c("<-1.96","-1","1","1.96",">1.96"), fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"), horiz=TRUE, cex=0.8)
        rm(getisgrid)
        
      } # end condition on number of reads
      #sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
      #sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    } # end looping all genus with specific grid size
    # write the hotspot data for the spcific grid size
    sheetName = paste("gridSize=",grid_size,sep="")
    addWorksheet(XSLX_Obj, sheetName)
    writeDataTable(x=hotspost_df, wb = XSLX_Obj,sheet = sheetName,rowNames = TRUE)
    
    # summary stats per grid size
    total_Bacteria_hotspots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96)
    total_Fungi_hotspots=sum(hotspost_df$HOTSPOT_all_Fungi>=1.96)
    total_Bacteria_coldspots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96)
    total_Fungi_coldspots=sum(hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_Bacteria_Fungi_HotSpots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
    overlap_Bacteria_Fungi_ColdSpots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_HotSpot_Bacteria_ColdSpot_Fungi=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_ColdSpot_Bacteria_HotSpot_Fungi=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
    global_moran_all_Bacteria_stat=global_moran_all_Bacteria[["estimate"]][["Moran I statistic"]]
    global_moran_all_Bacteria_p=global_moran_all_Bacteria[["p.value"]]
    global_moran_all_Fungi_stat=global_moran_all_Fungi[["estimate"]][["Moran I statistic"]]
    global_moran_all_Fungi_p=global_moran_all_Fungi[["p.value"]]
    record=data.frame(grid_size=grid_size,
                      total_Bacteria_hotspots=total_Bacteria_hotspots,
                      total_Fungi_hotspots=total_Fungi_hotspots,
                      total_Bacteria_coldspots=total_Bacteria_coldspots,
                      total_Fungi_coldspots=total_Fungi_coldspots,
                      overlap_Bacteria_Fungi_HotSpots=overlap_Bacteria_Fungi_HotSpots,
                      overlap_Bacteria_Fungi_ColdSpots=overlap_Bacteria_Fungi_ColdSpots,
                      overlap_HotSpot_Bacteria_ColdSpot_Fungi=overlap_HotSpot_Bacteria_ColdSpot_Fungi,
                      overlap_ColdSpot_Bacteria_HotSpot_Fungi=overlap_ColdSpot_Bacteria_HotSpot_Fungi,
                      total_pixls=nrow(hotspost_df),
                      global_moran_all_Bacteria_stat=global_moran_all_Bacteria_stat,
                      global_moran_all_Bacteria_p=global_moran_all_Bacteria_p,
                      global_moran_all_Fungi_stat=global_moran_all_Fungi_stat,
                      global_moran_all_Fungi_p=global_moran_all_Fungi_p
    )
    hotspots_overlaps_per_grid_size=rbind(hotspots_overlaps_per_grid_size,record)
    rm (hotspost_df,"getisgrid.df","getisraster",
        "overlap_Bacteria_Fungi_ColdSpots","overlap_Bacteria_Fungi_HotSpots",
        "record","sheetName","neighbors","template","weighted_neighbors")
    
    rm ("box","breaks","col")
    
  } # end loop of different grid sizes
  # write the hotspots_overlaps_per_grid_size data
  write.table(file = out_summary_overlaps,x = hotspots_overlaps_per_grid_size,quote = F,sep = ";",row.names = F)
  saveWorkbook(XSLX_Obj, file = out_xls, overwrite = TRUE)
  
  
  rm ("data" ,"data_file","data_with_xy","genus_name","grid_size","hotspots_overlaps_per_grid_size","i",
      "out_pdfs","out_summary_overlaps","out_xls","palette","pixelsize","spatial_data",
      "XSLX_Obj","xy","xy_list","total")
  
  dev.off()
}
# for easy comparision plot again now for each grid size
for (smpl in c("A1","A2","B1","B2","C1","C2","D1"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots_compact.pdf",sep="")
  pdf(out_pdfs,width = 14,height = 14)
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  
  i=0
  for (genus_name in names(spatial_data)) # loop all genus of a sample
  {
    i=i+1
    type="Bacteria"
    if (i>50&i<100) {type="Fungi"}
    if (i>100) {type=""}
    total=sum(spatial_data[[genus_name]])
    if(total>5)
    {
      par(mfrow = c(1, 1))
      p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name," - ",type," ",i,"-","total reads",total),cex=1))
      plot(p)
      getisgrid_2=Getis_Ord_GI_per_grid(grid_size = 2,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_3=Getis_Ord_GI_per_grid(grid_size = 3,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_6=Getis_Ord_GI_per_grid(grid_size = 6,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_12=Getis_Ord_GI_per_grid(grid_size = 12,spatial_data = spatial_data,genus_name = genus_name)
      breaks = c(-20, -1.96, -1, 1, 1.96, 20)
      palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(4, 4, 1, 1))
      grid_size=2
      col = palette[cut(getisgrid_2$HOTSPOT, breaks)]
      plot(getisgrid_2, col=col,main=paste(genus_name,"-",type," ",i,"-",grid_size,"x",grid_size))
      grid_size=3
      col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
      plot(getisgrid_3, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      grid_size=6
      col = palette[cut(getisgrid_6$HOTSPOT, breaks)]
      plot(getisgrid_6, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      grid_size=12
      col = palette[cut(getisgrid_12$HOTSPOT, breaks)]
      plot(getisgrid_12, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
             legend =c("<-1.96","-1","1","1.96",">1.96"), 
             fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
             xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    }
  }
  dev.off()
}

# for easy comparision plot again now for each grid size only all Bacteria all Fungi
for (smpl in c("A1","A2","B1","B2","C1","C2","D1"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.AllBacFun.UnderTissue.hotspots.pdf",sep="")
  pdf(out_pdfs,width = 14,height = 14)
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  
  i=0
  for (genus_name in names(spatial_data)[c(101,102)]) # loop all genus of a sample
  {
    i=i+1
    type="Bacteria"
    if (i>50&i<100) {type="Fungi"}
    if (i>100) {type=""}
    total=sum(spatial_data[[genus_name]])
    if(total>5)
    {
      par(mfrow = c(1, 1))
      p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name,"total reads",total),cex=1))
      plot(p)
      getisgrid_2=Getis_Ord_GI_per_grid(grid_size = 2,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_3=Getis_Ord_GI_per_grid(grid_size = 3,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_6=Getis_Ord_GI_per_grid(grid_size = 6,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_12=Getis_Ord_GI_per_grid(grid_size = 12,spatial_data = spatial_data,genus_name = genus_name)
      breaks = c(-20, -1.96, -1, 1, 1.96, 20)
      palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(4, 4, 1, 1))
      grid_size=2
      col = palette[cut(getisgrid_2$HOTSPOT, breaks)]
      plot(getisgrid_2, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_2),"regions"))
      grid_size=3
      col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
      plot(getisgrid_3, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_3),"regions"))
      grid_size=6
      col = palette[cut(getisgrid_6$HOTSPOT, breaks)]
      plot(getisgrid_6, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_6),"regions"))
      grid_size=12
      col = palette[cut(getisgrid_12$HOTSPOT, breaks)]
      plot(getisgrid_12, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_12),"regions"))
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
             legend =c("<-1.96","-1","1","1.96",">1.96"), 
             fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
             xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    }
  }
  dev.off()
}

# just all_Bacteria all_Fungi grid size 3 hotspots and observed
for (smpl in c("A1","A2","B1","B2","C1","C2","D1"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI13/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI13_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.AllBacFun.UnderTissue.hotspots3x3.pdf",sep="")
  pdf(out_pdfs,width = 7,height = 7)
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  
  i=0
  for (genus_name in names(spatial_data)[c(101,102)]) # loop all genus of a sample
  {
    i=i+1
    type="Bacteria"
    if (i>50&i<100) {type="Fungi"}
    if (i>100) {type=""}
    total=sum(spatial_data[[genus_name]])
    if(total>5)
    {
      par(mfrow = c(1, 1))
      p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name,"total reads",total),cex=1))
      plot(p)
      getisgrid_3=Getis_Ord_GI_per_grid(grid_size = 3,spatial_data = spatial_data,genus_name = genus_name)
      breaks = c(-20, -1.96, -1, 1, 1.96, 20)
      palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
      grid_size=3
      col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
      plot(getisgrid_3, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_3),"regions"))
      legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
             legend =c("<-1.96","-1","1","1.96",">1.96"), 
             fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
             xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    }
  }
  dev.off()
}



############# START EXPRESSION AND BACTERIAL AND FUNGI #######
#####

hotspots_overlaps_per_grid_size=data.frame()
exp="OMNI12"
exp="OMNI13"

smpl_list=c()
if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}

# one PDF for all PDFs of the all_BACTERIAL all_FUNGI
# out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.AllBacFun.UnderTissue.all_sections.grid2x2.hotspots.pdf",sep="")

for (smpl in smpl_list)
{
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.AllBacFun.UnderTissue.AND_Genes_grt10.grid2x2.hotspots.pdf",sep="")
  pdf(out_pdfs,width = 7,height = 7)
  message (paste("-- ",smpl))
  hotspost_df=data.frame()
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
  # subset genes
  expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=10])]
  # expression_data_with_xy_subset=expression_data_with_xy
  expression_data_with_xy_subset$all_genes=rowSums(expression_data_with_xy_subset[,3:NCOL(expression_data_with_xy_subset)])

  # Bacterial and Fungi
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

  # create the dataframe with the analyzed data -> would be used for creating the spatial object
  house_keeping_genes=c("AT3G18780",
                        "AT1G75780",
                        "AT1G49240",
                        "AT5G62690",
                        "AT2G28390")
  Expression_Fungi_Bacteria_df=expression_data_with_xy_subset[,grepl(paste0(c("x","y",house_keeping_genes,"all_genes"),collapse = "|"),names(expression_data_with_xy_subset))]
  # Expression_Fungi_Bacteria_df=expression_data_with_xy_subset[,grepl(paste0(c("x","y","all_genes"),collapse = "|"),names(expression_data_with_xy_subset))]
  Expression_Fungi_Bacteria_df=merge(Expression_Fungi_Bacteria_df,data_with_xy[,c("x","y","all_Bacteria","all_Fungi")],by=c("x","y"),all=T)
  # remove positions with no genes according to the previous criteria ->> need to figure this out!
  sum(is.na(Expression_Fungi_Bacteria_df$all_genes))
  # Expression_Fungi_Bacteria_df=Expression_Fungi_Bacteria_df[!is.na(Expression_Fungi_Bacteria_df$all_genes),]
  # View(Expression_Fungi_Bacteria_df[is.na(Expression_Fungi_Bacteria_df$all_genes),])

  # create one spatial object only with the (1) expression of house keeping, (2) all expression sum > 10 reads, (3) all Bacterial (4) all fungi
  spatial_data=SpatialPointsDataFrame(coords = Expression_Fungi_Bacteria_df[,c(1:2)],data = Expression_Fungi_Bacteria_df[,3:NCOL(Expression_Fungi_Bacteria_df)])
  
  # par(oma = c(4,1,1,1), mfrow = c(6, 2), mar = c(4, 4, 1, 1))
  for (layer in names(spatial_data)) {
    p=spplot(spatial_data,col='transparent',zcol=layer,main=list(label=layer),cex=1)
    plot(p)
    
    # no grid - need to take special care for NA - hot implemented yet
    # grid_size=1
    # getisgrid=Getis_Ord_GI_per_spot(spatial_data = spatial_data,layer = layer)
    # breaks = c(-20, -1.96, -1, 1, 1.96, 20)
    # palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
    # col = palette[cut(getisgrid$HOTSPOT, breaks)]
    # plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
    
    grid_size=2
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
    brewer.pal(n = 3, name = "RdBu")
    breaks = c(-20, -1.96, -1, 1, 1.96, 20)
    palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
    col = palette[cut(getisgrid$HOTSPOT, breaks)]
    plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
    legend("bottom", inset=.02, title="Z score (local spatial G[i] statistic)",
           legend =c("<-1.96","-1","1","1.96",">1.96"), 
           fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
           xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    
    pmap=tm_shape(getisgrid) + 
      tm_fill("HOTSPOT", 
              palette = "-RdBu",
              style = "pretty", title="G stat",n=5,legend.reverse=T) +
      tm_borders(alpha=.4) + tm_legend(legend.text.size=0.5,legend.title.size=0.6)# +
    print (pmap)
    #tm_layout(legend.position = c("right", "top"), 
    #          title= '% of Population of Black Race', 
    #          title.position = c('center', 'top'))
    
    getisgrid.df=as.data.frame(getisgrid)
    names(getisgrid.df)=paste(names(getisgrid.df),layer,sep="_")
    getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
    names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
    getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
    getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
    row.names(getisgrid.df)=getisgrid.df$xy
    getisgrid.df=getisgrid.df[,-1]        
    if (layer=="all_Bacteria") {
      global_moran_all_Bacteria=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
    }
    if (layer=="all_Fungi") {
      global_moran_all_Fungi=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
    }
    if (layer=="all_genes") {
      global_moran_all_genes=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
    }
    if (nrow(hotspost_df)==0)
    {
      hotspost_df=getisgrid.df
    } else {
      hotspost_df=merge(hotspost_df,getisgrid.df[,1:4],by="row.names",all=T)
      row.names(hotspost_df)=hotspost_df$Row.names
      hotspost_df=hotspost_df[,-1]
    }
    
  } # finish loop over all genes per grid size
  
  
  # summary stats per grid size
  total_Bacteria_hotspots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)>0)
  total_Fungi_hotspots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)>0)
  total_Genes_hotspots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_genes<0.05 & sign (hotspost_df$HOTSPOT_all_genes)>0)
  
  total_Bacteria_coldspots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)<0)
  total_Fungi_coldspots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)<0)
  total_Genes_coldpots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_genes<0.05 & sign (hotspost_df$HOTSPOT_all_genes)<0)
  
  overlap_Bacteria_Fungi_HotSpots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)>0 &
                                      hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)>0)
    
  overlap_Bacteria_Fungi_ColdSpots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)<0 &
                                       hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)<0)
  
  overlap_HotSpot_Bacteria_ColdSpot_Fungi=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)>0 & 
                                              hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)<0)
  
  overlap_ColdSpot_Bacteria_HotSpot_Fungi=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)<0 & 
                                              hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)>0)
  
  overlap_Bacteria_Genes_HotSpots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)>0 &
                                      hotspost_df$HOTSPOT.p.SP_FDR_all_genes<0.05 & sign (hotspost_df$HOTSPOT_all_genes)>0)
                              
  overlap_Bacteria_Gense_ColdSpots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)<0 &
                                       hotspost_df$HOTSPOT.p.SP_FDR_all_genes<0.05 & sign (hotspost_df$HOTSPOT_all_genes)<0)
  
  overlap_Fungi_Genes_HotSpots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)>0 &
                                   hotspost_df$HOTSPOT.p.SP_FDR_all_genes<0.05 & sign (hotspost_df$HOTSPOT_all_genes)>0)
    
  overlap_Fungi_Gense_ColdSpots=sum(hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)<0 &
                                    hotspost_df$HOTSPOT.p.SP_FDR_all_genes<0.05 & sign (hotspost_df$HOTSPOT_all_genes)<0)
  
  global_moran_all_Bacteria_stat=global_moran_all_Bacteria[["estimate"]][["Moran I statistic"]]
  global_moran_all_Bacteria_p=global_moran_all_Bacteria[["p.value"]]
  global_moran_all_Fungi_stat=global_moran_all_Fungi[["estimate"]][["Moran I statistic"]]
  global_moran_all_Fungi_p=global_moran_all_Fungi[["p.value"]]
  
  global_moran_all_Genes_stat=global_moran_all_genes[["estimate"]][["Moran I statistic"]]
  global_moran_all_Genes_p=global_moran_all_genes[["p.value"]]
  
  record=data.frame(sample=smpl,
                    grid_size=grid_size,
                    total_Bacteria_hotspots=total_Bacteria_hotspots,
                    total_Fungi_hotspots=total_Fungi_hotspots,
                    total_Bacteria_coldspots=total_Bacteria_coldspots,
                    total_Fungi_coldspots=total_Fungi_coldspots,
                    overlap_Bacteria_Fungi_HotSpots=overlap_Bacteria_Fungi_HotSpots,
                    overlap_Bacteria_Fungi_ColdSpots=overlap_Bacteria_Fungi_ColdSpots,
                    overlap_HotSpot_Bacteria_ColdSpot_Fungi=overlap_HotSpot_Bacteria_ColdSpot_Fungi,
                    overlap_ColdSpot_Bacteria_HotSpot_Fungi=overlap_ColdSpot_Bacteria_HotSpot_Fungi,
                    overlap_Bacteria_Genes_HotSpots=overlap_Bacteria_Genes_HotSpots,
                    overlap_Bacteria_Gense_ColdSpots=overlap_Bacteria_Gense_ColdSpots,
                    overlap_Fungi_Genes_HotSpots=overlap_Fungi_Genes_HotSpots,
                    overlap_Fungi_Gense_ColdSpots=overlap_Fungi_Gense_ColdSpots,
                    
                    total_pixls=nrow(hotspost_df),
                    global_moran_all_Bacteria_stat=global_moran_all_Bacteria_stat,
                    global_moran_all_Bacteria_p=global_moran_all_Bacteria_p,
                    global_moran_all_Fungi_stat=global_moran_all_Fungi_stat,
                    global_moran_all_Fungi_p=global_moran_all_Fungi_p,
                    global_moran_all_Genes_stat=global_moran_all_Genes_stat,
                    global_moran_all_Genes_p=global_moran_all_Genes_p
  )
  hotspots_overlaps_per_grid_size=rbind(hotspots_overlaps_per_grid_size,record)

  # UpSetR
  library(ComplexHeatmap)
  hotSpots_list=list(Fungi=row.names(hotspost_df)[hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)>0],
                     Bacteria=row.names(hotspost_df)[hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)>0],
                     Genes=row.names(hotspost_df)[hotspost_df$HOTSPOT.p.SP_FDR_all_genes<0.05 & sign (hotspost_df$HOTSPOT_all_genes)>0])
  hotSpots_PA=as.data.frame(list_to_matrix(hotSpots_list))                

  #pdf(file = paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/OMNI12/data/multimaped_filtered/MMSEQS2/spatial_info/OMNI12.All_Bacterial_ITS_Probs.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.Top50.","corr.upsetR.pdf",sep=""),
  #    height=5, width=20
  #    
  # )
  
  # upset(fromList(nodes_list),empty.intersections = "on")#, order.by = "freq"
  # movies <- read.csv( system.file("extdata", "movies.csv", package = "UpSetR"), header=T, sep=";" )
  # upset(as.data.frame(test),empty.intersections = "on")
  # with complexHeatMap: https://jokergoo.github.io/ComplexHeatmap-reference/book/upset-plot.html
  m = make_comb_mat(hotSpots_list) # https://support.bioconductor.org/p/118557/
  cs = comb_size(m)
  row_size = set_size(m)
  
  ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))))
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
  library(eulerr)
  
  intersects_hotSpots=c("Fungi" = 0, "Bacteria" = 0, "Genes" = 0,
               "Fungi&Bacteria&Genes"=0,  
               "Fungi&Bacteria"=0,
               "Bacteria&Genes"=0,
               "Fungi&Genes"=0)
  for (i in names(cs)) { # assign
    if (i=="100") {intersects_hotSpots[["Fungi"]]=cs[[i]]}
    if (i=="010") {intersects_hotSpots[["Bacteria"]]=cs[[i]]}
    if (i=="001") {intersects_hotSpots[["Genes"]]=cs[[i]]}
    if (i=="111") {intersects_hotSpots[["Fungi&Bacteria&Genes"]]=cs[[i]]}
    if (i=="110") {intersects_hotSpots[["Fungi&Bacteria"]]=cs[[i]]}
    if (i=="101") {intersects_hotSpots[["Fungi&Genes"]]=cs[[i]]}
    if (i=="011") {intersects_hotSpots[["Bacteria&Genes"]]=cs[[i]]}
  }
  
  fit2 <- euler(intersects_hotSpots)
  plot(fit2,quantities = TRUE,main=paste("HotSpots ",smpl,sep=""))
  
  # For cold spots
  coldSpots_list=list(Fungi=row.names(hotspost_df)[hotspost_df$HOTSPOT.p.SP_FDR_all_Fungi<0.05 & sign (hotspost_df$HOTSPOT_all_Fungi)<0],
                      Bacteria=row.names(hotspost_df)[hotspost_df$HOTSPOT.p.SP_FDR_all_Bacteria<0.05 & sign (hotspost_df$HOTSPOT_all_Bacteria)<0],
                      Genes=row.names(hotspost_df)[hotspost_df$HOTSPOT.p.SP_FDR_all_genes<0.05 & sign (hotspost_df$HOTSPOT_all_genes)<0])
  m = make_comb_mat(coldSpots_list) # https://support.bioconductor.org/p/118557/
  cs = comb_size(m)
  library(eulerr)
  
  intersects_coldSpots=c("Fungi" = 0, "Bacteria" = 0, "Genes" = 0,
                        "Fungi&Bacteria&Genes"=0,  
                        "Fungi&Bacteria"=0,
                        "Bacteria&Genes"=0,
                        "Fungi&Genes"=0)
  for (i in names(cs)) { # assign
    if (i=="100") {intersects_coldSpots[["Fungi"]]=cs[[i]]}
    if (i=="010") {intersects_coldSpots[["Bacteria"]]=cs[[i]]}
    if (i=="001") {intersects_coldSpots[["Genes"]]=cs[[i]]}
    if (i=="111") {intersects_coldSpots[["Fungi&Bacteria&Genes"]]=cs[[i]]}
    if (i=="110") {intersects_coldSpots[["Fungi&Bacteria"]]=cs[[i]]}
    if (i=="101") {intersects_coldSpots[["Fungi&Genes"]]=cs[[i]]}
    if (i=="011") {intersects_coldSpots[["Bacteria&Genes"]]=cs[[i]]}
  }
  fit2 <- euler(intersects_coldSpots)
  plot(fit2,quantities = TRUE,main=paste("ColdSpots ",smpl,sep=""))
  dev.off()
}
# write the summing table
write.table(file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,".Top50_ITS_16S_Probs_and_UNKNOWN.AllBacFun.UnderTissue.AND_Genes_grt10.grid2x2.hotspots.csv",sep=""),
            x=hotspots_overlaps_per_grid_size,quote = F,sep = ";",row.names = F,col.names = T)


############# DO SOME PLOTING OF FEW MARKER GENES AND TOTAL MICROBIAL ############# 
############# CLEAN START
# ploting functions
library(RColorBrewer)
library(tmap)

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

# read the data of all samples into list of sp objects
# new test based on: https://michaelminn.net/tutorials/r-point-analysis/index.html
library(sp)
# read all data files to a list of spatial object indexed by their sample names
all_samples_spatial_data=list()
message ("## Read all data into spatial objects")
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (smpl in samples)
  {
    sample_name=paste(Exp,smpl,sep="_")
    message(paste("--",sample_name))
    data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
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
    expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",Exp,".",smpl,".tsv",sep="")
    expression_data=read.delim(expression_data_file,sep="\t",stringsAsFactors = F)
    expression_data=t(expression_data) # row = position; col = gene
    expression_data=as.data.frame(expression_data)
    expression_data$all_Expression=rowSums(expression_data)
    
    
    # expression_sum=data.frame(all_Expression=rowSums(expression_data))
    xy_list=strsplit(row.names(expression_data), "x")
    xy=as.data.frame(t(as.data.frame(xy_list)))
    xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
    names(xy)=c("x","y")
    expression_data_with_xy=cbind(xy,expression_data)
    expression_data_with_xy$x=as.numeric(expression_data_with_xy$x)
    expression_data_with_xy$y=as.numeric(expression_data_with_xy$y)
    
    joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
    
    # save the option to have all expression for the later use...    
    # # Only genes with more than overall 10 raeds
    # gene_sum=data.frame(count=colSums(expression_data))
    # gene_sum$percent_of_data=gene_sum$count/(sum(gene_sum$count))
    ## sum(gene_sum$percent[gene_sum$count>=10])
    ## row.names(gene_sum)
    # gene_sum_ordered=gene_sum[order(gene_sum$count,decreasing = T),]
    
    ## create the spatial object
    # xy_list=strsplit(row.names(expression_data), "x")
    # xy=as.data.frame(t(as.data.frame(xy_list)))
    # xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
    # names(xy)=c("x","y")
    # expression_data_with_xy=cbind(xy,expression_data)
    # expression_data_with_xy$x=as.numeric(expression_data_with_xy$x)
    # expression_data_with_xy$y=as.numeric(expression_data_with_xy$y)
    # 
    # rownames(expression_data_with_xy)=rownames(expression_data)
    ## subset genes
    # RNA_cutoff=100 # Only genes with more than 100 reads are accounted
    # expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
    
    ## Join the expression and all_Fungi all_Bacteria
    # joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
    
    
    # create the spatial object
    spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
    
    #store it in the a list
    all_samples_spatial_data[[sample_name]]=spatial_data
    
    # clean up
    rm (sample_name,data_file,data,xy_list,xy,data_with_xy,
        expression_data_file,expression_data,
        joint_Expression_Bacteria_Fungi_data,spatial_data) 
    # gene_sum,gene_sum_ordered,expression_data_with_xy,RNA_cutoff,expression_data_with_xy_subset
  }
}
rm (Exp,smpl) 

sample_name="OMNI12_A1"
spatial_data=all_samples_spatial_data[[sample_name]]
house_keeping_genes=c("AT3G18780",
                      "AT1G75780",
                      "AT1G49240",
                      "AT5G62690",
                      "AT2G28390")
low_Moran_examples=c("AT3G46530","AT1G80600",
                     "AT5G42740","AT3G23820",
                     "AT2G18790","AT1G64550")
sum_layers=c("all_Expression","all_Bacteria","all_Fungi")
layers_to_plot_list=c(house_keeping_genes,low_Moran_examples,sum_layers)

pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/OMNI12_A1_OMNI13_B1_marker_genes_LowMoran_and_housekeepenig.pdf"))
pdf(pdf_file,width = 42, height = 14)

samples_to_plot=c("OMNI12_A1","OMNI13_B1")
for (sample_to_plot in samples_to_plot) #names(all_samples_spatial_data))
{
  plots_list=list()
  for (type in c("localG")) #
  { 
    for (layer_to_plot in layers_to_plot_list)
    {
      plot_name=paste(sample_to_plot,layer_to_plot,type,sep="_")
      if (type=="localG")
      {
        plots_list[[plot_name]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,sep=""))
      }
      if (type=="reads") 
      {
        plots_list[[plot_name]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(plot_name,sep=""))
      }
    }
  }
  # plots_list 
  # https://github.com/r-tmap/tmap/issues/511
  
  # current.mode <- tmap_mode("plot")
  ## tmap_arrange(plots_list, widths = c(.75, .75))
  
  # plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
  # plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
  # plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
  # plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
  
  # for (i in 3:4) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
  print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 7, nrow = 2))
}
dev.off()


#### END PLOTING MARKER GENES
#### For each section plot the Bacteria, Fungi, Expression localG with tmap
#### START VERY CLEAN --> clean the entire workspace
#### LAST UPDATE AND RUN: 20Apr2022
#### THIS IS THE FINAL VERSION USED TO CREATE DATA FOR FIGURE 3 and 
#### 
# ploting functions
library(RColorBrewer)
library(tmap)

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

get_hotspots_maps=function (spatial_data,layer,grid_size=2,title=NA,only_significant="NO")
{
  if (is.na(title)) {title=layer}
  # pdf
  # pdf (file = out_pdf,height = 7,width = 7)
  
  # calculate the localG stas
  if (length(spatial_data[[layer]])>0)
  {
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = layer)
    if (toupper(only_significant)=="YES") { # make the non significant hotspots value as NA
      getisgrid$HOTSPOT[getisgrid$HOTSPOT.p.SP_FDR>0.05]=NA
    }
    pmap=tm_shape(getisgrid) + 
      tm_fill("HOTSPOT", 
              palette = "-RdBu",
              style = "pretty", title="G stat",n=5,legend.reverse=T, midpoint=0,textNA = "NS") +
      tm_borders(alpha=.4) + tm_legend(legend.text.size=0.5,legend.title.size=0.6) +
      tm_layout(title = title,title.size=0.8)
  } else {
    pmap=NULL
  }
  return (pmap)
}
get_hotspots_maps_fixed=function (spatial_data,layer,grid_size=2,title=NA,only_significant="NO")
{
  if (is.na(title)) {title=layer}
  # pdf
  # pdf (file = out_pdf,height = 7,width = 7)
  
  # calculate the localG stas
  if (length(spatial_data[[layer]])>0)
  {
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = layer)
    if (toupper(only_significant)=="YES") { # make the non significant hotspots value as NA
      getisgrid$HOTSPOT[getisgrid$HOTSPOT.p.SP_FDR>0.05]=NA
    }
    # https://geocompr.github.io/post/2019/tmap-color-scales/
    message(paste("== ",title,"min:",min(getisgrid$HOTSPOT,na.rm = T),"max:",max(getisgrid$HOTSPOT,na.rm=T)))
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

# read the data of all samples into list of sp objects
# new test based on: https://michaelminn.net/tutorials/r-point-analysis/index.html
library(sp)
# read all data files to a list of spatial object indexed by their sample names
all_samples_spatial_data=list()
all_samples_GeneExpression_spatial_data=list() # Spatial object for each of the genes
all_samples_GeneExpression_and_microbial_abundance_spatial_data=list() # Spatial object for each of the genes and total Top50 Bacteria/Fungi genera
all_samples_total_gene_expression=data.frame() # for each gene in each section hold the total expresion across the section pixels
message ("## Read all data into spatial objects")
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (smpl in samples)
  {
    sample_name=paste(Exp,smpl,sep="_")
    message(paste("--",sample_name))
    ### Create new under the tissue matrix based on the new expression data
    all_Bacterial_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Bacterial_and_UNKNOWN/",Exp,"_",smpl,".All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos.csv",sep="")
    bacteria_data=read.delim(all_Bacterial_data_file,sep=";",stringsAsFactors = F)
    row.names(bacteria_data)=bacteria_data$tax
    bacteria_data=bacteria_data[,-1]
    t_bacteria_data=t(bacteria_data)
    bacteria_top_50_names=colnames(t_bacteria_data)[1:50]
    
    message(paste("[INFO] Total number of pixels in Bacteria data '",all_Bacterial_data_file,"': ",NROW(t_bacteria_data),sep=""))
    
    all_Fungi_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Fungi_and_UNKNOWN/",Exp,"_",smpl,".All_ITS_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos.csv",sep="")
    fungi_data=read.delim(all_Fungi_data_file,sep=";",stringsAsFactors = F)
    row.names(fungi_data)=fungi_data$tax
    fungi_data=fungi_data[,-1]
    t_fungi_data=t(fungi_data)
    fungi_top_50_names=colnames(t_fungi_data)[1:50]
    message(paste("[INFO] Total number of pixels in Fungi data '",all_Fungi_data_file,"': ",NROW(t_fungi_data),sep=""))
    # join bacteria_fungi
    bacteria_fungi_joined=merge(t_bacteria_data,t_fungi_data,by="row.names",all=T)
    row.names(bacteria_fungi_joined)=bacteria_fungi_joined$Row.names
    bacteria_fungi_joined=bacteria_fungi_joined[,-1]
    
    message(paste("[INFO] Total number of pixels in joined data: ",NROW(bacteria_fungi_joined),sep=""))
    # info the number of rows (pixels) with NA
    message(paste("[INFO] Total number of pixels with missing data in joined data (will be transfered to 0): ",sum(is.na(rowSums(bacteria_fungi_joined)))," with total of: ",sum(is.na(bacteria_fungi_joined))," cells",sep=""))
    bacteria_fungi_joined[is.na(bacteria_fungi_joined)]=0
    message(paste("[INFO] Total number of pixels with total microbial abundance 0: ",sum(rowSums(bacteria_fungi_joined)==0),sep=""))
    if (sum(rowSums(bacteria_fungi_joined)==0)>0)
    {
      bacteria_fungi_joined=bacteria_fungi_joined[rowSums(bacteria_fungi_joined)>0,]
    }
    
    # Take the sum of top 50 Bacteria and top 50 Fungi
    xy_list=strsplit(row.names(bacteria_fungi_joined), "x")
    xy=as.data.frame(t(as.data.frame(xy_list)))
    xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
    names(xy)=c("x","y")
    data_with_xy=data.frame(row.names=row.names(bacteria_fungi_joined),xy,all_Bacteria=rowSums(bacteria_fungi_joined[,bacteria_top_50_names]),all_Fungi=rowSums(bacteria_fungi_joined[,fungi_top_50_names]))
    
    # old version
    # data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
    # data=read.delim(data_file,sep=";",stringsAsFactors = F)
    # data$all_Bacteria=rowSums(data[,1:50])
    # data$all_Fungi=rowSums(data[,51:100])
    
    # xy_list=strsplit(row.names(data), "x")
    # xy=as.data.frame(t(as.data.frame(xy_list)))
    # xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
    # names(xy)=c("x","y")
    # data_with_xy=cbind(xy,data)
    # rownames(data_with_xy)=rownames(data)
    # data_with_xy$x=as.numeric(data_with_xy$x)
    # data_with_xy$y=as.numeric(data_with_xy$y)
    
    # Now with the expression
    # expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",Exp,".",smpl,".tsv",sep="")
    
    # Updated expression data - Apr2022
    # local version - Apr2022
    expression_data_file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/RNA_data/NewestData/220405_semiwild_dataset_rawcounts_filtered.",Exp,".",smpl,".tsv",sep="") # the latest
    # cluster version - Apr2022
    # expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/Apr2022/220405_semiwild_dataset_rawcounts_filtered.",Exp,".",smpl,".tsv",sep="") # the latest
    message(paste("[INFO] Expression data: '",expression_data_file,"'",sep=""))
    expression_data=read.delim(expression_data_file,sep="\t",stringsAsFactors = F)
    expression_data=t(expression_data) # row = position; col = gene
    
    expression_sum=data.frame(all_Expression=rowSums(expression_data))
    xy_list=strsplit(row.names(expression_data), "x")
    xy=as.data.frame(t(as.data.frame(xy_list)))
    xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
    names(xy)=c("x","y")
    expression_sum_with_xy=cbind(xy,expression_sum)
    expression_sum_with_xy$x=as.numeric(expression_sum_with_xy$x)
    expression_sum_with_xy$y=as.numeric(expression_sum_with_xy$y)
    
    expression_data_with_xy=data.frame(xy,expression_data)
    
    message (paste ("[INFO] Total position with expression values:",NROW(expression_sum_with_xy)))
    
    # validate the expression sum is greater than 0
    message (paste ("[INFO] Total position with 0 expression (would be filtered):",sum(expression_sum_with_xy$all_Expression==0)))
    if (sum(expression_sum_with_xy$all_Expression==0)>0)
    {
      expression_sum_with_xy=expression_sum_with_xy[expression_sum_with_xy$all_Expression>0,]
    }
    joint_Expression_Bacteria_Fungi_data=merge(expression_sum_with_xy,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"),all.x=T)
    joint_Expression_Bacteria_Fungi_data[is.na(joint_Expression_Bacteria_Fungi_data)]=0
    message (paste ("[INFO] Final total position under the tissue:",NROW(joint_Expression_Bacteria_Fungi_data)))
    
    # Create an object with sum of Top50 Bacteria, Top50 Fungi and the full expression pattern for all genes
    joint_full_Expreesion_pattern_and_total_Microbial_abundance=merge(expression_data_with_xy,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"),all.x=T)
    joint_full_Expreesion_pattern_and_total_Microbial_abundance$x=as.numeric(joint_full_Expreesion_pattern_and_total_Microbial_abundance$x)
    joint_full_Expreesion_pattern_and_total_Microbial_abundance$y=as.numeric(joint_full_Expreesion_pattern_and_total_Microbial_abundance$y)
    joint_full_Expreesion_pattern_and_total_Microbial_abundance[is.na(joint_full_Expreesion_pattern_and_total_Microbial_abundance)]=0
    message (paste ("[INFO] Final total position under the tissue for alll genes:",NROW(joint_full_Expreesion_pattern_and_total_Microbial_abundance)))
    
    # save the spatial object for all the genes expression (gene by gene) and total microbial abundance 
    spatial_gene_expression_and_microbial_abundance_data=SpatialPointsDataFrame(coords = joint_full_Expreesion_pattern_and_total_Microbial_abundance[,c(1:2)],data = joint_full_Expreesion_pattern_and_total_Microbial_abundance[,3:NCOL(joint_full_Expreesion_pattern_and_total_Microbial_abundance)])
    
    #store it in the a list
    all_samples_GeneExpression_and_microbial_abundance_spatial_data[[sample_name]]=spatial_gene_expression_and_microbial_abundance_data
    
    
    
    # double make sure no position with row sums 0
    joint_Expression_Bacteria_Fungi_data.qa=rowSums(joint_Expression_Bacteria_Fungi_data[,3:5],na.rm = T)
    message (paste ("[INFO] QA2, Positions left with total microbial+expression zero:",sum(joint_Expression_Bacteria_Fungi_data.qa==0)))
    rm(joint_Expression_Bacteria_Fungi_data.qa)
    
    # save the final list of pixels for later use and substing
    xy_under_the_tissue=data.frame(x=joint_Expression_Bacteria_Fungi_data$x,y=joint_Expression_Bacteria_Fungi_data$y,xy=paste(joint_Expression_Bacteria_Fungi_data$x,joint_Expression_Bacteria_Fungi_data$y,sep="x"))
    under_the_tissue_pixles_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/",Exp,"_",smpl,".Bacteria_Fungi_expression_220405.under_the_tissue_pixels.txt",sep="")
    write.table(file=under_the_tissue_pixles_file,x=xy_under_the_tissue,sep="\t",quote = F,row.names = F)
    message(paste ("[INFO] all positions under the tissue were written to: ",under_the_tissue_pixles_file,sep=""))
    
    # Write also the sum of Bacterial/Fungi data to be used later for the correlations etc.
    bacteria_under_the_tissue_Top50_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Bacterial_and_UNKNOWN/Under_tissue/Apr22_expression_220405/",Exp,"_",smpl,".All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top_50_genus.under_tissue.based_on_joined_Bacteria_Fungi_Expression_220405.spatial_pos.csv",sep="")
    
    bacteria_sum_for_writing=as.data.frame(t(data.frame(row.names=paste(joint_Expression_Bacteria_Fungi_data$x,joint_Expression_Bacteria_Fungi_data$y,sep="x"),
                                                        Top50_Bacteria_genus=joint_Expression_Bacteria_Fungi_data$all_Bacteria)))
    bacteria_sum_for_writing=data.frame(tax="Top50_Bacteria_genus",bacteria_sum_for_writing)
    # str(bacteria_sum_for_writing)
    write.table(file=bacteria_under_the_tissue_Top50_file,x=bacteria_sum_for_writing,quote = F,row.names = F,sep = ";")
    message(paste("[INFO] The sum of top 50 Bacterial genus under the tissue was written to:",bacteria_under_the_tissue_Top50_file))
    
    fungi_under_the_tissue_Top50_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Fungi_and_UNKNOWN/Under_tissue/Apr22_expression_220405/",Exp,"_",smpl,".All_ITS_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top_50_genus.under_tissue.based_on_joined_Bacteria_Fungi_Expression_220405.spatial_pos.csv",sep="")
    fungi_sum_for_writing=as.data.frame(t(data.frame(row.names=paste(joint_Expression_Bacteria_Fungi_data$x,joint_Expression_Bacteria_Fungi_data$y,sep="x"),
                                                        Top50_ITS_genus=joint_Expression_Bacteria_Fungi_data$all_Fungi)))
    fungi_sum_for_writing=data.frame(tax="Top50_ITS_genus",fungi_sum_for_writing)
    # str(fungi_sum_for_writing)
    write.table(file=fungi_under_the_tissue_Top50_file,x=fungi_sum_for_writing,quote = F,row.names = F,sep = ";")
    message(paste("[INFO] The sum of top 50 ITS genus under the tissue was written to:",fungi_under_the_tissue_Top50_file))
    
    # make sure the expression data (gene-by-gene) positions is the same as the joined object
    # subste the expression data to be with the same positions as the joined Microbial expression position (should be the same)
    message(paste("[INFO] Total position on the original expression file was:",NROW(expression_data_with_xy),"- if needed subset it to be the same final pixels having Expresion and Microbial data, see above (that is total of:",NROW(xy_under_the_tissue),"pixels)",sep=" "))
    expression_data_with_xy=expression_data_with_xy[paste(expression_data_with_xy$x,expression_data_with_xy$y,sep="x") %in% xy_under_the_tissue$xy,]
    message(paste("[INFO] QA - Total pixels after subseting:",NROW(expression_data_with_xy),"- of which with NA values:",sum(is.na(xy_under_the_tissue)),sep=" "))
    # save the spatial object for all the genes expression (gene by gene)
    expression_data_with_xy$x=as.numeric(expression_data_with_xy$x)
    expression_data_with_xy$y=as.numeric(expression_data_with_xy$y)
    
    spatial_gene_expression_data=SpatialPointsDataFrame(coords = expression_data_with_xy[,c(1:2)],data = expression_data_with_xy[,3:NCOL(expression_data_with_xy)])
    
    # save it for furter usage
    Expression_under_the_tissue_all_genes_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/Apr2022/Under_tissue/","220405_semiwild_dataset_rawcounts_filtered.under_tissue.based_on_joined_Bacteria_Fungi_Expression_220405.spatial_pos.",Exp,".",smpl,".tsv",sep="")
    Expression_under_the_tissue_all_genes_for_writing=data.frame(row.names=paste("X",paste(expression_data_with_xy$x,expression_data_with_xy$y,sep="x"),sep=""),
                                                                 expression_data_with_xy[,3:NCOL(expression_data_with_xy)])
    Expression_under_the_tissue_all_genes_for_writing=t(Expression_under_the_tissue_all_genes_for_writing)  
    write.table(file=Expression_under_the_tissue_all_genes_file,x=Expression_under_the_tissue_all_genes_for_writing,quote = F,row.names = T,sep = "\t")
    message(paste("[INFO] The Expression under the tissue (all genes) was written to:",Expression_under_the_tissue_all_genes_file))
    
    #store it in the a list
    all_samples_GeneExpression_spatial_data[[sample_name]]=spatial_gene_expression_data
    
    # save the total gene expression across the section
    gene_expression=data.frame(gene=colnames(expression_data))
    gene_expression[[sample_name]]=colSums(expression_data)
    if (NROW(all_samples_total_gene_expression)==0)
    {
      all_samples_total_gene_expression=gene_expression
    } else {
      all_samples_total_gene_expression=merge(all_samples_total_gene_expression,gene_expression,by="gene",all=T)
    }
    
    # save the option to have all expression for the later use...    
    # # Only genes with more than overall 10 raeds
    # gene_sum=data.frame(count=colSums(expression_data))
    # gene_sum$percent_of_data=gene_sum$count/(sum(gene_sum$count))
    ## sum(gene_sum$percent[gene_sum$count>=10])
    ## row.names(gene_sum)
    # gene_sum_ordered=gene_sum[order(gene_sum$count,decreasing = T),]
    
    ## create the spatial object
    # xy_list=strsplit(row.names(expression_data), "x")
    # xy=as.data.frame(t(as.data.frame(xy_list)))
    # xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
    # names(xy)=c("x","y")
    # expression_data_with_xy=cbind(xy,expression_data)
    # expression_data_with_xy$x=as.numeric(expression_data_with_xy$x)
    # expression_data_with_xy$y=as.numeric(expression_data_with_xy$y)
    
    # rownames(expression_data_with_xy)=rownames(expression_data)
    ## subset genes
    # RNA_cutoff=100 # Only genes with more than 100 reads are accounted
    # expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
    
    ## Join the expression and all_Fungi all_Bacteria
    # joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
    

    # create the spatial object for sum of genes and sum of Top50 microbial Genus
    spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
    
    #store it in the a list
    all_samples_spatial_data[[sample_name]]=spatial_data
    
    # clean up
    rm (sample_name,all_Bacterial_data_file,bacteria_data,bacteria_top_50_names,t_bacteria_data,bacteria_sum_for_writing,bacteria_under_the_tissue_Top50_file,
        all_Fungi_data_file,fungi_data,fungi_top_50_names,t_fungi_data,fungi_sum_for_writing,fungi_under_the_tissue_Top50_file,                   
        bacteria_fungi_joined,xy_list,xy,data_with_xy,
        expression_data_file,expression_data,
        joint_Expression_Bacteria_Fungi_data,spatial_data,expression_sum,expression_sum_with_xy,
        xy_under_the_tissue,under_the_tissue_pixles_file,
        spatial_gene_expression_data,gene_expression,expression_data_with_xy,
        Expression_under_the_tissue_all_genes_file, Expression_under_the_tissue_all_genes_for_writing) 
    # data,data_file
    # gene_sum,gene_sum_ordered,expression_data_with_xy,RNA_cutoff,expression_data_with_xy_subset
  }
}
rm (Exp,smpl) 

total_gene_expression_file="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/220405_semiwild_dataset_rawcounts_filtered.Apr2022.OMNI12_OMNI13.Total_gene_expressions_alll_sections.tsv"
message(paste("[INFO] Writing the total gene expression profiles in all OMNI12 OMNI13 sections to:",total_gene_expression_file))
write.table(x=all_samples_total_gene_expression,file=total_gene_expression_file,row.names = F,quote = F,sep = "\t")

# Save the list with all the spatial objects
message(paste("[INFO] Save the spatial objects list for the Expression and Microbial sums to file: ~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Top50_Bacterial_Genus_Top50_ITS_Genus_All_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.SUMS.spatial_objects_list.rds"))
saveRDS(all_samples_spatial_data, file = "~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Top50_Bacterial_Genus_Top50_ITS_Genus_All_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.SUMS.spatial_objects_list.rds")
message(paste("[INFO] Save the spatial objects list for all the gene expression profiles to file: ~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.spatial_objects_list.rds"))
saveRDS(all_samples_GeneExpression_spatial_data, file = "~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.spatial_objects_list.rds")
message(paste("[INFO] Save the spatial objects list for all the gene expression profiles and the total top50 Bacteria/Fungi genera to file: ~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.and_Top50_Bacteria_Fungi_Sum.spatial_objects_list.rds"))
saveRDS(all_samples_GeneExpression_and_microbial_abundance_spatial_data, file = "~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.and_Top50_Bacteria_Fungi_Sum.spatial_objects_list.rds")


# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/all_selected_gene.pdf"))

# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/OMNI12_OMNI13_All_sections_total_Expression_Bacteria_Fungi_Hotspots.pdf"))
# Apr2022 expression data
for (Fixed_sacle in c(T,F)) {
  pdf_file=NA
  if (Fixed_sacle) {
    pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/OMNI12_OMNI13_All_sections_total_Expression_220405_Bacteria_Fungi_Hotspots.Fixed_Sacle.pdf"))
  } else {
    pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/OMNI12_OMNI13_All_sections_total_Expression_220405_Bacteria_Fungi_Hotspots.pdf"))
  }
  pdf(pdf_file,width = 21,height = 14)
  
  for (sample_to_plot in names(all_samples_spatial_data))
  {
    plots_list=list()
    for (type in c("localG","reads"))
    {
      for (layer_to_plot in c("all_Bacteria","all_Fungi","all_Expression"))
      {
        plot_name=paste(sample_to_plot,layer_to_plot,type,sep="_")
        if (type=="localG")
        {
          if (Fixed_sacle) {
            plots_list[[plot_name]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,sep=""))
          }else {
            plots_list[[plot_name]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,sep=""))
          }
        }
        if (type=="reads") 
        {
          plots_list[[plot_name]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(plot_name,sep=""))
        }
      }
    }
    # plots_list 
    # https://github.com/r-tmap/tmap/issues/511
    
    # current.mode <- tmap_mode("plot")
    ## tmap_arrange(plots_list, widths = c(.75, .75))
    
    # plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
    # plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
    # plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
    # plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
    
    # for (i in 3:4) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
    print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 3, nrow = 2))
  }
  dev.off()
}


### PLOT ONLY THE SIGNIFICANT HOTSPOTS
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/OMNI12_OMNI13_All_sections_total_Expression_Bacteria_Fungi_Hotspots.ONLY_SIGNIFICANT.pdf"))

# Apr2022 version
for (Fixed_sacle in c(T,F)) {
  pdf_file=NA
  if (Fixed_sacle) {
    pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/OMNI12_OMNI13_All_sections_total_Expression_220405_Bacteria_Fungi_Hotspots.ONLY_SIGNIFICANT.Fixed_Sacle.pdf"))
  } else {
    pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/OMNI12_OMNI13_All_sections_total_Expression_220405_Bacteria_Fungi_Hotspots.ONLY_SIGNIFICANT.pdf"))
  }

  pdf(pdf_file,width = 21,height = 14)
  
  for (sample_to_plot in names(all_samples_spatial_data))
  {
    plots_list=list()
    for (type in c("localG","reads"))
    {
      for (layer_to_plot in c("all_Bacteria","all_Fungi","all_Expression"))
      {
        plot_name=paste(sample_to_plot,layer_to_plot,type,sep="_")
        if (type=="localG")
        {
          if (Fixed_sacle) {
            plots_list[[plot_name]]=get_hotspots_maps_fixed(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,sep=""),only_significant = "YES")
          } else {
            plots_list[[plot_name]]=get_hotspots_maps(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,grid_size = 2,title = paste(plot_name,sep=""),only_significant = "YES")
          }
        }
        if (type=="reads") 
        {
          plots_list[[plot_name]]=get_dot_map(spatial_data = all_samples_spatial_data[[sample_to_plot]],layer = layer_to_plot,title = paste(plot_name,sep=""))
        }
      }
    }
    # plots_list 
    # https://github.com/r-tmap/tmap/issues/511
    
    # current.mode <- tmap_mode("plot")
    ## tmap_arrange(plots_list, widths = c(.75, .75))
    
    # plots_list[[1]]=plots_list[[1]]+tm_layout(main.title = paste (gene_id))
    # plots_list[[2]]=plots_list[[2]]+tm_layout(main.title = gene_name)
    # plots_list[[3]]=plots_list[[3]]+tm_layout(main.title = paste ("nBac",nBac,"nFun",nFungi))
    # plots_list[[4]]=plots_list[[4]]+tm_layout(main.title = response_type)
    
    # for (i in 3:4) {plots_list[[i]]=plots_list[[i]]+tm_layout(main.title = " ")}
    print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 3, nrow = 2))
  }
  dev.off()
}

# Calculate overlaps Hotspots/Bacteria/Fungi
# Get hotspots.df for each section Bacteri/Fungi/Expression and intersect

# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/OMNI12_OMNI13_All_sections_intersect_Expression_Bacteria_Fungi_Hotspots.pdf"))
# Apr2022
pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/OMNI12_OMNI13_All_sections_intersect_Expression_220405_Bacteria_Fungi_Hotspots.pdf"))
pdf(pdf_file)

hotsposts_intersect_Bacteria_Fungi_Expression_stats=data.frame()
grid_size=2
library(eulerr)
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (smpl in samples)
  {
    hotspots_df=data.frame()
    section_localG_stats=data.frame() # save also the raw results data for each of the layers 
    sample_name=paste(Exp,smpl,sep="_")
    message(paste("--",sample_name))
    spatial_data=all_samples_spatial_data[[sample_name]]
    for (marker in c("all_Bacteria","all_Fungi","all_Expression"))
    {
      getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = marker)
      if (sum(is.na(getisgrid$HOTSPOT))>0) # NA in the Hotspots data -- report and remove the pixel from download analyses
      {
        message(paste("[ERROR] ",sum(is.na(getisgrid$HOTSPOT))," grid pixels were NA for layer: ",marker," (pixels: ",
                      paste(paste(coordinates(getisgrid[is.na(getisgrid$HOTSPOT),])[,1],coordinates(getisgrid[is.na(getisgrid$HOTSPOT),])[,2],sep="x"),collapse = ' '),
                      ") and would be filtered!",sep="")
                )
        getisgrid=getisgrid[!is.na(getisgrid$HOTSPOT),]
      }
      ## save the data matrix for later usage (e.g the BORUTA)
      localG_matrix_file="NA"
      if (marker=="all_Bacteria")
      {
        localG_matrix_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Bacterial_and_UNKNOWN/Under_tissue/Apr22_expression_220405/hotspots/",Exp,"_",smpl,".All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top_50_genus.",grid_size,"x",grid_size,"_Getis_Ord_localG_hotspots.under_tissue_based_on_joined_Bacteria_Fungi_Expression_220405.spatial_pos.csv",sep="")  
      }
      if (marker=="all_Fungi")
      {
        localG_matrix_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Fungi_and_UNKNOWN/Under_tissue/Apr22_expression_220405/hotspots/",Exp,"_",smpl,".All_ITS_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top_50_genus.",grid_size,"x",grid_size,"_Getis_Ord_localG_hotspots.under_tissue_based_on_joined_Bacteria_Fungi_Expression_220405.spatial_pos.csv",sep="")  
      }
      if (marker=="all_Expression")
      {
        localG_matrix_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/Apr2022/hotspots/",Exp,"_",smpl,".220405_semiwild_dataset_rawcounts_filtered.SUM_ALL_EXPRESSION.",grid_size,"x",grid_size,"_Getis_Ord_localG_hotspots.under_tissue_based_on_joined_Bacteria_Fungi_Expression_220405.spatial_pos.csv",sep="")  
      }
      localG_data_for_writing=data.frame(row.names=paste(coordinates(getisgrid)[,1],coordinates(getisgrid)[,2],sep="x"))
      layer_name_for_writing=paste(marker,"_LocalG",sep="")
      localG_data_for_writing[[layer_name_for_writing]]=getisgrid$HOTSPOT
      localG_matrix_for_writing=as.data.frame(t(localG_data_for_writing))
      localG_matrix_for_writing=data.frame(tax=layer_name_for_writing,localG_matrix_for_writing)
      # str(bacteria_sum_for_writing)
      write.table(file=localG_matrix_file,x=localG_matrix_for_writing,quote = F,row.names = F,sep = ";")
      message(paste("[INFO] The localG",grid_size,"x",grid_size,"for",marker,"under the tissue was written to:",localG_matrix_file))
      
      ## also save in the agregate per section
      
      ## TODO: Should we also extract here the grix actual pixels?
      
      # add the extra data about the pValue and significance
      layer_name_for_writing.reads=paste(marker,"_reads_on_grid")
      layer_name_for_writing.p=paste(layer_name_for_writing,".p",sep="")
      layer_name_for_writing.p.SP_FDR=paste(layer_name_for_writing,".p.SP_FDR",sep="")
      localG_data_for_writing[[layer_name_for_writing.reads]]=getisgrid$layer
      localG_data_for_writing[[layer_name_for_writing.p]]=getisgrid$HOTSPOT.p
      localG_data_for_writing[[layer_name_for_writing.p.SP_FDR]]=getisgrid$HOTSPOT.p.SP_FDR
      rm (layer_name_for_writing.p,layer_name_for_writing.p.SP_FDR)
      if(NROW(section_localG_stats)==0)
      {
        section_localG_stats=data.frame(xy=row.names(localG_data_for_writing),localG_data_for_writing)
      } else {
        section_localG_stats=merge(section_localG_stats,
                                   data.frame(xy=row.names(localG_data_for_writing),localG_data_for_writing),
                                   by="xy",all=T)
      }
      # for all
      hotspots_data=data.frame(x=coordinates(getisgrid)[,1],
                               y=coordinates(getisgrid)[,2])
      hotspots_data[[marker]]=0
      hotspots_data[[marker]][getisgrid$HOTSPOT>0&getisgrid$HOTSPOT.p.SP_FDR<=0.05]=1 # hotspot
      hotspots_data[[marker]][getisgrid$HOTSPOT<0&getisgrid$HOTSPOT.p.SP_FDR<=0.05]=-1 # coldspot
      if (NROW(hotspots_df)==0)
      {
        hotspots_df=hotspots_data
      } else {
        hotspots_df=merge(hotspots_df,hotspots_data,by=c("x","y"))
      }
      rm (getisgrid,hotspots_data)
    }
    
    # write the full stats of all layers
    all_layers_stats_file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/raw_data_localG/",Exp,".",smpl,".localG_ststs.Total_expression_Top50_Microbial_Genus.tsv",sep="")
    write.table(file=all_layers_stats_file,x=section_localG_stats,row.names = F,quote = F,sep = "\t")
    rm (all_layers_stats_file,section_localG_stats)
    
    hotspots_df$xy=paste(hotspots_df$x,hotspots_df$y,sep="x")
    ## Do the Venn
    library(ComplexHeatmap)
    hotSpots_list=list(Fungi=hotspots_df$xy[hotspots_df$all_Fungi==1],
                       Bacteria=hotspots_df$xy[hotspots_df$all_Bacteria==1],
                       Genes=hotspots_df$xy[hotspots_df$all_Expression==1])
    m = make_comb_mat(hotSpots_list) # https://support.bioconductor.org/p/118557/
    cs = comb_size(m)
    # row_size = set_size(m)
    # 
    # ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))))
    # ht = draw(ht)
    # co = column_order(ht)
    # row_od = row_order(ht)
    
    # nc = ncol(m)
    
    # decorate_annotation("Intersection\nsize", {
    #  grid.text(cs[co], 
    #            x = 1:nc, 
    #            y = unit(cs[co], "native") + unit(1, "mm"), 
    #            gp = gpar(fontsize = 9), 
    #            just = "bottom",
    #            default.units = "native")
    # })
    # decorate_annotation("Set size", {
    #  grid.text(row_size[row_od], 
    #            unit(row_size[row_od], "native") + unit(1, "mm"), 
    #            rev(seq_len(length(row_size))), 
    #            default.units = "native", just = "bottom", rot = -90,
    #            gp = gpar(fontsize = 10))
    #})
    library(eulerr)
    
    intersects_hotSpots=c("Fungi" = 0, "Bacteria" = 0, "Genes" = 0,
                          "Fungi&Bacteria&Genes"=0,  
                          "Fungi&Bacteria"=0,
                          "Bacteria&Genes"=0,
                          "Fungi&Genes"=0)
    for (i in names(cs)) { # assign
      if (i=="100") {intersects_hotSpots[["Fungi"]]=cs[[i]]}
      if (i=="010") {intersects_hotSpots[["Bacteria"]]=cs[[i]]}
      if (i=="001") {intersects_hotSpots[["Genes"]]=cs[[i]]}
      if (i=="111") {intersects_hotSpots[["Fungi&Bacteria&Genes"]]=cs[[i]]}
      if (i=="110") {intersects_hotSpots[["Fungi&Bacteria"]]=cs[[i]]}
      if (i=="101") {intersects_hotSpots[["Fungi&Genes"]]=cs[[i]]}
      if (i=="011") {intersects_hotSpots[["Bacteria&Genes"]]=cs[[i]]}
    }
    
    # fit the Venn
    intersect_hotspots_fit <- euler(intersects_hotSpots)
    
    print(plot(intersect_hotspots_fit,
         quantities = list(type = c("counts", "percent"), font=8, round=2, cex=1) ,
         main=paste("Hotspots",sample_name,sep=" - "),
         fills = c("#8da0cb","#66c2a5","#fc8d62"))) # )
    

    # collect data for the intersections stats
    Sum_Hotspots_count=sum(intersects_hotSpots)
    
    record=data.frame(sample=sample_name,
                      pixels=nrow(spatial_data),
                      Sum_Bacterial_reads=sum(spatial_data$all_Bacteria),
                      Sum_Fungi_reads=sum(spatial_data$all_Fungi),
                      Sum_Expression_reads=sum(spatial_data$all_Expression),
                      Bacterial_Hotspots_count=sum(hotspots_df$all_Bacteria>0),
                      Bacterial_Coldspots_count=sum(hotspots_df$all_Bacteria<0),
                      Significant_Bacterial_count=sum(hotspots_df$all_Bacteria!=0),
                      
                      Fungal_Hotspots_count=sum(hotspots_df$all_Fungi>0),
                      Fungal_Coldspots_count=sum(hotspots_df$all_Fungi<0),
                      Significant_Fungal_count=sum(hotspots_df$all_Fungi!=0),
                      
                      Expression_Hotspots_count=sum(hotspots_df$all_Expression>0),
                      Expression_Coldspots_count=sum(hotspots_df$all_Expression<0),
                      Significant_Expression_count=sum(hotspots_df$all_Expression!=0),
                      
                      Sum_Hotspots_count=sum(intersects_hotSpots),
                      
                      Bacterial_Fungal_Expression_hotspots_intersect=intersects_hotSpots[["Fungi&Bacteria&Genes"]],
                      Bacteria_uniq_hotspots=intersects_hotSpots[["Bacteria"]],
                      Expression_uniq_hotspots=intersects_hotSpots[["Genes"]],
                      Fungi_uniq_hotspots=intersects_hotSpots[["Fungi"]],
                      Bacterial_Fungal_hotspots_intersect=intersects_hotSpots[["Fungi&Bacteria"]],
                      Bacterial_Expression_hotspots_intersect=intersects_hotSpots[["Bacteria&Genes"]],
                      Fungal_Expression_hotspots_intersect=intersects_hotSpots[["Fungi&Genes"]],
                      
                      Bacterial_Fungal_Expression_hotspots_intersect_frac=intersects_hotSpots[["Fungi&Bacteria&Genes"]]/Sum_Hotspots_count,
                      Bacteria_uniq_hotspots_frac=intersects_hotSpots[["Bacteria"]]/Sum_Hotspots_count,
                      Expression_uniq_hotspots_frac=intersects_hotSpots[["Genes"]]/Sum_Hotspots_count,
                      Fungi_uniq_hotspots_frac=intersects_hotSpots[["Fungi"]]/Sum_Hotspots_count,
                      Bacterial_Fungal_hotspots_intersect_frac=intersects_hotSpots[["Fungi&Bacteria"]]/Sum_Hotspots_count,
                      Bacterial_Expression_hotspots_intersect_frac=intersects_hotSpots[["Bacteria&Genes"]]/Sum_Hotspots_count,
                      Fungal_Expression_hotspots_intersect_frac=intersects_hotSpots[["Fungi&Genes"]]/Sum_Hotspots_count
                      )
    hotsposts_intersect_Bacteria_Fungi_Expression_stats=rbind(hotsposts_intersect_Bacteria_Fungi_Expression_stats,record)
    rm (spatial_data,marker,Sum_Hotspots_count)
    rm (intersect_hotspots_fit,intersects_hotSpots,hotSpots_list,m,cs,i)
  } # end sample
} # end experiment
dev.off()

# save teh intersections stast
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

hotsposts_intersect_Bacteria_Fungi_Expression_stats_and_label=merge(hotsposts_intersect_Bacteria_Fungi_Expression_stats,samples_label,by = "sample")

# save
# write.table(x=hotsposts_intersect_Bacteria_Fungi_Expression_stats_and_label,file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Compare_Hotspots_OMNI12_OMNI13_Expression_Fungi_Bacteria.csv",sep=";",quote = F,row.names = F)
# Apr2022
write.table(x=hotsposts_intersect_Bacteria_Fungi_Expression_stats_and_label,file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Compare_Hotspots_OMNI12_OMNI13_220405_Expression_Fungi_Bacteria.csv",sep=";",quote = F,row.names = F)

library(reshape)
mdata <- melt(hotsposts_intersect_Bacteria_Fungi_Expression_stats_and_label, id=c("sample"))
mdata_frac=mdata[grepl(x = mdata$variable,pattern = "_frac"),]
mdata_frac=merge(mdata_frac,samples_label)
mdata_frac$value=as.numeric(mdata_frac$value)
mdata_frac$value=signif(x = mdata_frac$value*100,3)
mdata_frac$variable=factor(mdata_frac$variable,levels=c("Bacterial_Fungal_Expression_hotspots_intersect_frac",
                                                 "Bacterial_Expression_hotspots_intersect_frac",
                                                 "Fungal_Expression_hotspots_intersect_frac",
                                                 "Bacterial_Fungal_hotspots_intersect_frac",
                                                 "Expression_uniq_hotspots_frac",
                                                 "Bacteria_uniq_hotspots_frac" ,
                                                 "Fungi_uniq_hotspots_frac"))

library(ggplot2)
library(RColorBrewer)
p=ggplot(mdata_frac, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
                                                        "Shared Bacterial-Fungal hotspots","Expression-unique","Bacterial-unique","Fungal-unique"),name="") +
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))

# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Fraction_of_shared_significant_Hotspots_OMNI12_OMNI13_Expression_Fungi_Bacteria.pdf",height = 7,width = 14)
# Apr2022
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Fraction_of_shared_significant_Hotspots_OMNI12_OMNI13_220405_Expression_Fungi_Bacteria.pdf",height = 7,width = 14)
print (p)
dev.off()

# OUT of the expressions hotspots
mdata_frac_expression=mdata_frac[(grepl(x=mdata_frac$variable,"Expression")),] # &!(grepl(x=mdata_frac$variable,"Expression_uniq"))
p=ggplot(mdata_frac_expression, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  # scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
  #                                                       "Shared Bacterial-Fungal hotspots","Expression-unique","Bacterial-unique","Fungal-unique"),name="") +
  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
                                                        "Expression-unique"),name="") +
  
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Fraction_of_shared_significant_Hotspots_OMNI12_OMNI13_220405_Expression_Fungi_Bacteria.shared_fraction.pdf",height = 7,width = 14)
print (p)
dev.off()

# all genes
tmp.data=hotsposts_intersect_Bacteria_Fungi_Expression_stats_and_label[,c("sample","Bacterial_Fungal_Expression_hotspots_intersect","Bacterial_Expression_hotspots_intersect",
                                                                          "Fungal_Expression_hotspots_intersect",
                                                                          "Expression_uniq_hotspots"),]

tmp.data_frac=tmp.data[,2:5]/rowSums(tmp.data[,2:5])
tmp.data_frac$sample=tmp.data$sample
mdata_expression_hotspots=melt(tmp.data_frac, id=c("sample"))
mdata_expression_hotspots=merge(mdata_expression_hotspots,samples_label)
# mdata_selected_expression_hotspots$value=as.numeric(mdata_selected_expression_hotspots$value)
mdata_expression_hotspots$value=signif(x = mdata_expression_hotspots$value*100,3)
#mdata_selected_expression_hotspots$variable=factor(mdata_selected_expression_hotspots$variable,levels=c("Bacterial_Fungal_Expression_hotspots_intersect_frac",
#                                                                          "Bacterial_Expression_hotspots_intersect_frac",
#                                                                          "Fungal_Expression_hotspots_intersect_frac",
#                                                                          "Bacterial_Fungal_hotspots_intersect_frac",
#                                                                          "Expression_uniq_hotspots_frac",
#                                                                          "Bacteria_uniq_hotspots_frac" ,
#                                                                          "Fungi_uniq_hotspots_frac"))




p=ggplot(mdata_expression_hotspots, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  #  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
  #                                                        "Shared Bacterial-Fungal hotspots","Expression-unique","Bacterial-unique","Fungal-unique"),name="") +
  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
                                                        "Expression-unique"),name="") +
  labs(y="Expression hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))

pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Fraction_of_shared_significant_Hotspots_OMNI12_OMNI13_220405_Expression_Fungi_Bacteria.shared_fraction.compare_out_of_expression_HS.pdf",height = 5,width =10)
print (p)
dev.off()

### CALCULATE PER GENE HOTSPOTS AND CREATE THEIR MATRIX (BASED ON DIFFERENT RNACutoff)
# FOR QA
# all_samples_total_gene_expression=read.delim(file="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/220405_semiwild_dataset_rawcounts_filtered.Apr2022.OMNI12_OMNI13.Total_gene_expressions_alll_sections.tsv",sep="\t",stringsAsFactors = F)
# all_samples_GeneExpression_spatial_data=readRDS("~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.spatial_objects_list.rds")
grid_size=2
for (RNA_Sum_Cutoff in c(100,10))
{
  for (Exp in c("OMNI12","OMNI13"))
  {
    samples=c("A1","A2","B1","B2","C1","C2")
    if (Exp=="OMNI13") {samples=c(samples,"D1")}
    for (smpl in samples)
    {
      sum_stat=data.frame()
      sample_name=paste(Exp,smpl,sep="_")
      message(paste("--",sample_name))
      spatial_data=all_samples_GeneExpression_spatial_data[[sample_name]]
      
      genes_above_cutoff_list=all_samples_total_gene_expression$gene[all_samples_total_gene_expression[[sample_name]]>RNA_Sum_Cutoff]
      message(paste("[INFO] Total genes above RNASum cutoof (>",RNA_Sum_Cutoff,"): ",length(genes_above_cutoff_list),sep=""))
      genes_above_cutoff_list.non_amb=genes_above_cutoff_list[!grepl(x=genes_above_cutoff_list,pattern = "ambiguous")]
      message(paste("[INFO] Total non ambigues genes above RNASum cutoof (>",RNA_Sum_Cutoff,"): ",length(genes_above_cutoff_list.non_amb),sep=""))
      sample_localG_data=data.frame()
      for (marker in genes_above_cutoff_list.non_amb)
      {
        getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = marker)
        layer_localG_data=data.frame(x=coordinates(getisgrid)[,1],y=coordinates(getisgrid)[,2])
        layer_localG_data[[paste(marker,"_localG",sep="")]]=getisgrid$HOTSPOT
       
        if (NROW(sample_localG_data)==0)
        {
          sample_localG_data=layer_localG_data
        } else {
          sample_localG_data=merge(sample_localG_data,layer_localG_data,by=c("x","y"),all=T)
        }
        
        # Summary stats
        ## total reads
        total_reads=all_samples_total_gene_expression[all_samples_total_gene_expression$gene==marker,sample_name]
        ## total significant hotspots
        total_significant_HS=sum(getisgrid$HOTSPOT>0&getisgrid$HOTSPOT.p.SP_FDR<=0.05)
        total_significant_CS=sum(getisgrid$HOTSPOT<0&getisgrid$HOTSPOT.p.SP_FDR<=0.05)
        rec=data.frame(gene=marker,total_reads=total_reads,total_significant_HS=total_significant_HS,total_significant_CS=total_significant_CS)
        sum_stat=rbind(sum_stat,rec)
        # Add Moran later?
        
        rm (getisgrid,layer_localG_data,
            total_reads,total_significant_HS,total_significant_CS,rec)
      }
      # sum_stat
      summary_statsitics_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/Apr2022/hotspots/RNASumCutoff",RNA_Sum_Cutoff,"/","220405_semiwild_dataset_rawcounts_filtered.PerGene_localG_",grid_size,"x",grid_size,".Genes_with_total_sum_grt",RNA_Sum_Cutoff,".summary_stats.",Exp,".",smpl,".tsv",sep="")
      message (paste("[INFO] Writing the summary statistics per gene with total reads above ",RNA_Sum_Cutoff," to: ",summary_statsitics_file,sep=""))
      write.table(file=summary_statsitics_file,x=sum_stat,row.names = F,quote = F,sep = "\t")
      
      # sample_localG_data
      # tab delim, rows=gene (as row.names), col=grid pixels
      # /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.A2.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.for_cor.tsv
      all_genes_localG_matrix_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/Apr2022/hotspots/RNASumCutoff",RNA_Sum_Cutoff,"/","220405_semiwild_dataset_rawcounts_filtered.PerGene_localG_",grid_size,"x",grid_size,".Genes_with_total_sum_grt",RNA_Sum_Cutoff,".",Exp,".",smpl,".tsv",sep="")
      sample_localG_matrix=as.data.frame(t(data.frame(row.names=paste(sample_localG_data$x,sample_localG_data$y,sep="x"),sample_localG_data[,3:NCOL(sample_localG_data)])))
      message (paste("[INFO] Writing the localG matrix for all the genes with total reads above ",RNA_Sum_Cutoff," to: ",all_genes_localG_matrix_file,sep=""))
      write.table(file=all_genes_localG_matrix_file,x=sample_localG_matrix,sep = "\t",quote = F,row.names = T)
      
      rm (sample_localG_matrix,all_genes_localG_matrix_file,sample_localG_data,spatial_data,genes_above_cutoff_list,marker,sample_name,
          summary_statsitics_file,sum_stat,genes_above_cutoff_list.non_amb)
    }
  }
}


### CLEAN START ####
######### Apr2022-May2022 Figure 3 DATA ########
############################################################################

## Read the Hotspots information
all_samples_spatial_data=readRDS(file = "~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Top50_Bacterial_Genus_Top50_ITS_Genus_All_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.SUMS.spatial_objects_list.rds")

#### ALL THE ANALYSES OF THIS SECTION WERE USED FOR FIGURE 3 OF THE MANUSCRIPT
#### Compare the FDR significant Hotspots for Bacterial+Fungi
experiments=c("OMNI12","OMNI13")

hotsposts_intersect_stats=data.frame()
grid_size=2
RNA_cutoff=100
# pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Bacteria_to_Fungi_Ratio_in_Hotspots_type.log2.pdf",width = 8,height = 7)
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Bacteria_to_Fungi_Ratio_in_Hotspots_type.log2.Apr2022.pdf",width = 8,height = 7)

for (exp in experiments) {
  smpl_list=c()
  if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
  if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}
  for (smpl in smpl_list)
  {
    message(paste("--",smpl))
    
    # hotspots_and_grid_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,".",grid_size,"x",grid_size,".Getis_Ord_hotspots_and_layer_data.","RNA_cutoff_",RNA_cutoff,".joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv",sep="")
    # hotspots_and_grid_data=read.delim(file = hotspots_and_grid_data_file,sep = ";",stringsAsFactors = F)
    
    hotspots_and_grid_data_file=paste("/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/raw_data_localG/",exp,".",smpl,".localG_ststs.Total_expression_Top50_Microbial_Genus.tsv",sep="")
    hotspots_and_grid_data=read.delim(file = hotspots_and_grid_data_file,sep = "\t",stringsAsFactors = F)
    
    Bacteria_Fungi_data=hotspots_and_grid_data[,grepl(pattern = "Bacteria|Fungi",x = names(hotspots_and_grid_data),ignore.case = T)]
 
    Sum_Fungi_reads=sum(Bacteria_Fungi_data$all_Fungi._reads_on_grid)
    Sum_Bacterial_reads=sum(Bacteria_Fungi_data$all_Bacteria._reads_on_grid)
    
    significant_Bacteria=Bacteria_Fungi_data[Bacteria_Fungi_data$all_Bacteria_LocalG.p.SP_FDR<=0.05,]
    significant_Fungi=Bacteria_Fungi_data[Bacteria_Fungi_data$all_Fungi_LocalG.p.SP_FDR<=0.05,]
    
    Bacterial_Hotspots_count=sum(significant_Bacteria$all_Bacteria_LocalG>0)
    Bacterial_Hotspots_grid_pixels_names=row.names(significant_Bacteria)[significant_Bacteria$all_Bacteria_LocalG>0]
    
    Bacterial_Coldspots_count=sum(significant_Bacteria$all_Bacteria_LocalG<0)
    
    Fungal_Hotspots_count=sum(significant_Fungi$all_Fungi_LocalG>0)
    Fungal_Hotspots_grid_pixels_names=row.names(significant_Fungi)[significant_Fungi$all_Fungi_LocalG>0]
    
    Fungal_Coldspots_count=sum(significant_Fungi$all_Fungi_LocalG<0)
    
    
    Bacterial_Fungal_hotspots_intersect=length(intersect(row.names(significant_Bacteria[significant_Bacteria$all_Bacteria_LocalG>0,]),row.names(significant_Fungi[significant_Fungi$all_Fungi_LocalG>0,])))
    Bacterial_Fungal_coldspots_intersect=length(intersect(row.names(significant_Bacteria[significant_Bacteria$all_Bacteria_LocalG<0,]),row.names(significant_Fungi[significant_Fungi$all_Fungi_LocalG<0,])))
    
    
    # get the Grid-pixels of the Bacteria-Fungi-shared-significant Hotstpots
    
    Bacteria_Fungi_Hotspots_intersect_data=merge(significant_Bacteria[significant_Bacteria$all_Bacteria_LocalG>0,],significant_Fungi[significant_Fungi$all_Fungi_LocalG>0,],by="row.names")
    
    Bacteria_Hotspots_unique_grid_pixels_names=Bacterial_Hotspots_grid_pixels_names[!(Bacterial_Hotspots_grid_pixels_names %in% Fungal_Hotspots_grid_pixels_names)]
    Bacteria_Hotspots_unique_grid_data=significant_Bacteria[Bacteria_Hotspots_unique_grid_pixels_names,]
    
    Fungal_Hotspots_unique_grid_pixels_names=Fungal_Hotspots_grid_pixels_names[!(Fungal_Hotspots_grid_pixels_names %in% Bacterial_Hotspots_grid_pixels_names)]
    Fungal_Hotspots_unique_grid_data=significant_Fungi[Fungal_Hotspots_unique_grid_pixels_names,]
    
    
    Bacterial_Fungal_ratios=data.frame()
    Bacterial_Fungal_ratios=rbind(Bacterial_Fungal_ratios,data.frame(Bacterial_Fungi_ratio=Bacteria_Fungi_Hotspots_intersect_data$all_Bacteria._reads_on_grid.x/Bacteria_Fungi_Hotspots_intersect_data$all_Fungi._reads_on_grid.x,type="Shared"))
    Bacterial_Fungal_ratios=rbind(Bacterial_Fungal_ratios,data.frame(Bacterial_Fungi_ratio=Fungal_Hotspots_unique_grid_data$all_Bacteria._reads_on_grid/Fungal_Hotspots_unique_grid_data$all_Fungi._reads_on_grid,type="Fungi_unique"))
    Bacterial_Fungal_ratios=rbind(Bacterial_Fungal_ratios,data.frame(Bacterial_Fungi_ratio=Bacteria_Hotspots_unique_grid_data$all_Bacteria._reads_on_grid/Bacteria_Hotspots_unique_grid_data$all_Fungi._reads_on_grid,type="Bacteria_unique"))
    
    Bacterial_Fungal_by_hotspot_type=data.frame()
    Bacterial_Fungal_by_hotspot_type=rbind(Bacterial_Fungal_by_hotspot_type,data.frame(Bacterial_count=Bacteria_Fungi_Hotspots_intersect_data$all_Bacteria._reads_on_grid.x,Fungal_count=Bacteria_Fungi_Hotspots_intersect_data$all_Fungi._reads_on_grid.x,type="Shared"))
    Bacterial_Fungal_by_hotspot_type=rbind(Bacterial_Fungal_by_hotspot_type,data.frame(Bacterial_count=Fungal_Hotspots_unique_grid_data$all_Bacteria._reads_on_grid,Fungal_count=Fungal_Hotspots_unique_grid_data$all_Fungi._reads_on_grid,type="Fungi_unique"))
    Bacterial_Fungal_by_hotspot_type=rbind(Bacterial_Fungal_by_hotspot_type,data.frame(Bacterial_count=Bacteria_Hotspots_unique_grid_data$all_Bacteria._reads_on_grid,Fungal_count=Bacteria_Hotspots_unique_grid_data$all_Fungi._reads_on_grid,type="Bacteria_unique"))
    Bacterial_Fungal_by_hotspot_type$Bacterial_log2=log(Bacterial_Fungal_by_hotspot_type$Bacterial_count,2)
    Bacterial_Fungal_by_hotspot_type$Fungal_log2=log(Bacterial_Fungal_by_hotspot_type$Fungal_count,2)
    
    
    points_count=data.frame(type=c("Shared","Fungi_unique","Bacteria_unique"),
                            num_of_hotspots=c(NROW(Bacteria_Fungi_Hotspots_intersect_data),NROW(Fungal_Hotspots_unique_grid_data),NROW(Bacteria_Hotspots_unique_grid_data)),
                            ypos=c(max(Bacteria_Fungi_Hotspots_intersect_data$all_Bacteria._reads_on_grid.x/Bacteria_Fungi_Hotspots_intersect_data$all_Fungi._reads_on_grid.x)+0.5,
                                   max(Fungal_Hotspots_unique_grid_data$all_Bacteria._reads_on_grid/Fungal_Hotspots_unique_grid_data$all_Fungi._reads_on_grid)+0.5,
                                   max(Bacteria_Hotspots_unique_grid_data$all_Bacteria._reads_on_grid/Bacteria_Hotspots_unique_grid_data$all_Fungi._reads_on_grid)+0.5)
    )
    
    # Message if there are significant HS with 0 reads...
    if (sum(Fungal_Hotspots_unique_grid_data$all_Fungi._reads_on_grid==0)>0) {message (paste("[WARNING]: in",hotspots_and_grid_data_file,"there is",sum(Fungal_Hotspots_unique_grid_data$all_Fungi._reads_on_grid==0),"Fungal hotspots with 0 Fungi reads on it"))}
    if (sum(Bacteria_Hotspots_unique_grid_data$all_Bacteria._reads_on_grid==0)>0) {message (paste("[WARNING]: in",hotspots_and_grid_data_file,"there is",sum(Bacteria_Hotspots_unique_grid_data$all_Bacteria._reads_on_grid==0),"Bacterial hotspots with 0 Bacterial reads on it"))}
    if (sum(Bacteria_Fungi_Hotspots_intersect_data$all_Fungi._reads_on_grid.x==0)>0) {message (paste("[WARNING]: in",hotspots_and_grid_data_file,"there is",sum(Bacteria_Fungi_Hotspots_intersect_data$all_Fungi._reads_on_grid.x==0),"Bacterial-Fungi hotspots with 0 Fungi reads on it"))}
    
    
    
    p=ggplot(data = Bacterial_Fungal_ratios, aes(x=type, y=Bacterial_Fungi_ratio, fill=type)) +
      geom_boxplot() +
      scale_fill_manual(values=c("#66c2a5","#8da0cb","#fc8d62")) +
      geom_jitter(color="black", size=0.4, alpha=0.9) +
      theme_bw() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      labs(x="Hotspot type",y="Bacteria to Fungi ratio",title=paste(exp,smpl,sep=".")) +
      geom_text(data = points_count, aes(label = paste("n=",num_of_hotspots," (",signif(num_of_hotspots/sum(num_of_hotspots),2)*100,"%)",sep=""), y = max(ypos)), 
                position = position_dodge(width = .75), 
                show.legend = FALSE ) +
      scale_x_discrete(labels=c("Bacteria-unique","Fungi-unique","Shared Bacteria-Fungi"))
    print(p)
    
    # as scatter
    p1=ggplot(Bacterial_Fungal_by_hotspot_type, aes(x=Fungal_log2, y=Bacterial_log2, color=type)) +
      geom_point() + 
      geom_smooth(method=glm, aes(fill=type,color=type)) +
      theme_bw()+
      theme(
        plot.title = element_text(size=11)
      ) +
      labs(x="log2 Fungal count",y="log2 Bacterial count",title=paste(exp,smpl,sep=".")) +
      scale_color_manual(values=c("#66c2a5","#8da0cb","#fc8d62"),labels=c("Bacteria-unique","Fungi-unique","Shared Bacteria-Fungi"),name="Hotspot type")+
      scale_fill_manual(values=c("#66c2a5","#8da0cb","#fc8d62"),labels=c("Bacteria-unique","Fungi-unique","Shared Bacteria-Fungi"),name="Hotspot type") 
    print(p1)  
    
    # Do e have points with Bacteria only data
    #summary()
    #wilcox.test(x = Bacteria_Hotspots_unique_grid_data$layer_all_Bacteria/Bacteria_Hotspots_unique_grid_data$layer_all_Fungi,y=Bacteria_Fungi_Hotspots_intersect_data$layer_all_Bacteria.x/Bacteria_Fungi_Hotspots_intersect_data$layer_all_Fungi.x)
    #wilcox.test(x = Fungal_Hotspots_unique_grid_data$layer_all_Bacteria/Fungal_Hotspots_unique_grid_data$layer_all_Fungi,y=Bacteria_Fungi_Hotspots_intersect_data$layer_all_Bacteria.x/Bacteria_Fungi_Hotspots_intersect_data$layer_all_Fungi.x)
    
    # get the Grid-pixels of the Fungi-unique-significant Hotstpots
    
    # get the Grid-pixels of the Bacteria-unique-significant Hotstpots
    
    record=data.frame(sample=paste(exp,smpl,sep="."),
                      pixels=nrow(hotspots_and_grid_data),
                      Sum_Bacterial_reads=Sum_Bacterial_reads,
                      Sum_Fungi_reads=Sum_Fungi_reads,
                      Significant_Bacteria_count=nrow(significant_Bacteria),
                      Bacterial_Hotspots_count=Bacterial_Hotspots_count,
                      Bacterial_Coldspots_count=Bacterial_Coldspots_count,
                      
                      Significant_Fungal_count=nrow(significant_Fungi),
                      Fungal_Hotspots_count=Fungal_Hotspots_count,
                      Fungal_Coldspots_count=Fungal_Coldspots_count,
                      
                      Bacterial_Fungal_hotspots_intersect=Bacterial_Fungal_hotspots_intersect,
                      Bacterial_Fungal_coldspots_intersect=Bacterial_Fungal_coldspots_intersect,
                      Bacterial_Fungal_hotspots_intersect_frac=Bacterial_Fungal_hotspots_intersect/(Fungal_Hotspots_count+Bacterial_Hotspots_count-Bacterial_Fungal_hotspots_intersect),
                      Bacterial_Fungal_coldspots_intersect_frac=Bacterial_Fungal_coldspots_intersect/(Fungal_Coldspots_count+Bacterial_Coldspots_count-Bacterial_Fungal_coldspots_intersect),
                      
                      # unique bacteria fraction
                      Bacterial_hotspots_unique_frac=(Bacterial_Hotspots_count-Bacterial_Fungal_hotspots_intersect)/(Fungal_Hotspots_count+Bacterial_Hotspots_count-Bacterial_Fungal_hotspots_intersect),
                      Bacterial_coldspots_unique_frac=(Bacterial_Coldspots_count-Bacterial_Fungal_coldspots_intersect)/(Fungal_Coldspots_count+Bacterial_Coldspots_count-Bacterial_Fungal_coldspots_intersect),
                      
                      # unique fungi fraction
                      Fungal_hotspots_unique_frac=(Fungal_Hotspots_count-Bacterial_Fungal_hotspots_intersect)/(Fungal_Hotspots_count+Bacterial_Hotspots_count-Bacterial_Fungal_hotspots_intersect),
                      Fungal_coldspots_unique_frac=(Fungal_Coldspots_count-Bacterial_Fungal_coldspots_intersect)/(Fungal_Coldspots_count+Bacterial_Coldspots_count-Bacterial_Fungal_coldspots_intersect),
                      
                      file=hotspots_and_grid_data_file)
    hotsposts_intersect_stats=rbind(hotsposts_intersect_stats,record)
    rm ("Bacteria_Fungi_data","Bacteria_Fungi_Hotspots_intersect_data","Bacteria_Hotspots_unique_grid_data","Bacteria_Hotspots_unique_grid_pixels_names",
        "Bacterial_Coldspots_count","Bacterial_Fungal_by_hotspot_type","Bacterial_Fungal_coldspots_intersect","Bacterial_Fungal_hotspots_intersect",
        "Bacterial_Fungal_ratios","Bacterial_Hotspots_count","Bacterial_Hotspots_grid_pixels_names",
        "Fungal_Coldspots_count","Fungal_Hotspots_count","Fungal_Hotspots_grid_pixels_names","Fungal_Hotspots_unique_grid_data",
        "Fungal_Hotspots_unique_grid_pixels_names","hotspots_and_grid_data","hotspots_and_grid_data_file","p","p1","points_count",                              
        "record","significant_Bacteria","significant_Fungi","Sum_Bacterial_reads","Sum_Fungi_reads")
  }
}
dev.off()

# percent of the array
hotsposts_intersect_stats$Bacterial_Coldspots_percent=hotsposts_intersect_stats$Bacterial_Coldspots_count/hotsposts_intersect_stats$pixels*100
hotsposts_intersect_stats$Bacterial_Hotspots_percent=hotsposts_intersect_stats$Bacterial_Hotspots_count/hotsposts_intersect_stats$pixels*100

hotsposts_intersect_stats$Fungal_Coldspots_percent=hotsposts_intersect_stats$Fungal_Coldspots_count/hotsposts_intersect_stats$pixels*100
hotsposts_intersect_stats$Fungal_Hotspots_percent=hotsposts_intersect_stats$Fungal_Hotspots_count/hotsposts_intersect_stats$pixels*100




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

hotsposts_intersect_stats=merge(hotsposts_intersect_stats,samples_label,by = "sample")


# save
# write.table(x=hotsposts_intersect_stats,file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/data/Compare_Hotspots_OMNI12_OMNI13.csv",sep=";",quote = F,row.names = F)
write.table(x=hotsposts_intersect_stats,file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Compare_Hotspots_OMNI12_OMNI13_Apr2022.csv",sep=";",quote = F,row.names = F)

library(reshape)
mdata <- melt(hotsposts_intersect_stats, id=c("sample"))

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


mdata_intersect=mdata[grepl(x = mdata$variable,pattern = "_intersect_frac"),]
mdata_intersect=merge(mdata_intersect,samples_label)
mdata_intersect$value=as.numeric(mdata_intersect$value)
mdata_intersect$value=signif(x = mdata_intersect$value*100,3)
library(ggplot2)
library(RColorBrewer)
p=ggplot(mdata_intersect, aes(fill=variable, y=value, x=lable)) + 
  geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c("#ca0020", "#0571b0"),labels=c("Hotspots","Coldspots"),name="") + 
  labs(y="Shared (%)",x="")+
  scale_y_continuous(breaks=seq(0,70,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Compare_Significant_Hotspots_OMNI12_OMNI13.pdf")
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Compare_Significant_Hotspots_OMNI12_OMNI13_Apr2022.pdf")
print (p)
dev.off()

# stacked barplot with the % of shared, unique bacteria unique Fungi
mdata_intersect_unique=mdata[grepl(x = mdata$variable,pattern = "_frac"),]
mdata_intersect_unique=merge(mdata_intersect_unique,samples_label)
mdata_intersect_unique$value=as.numeric(mdata_intersect_unique$value)
mdata_intersect_unique$value=signif(x = mdata_intersect_unique$value*100,3)

# separate into hotspots and coldspots
mdata_intersect_unique_hotspots=mdata_intersect_unique[grepl(x = mdata_intersect_unique$variable,pattern = "_hotspots_"),]
mdata_intersect_unique_coldspots=mdata_intersect_unique[grepl(x = mdata_intersect_unique$variable,pattern = "_coldspots_"),]

p=ggplot(mdata_intersect_unique_hotspots, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  
  #scale_fill_manual(labels=c("Hotspots","Coldspots"),name="") + 
  scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"),labels=c("Shared Bacterial-Fungal hotspots", "Bacterial-unique","Fungal-unique"),name="") +  # 
  
  #  scale_fill_discrete("", 
  #                      labels=c("Shared Bacterial-Fungal hotspots", "Bacterial-unique","Fungal-unique")) +
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Percent_of_Significant_shared_and_unique_Hotspots_OMNI12_OMNI13.pdf",height = 7,width = 10)
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Percent_of_Significant_shared_and_unique_Hotspots_OMNI12_OMNI13.Apr2022.pdf",height = 7,width = 10)

print (p)
dev.off()

mdata_intersect_unique_hotspots_and_coldpots=rbind(data.frame(mdata_intersect_unique_hotspots,type="hotspot"),data.frame(mdata_intersect_unique_coldspots,type="coldspot"))
p=ggplot(mdata_intersect_unique_hotspots_and_coldpots, aes(fill=variable, y=value, x=lable)) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  #scale_fill_manual(values = c("#ca0020", "#0571b0"),labels=c("Hotspots","Coldspots"),name="") + 
  labs(y="Shared (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  facet_grid( ~ type)
p

# create the flat df for hotspots of differnt types
mdatat_hotspots_percent=data.frame()
for (i in 1:nrow(hotsposts_intersect_stats))
{
  # each row will become 4 in the melted df
  record_fun=data.frame(sample=rep(hotsposts_intersect_stats$sample[i],2),
                        tax=rep("Fungi",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Fungal_Hotspots_percent[i],hotsposts_intersect_stats$Fungal_Coldspots_percent[i]))
  record_bac=data.frame(sample=rep(hotsposts_intersect_stats$sample[i],2),
                        tax=rep("Bacteria",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Bacterial_Hotspots_percent[i],hotsposts_intersect_stats$Bacterial_Coldspots_percent[i]))
  mdatat_hotspots_percent=rbind(mdatat_hotspots_percent,record_fun,record_bac)
}
mdatat_hotspots_percent=merge(mdatat_hotspots_percent,samples_label)

p=ggplot(mdatat_hotspots_percent, aes(fill=type, y=value, x=tax)) + 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c("#0571b0","#ca0020" ),name="") +  # labels=c("Bacteria","Fungi")
  labs(y="% of array",x="")+
  scale_y_continuous(breaks=seq(0,40,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_grid( ~ lable)

# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Percent_of_Significant_Hotspots_OMNI12_OMNI13.pdf",height = 7,width = 10)
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Percent_of_Significant_Hotspots_OMNI12_OMNI13.Apr2022.pdf",height = 7,width = 10)
print (p)
dev.off()

# total number of reads of each type
mdata_sum_reads=mdata[grepl(x = mdata$variable,pattern = "Sum_"),]
mdata_sum_reads=merge(mdata_sum_reads,samples_label)
mdata_sum_reads$value=as.numeric(mdata_sum_reads$value)
library(ggplot2)
library(RColorBrewer)
library(scales)
p=ggplot(mdata_sum_reads, aes(fill=variable, y=value, x=lable)) + 
  geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c( "#018571","#a6611a"),name="",labels=c("Total Bacterial reads","Total Fungi reads")) + #
  labs(y="Number of reads",x="") +
  scale_y_continuous(label=comma) +
  # scale_y_continuous(breaks=seq(0,70,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/NumberOfReads_OMNI12_OMNI13.pdf")
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/NumberOfReads_OMNI12_OMNI13.Apr2022.pdf")
print (p)
dev.off()


# interkigdom as function of shared hotspots
# based on "~/Dropbox/PostDoc/Projects/16S_array/compare_ITS_16S_networks.R"
k1_k2_counts=data.frame()
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (smpl in samples) 
  {
    net_name=paste(Exp,smpl,sep=".")
    message(paste("--",net_name))
    spatial_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Microbial_Networks/joint_Bacteria_Fungi_and_UNKNOWN_transposed/permatfull_col_shuffle/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.spearman_corr_and_permatfull_col_shuffle.1000.empP.csv",sep="")
    # hotspots_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.AND.210708_raw_counts_mtrbcpfiltered_RNA_cutoff_",RNA_cutoff,".",grid_size,"x",grid_size,".Getis_Ord_hotspots",".spearman_network.microbial_only.csv",sep="")
    # spatial_and_hotspots_significant_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.ABUNDANCE_and_HOTSPOTS_significant_agreement_network.csv",sep="")
    spatial_net_raw=read.delim(file=spatial_net_file,sep=";",stringsAsFactors = F)
    # hotspots_net_raw=read.delim(file = hotspots_net_file, sep=";",stringsAsFactors = F)
    # hotspots_net_raw=unique(hotspots_net_raw)
    # hotspots_net_raw$s1_s2=paste(hotspots_net_raw$row,hotspots_net_raw$column,sep="_")
    spatial_net_signif=spatial_net_raw[spatial_net_raw$p.BH<0.05&spatial_net_raw$empP.BH<0.05&!is.na(spatial_net_raw$empP.BH),]
    # hotspots_net_signif=hotspots_net_raw[hotspots_net_raw$p.BH<0.05&!is.na(hotspots_net_raw$p.BH),]
    spatial_net_signif$k1_k2=paste(spatial_net_signif$s1_domain,spatial_net_signif$s2_domain,sep="_")
    count_t=table(spatial_net_signif$k1_k2)
    message (paste ("[INFO]",net_name,"has the following type of interactions:",paste(names(count_t),collapse = " | ")))# QA
    count_df=data.frame(sample=net_name,Bacteria_Bacteria=count_t[["Bacteria_Bacteria"]],Bacteria_Fungi=count_t[["Bacteria_Fungi"]],Fungi_Fungi=count_t[["Fungi_Fungi"]],net_file=spatial_net_file)
    k1_k2_counts=rbind(k1_k2_counts,count_df)
    rm (net_name,spatial_net_file,spatial_net_raw,spatial_net_signif,count_t,count_df)
  }
}

k1_k2_counts$Bacteria_Bacteria_frac=k1_k2_counts$Bacteria_Bacteria/rowSums(k1_k2_counts[,c(2:4)])
k1_k2_counts$Bacteria_Fungi_frac=k1_k2_counts$Bacteria_Fungi/rowSums(k1_k2_counts[,c(2:4)])
k1_k2_counts$Fungi_Fungi_frac=k1_k2_counts$Fungi_Fungi/rowSums(k1_k2_counts[,c(2:4)])

# plot the ratio of interaction type per section
# create the flat df for hotspots of differnt types
library(reshape)
m_k1_k2_counts_all <- melt(k1_k2_counts, id=c("sample"))

m_k1_k2_frac=m_k1_k2_counts_all[grepl(x=m_k1_k2_counts_all$variable,pattern = "_frac"),]
m_k1_k2_count=m_k1_k2_counts_all[!grepl(x=m_k1_k2_counts_all$variable,pattern = "_frac|net"),]

m_k1_k2_count_and_frac=m_k1_k2_frac
m_k1_k2_count_and_frac$variable=gsub(m_k1_k2_count_and_frac$variable,pattern = "_frac",replacement = "")
names(m_k1_k2_count_and_frac)[3]="frac"

m_k1_k2_count_and_frac=merge(m_k1_k2_count_and_frac,m_k1_k2_count,by=c("sample","variable"))
m_k1_k2_count_and_frac$frac=as.numeric(m_k1_k2_count_and_frac$frac)
m_k1_k2_count_and_frac$value=as.numeric(m_k1_k2_count_and_frac$value)

m_k1_k2_count_and_frac=merge(m_k1_k2_count_and_frac,samples_label,by="sample")
# With total count
p=ggplot(m_k1_k2_count_and_frac, aes(fill=variable, y=frac*100, x=lable,label=paste(value,"\n",signif(x = frac*100,digits = 2),"%",sep=""))) + 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+ 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb"),labels=c("Bacteria-Bacteria","Bacteria-Fungi","Fungi-Fungi"),name="") +  # 
  labs(y="Microbial interactions",x="")+
  scale_y_continuous(breaks=seq(0,100,10)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))  
# facet_grid( ~ lable)
# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Significant_reads_based_microbial_interactions_by_type_OMNI12_OMNI13.pdf",height = 7,width = 10)
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Significant_reads_based_microbial_interactions_by_type_OMNI12_OMNI13.Apr2022_Hotspots.pdf",height = 7,width = 10)
print (p)
dev.off()


# with percent
mdatat_hotspots_percent=data.frame()
for (i in 1:nrow(k1_k2_counts))
{
  # each row will become 4 in the melted df
  record_fun=data.frame(sample=rep(k1_k2_counts$sample[i],2),
                        tax=rep("Fungi",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Fungal_Hotspots_percent[i],hotsposts_intersect_stats$Fungal_Coldspots_percent[i]))
  record_bac=data.frame(sample=rep(hotsposts_intersect_stats$sample[i],2),
                        tax=rep("Bacteria",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Bacterial_Hotspots_percent[i],hotsposts_intersect_stats$Bacterial_Coldspots_percent[i]))
  mdatat_hotspots_percent=rbind(mdatat_hotspots_percent,record_fun,record_bac)
  rm(record_fun,record_bac)
}
mdatat_hotspots_percent=merge(mdatat_hotspots_percent,samples_label)

p=ggplot(mdatat_hotspots_percent, aes(fill=type, y=value, x=tax)) + 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c("#0571b0","#ca0020" ),name="") +  # labels=c("Bacteria","Fungi")
  labs(y="% of array",x="")+
  scale_y_continuous(breaks=seq(0,40,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_grid( ~ lable)

# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Percent_of_Significant_Hotspots_OMNI12_OMNI13.pdf",height = 7,width = 10)
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Percent_of_Significant_Hotspots_OMNI12_OMNI13.Apr2022_Hotspots.pdf",height = 7,width = 10)
print (p)
dev.off()



hotsposts_intersect_stats_and_microbial_net_stats=merge(hotsposts_intersect_stats,k1_k2_counts,by="sample")

# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
# https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
lm_eqn <- function(df,x,y){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                   list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), 
                        pvalue = format(summary(m)$coefficients[2,4], digits = 2)))
  as.character(as.expression(eq));
}

# lm model for the Bacteria-Fungi interaction
m <- lm(Bacterial_Fungal_hotspots_intersect_frac ~ Bacteria_Fungi_frac, hotsposts_intersect_stats_and_microbial_net_stats);
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), 
                      pvalue = format(summary(m)$coefficients[2,4], digits = 2)))
hotsposts_intersect_stats_and_microbial_net_stats$sum_all_reads=hotsposts_intersect_stats_and_microbial_net_stats$Sum_Bacterial_reads+hotsposts_intersect_stats_and_microbial_net_stats$Sum_Fungi_reads
summary(hotsposts_intersect_stats_and_microbial_net_stats$sum_all_read)
p=ggplot(hotsposts_intersect_stats_and_microbial_net_stats, aes(x=Bacteria_Fungi_frac*100, y=Bacterial_Fungal_hotspots_intersect_frac*100)) + 
  geom_point(aes(size = sum_all_reads, colour = sum_all_reads))+ # colour = sum_all_reads,
  scale_size_binned() +
  # scale_fill_gradient(low = "grey", high = "red") +
  geom_text(aes(label=paste(lable," (",format(sum_all_reads,big.mark = ",",trim = TRUE),")",sep="")),hjust=0, vjust=0) +
  labs(y="Shared Bacterial-Fungal hotspots (%)",x="Bacteria-Fungi interactions (%)",size="Number of reads",color="Number of reads") +
  geom_smooth(method=lm) + 
  geom_text(x = 31, y = 64, 
            label = as.character(as.expression(eq)), parse = TRUE) +
  theme_bw()
# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI12_OMNI13_shared_Hotspots_vs_Bacteria_Fungi_interactions.pdf")
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/OMNI12_OMNI13_shared_Hotspots_vs_Bacteria_Fungi_interactions.Apr2022_Hotspots.pdf")
print (p)
dev.off()
# lm model
m <- lm(Bacterial_hotspots_unique_frac ~ Bacteria_Bacteria_frac, hotsposts_intersect_stats_and_microbial_net_stats);
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), 
                      pvalue = format(summary(m)$coefficients[2,4], digits = 2)))

p=ggplot(hotsposts_intersect_stats_and_microbial_net_stats, aes(x=Bacteria_Bacteria_frac*100, y=Bacterial_hotspots_unique_frac*100)) + 
  geom_point(aes(size = Sum_Bacterial_reads, colour = Sum_Bacterial_reads))+ # colour = sum_all_reads,
  scale_size_binned() +
  # scale_fill_gradient(low = "grey", high = "red") +
  geom_text(aes(label=paste(lable," (",format(Sum_Bacterial_reads,big.mark = ",",trim = TRUE),")",sep="")),hjust=0, vjust=0) +
  labs(y="Bacterial-unique hotspots (%)",x="Bacteria-Bacteria interactions (%)",size="Number of Bacterial reads",color="Number of Bacterial reads") +
  geom_smooth(method=lm) + 
  geom_text(x = 42, y = 78, 
            label = as.character(as.expression(eq)), parse = TRUE) +
  theme_bw()
# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI12_OMNI13_Bacterial_unique_Hotspots_vs_Bacteria_Bacteria_interactions.pdf")
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/OMNI12_OMNI13_Bacterial_unique_Hotspots_vs_Bacteria_Bacteria_interactions.Apr2022_Hotspots.pdf")
print (p)
dev.off()
# lm model
m <- lm(Fungal_hotspots_unique_frac ~ Fungi_Fungi_frac, hotsposts_intersect_stats_and_microbial_net_stats);
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), 
                      pvalue = format(summary(m)$coefficients[2,4], digits = 2)))

p=ggplot(hotsposts_intersect_stats_and_microbial_net_stats, aes(x=Fungi_Fungi_frac*100, y=Fungal_hotspots_unique_frac*100)) + 
  geom_point(aes(size = Sum_Fungi_reads, colour = Sum_Fungi_reads))+ # colour = sum_all_reads,
  scale_size_binned() +
  # scale_fill_gradient(low = "grey", high = "red") +
  geom_text(aes(label=paste(lable," (",format(Sum_Fungi_reads,big.mark = ",",trim = TRUE),")",sep="")),hjust=0, vjust=0) +
  
  labs(y="Fungal-unique hotspots (%)",x="Fungi-Fungi interactions (%)",size="Number of Fungal reads",color="Number of Fungal reads") +
  geom_smooth(method=lm) + 
  geom_text(x = 10, y = 36.5, 
            label = as.character(as.expression(eq)), parse = TRUE) +
  theme_bw()
# pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI12_OMNI13_Fungi_unique_Hotspots_vs_Fungi_Fungi_interactions.pdf")
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/OMNI12_OMNI13_Fungi_unique_Hotspots_vs_Fungi_Fungi_interactions.Apr2022_Hotspots.pdf")
print (p)
dev.off()

cor.test(hotsposts_intersect_stats_and_microbial_net_stats$Fungi_Fungi_frac,hotsposts_intersect_stats_and_microbial_net_stats$Fungal_hotspots_unique_frac)

# correlation between the number of reads and the proprtion of significan hot/coldspots
cor.test(hotsposts_intersect_stats_and_microbial_net_stats$Sum_Bacterial_reads,hotsposts_intersect_stats_and_microbial_net_stats$Bacterial_Hotspots_count,method = "spearman")
cor.test(hotsposts_intersect_stats_and_microbial_net_stats$Sum_Fungi_reads,hotsposts_intersect_stats_and_microbial_net_stats$Fungal_Hotspots_count,method = "spearman")


# save the data
write.table(x=hotsposts_intersect_stats_and_microbial_net_stats,file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/hotsposts_intersect_stats_and_microbial_net_stats.txt",sep = "\t",quote = F,row.names = F)
write.table(x=k1_k2_counts,file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Apr2022/Top50_spearman_corr_and_permatfull_col_shuffle.1000.empPmicrobial_net_stats.txt",sep = "\t",quote = F,row.names = F)

### CONTINUE UPDATE HERE!
### END FIG3 PLOTS ####



### MICROBIAL NETWORK PLOTS - MAYBE THIS SHOULD BE SOMEWHERE ELES ####
### START CLEAN ###
microbial_networks=data.frame()
species_domians_info=data.frame()
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (smpl in samples) 
  {
    net_name=paste(Exp,smpl,sep="_")
    message(paste("--",net_name))
    spatial_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Microbial_Networks/joint_Bacteria_Fungi_and_UNKNOWN_transposed/permatfull_col_shuffle/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.spearman_corr_and_permatfull_col_shuffle.1000.empP.csv",sep="")
    # hotspots_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.AND.210708_raw_counts_mtrbcpfiltered_RNA_cutoff_",RNA_cutoff,".",grid_size,"x",grid_size,".Getis_Ord_hotspots",".spearman_network.microbial_only.csv",sep="")
    # spatial_and_hotspots_significant_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.ABUNDANCE_and_HOTSPOTS_significant_agreement_network.csv",sep="")
    spatial_net_raw=read.delim(file=spatial_net_file,sep=";",stringsAsFactors = F)
    # hotspots_net_raw=read.delim(file = hotspots_net_file, sep=";",stringsAsFactors = F)
    # hotspots_net_raw=unique(hotspots_net_raw)
    # hotspots_net_raw$s1_s2=paste(hotspots_net_raw$row,hotspots_net_raw$column,sep="_")
    spatial_net_signif=spatial_net_raw[spatial_net_raw$p.BH<0.05&spatial_net_raw$empP.BH<0.05&!is.na(spatial_net_raw$empP.BH),]
    # hotspots_net_signif=hotspots_net_raw[hotspots_net_raw$p.BH<0.05&!is.na(hotspots_net_raw$p.BH),]
    spatial_net_signif$k1_k2=paste(spatial_net_signif$s1_domain,spatial_net_signif$s2_domain,sep="_")
    # compact_network=spatial_net_signif[,c("s1","s2","corr","p","p.BH","empP","empP.BH","sum00","s1_rank","s2_rank")]
    compact_network=spatial_net_signif[,c("s1","s2","corr","s1_rank","s2_rank")]
    names(compact_network)[3:5]=paste(net_name,names(compact_network)[3:5],sep=".")
    if (NROW(microbial_networks)==0) {microbial_networks=compact_network}
    else {microbial_networks=merge(microbial_networks,compact_network,by=c("s1","s2"),all=T)}
    species_domians_info=rbind(species_domians_info,spatial_net_signif[,c("s1_s2","k1_k2")])
  }
}
 microbial_networks=unique(microbial_networks)
 
# seperate col by sep https://stackoverflow.com/questions/7069076/split-column-at-delimiter-in-data-frame
library(tidyr)
species_domians_info_per_s=separate(data = species_domians_info, col = s1_s2, into = c("s1", "s2"), sep = "_")
species_domians_info_per_s=separate(data = species_domians_info_per_s, col = k1_k2, into = c("k1", "k2"), sep = "_")
species_domians_info_per_s=data.frame(s=c(species_domians_info_per_s$s1,species_domians_info_per_s$s2),k=c(species_domians_info_per_s$k1,species_domians_info_per_s$k2))
species_domians_info_per_s=unique(species_domians_info_per_s)

# add average_r and median
library(matrixStats)
microbial_networks$average_r=rowMeans(microbial_networks[,grepl(x=names(microbial_networks),pattern = ".corr",fixed = T)],na.rm = T)
microbial_networks$median_r=rowMedians(as.matrix(microbial_networks[,grepl(x=names(microbial_networks),pattern = ".corr",fixed = T)]),na.rm = T)
microbial_networks$section_count=rowSums(!is.na(microbial_networks[,grepl(x=names(microbial_networks),pattern = ".corr",fixed = T)]))

microbial_networks$s1_s2=paste(microbial_networks$s1,microbial_networks$s2,sep="_")
microbial_networks=merge(microbial_networks,species_domians_info,by="s1_s2",all.x=T)

all_setions_sub_net=microbial_networks[microbial_networks$section_count==13,]

# calciulate the rank for the taxa interacting in all sections
species_rank=data.frame(s=microbial_networks$s1,average_rnak=rowMeans(as.matrix(microbial_networks[,grepl(x=names(microbial_networks),pattern = ".s1_rank",fixed = T)])))
species_rank=species_rank[!is.na(species_rank$average_rnak),]
species_rank=unique(species_rank)
species_rank=rbind(species_rank,data.frame(s=microbial_networks$s2,average_rnak=rowMeans(microbial_networks[,grepl(x=names(microbial_networks),pattern = ".s2_rank",fixed = T)])))
species_rank=unique(species_rank)
species_rank=species_rank[!is.na(species_rank$average_rnak),]

# create the network for the nodes repeatedly associated
# https://www.jessesadler.com/post/network-analysis-with-r/
# https://slcladal.github.io/net.html # seems as nicer visualiztion

# kingdom_for_all_species

library(igraph)
all_setions_sub_net_edges=data.frame(from=all_setions_sub_net$s1,to=all_setions_sub_net$s2,weight=all_setions_sub_net$average_r)
all_setions_sub_net_nodes=data.frame(label=unique(c(all_setions_sub_net$s1,all_setions_sub_net$s2)))
all_setions_sub_net_nodes=merge(all_setions_sub_net_nodes,species_domians_info_per_s,by.x="label",by.y="s",all.x=T)
all_setions_sub_net_igraph=graph_from_data_frame(d = all_setions_sub_net_edges, vertices = all_setions_sub_net_nodes, directed = FALSE)

plot(all_setions_sub_net_igraph)

# very usefull! 
# https://ona-book.org/community.html
# https://mr.schochastics.net/material/netvizr/#layout
library(ggraph)

V(all_setions_sub_net_igraph)
set.seed(123)
p=ggraph(all_setions_sub_net_igraph) + 
  geom_edge_link(color = "grey") + # arrow = arrow(length = unit(0.2, "cm"))
  #geom_node_point(size = 2, color = "blue") +
   geom_node_point(aes(color = as.factor(k)), size = 5, 
                  show.legend = TRUE) +
  geom_node_text(aes(label = name, size = degree(all_setions_sub_net_igraph)),
    family = "serif", repel = TRUE,show.legend = FALSE)+
   labs(color = "") +
  theme_void()

network_degrees=as.data.frame(degree(all_setions_sub_net_igraph))
network_degrees=merge(network_degrees,species_domians_info_per_s,by.x="row.names",by.y="s")
names(network_degrees)=c("genus","degree","kingdom")
network_degrees=network_degrees[order(network_degrees$degree,decreasing = T),]
network_degrees=merge(network_degrees,species_rank,all.x=T,by.x="genus",by.y="s")
cor.test(network_degrees$degree,network_degrees$average_rnak,method = "spearman")
network_degrees=network_degrees[order(network_degrees$degree,decreasing = T),]
names(network_degrees)[4]="genus_average_rank"
hub_data=hub_score(all_setions_sub_net_igraph) # calculate the hub score for each node
network_degrees_and_hub_score=merge(network_degrees,as.data.frame(hub_data$vector),by.x="genus",by.y="row.names",all.x=T)
names(network_degrees_and_hub_score)[5]="hub_score"
network_degrees_and_hub_score=network_degrees_and_hub_score[order(network_degrees_and_hub_score$hub_score,decreasing = T),]

cor.test(network_degrees_and_hub_score$hub_score,network_degrees_and_hub_score$genus_average_rank,method = "spearman")

write.table(file="~/Dropbox/PostDoc/Projects/16S_array/Microbial_networks/sub_network_in_all_13_sections.degrees.txt",x=network_degrees,sep="\t",row.names = F,quote = F)
write.table(file="~/Dropbox/PostDoc/Projects/16S_array/Microbial_networks/sub_network_in_all_13_sections.degrees_and_hub_score.txt",x=network_degrees_and_hub_score,sep="\t",row.names = F,quote = F)

pdf ("~/Dropbox/PostDoc/Projects/16S_array/Microbial_networks/sub_network_in_all_13_sections.pdf",width = 9)
print(p)
dev.off()

# detect communities using Louvain
communities <- cluster_louvain(all_setions_sub_net_igraph)

# assign as a vertex property
V(all_setions_sub_net_igraph)$community <- membership(communities)
membership(communities)
sizes(communities)

# subset strong correlations
table(all_setions_sub_net[all_setions_sub_net$average_r>=0.35,"k1_k2"])
all_setions_sub_net_grt0.35=all_setions_sub_net[all_setions_sub_net$average_r>=0.35,]
all_setions_sub_net_grt0.35_edges=data.frame(from=all_setions_sub_net_grt0.35$s1,to=all_setions_sub_net_grt0.35$s2,weight=all_setions_sub_net_grt0.35$average_r)
all_setions_sub_net_grt0.35_nodes=data.frame(label=unique(c(all_setions_sub_net_grt0.35$s1,all_setions_sub_net_grt0.35$s2)))
all_setions_sub_net_grt0.35_nodes=merge(all_setions_sub_net_grt0.35_nodes,species_domians_info_per_s,by.x="label",by.y="s",all.x=T)
all_setions_sub_net_grt0.35_igraph=graph_from_data_frame(d = all_setions_sub_net_grt0.35_edges, vertices = all_setions_sub_net_grt0.35_nodes, directed = FALSE)

# plot
ggraph(all_setions_sub_net_grt0.35_igraph) + 
  geom_edge_link(color = "grey") + # arrow = arrow(length = unit(0.2, "cm"))
  #geom_node_point(size = 2, color = "blue") +
   geom_node_point(aes(color = as.factor(k)), size = 5, 
                  show.legend = TRUE) +
  geom_node_text(aes(label = name), # size = degree(all_setions_sub_net_grt0.35_igraph)
    family = "serif", repel = TRUE,show.legend = FALSE)+
   labs(color = "") +
  theme_void()

# plot the hotspots of the network members in OMNI13.B1
## Read the Full Microbial and Fungi data
library(sp)
get_spatial_object_for_all_Microbial=function(experiment,section)
{
  # read the list of under the tissue
  
  # read all Bacteria data
  all_Bacterial_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",experiment,"/spatial_tables/Bacterial_and_UNKNOWN/",experiment,"_",section,".All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos.csv",sep="")
  bacteria_data=read.delim(all_Bacterial_data_file,sep=";",stringsAsFactors = F)
  row.names(bacteria_data)=bacteria_data$tax
  bacteria_data=bacteria_data[,-1]
  t_bacteria_data=t(bacteria_data)
  bacteria_top_50_names=colnames(t_bacteria_data)[1:50]
    
  message(paste("[INFO] Total number of pixels in Bacteria data '",all_Bacterial_data_file,"': ",NROW(t_bacteria_data),sep=""))
  
  # read Fungi data  
  all_Fungi_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",experiment,"/spatial_tables/Fungi_and_UNKNOWN/",experiment,"_",section,".All_ITS_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos.csv",sep="")
  fungi_data=read.delim(all_Fungi_data_file,sep=";",stringsAsFactors = F)
  row.names(fungi_data)=fungi_data$tax
  fungi_data=fungi_data[,-1]
  t_fungi_data=t(fungi_data)
  fungi_top_50_names=colnames(t_fungi_data)[1:50]
  message(paste("[INFO] Total number of pixels in Fungi data '",all_Fungi_data_file,"': ",NROW(t_fungi_data),sep=""))
  
  # join bacteria_fungi
  bacteria_fungi_joined=merge(t_bacteria_data,t_fungi_data,by="row.names",all=T)
  row.names(bacteria_fungi_joined)=bacteria_fungi_joined$Row.names
  bacteria_fungi_joined=bacteria_fungi_joined[,-1]
  
  message(paste("[INFO] Total number of pixels in joined data: ",NROW(bacteria_fungi_joined),sep=""))
  # info the number of rows (pixels) with NA
  message(paste("[INFO] Total number of pixels with missing data in joined data (will be transfered to 0): ",sum(is.na(rowSums(bacteria_fungi_joined)))," with total of: ",sum(is.na(bacteria_fungi_joined))," cells",sep=""))
  bacteria_fungi_joined[is.na(bacteria_fungi_joined)]=0
  message(paste("[INFO] Total number of pixels with total microbial abundance 0: ",sum(rowSums(bacteria_fungi_joined)==0),sep=""))
  if (sum(rowSums(bacteria_fungi_joined)==0)>0)
  {
    bacteria_fungi_joined=bacteria_fungi_joined[rowSums(bacteria_fungi_joined)>0,]
  }
    
  # create a dataframe with the Top50 Bacteria and Fungi
  xy_list=strsplit(row.names(bacteria_fungi_joined), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=data.frame(row.names=row.names(bacteria_fungi_joined),
                          xy,
                          bacteria_fungi_joined[,c(bacteria_top_50_names,fungi_top_50_names)],
                          all_Bacteria=rowSums(bacteria_fungi_joined[,bacteria_top_50_names]),
                          all_Fungi=rowSums(bacteria_fungi_joined[,fungi_top_50_names]))
  
  # extraact only the pixels under the tissue
  under_the_tissue_pixles_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",experiment,"/spatial_tables/",experiment,"_",section,".Bacteria_Fungi_expression_220405.under_the_tissue_pixels.txt",sep="")
  message(paste ("[INFO] read the info about positions under the tissue from: ",under_the_tissue_pixles_file,sep=""))
  xy_under_the_tissue=read.delim(file=under_the_tissue_pixles_file,sep="\t",stringsAsFactors = F)
  
  data_under_the_tissue=merge(xy_under_the_tissue,data_with_xy,all.x=T,by=c("x","y"))
  message (paste ("[INFO] Final total position under the tissue for all Microbes:",NROW(xy_under_the_tissue)))
    
  # create the spatial object
  data_under_the_tissue.SP=SpatialPointsDataFrame(coords = data_under_the_tissue[,c(1:2)],data = data_under_the_tissue[,3:NCOL(data_under_the_tissue)])

  # clean up
  rm (all_Bacterial_data_file,bacteria_data,bacteria_top_50_names,t_bacteria_data,
        all_Fungi_data_file,fungi_data,fungi_top_50_names,t_fungi_data,                   
        bacteria_fungi_joined,xy_list,xy,data_with_xy,
        data_under_the_tissue,
        xy_under_the_tissue,under_the_tissue_pixles_file) 
    # data,data_file
    # gene_sum,gene_sum_ordered,expression_data_with_xy,RNA_cutoff,expression_data_with_xy_subset
  return (data_under_the_tissue.SP)
}

# spatial objects plotting
# ploting functions
library(RColorBrewer)
library(tmap)

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

get_hotspots_maps=function (spatial_data,layer,grid_size=2,title=NA,only_significant="NO")
{
  if (is.na(title)) {title=layer}
  # pdf
  # pdf (file = out_pdf,height = 7,width = 7)
  
  # calculate the localG stas
  if (length(spatial_data[[layer]])>0)
  {
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = layer)
    if (toupper(only_significant)=="YES") { # make the non significant hotspots value as NA
      getisgrid$HOTSPOT[getisgrid$HOTSPOT.p.SP_FDR>0.05]=NA
    }
    pmap=tm_shape(getisgrid) + 
      tm_fill("HOTSPOT", 
              palette = "-RdBu",
              style = "pretty", title="G stat",n=5,legend.reverse=T, midpoint=0,textNA = "NS") +
      tm_borders(alpha=.4) + tm_legend(legend.text.size=0.5,legend.title.size=0.6) +
      tm_layout(title = title,title.size=0.8)
  } else {
    pmap=NULL
  }
  return (pmap)
}
get_hotspots_maps_fixed=function (spatial_data,layer,grid_size=2,title=NA,only_significant="NO")
{
  if (is.na(title)) {title=layer}
  # pdf
  # pdf (file = out_pdf,height = 7,width = 7)
  
  # calculate the localG stas
  if (length(spatial_data[[layer]])>0)
  {
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = layer)
    if (toupper(only_significant)=="YES") { # make the non significant hotspots value as NA
      getisgrid$HOTSPOT[getisgrid$HOTSPOT.p.SP_FDR>0.05]=NA
    }
    # https://geocompr.github.io/post/2019/tmap-color-scales/
    message(paste("== ",title,"min:",min(getisgrid$HOTSPOT,na.rm = T),"max:",max(getisgrid$HOTSPOT,na.rm=T)))
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

experiment = "OMNI13"
section = "B1"
spatial_data_obj=get_spatial_object_for_all_Microbial(experiment = experiment,section = section)
plots_list=list()
for (tax in c(all_setions_sub_net_grt0.35_nodes$label,"all_Bacteria","all_Fungi"))
{
  plots_list[[tax]]=get_hotspots_maps_fixed(spatial_data = spatial_data_obj,layer = tax,grid_size = 2,title = tax,only_significant = "YES")
}
pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Microbial_networks/HotSpots_of_nodes.all_setions_sub_net_grt0.35_nodes.",experiment,"_",section,".ONLY_SIGNIFICANT.Fixed_Sacle.pdf",sep=""))
pdf(pdf_file,width = 28,height = 28)
print(tmap_arrange(plots_list, widths = c(.75, .75),ncol = 4, nrow = 4))  
dev.off()

# plot the network itself
p=ggraph(all_setions_sub_net_grt0.35_igraph) + 
  geom_edge_link(color = "grey") + # arrow = arrow(length = unit(0.2, "cm"))
  #geom_node_point(size = 2, color = "blue") +
   geom_node_point(aes(color = as.factor(k)), size = 5, 
                  show.legend = TRUE) +
  geom_node_text(aes(label = name), # size = degree(all_setions_sub_net_grt0.35_igraph)
    family = "serif", repel = TRUE,show.legend = FALSE)+
   labs(color = "") +
  theme_void()
pdf(paste("~/Dropbox/PostDoc/Projects/16S_array/Microbial_networks/all_setions_sub_net_grt0.35.pdf",sep=""))
print(p)
dev.off()

p1=ggraph(all_setions_sub_net_grt0.35_igraph,layout = "circle") + 
  geom_edge_link(color = "grey") + # arrow = arrow(length = unit(0.2, "cm"))
  #geom_node_point(size = 2, color = "blue") +
   geom_node_point(aes(color = as.factor(k)), size = 5, 
                  show.legend = TRUE) +
  geom_node_text(aes(label = name), # size = degree(all_setions_sub_net_grt0.35_igraph)
    family = "serif", repel = TRUE,show.legend = FALSE)+
   labs(color = "") +
  theme_void()
pdf(paste("~/Dropbox/PostDoc/Projects/16S_array/Microbial_networks/all_setions_sub_net_grt0.35.circle.pdf",sep=""),width = 5,height = 5)
print(p1)
dev.off()

# all_samples_spatial_data=readRDS(file = "~/Dropbox/PostDoc/Projects/16S_array/spatial_objects/OMNI12_OMNI13_Top50_Bacterial_Genus_Top50_ITS_Genus_All_Expression_Apr2022_220405_semiwild_dataset_rawcounts_filtered.SUMS.spatial_objects_list.rds")
str(spatial_data)


set.seed(123)
ggraph(all_setions_sub_net_igraph, layout = "fr") +
  geom_edge_link(color =  "grey") +
  geom_node_point(aes(color = as.factor(community)), # size = as.factor(leader)
                  show.legend = FALSE) +
  theme_void()

### CONTINUE RUN HERE! ####

### START HERE STILL TO UPDATE AFTER NEW BORUTA !!!  ####
# Do the same analyzes with the gens selected by BUROTA

# read the data and sum separately the selected genes
selected_genes=read.delim(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Host_response/FINAL_selected_by_Boruta_AND_significant_0.01_spearman_all_sections_ell_exp_Bacteria_Fungi_WO_HS_separation.detailed.txt.gene_ids",stringsAsFactors = F,header = F)
# read all data files to a list of spatial object indexed by their sample names
all_samples_spatial_data_with_selected=list()
expression_sums_stast=data.frame()
message ("## Read all data into spatial objects")
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (smpl in samples)
  {
    sample_name=paste(Exp,smpl,sep="_")
    message(paste("--",sample_name))
    data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
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
    expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",Exp,"/Filtered_CpMtRb/210708_raw_counts_mtrbcpfiltered.",Exp,".",smpl,".tsv",sep="")
    expression_data=read.delim(expression_data_file,sep="\t",stringsAsFactors = F)
    expression_data=t(expression_data) # row = position; col = gene
    
    expression_sum=data.frame(all_Expression=rowSums(expression_data))
    xy_list=strsplit(row.names(expression_data), "x")
    xy=as.data.frame(t(as.data.frame(xy_list)))
    xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
    names(xy)=c("x","y")
    expression_sum_with_xy=cbind(xy,expression_sum)
    expression_sum_with_xy$x=as.numeric(expression_sum_with_xy$x)
    expression_sum_with_xy$y=as.numeric(expression_sum_with_xy$y)
    
    selected_genes_expressions=expression_data[,colnames(expression_data) %in% selected_genes$V1]
    selected_genes_expression_sum=rowSums(selected_genes_expressions)
    expression_sum_with_xy$selected_genes_expression_sum=selected_genes_expression_sum
      
    joint_Expression_Bacteria_Fungi_data=merge(expression_sum_with_xy,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
    
    expression_stat_rec=data.frame(sample=sample_name,
                                 selected_genes_sum=sum(selected_genes_expression_sum),seleted_genes_expressed=sum(colSums(selected_genes_expressions)>0),
                                 all_expression_sum=sum(expression_sum),all_expresion_gene_expresed=sum(colSums(expression_data)>0))
    
    # save the option to have all expression for the later use...    
    # # Only genes with more than overall 10 raeds
    # gene_sum=data.frame(count=colSums(expression_data))
    # gene_sum$percent_of_data=gene_sum$count/(sum(gene_sum$count))
    ## sum(gene_sum$percent[gene_sum$count>=10])
    ## row.names(gene_sum)
    # gene_sum_ordered=gene_sum[order(gene_sum$count,decreasing = T),]
    
    ## create the spatial object
    # xy_list=strsplit(row.names(expression_data), "x")
    # xy=as.data.frame(t(as.data.frame(xy_list)))
    # xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
    # names(xy)=c("x","y")
    # expression_data_with_xy=cbind(xy,expression_data)
    # expression_data_with_xy$x=as.numeric(expression_data_with_xy$x)
    # expression_data_with_xy$y=as.numeric(expression_data_with_xy$y)
    
    # rownames(expression_data_with_xy)=rownames(expression_data)
    ## subset genes
    # RNA_cutoff=100 # Only genes with more than 100 reads are accounted
    # expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
    
    ## Join the expression and all_Fungi all_Bacteria
    # joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
    
    
    # create the spatial object
    spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
    
    #store it in the a list
    all_samples_spatial_data_with_selected[[sample_name]]=spatial_data
    
    expression_sums_stast=rbind(expression_sums_stast,expression_stat_rec)
    # clean up
    rm (sample_name,data_file,data,xy_list,xy,data_with_xy,
        expression_data_file,expression_data,expression_stat_rec,
        joint_Expression_Bacteria_Fungi_data,spatial_data) 
    # gene_sum,gene_sum_ordered,expression_data_with_xy,RNA_cutoff,expression_data_with_xy_subset
  }
}
rm (Exp,smpl)


hotsposts_intersect_Bacteria_Fungi_selected_Expression_stats=data.frame()
grid_size=2
library(eulerr)
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (smpl in samples)
  {
    hotspots_df=data.frame()
    sample_name=paste(Exp,smpl,sep="_")
    message(paste("--",sample_name))
    spatial_data=all_samples_spatial_data_with_selected[[sample_name]]
    for (marker in c("all_Bacteria","all_Fungi","selected_genes_expression_sum")) # all_Expression
    {
      getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,layer_name = marker)
      hotspots_data=data.frame(x=coordinates(getisgrid)[,1],
                               y=coordinates(getisgrid)[,2])
      hotspots_data[[marker]]=0
      hotspots_data[[marker]][getisgrid$HOTSPOT>0&getisgrid$HOTSPOT.p.SP_FDR<=0.05]=1 # hotspot
      hotspots_data[[marker]][getisgrid$HOTSPOT<0&getisgrid$HOTSPOT.p.SP_FDR<=0.05]=-1 # coldspot
      if (NROW(hotspots_df)==0)
      {
        hotspots_df=hotspots_data
      } else {
        hotspots_df=merge(hotspots_df,hotspots_data,by=c("x","y"))
      }
      rm (getisgrid,hotspots_data)
    }
    
    hotspots_df$xy=paste(hotspots_df$x,hotspots_df$y,sep="x")
    ## Do the Venn
    library(ComplexHeatmap)
    hotSpots_list=list(Fungi=hotspots_df$xy[hotspots_df$all_Fungi==1],
                       Bacteria=hotspots_df$xy[hotspots_df$all_Bacteria==1],
                       Genes=hotspots_df$xy[hotspots_df$selected_genes_expression_sum==1])
    m = make_comb_mat(hotSpots_list) # https://support.bioconductor.org/p/118557/
    cs = comb_size(m)
    # row_size = set_size(m)
    # 
    # ht = UpSet(m, top_annotation = upset_top_annotation(m, ylim = c(0, 1.1*max(cs))))
    # ht = draw(ht)
    # co = column_order(ht)
    # row_od = row_order(ht)
    
    # nc = ncol(m)
    
    # decorate_annotation("Intersection\nsize", {
    #  grid.text(cs[co], 
    #            x = 1:nc, 
    #            y = unit(cs[co], "native") + unit(1, "mm"), 
    #            gp = gpar(fontsize = 9), 
    #            just = "bottom",
    #            default.units = "native")
    # })
    # decorate_annotation("Set size", {
    #  grid.text(row_size[row_od], 
    #            unit(row_size[row_od], "native") + unit(1, "mm"), 
    #            rev(seq_len(length(row_size))), 
    #            default.units = "native", just = "bottom", rot = -90,
    #            gp = gpar(fontsize = 10))
    #})
    library(eulerr)
    
    intersects_hotSpots=c("Fungi" = 0, "Bacteria" = 0, "Genes" = 0,
                          "Fungi&Bacteria&Genes"=0,  
                          "Fungi&Bacteria"=0,
                          "Bacteria&Genes"=0,
                          "Fungi&Genes"=0)
    for (i in names(cs)) { # assign
      if (i=="100") {intersects_hotSpots[["Fungi"]]=cs[[i]]}
      if (i=="010") {intersects_hotSpots[["Bacteria"]]=cs[[i]]}
      if (i=="001") {intersects_hotSpots[["Genes"]]=cs[[i]]}
      if (i=="111") {intersects_hotSpots[["Fungi&Bacteria&Genes"]]=cs[[i]]}
      if (i=="110") {intersects_hotSpots[["Fungi&Bacteria"]]=cs[[i]]}
      if (i=="101") {intersects_hotSpots[["Fungi&Genes"]]=cs[[i]]}
      if (i=="011") {intersects_hotSpots[["Bacteria&Genes"]]=cs[[i]]}
    }
    
    # fit the Venn
    intersect_hotspots_fit <- euler(intersects_hotSpots)
    
    print(plot(intersect_hotspots_fit,
               quantities = list(type = c("counts", "percent"), font=8, round=2, cex=1) ,
               main=paste("Hotspots",sample_name,sep=" - "),
               fills = c("#8da0cb","#66c2a5","#fc8d62"))) # )
    
    
    # collect data for the intersections stats
    Sum_Hotspots_count=sum(intersects_hotSpots)
    
    record=data.frame(sample=sample_name,
                      pixels=nrow(spatial_data),
                      Sum_Bacterial_reads=sum(spatial_data$all_Bacteria),
                      Sum_Fungi_reads=sum(spatial_data$all_Fungi),
                      Sum_Expression_reads=sum(spatial_data$all_Expression),
                      Sum_slected_Expression_reads=sum(spatial_data$selected_genes_expression_sum),
                      Bacterial_Hotspots_count=sum(hotspots_df$all_Bacteria>0),
                      Bacterial_Coldspots_count=sum(hotspots_df$all_Bacteria<0),
                      Significant_Bacterial_count=sum(hotspots_df$all_Bacteria!=0),
                      
                      Fungal_Hotspots_count=sum(hotspots_df$all_Fungi>0),
                      Fungal_Coldspots_count=sum(hotspots_df$all_Fungi<0),
                      Significant_Fungal_count=sum(hotspots_df$all_Fungi!=0),
                      
                      # Expression_Hotspots_count=sum(hotspots_df$all_Expression>0),
                      # Expression_Coldspots_count=sum(hotspots_df$all_Expression<0),
                      # Significant_Expression_count=sum(hotspots_df$all_Expression!=0),
                      
                      Expression_Hotspots_count=sum(hotspots_df$selected_genes_expression_sum>0),
                      Expression_Coldspots_count=sum(hotspots_df$selected_genes_expression_sum<0),
                      Significant_Expression_count=sum(hotspots_df$selected_genes_expression_sum!=0),
                      
                      Sum_Hotspots_count=sum(intersects_hotSpots),
                      
                      Bacterial_Fungal_Expression_hotspots_intersect=intersects_hotSpots[["Fungi&Bacteria&Genes"]],
                      Bacteria_uniq_hotspots=intersects_hotSpots[["Bacteria"]],
                      Expression_uniq_hotspots=intersects_hotSpots[["Genes"]],
                      Fungi_uniq_hotspots=intersects_hotSpots[["Fungi"]],
                      Bacterial_Fungal_hotspots_intersect=intersects_hotSpots[["Fungi&Bacteria"]],
                      Bacterial_Expression_hotspots_intersect=intersects_hotSpots[["Bacteria&Genes"]],
                      Fungal_Expression_hotspots_intersect=intersects_hotSpots[["Fungi&Genes"]],
                      
                      Bacterial_Fungal_Expression_hotspots_intersect_frac=intersects_hotSpots[["Fungi&Bacteria&Genes"]]/Sum_Hotspots_count,
                      Bacteria_uniq_hotspots_frac=intersects_hotSpots[["Bacteria"]]/Sum_Hotspots_count,
                      Expression_uniq_hotspots_frac=intersects_hotSpots[["Genes"]]/Sum_Hotspots_count,
                      Fungi_uniq_hotspots_frac=intersects_hotSpots[["Fungi"]]/Sum_Hotspots_count,
                      Bacterial_Fungal_hotspots_intersect_frac=intersects_hotSpots[["Fungi&Bacteria"]]/Sum_Hotspots_count,
                      Bacterial_Expression_hotspots_intersect_frac=intersects_hotSpots[["Bacteria&Genes"]]/Sum_Hotspots_count,
                      Fungal_Expression_hotspots_intersect_frac=intersects_hotSpots[["Fungi&Genes"]]/Sum_Hotspots_count
    )
    hotsposts_intersect_Bacteria_Fungi_selected_Expression_stats=rbind(hotsposts_intersect_Bacteria_Fungi_selected_Expression_stats,record)
    rm (spatial_data,marker,Sum_Hotspots_count)
    rm (intersect_hotspots_fit,intersects_hotSpots,hotSpots_list,m,cs,i)
  } # end sample
} # end experiment

# save teh intersections stast
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

hotsposts_intersect_Bacteria_Fungi_selected_Expression_stats_and_label=merge(hotsposts_intersect_Bacteria_Fungi_selected_Expression_stats,samples_label,by = "sample")


# save
write.table(x=hotsposts_intersect_Bacteria_Fungi_selected_Expression_stats_and_label,file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Compare_Hotspots_OMNI12_OMNI13_BORUTA_selected_Expression_Fungi_Bacteria.csv",sep=";",quote = F,row.names = F)


library(reshape)
mdata_selected <- melt(hotsposts_intersect_Bacteria_Fungi_selected_Expression_stats_and_label, id=c("sample"))
mdata_selected_frac=mdata_selected[grepl(x = mdata_selected$variable,pattern = "_frac"),]
mdata_selected_frac=merge(mdata_selected_frac,samples_label)
mdata_selected_frac$value=as.numeric(mdata_selected_frac$value)
mdata_selected_frac$value=signif(x = mdata_selected_frac$value*100,3)
mdata_selected_frac$variable=factor(mdata_selected_frac$variable,levels=c("Bacterial_Fungal_Expression_hotspots_intersect_frac",
                                                        "Bacterial_Expression_hotspots_intersect_frac",
                                                        "Fungal_Expression_hotspots_intersect_frac",
                                                        "Bacterial_Fungal_hotspots_intersect_frac",
                                                        "Expression_uniq_hotspots_frac",
                                                        "Bacteria_uniq_hotspots_frac" ,
                                                        "Fungi_uniq_hotspots_frac"))

library(ggplot2)
library(RColorBrewer)
p=ggplot(mdata_selected_frac, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
                                                        "Shared Bacterial-Fungal hotspots","Expression-unique","Bacterial-unique","Fungal-unique"),name="") +
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))

pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Fraction_of_shared_significant_Hotspots_OMNI12_OMNI13_BORUTA_selected_Expression_Fungi_Bacteria.pdf",height = 7,width = 14)
print (p)
dev.off()


# Compare only the fraction of shared expression
mdata_selected_frac_expression=mdata_selected_frac[(grepl(x=mdata_selected_frac$variable,"Expression")),] # &!(grepl(x=mdata_selected_frac$variable,"Expression_uniq"))
p=ggplot(mdata_selected_frac_expression, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
#  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
#                                                        "Shared Bacterial-Fungal hotspots","Expression-unique","Bacterial-unique","Fungal-unique"),name="") +
  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
                                                        "Expression-unique",name="")) +
  
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Fraction_of_shared_significant_Hotspots_OMNI12_OMNI13_BORUTA_selected_Expression_Fungi_Bacteria.shared_fractions.pdf",height = 7,width = 14)
print (p)
dev.off()


# Make the 100% all the expression hotspots and indicate what fraction is unique, shred with Bacteria, shared with Fungi
tmp.data=hotsposts_intersect_Bacteria_Fungi_selected_Expression_stats_and_label[,c("sample","Bacterial_Fungal_Expression_hotspots_intersect","Bacterial_Expression_hotspots_intersect",
                                                                                   "Fungal_Expression_hotspots_intersect",
                                                                                   "Expression_uniq_hotspots"),]

tmp.data_frac=tmp.data[,2:5]/rowSums(tmp.data[,2:5])
tmp.data_frac$sample=tmp.data$sample
mdata_selected_expression_hotspots=melt(tmp.data_frac, id=c("sample"))
mdata_selected_expression_hotspots=merge(mdata_selected_expression_hotspots,samples_label)
# mdata_selected_expression_hotspots$value=as.numeric(mdata_selected_expression_hotspots$value)
mdata_selected_expression_hotspots$value=signif(x = mdata_selected_expression_hotspots$value*100,3)
#mdata_selected_expression_hotspots$variable=factor(mdata_selected_expression_hotspots$variable,levels=c("Bacterial_Fungal_Expression_hotspots_intersect_frac",
#                                                                          "Bacterial_Expression_hotspots_intersect_frac",
#                                                                          "Fungal_Expression_hotspots_intersect_frac",
#                                                                          "Bacterial_Fungal_hotspots_intersect_frac",
#                                                                          "Expression_uniq_hotspots_frac",
#                                                                          "Bacteria_uniq_hotspots_frac" ,
#                                                                          "Fungi_uniq_hotspots_frac"))




p=ggplot(mdata_selected_expression_hotspots, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  #  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
  #                                                        "Shared Bacterial-Fungal hotspots","Expression-unique","Bacterial-unique","Fungal-unique"),name="") +
  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
                                                        "Expression-unique"),name="") +
  labs(y="Expression hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))


### HERE TO CONTINUE...

mdata_expression_hotspots$group="ALL"
mdata_selected_expression_hotspots$group="SELECTED"

mdata_expression_compare=rbind(mdata_expression_hotspots,mdata_selected_expression_hotspots)

p=ggplot(mdata_expression_compare, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  #  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
  #                                                        "Shared Bacterial-Fungal hotspots","Expression-unique","Bacterial-unique","Fungal-unique"),name="") +
  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
                                                        "Expression-unique"),name="") +
  
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + facet_wrap(~group)

pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Fraction_of_shared_significant_Hotspots_OMNI12_OMNI13_BORUTA_Selected_and_Expression_Fungi_Bacteria.shared_fraction.compare_out_of_expression_HS.pdf",height = 7,width = 18)
print (p)
dev.off()




mdata_frac_expression_compare=rbind(mdata_frac_expression,mdata_selected_frac_expression)

p=ggplot(mdata_frac_expression_compare, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
#  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
#                                                        "Shared Bacterial-Fungal hotspots","Expression-unique","Bacterial-unique","Fungal-unique"),name="") +
  scale_fill_brewer(type="qual",palette="Set1",labels=c("Shared Expression-Bacterial-Fungal hotspots", "Shared Expression-Bacterial hotspots","Shared Expression-Fungal hotspots",
                                                        "Expression-unique"),name="") +
                                                        
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + facet_wrap(~group)

pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/HotSpots/Fraction_of_shared_significant_Hotspots_OMNI12_OMNI13_BORUTA_Selected_and_Expression_Fungi_Bacteria.shared_fraction.compare.pdf",height = 7,width = 18)
print (p)
dev.off()



### DRAFTS
#type = 
#    type = "seq",
#    palette = 1,
#    direction = 1,
#    aesthetics = "fill"
#  )
#  7-class Paired

  #scale_fill_manual(labels=c("Hotspots","Coldspots"),name="") + 
  #scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"),labels=c("Shared Bacterial-Fungal hotspots", "Bacterial-unique","Fungal-unique"),name="") +  # 
  
  #  scale_fill_discrete("", 
  #                      labels=c("Shared Bacterial-Fungal hotspots", "Bacterial-unique","Fungal-unique")) +
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))





p=ggplot(mdata_intersect, aes(fill=variable, y=value, x=lable)) + 
  geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  # scale_fill_manual(values = c("#ca0020", "#0571b0"),labels=c("Hotspots","Coldspots"),name="") + 
  labs(y="Shared (%)",x="")+
  scale_y_continuous(breaks=seq(0,70,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))

### DRAFTS

library(reshape)
mdata <- melt(hotsposts_intersect_stats, id=c("sample"))

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


mdata_intersect=mdata[grepl(x = mdata$variable,pattern = "_intersect_frac"),]
mdata_intersect=merge(mdata_intersect,samples_label)
mdata_intersect$value=as.numeric(mdata_intersect$value)
mdata_intersect$value=signif(x = mdata_intersect$value*100,3)
library(ggplot2)
library(RColorBrewer)
p=ggplot(mdata_intersect, aes(fill=variable, y=value, x=lable)) + 
  geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c("#ca0020", "#0571b0"),labels=c("Hotspots","Coldspots"),name="") + 
  labs(y="Shared (%)",x="")+
  scale_y_continuous(breaks=seq(0,70,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Compare_Significant_Hotspots_OMNI12_OMNI13.pdf")
print (p)
dev.off()

# stacked barplot with the % of shared, unique bacteria unique Fungi
mdata_intersect_unique=mdata[grepl(x = mdata$variable,pattern = "_frac"),]
mdata_intersect_unique=merge(mdata_intersect_unique,samples_label)
mdata_intersect_unique$value=as.numeric(mdata_intersect_unique$value)
mdata_intersect_unique$value=signif(x = mdata_intersect_unique$value*100,3)

# separate into hotspots and coldspots
mdata_intersect_unique_hotspots=mdata_intersect_unique[grepl(x = mdata_intersect_unique$variable,pattern = "_hotspots_"),]
mdata_intersect_unique_coldspots=mdata_intersect_unique[grepl(x = mdata_intersect_unique$variable,pattern = "_coldspots_"),]

p=ggplot(mdata_intersect_unique_hotspots, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  
  #scale_fill_manual(labels=c("Hotspots","Coldspots"),name="") + 
  scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"),labels=c("Shared Bacterial-Fungal hotspots", "Bacterial-unique","Fungal-unique"),name="") +  # 
  
  #  scale_fill_discrete("", 
  #                      labels=c("Shared Bacterial-Fungal hotspots", "Bacterial-unique","Fungal-unique")) +
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Percent_of_Significant_shared_and_unique_Hotspots_OMNI12_OMNI13.pdf",height = 7,width = 10)
print (p)
dev.off()

mdata_intersect_unique_hotspots_and_coldpots=rbind(data.frame(mdata_intersect_unique_hotspots,type="hotspot"),data.frame(mdata_intersect_unique_coldspots,type="coldspot"))
p=ggplot(mdata_intersect_unique_hotspots_and_coldpots, aes(fill=variable, y=value, x=lable)) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  #scale_fill_manual(values = c("#ca0020", "#0571b0"),labels=c("Hotspots","Coldspots"),name="") + 
  labs(y="Shared (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  facet_grid( ~ type)
p

# create the flat df for hotspots of differnt types
mdatat_hotspots_percent=data.frame()
for (i in 1:nrow(hotsposts_intersect_stats))
{
  # each row will become 4 in the melted df
  record_fun=data.frame(sample=rep(hotsposts_intersect_stats$sample[i],2),
                        tax=rep("Fungi",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Fungal_Hotspots_percent[i],hotsposts_intersect_stats$Fungal_Coldspots_percent[i]))
  record_bac=data.frame(sample=rep(hotsposts_intersect_stats$sample[i],2),
                        tax=rep("Bacteria",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Bacterial_Hotspots_percent[i],hotsposts_intersect_stats$Bacterial_Coldspots_percent[i]))
  mdatat_hotspots_percent=rbind(mdatat_hotspots_percent,record_fun,record_bac)
}
mdatat_hotspots_percent=merge(mdatat_hotspots_percent,samples_label)

p=ggplot(mdatat_hotspots_percent, aes(fill=type, y=value, x=tax)) + 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c("#0571b0","#ca0020" ),name="") +  # labels=c("Bacteria","Fungi")
  labs(y="% of array",x="")+
  scale_y_continuous(breaks=seq(0,40,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_grid( ~ lable)

pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Percent_of_Significant_Hotspots_OMNI12_OMNI13.pdf",height = 7,width = 10)
print (p)
dev.off()






### Do it as stackplot

#### END PLOTTING ALL EXPREESION ALL FUNGI ALL BACTERIA HOTSPOTS

  
  
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
}
# defense related genes
# dataset_to_plot=final_selected_genes_from_all_sections_with_all_details_annotated[grepl(x=final_selected_genes_from_all_sections_with_all_details_annotated$Gene.ontology..biological.process.,pattern = "defense"),]
# dataset_to_plot=dataset_to_plot[order(dataset_to_plot$sum_sections_Bacteria+dataset_to_plot$sum_sections_Fungi,
#                                       pmax(dataset_to_plot$r_LocalG.Bacteria.max_value,dataset_to_plot$r_LocalG.Fungi.max_value),
#                                       decreasing = T),]
# pdf_file=(paste("~/Dropbox/PostDoc/Projects/16S_array/Host_response/defense_response_selected_gene.pdf"))
# pdf(pdf_file,width = 21,height = 14) 

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
  if (i%%10==0) {message(paste("Finished ploting",i,"out of",NROW(dataset_to_plot),sep=" "))}
  flush.console()
  # tmap_mode(current.mode)
}

dev.off()



####





### OMNI6 C2 example
hotspost_df=data.frame()

microbial_data_file="/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI6/spatial_tables/Bacterial/191028-062_C2_S6_Ar_multimap_bacteria_R2.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos.txt"
expression_data_file="/Users/hashkenazy/SeaFile/Seafile/.seafile-data/file-cache/210d1574-abc7-4d2b-a9c3-da265cfa59b9/Data/Lableaf_OMNI6/211111_lableaf5050_rawdata.csv"
expression_data=read.delim(expression_data_file,sep=",",stringsAsFactors = F)
microbial_data=read.delim(microbial_data_file,sep=";",stringsAsFactors = F)

microbial_data=t(microbial_data) # row = position; col = gene
colnames(microbial_data)=microbial_data[1,]
microbial_data=microbial_data[-1,]
microbial_data=as.data.frame(microbial_data)
microbial_data_Pseudomonas=data.frame(row.names = row.names(microbial_data),Pseudomonas=as.numeric(microbial_data[,1]))
xy_list_Pseudomonas=strsplit(row.names(microbial_data_Pseudomonas), "x")
xy_Pseudomonas=as.data.frame(t(as.data.frame(xy_list_Pseudomonas)))
xy_Pseudomonas$V1=gsub(pattern = "X",replacement = "",x = xy_Pseudomonas$V1)
names(xy_Pseudomonas)=c("x","y")
microbial_data_Pseudomonas=cbind(microbial_data_Pseudomonas,xy_Pseudomonas)
microbial_data_Pseudomonas$x=as.numeric(microbial_data_Pseudomonas$x)
microbial_data_Pseudomonas$y=as.numeric(microbial_data_Pseudomonas$y)

expression_data=t(expression_data) # row = position; col = gene
colnames(expression_data)=expression_data[1,]
expression_data=expression_data[-1,]
expression_data=as.data.frame(expression_data)
for (i in 1:ncol(expression_data))
{
  expression_data[,i]=as.numeric(expression_data[,i])
}
summary(expression_data[1:10,1:10])
# Add associated genes
genes=c("AT3G41768","AT2G14610","AT1G29930","AT4G02380","AT3G50970",'AT3G02750',"AT4G30280","AT5G49720",'AT5G13420',"AT3G18060","AT2G16660", 
        "AT4G26690","AT2G18440",'AT3G02740',"AT2G27080","AT3G57230","AT4G26710","AT5G56750","AT3G17390","AT1G12920")
features_df=data.frame(row.names=row.names(expression_data),nFeature=rowSums(expression_data>0),nCount=rowSums(expression_data),expression_data[,genes])
row.names(features_df)=gsub(pattern = "_1$",perl = T,replacement = "",x = row.names(features_df))
xy_list=strsplit(row.names(features_df), "x")
xy=as.data.frame(t(as.data.frame(xy_list)))
xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
names(xy)=c("x","y")
features_data_with_xy=cbind(xy,features_df)
features_data_with_xy$x=as.numeric(features_data_with_xy$x)
features_data_with_xy$y=as.numeric(features_data_with_xy$y)

rownames(features_data_with_xy)=rownames(features_df)
features_data_with_xy$log_nFeature=log(features_data_with_xy$nFeature+1)
features_data_with_xy$log_nCount=log(features_data_with_xy$nCount+1)

features_data_with_xy=merge(features_data_with_xy,microbial_data_Pseudomonas,by=c("x","y"))
features_data_with_xy$log_Pseudomonas=log(features_data_with_xy$Pseudomonas+1)

spatial_data=SpatialPointsDataFrame(coords = features_data_with_xy[,c(1:2)],data = features_data_with_xy[,3:NCOL(features_data_with_xy)])

out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI6/spatial_tables/Bacterial/191028-062_C2_S6_Ar_multimap_bacteria_R2.usearch_unique_vs_NT_Jan2021.UMI_filtered.Pseudomoans_and_Expression.hotspots.pdf",sep="")
out_hotsposts_df=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI6/spatial_tables/Bacterial/191028-062_C2_S6_Ar_multimap_bacteria_R2.usearch_unique_vs_NT_Jan2021.UMI_filtered.Pseudomoans_and_Expression.hotspots.csv",sep="")

pdf(out_pdfs,width = 7,height = 7)
for (layer in names (spatial_data))
{
  grid_size=2
  p=spplot(spatial_data,col='transparent',zcol=layer,main=list(label=layer,cex=1)) 
  print(p)
  getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
  
  brewer.pal(n = 3, name = "RdBu")
  breaks = c(-20, -1.96, -1, 1, 1.96, 20)
  palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
  col = palette[cut(getisgrid$HOTSPOT, breaks)]
  plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
  legend("bottom", inset=.02, title="Z score (local spatial G[i] statistic)",
         legend =c("<-1.96","-1","1","1.96",">1.96"), 
         fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
         xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
  
  pmap=tm_shape(getisgrid) + 
    tm_fill("HOTSPOT", 
            palette = "-RdBu",
            style = "pretty", title="G stat",n=5,legend.reverse=T) +
    tm_borders(alpha=.4) + tm_legend(legend.text.size=0.5,legend.title.size=0.6)# +
  print (pmap)
  
  getisgrid.df=as.data.frame(getisgrid)
  names(getisgrid.df)=paste(names(getisgrid.df),layer,sep="_")
  getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
  names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
  getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
  getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
  row.names(getisgrid.df)=getisgrid.df$xy
  if (nrow(hotspost_df)==0) {
    hotspost_df=getisgrid.df[,c(2:8)]
  } else {
    hotspost_df=merge(hotspost_df,getisgrid.df[,c(2:5)],by="row.names")
    row.names(hotspost_df)=hotspost_df$Row.names
    hotspost_df=hotspost_df[,-1]
  }
}
dev.off()
write.table(x=hotspost_df,file=out_hotsposts_df,quote = F,sep = ";",row.names = T)

#### DRAFTS.....



nFeature=as.data.frame(rowSums(expression_data>0))
features_df=merge(features_df,nFeature,by="row.names",all=T)

expression_data=sapply( expression_data, as.numeric )
as.numeric()


nFeature=rowsum(expression_data>0)

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
# subset genes
RNA_cutoff=100 # Only genes with more than 100 reads are accounted
expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
spatial_data=SpatialPointsDataFrame(coords = expression_data_with_xy_subset[,c(1:2)],data = expression_data_with_xy_subset[,3:NCOL(expression_data_with_xy_subset)])
for (i in 1:length(names(spatial_data)))
{
  layer=names(spatial_data)[i]
  if (i%%100==0) {
    message(paste("\t--",i,"out of",length(names(spatial_data))))
  }
  grid_size=2
  getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
  getisgrid.df=as.data.frame(getisgrid)
  names(getisgrid.df)=paste(names(getisgrid.df),layer,sep="_")
  getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
  names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
  getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
  getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
  row.names(getisgrid.df)=getisgrid.df$xy
  getisgrid.df=getisgrid.df[,-1]        
  
  if (nrow(hotspost_df)==0)
  {
    hotspost_df=getisgrid.df
  } else {
    hotspost_df=merge(hotspost_df,getisgrid.df[,1:2],by="row.names",all=T)
    row.names(hotspost_df)=hotspost_df$Row.names
    hotspost_df=hotspost_df[,-1]
  }
}
# To save the hotspots df
# Once full with the layers 
hotspots_and_grid_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,grid_size,"x",grid_size,".Getis_Ord_hotspots_and_layer_data.","RNA_cutoff_",RNA_cutoff,".tsv",sep="")
write.table(x=hotspost_df,file = hotspots_and_grid_data_file,row.names = T,quote = F,sep = ";")

# Once only as a mtrix of the localG values
hotspots_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,grid_size,"x",grid_size,".Getis_Ord_hotspots_data.","RNA_cutoff_",RNA_cutoff,".tsv",sep="")
hotspots_no_layer=hotspost_df[,!grepl(pattern = "layer",x = colnames(hotspost_df))]
write.table(x=hotspots_no_layer,file = hotspots_data_file,row.names = T,quote = F,sep = ";")




############# HOTSPOTS FOR ONLY EXPRESSION DATA
############# START EXPRESSION AND BACTERIAL AND FUNGI #######
#####

hotspots_overlaps_per_grid_size=data.frame()
exp="OMNI12"
exp="OMNI13"

smpl_list=c()
if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}

for (smpl in smpl_list)
{
  message(paste("== ",smpl))
  hotspost_df=data.frame()
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
  # subset genes
  RNA_cutoff=100 # Only genes with more than 100 reads are accounted
  expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
  spatial_data=SpatialPointsDataFrame(coords = expression_data_with_xy_subset[,c(1:2)],data = expression_data_with_xy_subset[,3:NCOL(expression_data_with_xy_subset)])
  for (i in 1:length(names(spatial_data)))
  {
    layer=names(spatial_data)[i]
    if (i%%100==0) {
      message(paste("\t--",i,"out of",length(names(spatial_data))))
    }
    grid_size=2
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
    getisgrid.df=as.data.frame(getisgrid)
    names(getisgrid.df)=paste(names(getisgrid.df),layer,sep="_")
    getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
    names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
    getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
    getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
    row.names(getisgrid.df)=getisgrid.df$xy
    getisgrid.df=getisgrid.df[,-1]        
  
    if (nrow(hotspost_df)==0)
    {
      hotspost_df=getisgrid.df
    } else {
      hotspost_df=merge(hotspost_df,getisgrid.df[,1:2],by="row.names",all=T)
      row.names(hotspost_df)=hotspost_df$Row.names
      hotspost_df=hotspost_df[,-1]
    }
  }
  # To save the hotspots df
  # Once full with the layers 
  hotspots_and_grid_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,grid_size,"x",grid_size,".Getis_Ord_hotspots_and_layer_data.","RNA_cutoff_",RNA_cutoff,".tsv",sep="")
  write.table(x=hotspost_df,file = hotspots_and_grid_data_file,row.names = T,quote = F,sep = ";")
  
  # Once only as a mtrix of the localG values
  hotspots_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,grid_size,"x",grid_size,".Getis_Ord_hotspots_data.","RNA_cutoff_",RNA_cutoff,".tsv",sep="")
  hotspots_no_layer=hotspost_df[,!grepl(pattern = "layer",x = colnames(hotspost_df))]
  write.table(x=hotspots_no_layer,file = hotspots_data_file,row.names = T,quote = F,sep = ";")
  
}

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

# For each section take the genes hotspots and calculate the correlation network
exp="OMNI12"
exp="OMNI13"
grid_size=2
RNA_cutoff=100

smpl_list=c()
if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}

for (smpl in smpl_list)
{
  message (paste ("--",smpl))
  hotspots_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,grid_size,"x",grid_size,".Getis_Ord_hotspots_data.","RNA_cutoff_",RNA_cutoff,".tsv",sep="")
  hotspots_no_layer_genes=read.delim(file = hotspots_data_file,sep = ";",stringsAsFactors = F)
  
  hotspots_no_layer_genes=hotspots_no_layer[,c(1,5:ncol(hotspots_no_layer))]
  hotspots_no_layer_genes_network=rcorr(x = as.matrix(hotspots_no_layer_genes),type = "spearman")
  hotspots_no_layer_genes_network_flat=flattenCorrMatrix(hotspots_no_layer_genes_network$r,hotspots_no_layer_genes_network$P)
  
  network_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,grid_size,"x",grid_size,".Getis_Ord_hotspots_data.","RNA_cutoff_",RNA_cutoff,".spearman_network.tsv",sep="")
  write.table(x=hotspots_no_layer_genes_network_flat,file = network_file,row.names = T,quote = F,sep = ";")
  message (network_file)
}

# For the all_Bacterial all_Fungi
exp="OMNI12"
exp="OMNI13"
smpl_list=c()
if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}

for (smpl in smpl_list)
{
  message(paste("--",smpl))
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

  # create the spatial object
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,3:NCOL(data_with_xy)])
  for (i in 1:length(names(spatial_data)))
  {
    layer=names(spatial_data)[i]
    if (i%%10==0) {
      message(paste("\t--",i,"out of",length(names(spatial_data))))
    }
    grid_size=2
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
    getisgrid.df=as.data.frame(getisgrid)
    names(getisgrid.df)=paste(names(getisgrid.df),layer,sep="_")
    getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
    names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
    getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
    getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
    row.names(getisgrid.df)=getisgrid.df$xy
    getisgrid.df=getisgrid.df[,-1]        
    
    if (nrow(hotspost_df)==0)
    {
      hotspost_df=getisgrid.df
    } else {
      hotspost_df=merge(hotspost_df,getisgrid.df[,1:2],by="row.names",all=T)
      row.names(hotspost_df)=hotspost_df$Row.names
      hotspost_df=hotspost_df[,-1]
    }
  }
  # To save the hotspots df
  # Once full with the layers 

  hotspots_and_grid_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.",grid_size,"x",grid_size,".Getis_Ord_hotspots_and_layer_data",".tsv",sep="")
  write.table(x=hotspost_df,file = hotspots_and_grid_data_file,row.names = T,quote = F,sep = ";")
  
  # Once only as a mtrix of the localG values
  hotspots_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.",grid_size,"x",grid_size,".Getis_Ord_hotspots",".tsv",sep="")
  hotspots_no_layer=hotspost_df[,!grepl(pattern = "layer",x = colnames(hotspost_df))]
  write.table(x=hotspots_no_layer,file = hotspots_data_file,row.names = T,quote = F,sep = ";")
  
  # Do the networks
  hotspots_no_layer_data=hotspots_no_layer[,c(1,5:ncol(hotspots_no_layer))]
  hotspots_no_layer_network=rcorr(x = as.matrix(hotspots_no_layer_data),type = "spearman")
  hotspots_no_layer_network_flat=flattenCorrMatrix(hotspots_no_layer_network$r,hotspots_no_layer_network$P)
  
  # add the total abundance and taxonomy of each
  total_sums=as.data.frame(colSums(data_with_xy))
  row.names(total_sums)=paste("HOTSPOT_",row.names(total_sums),sep="")
  names(total_sums)=("row_total_sum")
  
  domain=data.frame(genus=c(paste("HOTSPOT_",names(data_with_xy)[3:102],sep=""),"HOTSPOT_all_Bacteria","HOTSPOT_all_Fungi"),domain=c(rep("Bacteria",50),rep("Fungi",50),"Bacteria","Fungi"))
  for (x in unique(hotspots_no_layer_network_flat$row))
  {
    hotspots_no_layer_network_flat$row_total_sum[hotspots_no_layer_network_flat$row==x]=total_sums$row_total_sum[row.names(total_sums)==x]
    hotspots_no_layer_network_flat$row_domain[hotspots_no_layer_network_flat$row==x]=domain$domain[domain$genus==x]
  }
  for (x in unique(hotspots_no_layer_network_flat$column))
  {
    hotspots_no_layer_network_flat$column_total_sum[hotspots_no_layer_network_flat$column==x]=total_sums$row_total_sum[row.names(total_sums)==x]
    hotspots_no_layer_network_flat$column_domain[hotspots_no_layer_network_flat$column==x]=domain$domain[domain$genus==x]
  }

  network_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.",grid_size,"x",grid_size,".Getis_Ord_hotspots",".spearman_network.tsv",sep="")
  write.table(x=hotspots_no_layer_network_flat,file = network_file,row.names = T,quote = F,sep = ";")
  
}

##### CLEAN START HOTSPOTS FOR joined Bacteria+Fungi+Expression ####
##### LAST RUN: 8Feb2022 After the bug fixed
##### Once for joined Bacteria+Fungi+Expression -> To avoid problems later when we do the correlations
# For the all_Bacterial all_Fungi
#
# new test based on: https://michaelminn.net/tutorials/r-point-analysis/index.html
library(sp)
library(raster)
library(spdep)
library(classInt)

Getis_Ord_GI_per_grid = function (grid_size,spatial_data,genus_name)
{
  pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
  box = round(extent(spatial_data) / pixelsize) * pixelsize
  template = raster(box, crs = spatial_data,
                    nrows = (box@ymax - box@ymin) / pixelsize, 
                    ncols = (box@xmax - box@xmin) / pixelsize)
  getisraster = rasterize(spatial_data, template, field = genus_name, fun = sum)
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

for (exp in c("OMNI12","OMNI13"))
{
  # exp="OMNI12"
  # exp="OMNI13"
  smpl_list=c()
  if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
  if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}
  
  for (smpl in smpl_list)
  {
    message(paste("--",exp,"-",smpl))
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
    # subset genes
    RNA_cutoff=100 # Only genes with more than 100 reads are accounted
    expression_data_with_xy_subset=expression_data_with_xy[,c("x","y",row.names(gene_sum)[gene_sum$count>=RNA_cutoff])]
    
    # Join the expression and all_Fungi all_Bacteria
    joint_Expression_Bacteria_Fungi_data=merge(expression_data_with_xy_subset,data_with_xy[,c("x","y","all_Fungi","all_Bacteria")],by=c("x","y"))
    
    #QA
    # write.table(file = "test_OMNI13_C2.data_used_for_hotspots_analyses.csv",joint_Expression_Bacteria_Fungi_data,sep=";",quote = F,row.names = F)
    
    # create the spatial object
    spatial_data=SpatialPointsDataFrame(coords = joint_Expression_Bacteria_Fungi_data[,c(1:2)],data = joint_Expression_Bacteria_Fungi_data[,3:NCOL(joint_Expression_Bacteria_Fungi_data)])
    
    # calculate the localG stas
    for (i in 1:length(names(spatial_data)))
    {
      layer=names(spatial_data)[i]
      if (i%%10==0) {
        message(paste("\t--",i,"out of",length(names(spatial_data))))
      }
      grid_size=2
      getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
      getisgrid.df=as.data.frame(getisgrid)
      names(getisgrid.df)=paste(names(getisgrid.df),layer,sep="_")
      getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
      names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
      getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
      getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
      row.names(getisgrid.df)=getisgrid.df$xy
      getisgrid.df=getisgrid.df[,-1]        
      
      if (nrow(hotspost_df)==0)
      {
        hotspost_df=getisgrid.df
      } else {
        hotspost_df=merge(hotspost_df,getisgrid.df[,1:4],by="row.names",all=T)
        row.names(hotspost_df)=hotspost_df$Row.names
        hotspost_df=hotspost_df[,-1]
      }
    }
    
    #### Save
    ## Once full with the layers 
    hotspots_and_grid_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,".",grid_size,"x",grid_size,".Getis_Ord_hotspots_and_layer_data.","RNA_cutoff_",RNA_cutoff,".joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv",sep="")
    write.table(x=hotspost_df,file = hotspots_and_grid_data_file,row.names = T,quote = F,sep = ";")
    # all layers file should not have changed. to validate that we copy it content before the rerun after 8Feb bug fix
    # mkdir /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/Before_Bug_Fix
    # ls -ltr /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.*.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350  142324434 Nov 24 17:38 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.A1.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350   78625829 Nov 24 18:03 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.A2.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350   36737446 Nov 24 18:13 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.B1.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350   44283602 Nov 24 18:26 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.B2.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350   61670776 Nov 24 18:39 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.C1.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350   34504504 Nov 24 18:47 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.C2.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # mv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.*.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/Before_Bug_Fix/
    
    
    # mkdir /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/Before_Bug_Fix
    # ls -ltr /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.*.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350  12707077 Nov 24 19:41 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.A1.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350  20173081 Nov 24 19:46 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.A2.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350  15302666 Nov 24 19:50 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.B1.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350   5302166 Nov 24 19:51 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.B2.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350  46667483 Nov 24 20:01 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.C1.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350   8830961 Nov 24 20:04 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.C2.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # -rw-r--r--  1 hashkenazy  350  49892398 Nov 24 20:16 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.D1.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv
    # mv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.*.2x2.Getis_Ord_hotspots_and_layer_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/Before_Bug_Fix/
    
    
    
    
    ## Save the hotspots df for correaltions purposes
    hotspots_no_layer=hotspost_df[,!grepl(pattern = "layer|HOTSPOT.p_|HOTSPOT.p.SP_FDR_",x = colnames(hotspost_df))]
    
    # The gene expresion layers localG values
    hotspots_expression=hotspots_no_layer[,!grepl(pattern = "all_Bacteria|all_Fungi",x=colnames(hotspots_no_layer))]
    # For the table: row.names=gene_id; values: expression on each spot
    hotspots_expression_table=as.data.frame(t(hotspots_expression))
    row.names(hotspots_expression_table)=gsub(pattern = "HOTSPOT_",replacement = "",fixed = T,x=row.names(hotspots_expression_table))
    hotspots_expression_table=hotspots_expression_table[!grepl(pattern = "X_Pcent|Y_Pcent|xy",x = row.names(hotspots_expression_table),ignore.case = T),]
    hotspots_expression_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,".",grid_size,"x",grid_size,".Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv",sep="")  
    write.table(x=hotspots_expression_table,file = hotspots_expression_data_file,row.names = T,quote = F,sep = "\t")
    
    # RNA layer should not have changhed we copy the data before rerun after bug fix 8Feb2022
    # mkdir /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/BeforeBugFix_8Feb2022
    
    # ls -ltr /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.*.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350  16652203 Nov 24 21:26 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.C1.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350  12091453 Nov 24 21:26 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.B2.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350   9944015 Nov 24 21:26 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.B1.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350  21279857 Nov 24 21:26 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.A2.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350  38908000 Nov 24 21:26 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.A1.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350   9454589 Nov 24 21:26 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.C2.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # mv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.*.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/BeforeBugFix_8Feb2022/
    # mv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI12.*.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv*.OK.tsv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI12/Filtered_CpMtRb/hotspots/BeforeBugFix_8Feb2022/  
    
    # mkdir /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/BeforeBugFix_8Feb2022
    # ls -ltr /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.*.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350  13507565 Nov 25 11:09 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.D1.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350   2372672 Nov 25 11:09 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.C2.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350  12561500 Nov 25 11:09 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.C1.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350   1389711 Nov 25 11:09 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.B2.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350   4187754 Nov 25 11:09 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.B1.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350   5542629 Nov 25 11:09 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.A2.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # -rw-r--r--  1 hashkenazy  350   3494065 Nov 25 11:09 /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.A1.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv
    # mv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.*.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/BeforeBugFix_8Feb2022/
    # mv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.OMNI13.*.2x2.Getis_Ord_hotspots_data.RNA_cutoff_100.joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.for_cor.tsv*OK.tsv /Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/OMNI13/Filtered_CpMtRb/hotspots/BeforeBugFix_8Feb2022/    
    
    
    # Once for the ITS data
    # For the table: row.names=NA, first col is tax; values: abundance on each spot
    hotspots_ITS=data.frame(row.names=row.names(hotspots_no_layer),"Getis_Ord_2x2_G"=hotspots_no_layer$HOTSPOT_all_Fungi)
    hotspots_ITS_table=as.data.frame(t(hotspots_ITS))
    hotspots_ITS_table$tax="Getis_Ord_2x2_G"
    hotspots_ITS_table=data.frame(tax=hotspots_ITS_table$tax,hotspots_ITS_table[,!grepl(x=names(hotspots_ITS_table),pattern = "tax")])
    out_file_ITS=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".All_ITS_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top50_genus.spatial_pos_UnderTissue_joint_RNA.",grid_size,"x",grid_size,".Getis_Ord_hotspots.csv",sep="")  
    write.table(file = out_file_ITS,x = hotspots_ITS_table,quote =  F,row.names = F,sep = ";")
    
    # Once for the Bacterial data ----> BUG!!!
    # hotspots_Bacteria=data.frame(row.names=row.names(hotspots_no_layer),"Getis_Ord_2x2_G"=hotspots_no_layer$HOTSPOT_all_Bacteria)
    # hotspots_Bacteria_table=as.data.frame(t(hotspots_ITS))
    # hotspots_Bacteria_table$tax="Getis_Ord_2x2_G"
    # hotspots_Bacteria_table=data.frame(tax=hotspots_ITS_table$tax,hotspots_ITS_table[,!grepl(x=names(hotspots_ITS_table),pattern = "tax")])
    # out_file_Bacteria=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top50_genus.spatial_pos_UnderTissue_joint_RNA.",grid_size,"x",grid_size,".Getis_Ord_hotspots.csv",sep="")  
    # write.table(file = out_file_Bacteria,x = hotspots_Bacteria_table,quote =  F,row.names = F,sep = ";")
    
    # Correct 8Feb2022
    hotspots_Bacteria=data.frame(row.names=row.names(hotspots_no_layer),"Getis_Ord_2x2_G"=hotspots_no_layer$HOTSPOT_all_Bacteria)
    hotspots_Bacteria_table=as.data.frame(t(hotspots_Bacteria))
    hotspots_Bacteria_table$tax="Getis_Ord_2x2_G"
    hotspots_Bacteria_table=data.frame(tax=hotspots_Bacteria_table$tax,hotspots_Bacteria_table[,!grepl(x=names(hotspots_Bacteria_table),pattern = "tax")])
    out_file_Bacteria=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".All_Bacterial_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.SUM_Top50_genus.spatial_pos_UnderTissue_joint_RNA.",grid_size,"x",grid_size,".Getis_Ord_hotspots.csv",sep="")  
    write.table(file = out_file_Bacteria,x = hotspots_Bacteria_table,quote =  F,row.names = F,sep = ";")
    
    rm("data","data_file","data_with_xy","expression_data",
       "expression_data_file","expression_data_with_xy","expression_data_with_xy_subset","gene_sum","gene_sum_ordered",
       "getisgrid","getisgrid.cooridinate_poligon_centers","getisgrid.df",
       "hotspost_df","hotspots_and_grid_data_file","hotspots_expression","hotspots_expression_table","hotspots_no_layer",
       "joint_Expression_Bacteria_Fungi_data","layer","spatial_data","xy","xy_list",
       hotspots_Bacteria,hotspots_Bacteria_table,hotspots_ITS,hotspots_ITS_table,
       out_file_Bacteria,out_file_ITS,i,hotspots_expression_data_file)
  }
}

######### CLEAN START ##############
######### OLD Figure 3 DATA ########
############################################################################

## Read 

#### FIG3 ---> STILL NEED TO MAKE SURE NONE OF THE ANALYSES ARE AFFECTED BY THE BUG CORRECTED ON 8Feb2022
#### ALL THE ANALYSES OF THIS SECTION WERE USED FOR FIGURE 3 OF THE MANUSCRIPT
#### Compare the FDR significant Hotspots for Bacterial+Fungi
experiments=c("OMNI12","OMNI13")

hotsposts_intersect_stats=data.frame()
grid_size=2
RNA_cutoff=100
pdf(file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Bacteria_to_Fungi_Ratio_in_Hotspots_type.log2.pdf",width = 8,height = 7)
for (exp in experiments) {
  smpl_list=c()
  if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
  if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}
  for (smpl in smpl_list)
  {
    message(paste("--",smpl))
    
    hotspots_and_grid_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,".",grid_size,"x",grid_size,".Getis_Ord_hotspots_and_layer_data.","RNA_cutoff_",RNA_cutoff,".joint_RNA_and_Sum_Top50_Fungi_Bacteria_pos.tsv",sep="")
    hotspots_and_grid_data=read.delim(file = hotspots_and_grid_data_file,sep = ";",stringsAsFactors = F)
    
    Bacteria_Fungi_data=hotspots_and_grid_data[,grepl(pattern = "Bacteria|Fungi",x = names(hotspots_and_grid_data),ignore.case = T)]
    Sum_Fungi_reads=sum(Bacteria_Fungi_data$layer_all_Fungi)
    Sum_Bacterial_reads=sum(Bacteria_Fungi_data$layer_all_Bacteria)
    
    significant_Bacteria=Bacteria_Fungi_data[Bacteria_Fungi_data$HOTSPOT.p.SP_FDR_all_Bacteria<0.05,]
    significant_Fungi=Bacteria_Fungi_data[Bacteria_Fungi_data$HOTSPOT.p.SP_FDR_all_Fungi<0.05,]
    
    Bacterial_Hotspots_count=sum(significant_Bacteria$HOTSPOT_all_Bacteria>0)
    Bacterial_Hotspots_grid_pixels_names=row.names(significant_Bacteria)[significant_Bacteria$HOTSPOT_all_Bacteria>0]
    
    Bacterial_Coldspots_count=sum(significant_Bacteria$HOTSPOT_all_Bacteria<0)
    
    Fungal_Hotspots_count=sum(significant_Fungi$HOTSPOT_all_Fungi>0)
    Fungal_Hotspots_grid_pixels_names=row.names(significant_Fungi)[significant_Fungi$HOTSPOT_all_Fungi>0]
    
    Fungal_Coldspots_count=sum(significant_Fungi$HOTSPOT_all_Fungi<0)
    
    
    Bacterial_Fungal_hotspots_intersect=length(intersect(row.names(significant_Bacteria[significant_Bacteria$HOTSPOT_all_Bacteria>0,]),row.names(significant_Fungi[significant_Fungi$HOTSPOT_all_Fungi>0,])))
    Bacterial_Fungal_coldspots_intersect=length(intersect(row.names(significant_Bacteria[significant_Bacteria$HOTSPOT_all_Bacteria<0,]),row.names(significant_Fungi[significant_Fungi$HOTSPOT_all_Fungi<0,])))
   
    
    # get the Grid-pixels of the Bacteria-Fungi-shared-significant Hotstpots
    
    Bacteria_Fungi_Hotspots_intersect_data=merge(significant_Bacteria[significant_Bacteria$HOTSPOT_all_Bacteria>0,],significant_Fungi[significant_Fungi$HOTSPOT_all_Fungi>0,],by="row.names")
    
    Bacteria_Hotspots_unique_grid_pixels_names=Bacterial_Hotspots_grid_pixels_names[!(Bacterial_Hotspots_grid_pixels_names %in% Fungal_Hotspots_grid_pixels_names)]
    Bacteria_Hotspots_unique_grid_data=significant_Bacteria[Bacteria_Hotspots_unique_grid_pixels_names,]
    
    Fungal_Hotspots_unique_grid_pixels_names=Fungal_Hotspots_grid_pixels_names[!(Fungal_Hotspots_grid_pixels_names %in% Bacterial_Hotspots_grid_pixels_names)]
    Fungal_Hotspots_unique_grid_data=significant_Fungi[Fungal_Hotspots_unique_grid_pixels_names,]
    
    
    Bacterial_Fungal_ratios=data.frame()
    Bacterial_Fungal_ratios=rbind(Bacterial_Fungal_ratios,data.frame(Bacterial_Fungi_ratio=Bacteria_Fungi_Hotspots_intersect_data$layer_all_Bacteria.x/Bacteria_Fungi_Hotspots_intersect_data$layer_all_Fungi.x,type="Shared"))
    Bacterial_Fungal_ratios=rbind(Bacterial_Fungal_ratios,data.frame(Bacterial_Fungi_ratio=Fungal_Hotspots_unique_grid_data$layer_all_Bacteria/Fungal_Hotspots_unique_grid_data$layer_all_Fungi,type="Fungi_unique"))
    Bacterial_Fungal_ratios=rbind(Bacterial_Fungal_ratios,data.frame(Bacterial_Fungi_ratio=Bacteria_Hotspots_unique_grid_data$layer_all_Bacteria/Bacteria_Hotspots_unique_grid_data$layer_all_Fungi,type="Bacteria_unique"))
    
    Bacterial_Fungal_by_hotspot_type=data.frame()
    Bacterial_Fungal_by_hotspot_type=rbind(Bacterial_Fungal_by_hotspot_type,data.frame(Bacterial_count=Bacteria_Fungi_Hotspots_intersect_data$layer_all_Bacteria.x,Fungal_count=Bacteria_Fungi_Hotspots_intersect_data$layer_all_Fungi.x,type="Shared"))
    Bacterial_Fungal_by_hotspot_type=rbind(Bacterial_Fungal_by_hotspot_type,data.frame(Bacterial_count=Fungal_Hotspots_unique_grid_data$layer_all_Bacteria,Fungal_count=Fungal_Hotspots_unique_grid_data$layer_all_Fungi,type="Fungi_unique"))
    Bacterial_Fungal_by_hotspot_type=rbind(Bacterial_Fungal_by_hotspot_type,data.frame(Bacterial_count=Bacteria_Hotspots_unique_grid_data$layer_all_Bacteria,Fungal_count=Bacteria_Hotspots_unique_grid_data$layer_all_Fungi,type="Bacteria_unique"))
    Bacterial_Fungal_by_hotspot_type$Bacterial_log2=log(Bacterial_Fungal_by_hotspot_type$Bacterial_count,2)
    Bacterial_Fungal_by_hotspot_type$Fungal_log2=log(Bacterial_Fungal_by_hotspot_type$Fungal_count,2)

    
    points_count=data.frame(type=c("Shared","Fungi_unique","Bacteria_unique"),
                            num_of_hotspots=c(NROW(Bacteria_Fungi_Hotspots_intersect_data),NROW(Fungal_Hotspots_unique_grid_data),NROW(Bacteria_Hotspots_unique_grid_data)),
                            ypos=c(max(Bacteria_Fungi_Hotspots_intersect_data$layer_all_Bacteria.x/Bacteria_Fungi_Hotspots_intersect_data$layer_all_Fungi.x)+0.5,
                                   max(Fungal_Hotspots_unique_grid_data$layer_all_Bacteria/Fungal_Hotspots_unique_grid_data$layer_all_Fungi)+0.5,
                                   max(Bacteria_Hotspots_unique_grid_data$layer_all_Bacteria/Bacteria_Hotspots_unique_grid_data$layer_all_Fungi)+0.5)
                            )

    p=ggplot(data = Bacterial_Fungal_ratios, aes(x=type, y=Bacterial_Fungi_ratio, fill=type)) +
      geom_boxplot() +
      scale_fill_manual(values=c("#66c2a5","#8da0cb","#fc8d62")) +
      geom_jitter(color="black", size=0.4, alpha=0.9) +
      theme_bw() +
      theme(
        legend.position="none",
        plot.title = element_text(size=11)
      ) +
      labs(x="Hotspot type",y="Bacteria to Fungi ratio",title=paste(exp,smpl,sep=".")) +
      geom_text(data = points_count, aes(label = paste("n=",num_of_hotspots," (",signif(num_of_hotspots/sum(num_of_hotspots),2)*100,"%)",sep=""), y = max(ypos)), 
                position = position_dodge(width = .75), 
                show.legend = FALSE ) +
      scale_x_discrete(labels=c("Bacteria-unique","Fungi-unique","Shared Bacteria-Fungi"))
    print(p)
    
    # as scatter
    p1=ggplot(Bacterial_Fungal_by_hotspot_type, aes(x=Fungal_log2, y=Bacterial_log2, color=type)) +
      geom_point() + 
      geom_smooth(method=glm, aes(fill=type,color=type)) +
      theme_bw()+
      theme(
        plot.title = element_text(size=11)
      ) +
      labs(x="log2 Fungal count",y="log2 Bacterial count",title=paste(exp,smpl,sep=".")) +
      scale_color_manual(values=c("#66c2a5","#8da0cb","#fc8d62"),labels=c("Bacteria-unique","Fungi-unique","Shared Bacteria-Fungi"),name="Hotspot type")+
      scale_fill_manual(values=c("#66c2a5","#8da0cb","#fc8d62"),labels=c("Bacteria-unique","Fungi-unique","Shared Bacteria-Fungi"),name="Hotspot type") 
    print(p1)  
    
    # Do e have points with Bacteria only data
    #summary()
    #wilcox.test(x = Bacteria_Hotspots_unique_grid_data$layer_all_Bacteria/Bacteria_Hotspots_unique_grid_data$layer_all_Fungi,y=Bacteria_Fungi_Hotspots_intersect_data$layer_all_Bacteria.x/Bacteria_Fungi_Hotspots_intersect_data$layer_all_Fungi.x)
    #wilcox.test(x = Fungal_Hotspots_unique_grid_data$layer_all_Bacteria/Fungal_Hotspots_unique_grid_data$layer_all_Fungi,y=Bacteria_Fungi_Hotspots_intersect_data$layer_all_Bacteria.x/Bacteria_Fungi_Hotspots_intersect_data$layer_all_Fungi.x)
    
    # get the Grid-pixels of the Fungi-unique-significant Hotstpots
    
    # get the Grid-pixels of the Bacteria-unique-significant Hotstpots
    
    record=data.frame(sample=paste(exp,smpl,sep="."),
                      pixels=nrow(hotspots_and_grid_data),
                      Sum_Bacterial_reads=Sum_Bacterial_reads,
                      Sum_Fungi_reads=Sum_Fungi_reads,
                      Significant_Bacteria_count=nrow(significant_Bacteria),
                      Bacterial_Hotspots_count=Bacterial_Hotspots_count,
                      Bacterial_Coldspots_count=Bacterial_Coldspots_count,
                      
                      Significant_Fungal_count=nrow(significant_Fungi),
                      Fungal_Hotspots_count=Fungal_Hotspots_count,
                      Fungal_Coldspots_count=Fungal_Coldspots_count,
                      
                      Bacterial_Fungal_hotspots_intersect=Bacterial_Fungal_hotspots_intersect,
                      Bacterial_Fungal_coldspots_intersect=Bacterial_Fungal_coldspots_intersect,
                      Bacterial_Fungal_hotspots_intersect_frac=Bacterial_Fungal_hotspots_intersect/(Fungal_Hotspots_count+Bacterial_Hotspots_count-Bacterial_Fungal_hotspots_intersect),
                      Bacterial_Fungal_coldspots_intersect_frac=Bacterial_Fungal_coldspots_intersect/(Fungal_Coldspots_count+Bacterial_Coldspots_count-Bacterial_Fungal_coldspots_intersect),
                      
                      # unique bacteria fraction
                      Bacterial_hotspots_unique_frac=(Bacterial_Hotspots_count-Bacterial_Fungal_hotspots_intersect)/(Fungal_Hotspots_count+Bacterial_Hotspots_count-Bacterial_Fungal_hotspots_intersect),
                      Bacterial_coldspots_unique_frac=(Bacterial_Coldspots_count-Bacterial_Fungal_coldspots_intersect)/(Fungal_Coldspots_count+Bacterial_Coldspots_count-Bacterial_Fungal_coldspots_intersect),
                      
                      # unique fungi fraction
                      Fungal_hotspots_unique_frac=(Fungal_Hotspots_count-Bacterial_Fungal_hotspots_intersect)/(Fungal_Hotspots_count+Bacterial_Hotspots_count-Bacterial_Fungal_hotspots_intersect),
                      Fungal_coldspots_unique_frac=(Fungal_Coldspots_count-Bacterial_Fungal_coldspots_intersect)/(Fungal_Coldspots_count+Bacterial_Coldspots_count-Bacterial_Fungal_coldspots_intersect),
                      
                      file=hotspots_and_grid_data_file)
    hotsposts_intersect_stats=rbind(hotsposts_intersect_stats,record)
    rm ("Bacteria_Fungi_data","Bacteria_Fungi_Hotspots_intersect_data","Bacteria_Hotspots_unique_grid_data","Bacteria_Hotspots_unique_grid_pixels_names",
        "Bacterial_Coldspots_count","Bacterial_Fungal_by_hotspot_type","Bacterial_Fungal_coldspots_intersect","Bacterial_Fungal_hotspots_intersect",
        "Bacterial_Fungal_ratios","Bacterial_Hotspots_count","Bacterial_Hotspots_grid_pixels_names",
        "Fungal_Coldspots_count","Fungal_Hotspots_count","Fungal_Hotspots_grid_pixels_names","Fungal_Hotspots_unique_grid_data",
        "Fungal_Hotspots_unique_grid_pixels_names","hotspots_and_grid_data","hotspots_and_grid_data_file","p","p1","points_count",                              
        "record","significant_Bacteria","significant_Fungi","Sum_Bacterial_reads","Sum_Fungi_reads")
  }
}
dev.off()
# percent of the array
hotsposts_intersect_stats$Bacterial_Coldspots_percent=hotsposts_intersect_stats$Bacterial_Coldspots_count/hotsposts_intersect_stats$pixels*100
hotsposts_intersect_stats$Bacterial_Hotspots_percent=hotsposts_intersect_stats$Bacterial_Hotspots_count/hotsposts_intersect_stats$pixels*100

hotsposts_intersect_stats$Fungal_Coldspots_percent=hotsposts_intersect_stats$Fungal_Coldspots_count/hotsposts_intersect_stats$pixels*100
hotsposts_intersect_stats$Fungal_Hotspots_percent=hotsposts_intersect_stats$Fungal_Hotspots_count/hotsposts_intersect_stats$pixels*100




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

hotsposts_intersect_stats=merge(hotsposts_intersect_stats,samples_label,by = "sample")


# save
write.table(x=hotsposts_intersect_stats,file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/data/Compare_Hotspots_OMNI12_OMNI13.csv",sep=";",quote = F,row.names = F)

library(reshape)
mdata <- melt(hotsposts_intersect_stats, id=c("sample"))

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
  

mdata_intersect=mdata[grepl(x = mdata$variable,pattern = "_intersect_frac"),]
mdata_intersect=merge(mdata_intersect,samples_label)
mdata_intersect$value=as.numeric(mdata_intersect$value)
mdata_intersect$value=signif(x = mdata_intersect$value*100,3)
library(ggplot2)
library(RColorBrewer)
p=ggplot(mdata_intersect, aes(fill=variable, y=value, x=lable)) + 
  geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c("#ca0020", "#0571b0"),labels=c("Hotspots","Coldspots"),name="") + 
  labs(y="Shared (%)",x="")+
  scale_y_continuous(breaks=seq(0,70,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Compare_Significant_Hotspots_OMNI12_OMNI13.pdf")
print (p)
dev.off()

# stacked barplot with the % of shared, unique bacteria unique Fungi
mdata_intersect_unique=mdata[grepl(x = mdata$variable,pattern = "_frac"),]
mdata_intersect_unique=merge(mdata_intersect_unique,samples_label)
mdata_intersect_unique$value=as.numeric(mdata_intersect_unique$value)
mdata_intersect_unique$value=signif(x = mdata_intersect_unique$value*100,3)

# separate into hotspots and coldspots
mdata_intersect_unique_hotspots=mdata_intersect_unique[grepl(x = mdata_intersect_unique$variable,pattern = "_hotspots_"),]
mdata_intersect_unique_coldspots=mdata_intersect_unique[grepl(x = mdata_intersect_unique$variable,pattern = "_coldspots_"),]

p=ggplot(mdata_intersect_unique_hotspots, aes(fill=variable, y=value, x=lable,label=paste(value,"%",sep=""))) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  
  #scale_fill_manual(labels=c("Hotspots","Coldspots"),name="") + 
  scale_fill_manual(values = c("#fc8d62","#66c2a5","#8da0cb"),labels=c("Shared Bacterial-Fungal hotspots", "Bacterial-unique","Fungal-unique"),name="") +  # 
  
#  scale_fill_discrete("", 
#                      labels=c("Shared Bacterial-Fungal hotspots", "Bacterial-unique","Fungal-unique")) +
  labs(y="Hotspots (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Percent_of_Significant_shared_and_unique_Hotspots_OMNI12_OMNI13.pdf",height = 7,width = 10)
print (p)
dev.off()

mdata_intersect_unique_hotspots_and_coldpots=rbind(data.frame(mdata_intersect_unique_hotspots,type="hotspot"),data.frame(mdata_intersect_unique_coldspots,type="coldspot"))
p=ggplot(mdata_intersect_unique_hotspots_and_coldpots, aes(fill=variable, y=value, x=lable)) + 
  # geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+
  #scale_fill_manual(values = c("#ca0020", "#0571b0"),labels=c("Hotspots","Coldspots"),name="") + 
  labs(y="Shared (%)",x="")+
  scale_y_continuous(breaks=seq(0,100,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) +
  facet_grid( ~ type)
p

# create the flat df for hotspots of differnt types
mdatat_hotspots_percent=data.frame()
for (i in 1:nrow(hotsposts_intersect_stats))
{
  # each row will become 4 in the melted df
  record_fun=data.frame(sample=rep(hotsposts_intersect_stats$sample[i],2),
                    tax=rep("Fungi",2),
                    type=c("Hotspot","Coldspot"),
                    value=c(hotsposts_intersect_stats$Fungal_Hotspots_percent[i],hotsposts_intersect_stats$Fungal_Coldspots_percent[i]))
  record_bac=data.frame(sample=rep(hotsposts_intersect_stats$sample[i],2),
                        tax=rep("Bacteria",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Bacterial_Hotspots_percent[i],hotsposts_intersect_stats$Bacterial_Coldspots_percent[i]))
  mdatat_hotspots_percent=rbind(mdatat_hotspots_percent,record_fun,record_bac)
}
mdatat_hotspots_percent=merge(mdatat_hotspots_percent,samples_label)

p=ggplot(mdatat_hotspots_percent, aes(fill=type, y=value, x=tax)) + 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c("#0571b0","#ca0020" ),name="") +  # labels=c("Bacteria","Fungi")
  labs(y="% of array",x="")+
  scale_y_continuous(breaks=seq(0,40,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_grid( ~ lable)

pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Percent_of_Significant_Hotspots_OMNI12_OMNI13.pdf",height = 7,width = 10)
print (p)
dev.off()


# total number of reads of each type
mdata_sum_reads=mdata[grepl(x = mdata$variable,pattern = "Sum_"),]
mdata_sum_reads=merge(mdata_sum_reads,samples_label)
mdata_sum_reads$value=as.numeric(mdata_sum_reads$value)
library(ggplot2)
library(RColorBrewer)
library(scales)
p=ggplot(mdata_sum_reads, aes(fill=variable, y=value, x=lable)) + 
  geom_bar(position="dodge", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c( "#018571","#a6611a"),name="",labels=c("Total Bacterial reads","Total Fungi reads")) + #
  labs(y="Number of reads",x="") +
  scale_y_continuous(label=comma) +
  # scale_y_continuous(breaks=seq(0,70,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/NumberOfReads_OMNI12_OMNI13.pdf")
print (p)
dev.off()



# interkigdom as function of shared hotspots
# based on "~/Dropbox/PostDoc/Projects/16S_array/compare_ITS_16S_networks.R"
k1_k2_counts=data.frame()
for (Exp in c("OMNI12","OMNI13"))
{
  samples=c("A1","A2","B1","B2","C1","C2")
  if (Exp=="OMNI13") {samples=c(samples,"D1")}
  for (smpl in samples) 
  {
    net_name=paste(Exp,smpl,sep=".")
    message(paste("--",net_name))
    spatial_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/Microbial_Networks/joint_Bacteria_Fungi_and_UNKNOWN_transposed/permatfull_col_shuffle/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.spearman_corr_and_permatfull_col_shuffle.1000.empP.csv",sep="")
    # hotspots_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.AND.210708_raw_counts_mtrbcpfiltered_RNA_cutoff_",RNA_cutoff,".",grid_size,"x",grid_size,".Getis_Ord_hotspots",".spearman_network.microbial_only.csv",sep="")
    # spatial_and_hotspots_significant_net_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",Exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",Exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.ABUNDANCE_and_HOTSPOTS_significant_agreement_network.csv",sep="")
    spatial_net_raw=read.delim(file=spatial_net_file,sep=";",stringsAsFactors = F)
    # hotspots_net_raw=read.delim(file = hotspots_net_file, sep=";",stringsAsFactors = F)
    # hotspots_net_raw=unique(hotspots_net_raw)
    # hotspots_net_raw$s1_s2=paste(hotspots_net_raw$row,hotspots_net_raw$column,sep="_")
    spatial_net_signif=spatial_net_raw[spatial_net_raw$p.BH<0.05&spatial_net_raw$empP.BH<0.05&!is.na(spatial_net_raw$empP.BH),]
    # hotspots_net_signif=hotspots_net_raw[hotspots_net_raw$p.BH<0.05&!is.na(hotspots_net_raw$p.BH),]
    spatial_net_signif$k1_k2=paste(spatial_net_signif$s1_domain,spatial_net_signif$s2_domain,sep="_")
    count_t=table(spatial_net_signif$k1_k2)
    count_df=data.frame(sample=net_name,Bacteria_Bacteria=count_t[["Bacteria_Bacteria"]],Bacteria_Fungi=count_t[["Bacteria_Fungi"]],Fungi_Fungi=count_t[["Fungi_Fungi"]],net_file=spatial_net_file)
    k1_k2_counts=rbind(k1_k2_counts,count_df)
  }
}

k1_k2_counts$Bacteria_Bacteria_frac=k1_k2_counts$Bacteria_Bacteria/rowSums(k1_k2_counts[,c(2:4)])
k1_k2_counts$Bacteria_Fungi_frac=k1_k2_counts$Bacteria_Fungi/rowSums(k1_k2_counts[,c(2:4)])
k1_k2_counts$Fungi_Fungi_frac=k1_k2_counts$Fungi_Fungi/rowSums(k1_k2_counts[,c(2:4)])

# plot the ratio of interaction type per section
# create the flat df for hotspots of differnt types
library(reshape)
m_k1_k2_counts_all <- melt(k1_k2_counts, id=c("sample"))

m_k1_k2_frac=m_k1_k2_counts_all[grepl(x=m_k1_k2_counts_all$variable,pattern = "_frac"),]
m_k1_k2_count=m_k1_k2_counts_all[!grepl(x=m_k1_k2_counts_all$variable,pattern = "_frac|net"),]

m_k1_k2_count_and_frac=m_k1_k2_frac
m_k1_k2_count_and_frac$variable=gsub(m_k1_k2_count_and_frac$variable,pattern = "_frac",replacement = "")
names(m_k1_k2_count_and_frac)[3]="frac"

m_k1_k2_count_and_frac=merge(m_k1_k2_count_and_frac,m_k1_k2_count,by=c("sample","variable"))
m_k1_k2_count_and_frac$frac=as.numeric(m_k1_k2_count_and_frac$frac)
m_k1_k2_count_and_frac$value=as.numeric(m_k1_k2_count_and_frac$value)

m_k1_k2_count_and_frac=merge(m_k1_k2_count_and_frac,samples_label,by="sample")
# With total count
p=ggplot(m_k1_k2_count_and_frac, aes(fill=variable, y=frac*100, x=lable,label=paste(value,"\n",signif(x = frac*100,digits = 2),"%",sep=""))) + 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+ 
  geom_text(size = 3, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("#66c2a5","#fc8d62","#8da0cb"),labels=c("Bacteria-Bacteria","Bacteria-Fungi","Fungi-Fungi"),name="") +  # 
  labs(y="Microbial interactions",x="")+
  scale_y_continuous(breaks=seq(0,100,10)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1))  
  # facet_grid( ~ lable)
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Significant_reads_based_microbial_interactions_by_type_OMNI12_OMNI13.pdf",height = 7,width = 10)
print (p)
dev.off()


# with percent
mdatat_hotspots_percent=data.frame()
for (i in 1:nrow(k1_k2_counts))
{
  # each row will become 4 in the melted df
  record_fun=data.frame(sample=rep(k1_k2_counts$sample[i],2),
                        tax=rep("Fungi",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Fungal_Hotspots_percent[i],hotsposts_intersect_stats$Fungal_Coldspots_percent[i]))
  record_bac=data.frame(sample=rep(hotsposts_intersect_stats$sample[i],2),
                        tax=rep("Bacteria",2),
                        type=c("Hotspot","Coldspot"),
                        value=c(hotsposts_intersect_stats$Bacterial_Hotspots_percent[i],hotsposts_intersect_stats$Bacterial_Coldspots_percent[i]))
  mdatat_hotspots_percent=rbind(mdatat_hotspots_percent,record_fun,record_bac)
}
mdatat_hotspots_percent=merge(mdatat_hotspots_percent,samples_label)

p=ggplot(mdatat_hotspots_percent, aes(fill=type, y=value, x=tax)) + 
  geom_bar(position="stack", stat="identity",colour = "black",width = 0.85)+ 
  scale_fill_manual(values = c("#0571b0","#ca0020" ),name="") +  # labels=c("Bacteria","Fungi")
  labs(y="% of array",x="")+
  scale_y_continuous(breaks=seq(0,40,5)) +
  theme_bw() + theme(axis.text.x = element_text(angle = 45,hjust=1)) + 
  facet_grid( ~ lable)

pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/Percent_of_Significant_Hotspots_OMNI12_OMNI13.pdf",height = 7,width = 10)
print (p)
dev.off()



hotsposts_intersect_stats_and_microbial_net_stats=merge(hotsposts_intersect_stats,k1_k2_counts,by="sample")

# GET EQUATION AND R-SQUARED AS STRING
# SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
# https://stackoverflow.com/questions/7549694/add-regression-line-equation-and-r2-on-graph
lm_eqn <- function(df,x,y){
  m <- lm(y ~ x, df);
  eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                   list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), 
                        pvalue = format(summary(m)$coefficients[2,4], digits = 2)))
  as.character(as.expression(eq));
}

# lm model for the Bacteria-Fungi interaction
m <- lm(Bacterial_Fungal_hotspots_intersect_frac ~ Bacteria_Fungi_frac, hotsposts_intersect_stats_and_microbial_net_stats);
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), 
                      pvalue = format(summary(m)$coefficients[2,4], digits = 2)))
hotsposts_intersect_stats_and_microbial_net_stats$sum_all_reads=hotsposts_intersect_stats_and_microbial_net_stats$Sum_Bacterial_reads+hotsposts_intersect_stats_and_microbial_net_stats$Sum_Fungi_reads
summary(hotsposts_intersect_stats_and_microbial_net_stats$sum_all_read)
p=ggplot(hotsposts_intersect_stats_and_microbial_net_stats, aes(x=Bacteria_Fungi_frac*100, y=Bacterial_Fungal_hotspots_intersect_frac*100)) + 
  geom_point(aes(size = sum_all_reads, colour = sum_all_reads))+ # colour = sum_all_reads,
  scale_size_binned() +
  # scale_fill_gradient(low = "grey", high = "red") +
  geom_text(aes(label=paste(lable," (",format(sum_all_reads,big.mark = ",",trim = TRUE),")",sep="")),hjust=0, vjust=0) +
  labs(y="Shared Bacterial-Fungal hotspots (%)",x="Bacteria-Fungi interactions (%)",size="Number of reads",color="Number of reads") +
  geom_smooth(method=lm) + 
  geom_text(x = 31, y = 64, 
            label = as.character(as.expression(eq)), parse = TRUE) +
  theme_bw()
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI12_OMNI13_shared_Hotspots_vs_Bacteria_Fungi_interactions.pdf")
print (p)
dev.off()
# lm model
m <- lm(Bacterial_hotspots_unique_frac ~ Bacteria_Bacteria_frac, hotsposts_intersect_stats_and_microbial_net_stats);
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), 
                      pvalue = format(summary(m)$coefficients[2,4], digits = 2)))

p=ggplot(hotsposts_intersect_stats_and_microbial_net_stats, aes(x=Bacteria_Bacteria_frac*100, y=Bacterial_hotspots_unique_frac*100)) + 
  geom_point(aes(size = Sum_Bacterial_reads, colour = Sum_Bacterial_reads))+ # colour = sum_all_reads,
  scale_size_binned() +
  # scale_fill_gradient(low = "grey", high = "red") +
  geom_text(aes(label=paste(lable," (",format(Sum_Bacterial_reads,big.mark = ",",trim = TRUE),")",sep="")),hjust=0, vjust=0) +
  labs(y="Bacterial-unique hotspots (%)",x="Bacteria-Bacteria interactions (%)",size="Number of Bacterial reads",color="Number of Bacterial reads") +
  geom_smooth(method=lm) + 
  geom_text(x = 42, y = 78, 
            label = as.character(as.expression(eq)), parse = TRUE) +
  theme_bw()
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI12_OMNI13_Bacterial_unique_Hotspots_vs_Bacteria_Bacteria_interactions.pdf")
print (p)
dev.off()
# lm model
m <- lm(Fungal_hotspots_unique_frac ~ Fungi_Fungi_frac, hotsposts_intersect_stats_and_microbial_net_stats);
eq <- substitute(italic(r)~"="~rvalue*","~italic(p)~"="~pvalue, 
                 list(rvalue = sprintf("%.2f",sign(coef(m)[2])*sqrt(summary(m)$r.squared)), 
                      pvalue = format(summary(m)$coefficients[2,4], digits = 2)))

p=ggplot(hotsposts_intersect_stats_and_microbial_net_stats, aes(x=Fungi_Fungi_frac*100, y=Fungal_hotspots_unique_frac*100)) + 
  geom_point(aes(size = Sum_Fungi_reads, colour = Sum_Fungi_reads))+ # colour = sum_all_reads,
  scale_size_binned() +
  # scale_fill_gradient(low = "grey", high = "red") +
  geom_text(aes(label=paste(lable," (",format(Sum_Fungi_reads,big.mark = ",",trim = TRUE),")",sep="")),hjust=0, vjust=0) +
  
  labs(y="Fungal-unique hotspots (%)",x="Fungi-Fungi interactions (%)",size="Number of Fungal reads",color="Number of Fungal reads") +
  geom_smooth(method=lm) + 
  geom_text(x = 10, y = 36.5, 
            label = as.character(as.expression(eq)), parse = TRUE) +
  theme_bw()
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI12_OMNI13_Fungi_unique_Hotspots_vs_Fungi_Fungi_interactions.pdf")
print (p)
dev.off()

cor.test(hotsposts_intersect_stats_and_microbial_net_stats$Fungi_Fungi_frac,hotsposts_intersect_stats_and_microbial_net_stats$Fungal_hotspots_unique_frac)

# correlation between the number of reads and the proprtion of significan hot/coldspots
cor.test(hotsposts_intersect_stats_and_microbial_net_stats$Sum_Bacterial_reads,hotsposts_intersect_stats_and_microbial_net_stats$Bacterial_Hotspots_count,method = "spearman")
cor.test(hotsposts_intersect_stats_and_microbial_net_stats$Sum_Fungi_reads,hotsposts_intersect_stats_and_microbial_net_stats$Fungal_Hotspots_count,method = "spearman")

# Venn Euler diagram For OMNI12 A2
library(eulerr)
# hotsposts_intersect_stats[hotsposts_intersect_stats$sample=="OMNI12.A2",]
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI12.A2.significant_hot_cold_spots.pdf")

OMNI12_A2_intersect_hotspots=c("Fungi" = hotsposts_intersect_stats$Fungal_Hotspots_count[hotsposts_intersect_stats$sample=="OMNI12.A2"]-hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A2"],
                               "Bacteria" = hotsposts_intersect_stats$Bacterial_Hotspots_count[hotsposts_intersect_stats$sample=="OMNI12.A2"]-hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A2"],
                               "Fungi&Bacteria"=hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A2"])
OMNI12_A2_intersect_hotspots_fit <- euler(OMNI12_A2_intersect_hotspots)
plot(OMNI12_A2_intersect_hotspots_fit,
     quantities = list(type = c("counts", "percent"), font=8, round=2, cex=1) ,
     main=paste("Hotspots P1.L1.2",sep=""),
     fills = c("#8da0cb","#66c2a5","#fc8d62")) # 

OMNI12_A2_intersect_coldspots=c("Fungi" = hotsposts_intersect_stats$Fungal_Coldspots_count[hotsposts_intersect_stats$sample=="OMNI12.A2"]-hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A2"],
                               "Bacteria" = hotsposts_intersect_stats$Bacterial_Coldspots_count[hotsposts_intersect_stats$sample=="OMNI12.A2"]-hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A2"],
                               "Fungi&Bacteria"=hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A2"])
OMNI12_A2_intersect_coldspots_fit <- euler(OMNI12_A2_intersect_coldspots)
plot(OMNI12_A2_intersect_coldspots_fit,
     quantities = list(type = c("counts", "percent"), font=8, round=2, cex=1) ,
     main=paste("Coldspots P1.L1.2",sep=""),
     fills = c("#8da0cb","#66c2a5","#fc8d62"))

dev.off()

# Venn Euler diagram For OMNI12 A1
library(eulerr)
# hotsposts_intersect_stats[hotsposts_intersect_stats$sample=="OMNI12.A1",]
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI12.A1.significant_hot_cold_spots.pdf")

OMNI12_A1_intersect_hotspots=c("Fungi" = hotsposts_intersect_stats$Fungal_Hotspots_count[hotsposts_intersect_stats$sample=="OMNI12.A1"]-hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A1"],
                               "Bacteria" = hotsposts_intersect_stats$Bacterial_Hotspots_count[hotsposts_intersect_stats$sample=="OMNI12.A1"]-hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A1"],
                               "Fungi&Bacteria"=hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A1"])
OMNI12_A1_intersect_hotspots_fit <- euler(OMNI12_A1_intersect_hotspots)
plot(OMNI12_A1_intersect_hotspots_fit,
     quantities = list(type = c("counts", "percent"), font=8, round=2, cex=1) ,
     main=paste("Hotspots P1.L1.1",sep=""),
     fills = c("#8da0cb","#66c2a5","#fc8d62")) # 

OMNI12_A1_intersect_coldspots=c("Fungi" = hotsposts_intersect_stats$Fungal_Coldspots_count[hotsposts_intersect_stats$sample=="OMNI12.A1"]-hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A1"],
                                "Bacteria" = hotsposts_intersect_stats$Bacterial_Coldspots_count[hotsposts_intersect_stats$sample=="OMNI12.A1"]-hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A1"],
                                "Fungi&Bacteria"=hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI12.A1"])
OMNI12_A1_intersect_coldspots_fit <- euler(OMNI12_A1_intersect_coldspots)
plot(OMNI12_A1_intersect_coldspots_fit,
     quantities = list(type = c("counts", "percent"), font=8, round=2, cex=1) ,
     main=paste("Coldspots P1.L1.1",sep=""),
     fills = c("#8da0cb","#66c2a5","#fc8d62"))

dev.off()

# Venn Euler diagram For OMNI13 B1
library(eulerr)
# hotsposts_intersect_stats[hotsposts_intersect_stats$sample=="OMNI12.A1",]
pdf (file="/Users/hashkenazy/Dropbox/PostDoc/Projects/16S_array/Figures_for_paper/OMNI13.B1.significant_hot_cold_spots.pdf")

OMNI13_B1_intersect_hotspots=c("Fungi" = hotsposts_intersect_stats$Fungal_Hotspots_count[hotsposts_intersect_stats$sample=="OMNI13.B1"]-hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI13.B1"],
                               "Bacteria" = hotsposts_intersect_stats$Bacterial_Hotspots_count[hotsposts_intersect_stats$sample=="OMNI13.B1"]-hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI13.B1"],
                               "Fungi&Bacteria"=hotsposts_intersect_stats$Bacterial_Fungal_hotspots_intersect[hotsposts_intersect_stats$sample=="OMNI13.B1"])
OMNI13_B1_intersect_hotspots_fit <- euler(OMNI13_B1_intersect_hotspots)
plot(OMNI13_B1_intersect_hotspots_fit,
     quantities = list(type = c("counts", "percent"), font=8, round=2, cex=1) ,
     main=paste("Hotspots P2.L1.4",sep=""),
     fills = c("#8da0cb","#66c2a5","#fc8d62")) # 

OMNI13_B1_intersect_coldspots=c("Fungi" = hotsposts_intersect_stats$Fungal_Coldspots_count[hotsposts_intersect_stats$sample=="OMNI13.B1"]-hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI13.B1"],
                                "Bacteria" = hotsposts_intersect_stats$Bacterial_Coldspots_count[hotsposts_intersect_stats$sample=="OMNI13.B1"]-hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI13.B1"],
                                "Fungi&Bacteria"=hotsposts_intersect_stats$Bacterial_Fungal_coldspots_intersect[hotsposts_intersect_stats$sample=="OMNI13.B1"])
OMNI13_B1_intersect_coldspots_fit <- euler(OMNI13_B1_intersect_coldspots)
plot(OMNI13_B1_intersect_coldspots_fit,
     quantities = list(type = c("counts", "percent"), font=8, round=2, cex=1) ,
     main=paste("Coldspots P2.L1.4",sep=""),
     fills = c("#8da0cb","#66c2a5","#fc8d62"))

dev.off()

##### END CLEAN START HOTSPOTS FOR joined Bacteria+Fungi+Expression ####

##### correlate microbial hotspots and RNA hotspots ####
# For the all_Bacterial all_Fungi
exp="OMNI12"
exp="OMNI13"
grid_size=2
RNA_cutoff=100
smpl_list=c()
if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}

for (smpl in smpl_list)
{
  message(paste("--",smpl))
  microbial_hotspots_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.",grid_size,"x",grid_size,".Getis_Ord_hotspots",".tsv",sep="")
  rna_hotspots_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/RNA/spatial_tables/",exp,"/Filtered_CpMtRb/hotspots/210708_raw_counts_mtrbcpfiltered.",exp,".",smpl,grid_size,"x",grid_size,".Getis_Ord_hotspots_data.","RNA_cutoff_",RNA_cutoff,".tsv",sep="")

  microbial_hotspots=read.delim(file = microbial_hotspots_file, stringsAsFactors = F, sep=";")
  microbial_hotspots=microbial_hotspots[,c(1,5:ncol(microbial_hotspots))]
  colnames(microbial_hotspots)[1:50]=paste("BACTERIA_",colnames(microbial_hotspots)[1:50],sep="")
  colnames(microbial_hotspots)[101]=paste("BACTERIA_",colnames(microbial_hotspots)[101],sep="")
  colnames(microbial_hotspots)[51:100]=paste("FUNGI_",colnames(microbial_hotspots)[51:100],sep="")
  colnames(microbial_hotspots)[102]=paste("FUNGI_",colnames(microbial_hotspots)[102],sep="")
  
  rna_hotspots=read.delim(file = rna_hotspots_file, stringsAsFactors = F, sep=";")
  rna_hotspots=rna_hotspots[,c(1,5:ncol(rna_hotspots))]
  colnames(rna_hotspots)=paste("RNA_",colnames(rna_hotspots),sep="")

  joined_hotspots=merge(rna_hotspots,microbial_hotspots,by="row.names",all=T)  
  row.names(joined_hotspots)=joined_hotspots$Row.names
  joined_hotspots=joined_hotspots[,-1]
  joind_data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.AND.210708_raw_counts_mtrbcpfiltered_RNA_cutoff_",RNA_cutoff,".",grid_size,"x",grid_size,".Getis_Ord_hotspots",".tsv",sep="")
  write.table(x=joined_hotspots,file = joind_data_file,row.names = T,quote = F,sep = ";")
  
  joined_hotspots_network=rcorr(x = as.matrix(joined_hotspots),type = "spearman")
  joined_hotspots_network_flat=flattenCorrMatrix(joined_hotspots_network$r,joined_hotspots_network$P)
  joined_hotspots_network_flat$p.BH=p.adjust(joined_hotspots_network_flat$p,method = "BH")
    
  network_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.AND.210708_raw_counts_mtrbcpfiltered_RNA_cutoff_",RNA_cutoff,".",grid_size,"x",grid_size,".Getis_Ord_hotspots",".spearman_network.tsv",sep="")
  
  write.table(x=joined_hotspots_network_flat,file = network_file,row.names = T,quote = F,sep = ";")
  message (network_file)
}



# extract the microbial part of the hotspot network
# For the all_Bacterial all_Fungi
library("stringr")

exp="OMNI12"
exp="OMNI13"
grid_size=2
RNA_cutoff=100
smpl_list=c()
if (exp=="OMNI12") {smpl_list=c("A1","A2","B1","B2","C1","C2")}
if (exp=="OMNI13") {smpl_list=c("A1","A2","B1","B2","C1","C2","D1")}

for (smpl in smpl_list)
{
  message(paste("--",smpl))
  
  network_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.AND.210708_raw_counts_mtrbcpfiltered_RNA_cutoff_",RNA_cutoff,".",grid_size,"x",grid_size,".Getis_Ord_hotspots",".spearman_network.tsv",sep="")
  joined_hotspots_network_flat=read.delim(file=network_file,sep=";",stringsAsFactors = F)
  
  microbial_network=joined_hotspots_network_flat[grepl(pattern = "BACTERIA|FUNGI",x = joined_hotspots_network_flat$row),]
  microbial_network=rbind(microbial_network,joined_hotspots_network_flat[grepl(pattern = "BACTERIA|FUNGI",x = joined_hotspots_network_flat$column),])
  microbial_network=microbial_network[!grepl(pattern = "RNA_HOTSPOT",x = microbial_network$column),]
  microbial_network=microbial_network[!grepl(pattern = "RNA_HOTSPOT",x = microbial_network$row),]
  microbial_network=microbial_network[!grepl(pattern = "_all_Bacteria|_all_Fungi",x = microbial_network$column),]
  microbial_network=microbial_network[!grepl(pattern = "_all_Bacteria|_all_Fungi",x = microbial_network$row),]
  microbial_network$p.BH=p.adjust(microbial_network$p,method="BH")
  sum(microbial_network$p.BH<0.05)
  microbial_network$s1_domain=str_split_fixed(microbial_network$row, "_", 2)
  microbial_network$s2_domain=str_split_fixed(microbial_network$column, "_", 2)
  
  microbial_network$row=gsub(pattern = "BACTERIA_HOTSPOT_",replacement = "",x = microbial_network$row)
  microbial_network$row=gsub(pattern = "FUNGI_HOTSPOT_",replacement = "",x = microbial_network$row)
  
  microbial_network$column=gsub(pattern = "BACTERIA_HOTSPOT_",replacement = "",x = microbial_network$column)
  microbial_network$column=gsub(pattern = "FUNGI_HOTSPOT_",replacement = "",x = microbial_network$column)
  
  # save
  microbial_network_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/",exp,"/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/",exp,"_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.AND.210708_raw_counts_mtrbcpfiltered_RNA_cutoff_",RNA_cutoff,".",grid_size,"x",grid_size,".Getis_Ord_hotspots",".spearman_network.microbial_only.csv",sep="")
  write.table(x=microbial_network,file = microbial_network_file,row.names = F,quote = F,sep = ";")
}
















  breaks = c(-20, -1.96, -1, 1, 1.96, 20)
  palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
  col = palette[cut(getisgrid$HOTSPOT, breaks)]
  plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
  
  
  getisgrid.df=as.data.frame(getisgrid)
  names(getisgrid.df)=paste(names(getisgrid.df),layer,sep="_")
  getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
  names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
  getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
  getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
  row.names(getisgrid.df)=getisgrid.df$xy
  getisgrid.df=getisgrid.df[,-1]        
  if (layer=="all_Bacteria") {
    global_moran_all_Bacteria=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
  }
  if (layer=="all_Fungi") {
    global_moran_all_Fungi=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
  }
  if (layer=="all_genes") {
    global_moran_all_genes=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
  }
  if (nrow(hotspost_df)==0)
  {
    hotspost_df=getisgrid.df
  } else {
    hotspost_df=merge(hotspost_df,getisgrid.df[,1:2],by="row.names",all=T)
    row.names(hotspost_df)=hotspost_df$Row.names
    hotspost_df=hotspost_df[,-1]
  }
  legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
         legend =c("<-1.96","-1","1","1.96",">1.96"), 
         fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
         xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
  
  # 
  
  # Bacterial and Fungi
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
  
  # create the dataframe with the analyzed data -> would be used for creating the spatial object
  house_keeping_genes=c("AT3G18780",
                        "AT1G75780",
                        "AT1G49240",
                        "AT5G62690",
                        "AT2G28390")
  Expression_Fungi_Bacteria_df=expression_data_with_xy_subset[,grepl(paste0(c("x","y",house_keeping_genes,"all_genes"),collapse = "|"),names(expression_data_with_xy_subset))]
  # Expression_Fungi_Bacteria_df=expression_data_with_xy_subset[,grepl(paste0(c("x","y","all_genes"),collapse = "|"),names(expression_data_with_xy_subset))]
  Expression_Fungi_Bacteria_df=merge(Expression_Fungi_Bacteria_df,data_with_xy[,c("x","y","all_Bacteria","all_Fungi")],by=c("x","y"),all=T)
  # remove positions with no genes according to the previous criteria ->> need to figure this out!
  sum(is.na(Expression_Fungi_Bacteria_df$all_genes))
  # Expression_Fungi_Bacteria_df=Expression_Fungi_Bacteria_df[!is.na(Expression_Fungi_Bacteria_df$all_genes),]
  # View(Expression_Fungi_Bacteria_df[is.na(Expression_Fungi_Bacteria_df$all_genes),])
  
  # create one spatial object only with the (1) expression of house keeping, (2) all expression sum > 10 reads, (3) all Bacterial (4) all fungi
  spatial_data=SpatialPointsDataFrame(coords = Expression_Fungi_Bacteria_df[,c(1:2)],data = Expression_Fungi_Bacteria_df[,3:NCOL(Expression_Fungi_Bacteria_df)])
  
  # par(oma = c(4,1,1,1), mfrow = c(6, 2), mar = c(4, 4, 1, 1))
  for (layer in names(spatial_data)) {
    p=spplot(spatial_data,col='transparent',zcol=layer,main=list(label=layer),cex=1)
    plot(p)
    
    # no grid
    grid_size=1
    getisgrid=Getis_Ord_GI_per_spot(spatial_data = spatial_data,layer = layer)
    breaks = c(-20, -1.96, -1, 1, 1.96, 20)
    palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
    col = palette[cut(getisgrid$HOTSPOT, breaks)]
    plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
    
    grid_size=2
    getisgrid=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
    breaks = c(-20, -1.96, -1, 1, 1.96, 20)
    palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
    col = palette[cut(getisgrid$HOTSPOT, breaks)]
    plot(getisgrid, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
    
    
    getisgrid.df=as.data.frame(getisgrid)
    names(getisgrid.df)=paste(names(getisgrid.df),layer,sep="_")
    getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
    names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
    getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
    getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
    row.names(getisgrid.df)=getisgrid.df$xy
    getisgrid.df=getisgrid.df[,-1]        
    if (layer=="all_Bacteria") {
      global_moran_all_Bacteria=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
    }
    if (layer=="all_Fungi") {
      global_moran_all_Fungi=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
    }
    if (layer=="all_genes") {
      global_moran_all_genes=globalMoran_per_grid(grid_size=grid_size,spatial_data = spatial_data,genus_name = layer)
    }
    if (nrow(hotspost_df)==0)
    {
      hotspost_df=getisgrid.df
    } else {
      hotspost_df=merge(hotspost_df,getisgrid.df[,1:2],by="row.names",all=T)
      row.names(hotspost_df)=hotspost_df$Row.names
      hotspost_df=hotspost_df[,-1]
    }
    legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
           legend =c("<-1.96","-1","1","1.96",">1.96"), 
           fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
           xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
  } # finish loop over all genes per grid size
  
  
  # summary stats per grid size
  total_Bacteria_hotspots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96)
  total_Fungi_hotspots=sum(hotspost_df$HOTSPOT_all_Fungi>=1.96)
  total_Bacteria_coldspots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96)
  total_Fungi_coldspots=sum(hotspost_df$HOTSPOT_all_Fungi<=-1.96)
  total_Genes_coldspots=sum(hotspost_df$HOTSPOT_all_genes<=-1.96)
  total_Genes_coldspots=sum(hotspost_df$HOTSPOT_all_genes<=-1.96)
  overlap_Bacteria_Fungi_HotSpots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
  overlap_Bacteria_Fungi_ColdSpots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
  overlap_HotSpot_Bacteria_ColdSpot_Fungi=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
  overlap_ColdSpot_Bacteria_HotSpot_Fungi=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
  
  overlap_Bacteria_Genes_HotSpots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_genes>=1.96)
  overlap_Bacteria_Gense_ColdSpots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_genes<=-1.96)
  
  overlap_Fungi_Genes_HotSpots=sum(hotspost_df$HOTSPOT_all_Fungi>=1.96&hotspost_df$HOTSPOT_all_genes>=1.96)
  overlap_Fungi_Gense_ColdSpots=sum(hotspost_df$HOTSPOT_all_Fungi<=-1.96&hotspost_df$HOTSPOT_all_genes<=-1.96)
  
  global_moran_all_Bacteria_stat=global_moran_all_Bacteria[["estimate"]][["Moran I statistic"]]
  global_moran_all_Bacteria_p=global_moran_all_Bacteria[["p.value"]]
  global_moran_all_Fungi_stat=global_moran_all_Fungi[["estimate"]][["Moran I statistic"]]
  global_moran_all_Fungi_p=global_moran_all_Fungi[["p.value"]]
  
  global_moran_all_Genes_stat=global_moran_all_genes[["estimate"]][["Moran I statistic"]]
  global_moran_all_Genes_p=global_moran_all_genes[["p.value"]]
  
  record=data.frame(sample=smpl,
                    grid_size=grid_size,
                    total_Bacteria_hotspots=total_Bacteria_hotspots,
                    total_Fungi_hotspots=total_Fungi_hotspots,
                    total_Bacteria_coldspots=total_Bacteria_coldspots,
                    total_Fungi_coldspots=total_Fungi_coldspots,
                    overlap_Bacteria_Fungi_HotSpots=overlap_Bacteria_Fungi_HotSpots,
                    overlap_Bacteria_Fungi_ColdSpots=overlap_Bacteria_Fungi_ColdSpots,
                    overlap_HotSpot_Bacteria_ColdSpot_Fungi=overlap_HotSpot_Bacteria_ColdSpot_Fungi,
                    overlap_ColdSpot_Bacteria_HotSpot_Fungi=overlap_ColdSpot_Bacteria_HotSpot_Fungi,
                    overlap_Bacteria_Genes_HotSpots=overlap_Bacteria_Genes_HotSpots,
                    overlap_Bacteria_Gense_ColdSpots=overlap_Bacteria_Gense_ColdSpots,
                    overlap_Fungi_Genes_HotSpots=overlap_Fungi_Genes_HotSpots,
                    overlap_Fungi_Gense_ColdSpots=overlap_Fungi_Gense_ColdSpots,
                    
                    total_pixls=nrow(hotspost_df),
                    global_moran_all_Bacteria_stat=global_moran_all_Bacteria_stat,
                    global_moran_all_Bacteria_p=global_moran_all_Bacteria_p,
                    global_moran_all_Fungi_stat=global_moran_all_Fungi_stat,
                    global_moran_all_Fungi_p=global_moran_all_Fungi_p,
                    global_moran_all_Genes_stat=global_moran_all_Genes_stat,
                    global_moran_all_Genes_p=global_moran_all_Genes_p
  )
  hotspots_overlaps_per_grid_size=rbind(hotspots_overlaps_per_grid_size,record)
  




####### DRAFTS  
  
  
  
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy_subset[,c(1,2)],data = data_with_xy_subset[,c(3:ncol(data_with_xy_subset))])
  
  
  p=spplot(spatial_data,col='transparent',zcol="all_genes") # main=list(label=paste(genus_name," - ",type," ",i),cex=1)
  plot(p)
  for (layer in c(row.names(gene_sum_ordered)[1:10],"all_genes")) {
    grid_size=3
    getisgrid_3=Getis_Ord_GI_per_grid(grid_size = grid_size,spatial_data = spatial_data,genus_name = layer)
    breaks = c(-20, -1.96, -1, 1, 1.96, 20)
    palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
    col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
    plot(getisgrid_3, col=col,main=paste(layer,"-",grid_size,"x",grid_size))
  }
  getisgrid_3=Getis_Ord_GI_per_grid(grid_size = 3,spatial_data = spatial_data,genus_name = "all_genes")
  breaks = c(-20, -1.96, -1, 1, 1.96, 20)
  palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
  col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
  plot(getisgrid_3, col=col) #,main=paste("all_genes","-",type," ",i,"-",grid_size,"x",grid_size)
  ### HERE
  
  breaks = c(-20, -1.96, -1, 1, 1.96, 20)
  palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
  par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(4, 4, 1, 1))
  grid_size=2
  col = palette[cut(getisgrid_2$HOTSPOT, breaks)]
  plot(getisgrid_2, col=col,main=paste(genus_name,"-",type," ",i,"-",grid_size,"x",grid_size))
  
  grid_size=3
  
  
  View(spatial_data[,"all_genes"])
  
  
  sum(test$count>20)
  hist(test)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.pdf",sep="")
  pdf(out_pdfs,width = 7,height = 7)
  out_xls=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.xlsx",sep="")
  XSLX_Obj=createWorkbook(out_xls)
  
  out_summary_overlaps=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots.Bacteria_Fungi_overlap.csv",sep="")
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  for (grid_size in c(2:15,24))
  {
    i=0
    hotspost_df=data.frame()
    global_moran_all_Bacteria=NA
    global_moran_all_Fungi=NA
    for (genus_name in names(spatial_data)) # loop all genus and specific grid size
    {
      i=i+1
      type="Bacteria"
      if (i>50&i<100) {type="Fungi"}
      if (i>100) {type=""}
      total=sum(spatial_data[[genus_name]])
      if(total>5)
      {
        if (grid_size==3) # plot only for the first time
        {
          # http://www.nickeubank.com/wp-content/uploads/2015/10/RGIS3_MakingMaps_part1_mappingVectorData.html
          # spatial_data_for_plot=spatial_data
          # spatial_data_for_plot[[genus_name]][spatial_data_for_plot[[genus_name]]==0]=NA
          # breaks.qt <- classIntervals(spatial_data_for_plot[[genus_name]], n = 5, style = "equal", intervalClosure = "right")
          # my.palette <- brewer.pal(n = 8, name = "PuRd")
          # p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name," - ",type," ",i),cex=1),cuts=9,at=c(-1,1,breaks.qt$brks[2:length(breaks.qt$brks)]))
          p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name," - ",type," ",i),cex=1))
          plot(p)
          rm(p)
        }
        
        # Create a regular grid of 1/8-mile-square crime areas via a raster
        # library(raster)
        message(paste("\t- genus ",genus_name," grid size=",grid_size,sep=""))
        pixelsize = grid_size # size of the grid box -> seems so by the plot(getisgrid)
        box = round(extent(spatial_data) / pixelsize) * pixelsize
        template = raster(box, crs = spatial_data,
                          nrows = (box@ymax - box@ymin) / pixelsize, 
                          ncols = (box@xmax - box@xmin) / pixelsize)
        getisraster = rasterize(spatial_data, template, field = genus_name, fun = sum)
        getisgrid = rasterToPolygons(getisraster)
        # plot(getisgrid)
        # Create the list of neighbors
        neighbors = poly2nb(getisgrid)
        weighted_neighbors = nb2listw(neighbors, zero.policy=T)
        
        # Perform the local G analysis (Getis-Ord GI*)
        getisgrid$HOTSPOT = as.vector(localG(getisgrid$layer, weighted_neighbors))
        # Color the grid cells based on the z-score
        # join with the others
        getisgrid.df=as.data.frame(getisgrid)
        names(getisgrid.df)=paste(names(getisgrid.df),genus_name,sep="_")
        getisgrid.cooridinate_poligon_centers=as.data.frame(coordinates(getisgrid))
        names(getisgrid.cooridinate_poligon_centers)=c("X_Pcent","Y_Pcent")
        getisgrid.df=merge(getisgrid.df,getisgrid.cooridinate_poligon_centers,by="row.names")
        getisgrid.df$xy=paste(getisgrid.df$X_Pcent,getisgrid.df$Y_Pcent,sep="x")
        row.names(getisgrid.df)=getisgrid.df$xy
        getisgrid.df=getisgrid.df[,-1]        
        if (genus_name=="all_Bacteria") {
          global_moran_all_Bacteria=moran.test(getisgrid$layer, weighted_neighbors)
        }
        if (genus_name=="all_Fungi") {
          global_moran_all_Fungi=moran.test(getisgrid$layer, weighted_neighbors)
        }
        if (nrow(hotspost_df)==0)
        {
          hotspost_df=getisgrid.df
        } else {
          hotspost_df=merge(hotspost_df,getisgrid.df[,1:2],by="row.names",all=T)
          row.names(hotspost_df)=hotspost_df$Row.names
          hotspost_df=hotspost_df[,-1]
        }
        breaks = c(-20, -1.96, -1, 1, 1.96, 20)
        palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
        col = palette[cut(getisgrid$HOTSPOT, breaks)]
        
        # Plot
        # spplot(getisgrid)
        # par(mfrow=c(2,1))
        # print (p)
        # par(mfrow=c(2,2))
        plot(getisgrid, col=col,main=paste("local G analysis",genus_name," - ",type," ",i,"\n",grid_size,"x",grid_size,"total read=",total))
        legend("bottomleft", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
               c("<-1.96","-1","1","1.96",">1.96"), fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"), horiz=TRUE, cex=0.8)
        rm(getisgrid)
        
      } # end condition on number of reads
      #sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
      #sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    } # end looping all genus with specific grid size
    # write the hotspot data for the spcific grid size
    sheetName = paste("gridSize=",grid_size,sep="")
    addWorksheet(XSLX_Obj, sheetName)
    writeDataTable(x=hotspost_df, wb = XSLX_Obj,sheet = sheetName,rowNames = TRUE)
    
    # summary stats per grid size
    total_Bacteria_hotspots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96)
    total_Fungi_hotspots=sum(hotspost_df$HOTSPOT_all_Fungi>=1.96)
    total_Bacteria_coldspots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96)
    total_Fungi_coldspots=sum(hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_Bacteria_Fungi_HotSpots=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
    overlap_Bacteria_Fungi_ColdSpots=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_HotSpot_Bacteria_ColdSpot_Fungi=sum(hotspost_df$HOTSPOT_all_Bacteria>=1.96&hotspost_df$HOTSPOT_all_Fungi<=-1.96)
    overlap_ColdSpot_Bacteria_HotSpot_Fungi=sum(hotspost_df$HOTSPOT_all_Bacteria<=-1.96&hotspost_df$HOTSPOT_all_Fungi>=1.96)
    global_moran_all_Bacteria_stat=global_moran_all_Bacteria[["estimate"]][["Moran I statistic"]]
    global_moran_all_Bacteria_p=global_moran_all_Bacteria[["p.value"]]
    global_moran_all_Fungi_stat=global_moran_all_Fungi[["estimate"]][["Moran I statistic"]]
    global_moran_all_Fungi_p=global_moran_all_Fungi[["p.value"]]
    record=data.frame(grid_size=grid_size,
                      total_Bacteria_hotspots=total_Bacteria_hotspots,
                      total_Fungi_hotspots=total_Fungi_hotspots,
                      total_Bacteria_coldspots=total_Bacteria_coldspots,
                      total_Fungi_coldspots=total_Fungi_coldspots,
                      overlap_Bacteria_Fungi_HotSpots=overlap_Bacteria_Fungi_HotSpots,
                      overlap_Bacteria_Fungi_ColdSpots=overlap_Bacteria_Fungi_ColdSpots,
                      overlap_HotSpot_Bacteria_ColdSpot_Fungi=overlap_HotSpot_Bacteria_ColdSpot_Fungi,
                      overlap_ColdSpot_Bacteria_HotSpot_Fungi=overlap_ColdSpot_Bacteria_HotSpot_Fungi,
                      total_pixls=nrow(hotspost_df),
                      global_moran_all_Bacteria_stat=global_moran_all_Bacteria_stat,
                      global_moran_all_Bacteria_p=global_moran_all_Bacteria_p,
                      global_moran_all_Fungi_stat=global_moran_all_Fungi_stat,
                      global_moran_all_Fungi_p=global_moran_all_Fungi_p
    )
    hotspots_overlaps_per_grid_size=rbind(hotspots_overlaps_per_grid_size,record)
    rm (hotspost_df,"getisgrid.df","getisraster",
        "overlap_Bacteria_Fungi_ColdSpots","overlap_Bacteria_Fungi_HotSpots",
        "record","sheetName","neighbors","template","weighted_neighbors")
    
    rm ("box","breaks","col")
    
  } # end loop of different grid sizes
  # write the hotspots_overlaps_per_grid_size data
  write.table(file = out_summary_overlaps,x = hotspots_overlaps_per_grid_size,quote = F,sep = ";",row.names = F)
  saveWorkbook(XSLX_Obj, file = out_xls, overwrite = TRUE)
  
  
  rm ("data" ,"data_file","data_with_xy","genus_name","grid_size","hotspots_overlaps_per_grid_size","i",
      "out_pdfs","out_summary_overlaps","out_xls","palette","pixelsize","spatial_data",
      "XSLX_Obj","xy","xy_list","total")
  
  dev.off()
}

# for easy comparision plot again now for each grid size
for (smpl in c("A1","A2","B1","B2","C1","C2"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.genus.UnderTissue.hotspots_compact.pdf",sep="")
  pdf(out_pdfs,width = 14,height = 14)
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  
  i=0
  for (genus_name in names(spatial_data)) # loop all genus of a sample
  {
    i=i+1
    type="Bacteria"
    if (i>50&i<100) {type="Fungi"}
    if (i>100) {type=""}
    total=sum(spatial_data[[genus_name]])
    if(total>5)
    {
      par(mfrow = c(1, 1))
      p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name," - ",type," ",i,"-","total reads",total),cex=1))
      plot(p)
      getisgrid_2=Getis_Ord_GI_per_grid(grid_size = 2,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_3=Getis_Ord_GI_per_grid(grid_size = 3,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_6=Getis_Ord_GI_per_grid(grid_size = 6,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_12=Getis_Ord_GI_per_grid(grid_size = 12,spatial_data = spatial_data,genus_name = genus_name)
      breaks = c(-20, -1.96, -1, 1, 1.96, 20)
      palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(4, 4, 1, 1))
      grid_size=2
      col = palette[cut(getisgrid_2$HOTSPOT, breaks)]
      plot(getisgrid_2, col=col,main=paste(genus_name,"-",type," ",i,"-",grid_size,"x",grid_size))
      grid_size=3
      col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
      plot(getisgrid_3, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      grid_size=6
      col = palette[cut(getisgrid_6$HOTSPOT, breaks)]
      plot(getisgrid_6, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      grid_size=12
      col = palette[cut(getisgrid_12$HOTSPOT, breaks)]
      plot(getisgrid_12, col=col,main=paste(genus_name,"-",type,i,"-",grid_size,"x",grid_size))
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
             legend =c("<-1.96","-1","1","1.96",">1.96"), 
             fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
             xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    }
  }
  dev.off()
}

# for easy comparision plot again now for each grid size only all Bacteria all Fungi
for (smpl in c("A1","A2","B1","B2","C1","C2"))
{
  message (paste("-- ",smpl))
  hotspots_overlaps_per_grid_size=data.frame()
  data_file=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.usearch_unique_vs_NT_Jan2021.UMI_filtered.genus.spatial_pos_UnderTissue.transposed.csv",sep="")
  data=read.delim(data_file,sep=";",stringsAsFactors = F)
  data$all_Bacteria=rowSums(data[,1:50])
  data$all_Fungi=rowSums(data[,51:100])
  
  out_pdfs=paste("/Volumes/spatial_array_metatranscriptomics/data/Spatial_transcriptomics/omni_array/Haim/MMSEQS2/OMNI12/spatial_tables/joint_Bacteria_Fungi_and_UNKNOWN_transposed/hotspots/OMNI12_",smpl,".Top50_ITS_16S_Probs_and_UNKNOWN.AllBacFun.UnderTissue.hotspots.pdf",sep="")
  pdf(out_pdfs,width = 14,height = 14)
  
  xy_list=strsplit(row.names(data), "x")
  xy=as.data.frame(t(as.data.frame(xy_list)))
  xy$V1=gsub(pattern = "X",replacement = "",x = xy$V1)
  names(xy)=c("x","y")
  data_with_xy=cbind(xy,data)
  rownames(data_with_xy)=rownames(data)
  data_with_xy$x=as.numeric(data_with_xy$x)
  data_with_xy$y=as.numeric(data_with_xy$y)
  spatial_data=SpatialPointsDataFrame(coords = data_with_xy[,c(1:2)],data = data_with_xy[,c(3:104)])
  #par(mfrow=c(2,1))
  
  i=0
  for (genus_name in names(spatial_data)[c(101,102)]) # loop all genus of a sample
  {
    i=i+1
    type="Bacteria"
    if (i>50&i<100) {type="Fungi"}
    if (i>100) {type=""}
    total=sum(spatial_data[[genus_name]])
    if(total>5)
    {
      par(mfrow = c(1, 1))
      p=spplot(spatial_data,col='transparent',zcol=genus_name,main=list(label=paste(genus_name,"total reads",total),cex=1))
      plot(p)
      getisgrid_2=Getis_Ord_GI_per_grid(grid_size = 2,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_3=Getis_Ord_GI_per_grid(grid_size = 3,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_6=Getis_Ord_GI_per_grid(grid_size = 6,spatial_data = spatial_data,genus_name = genus_name)
      getisgrid_12=Getis_Ord_GI_per_grid(grid_size = 12,spatial_data = spatial_data,genus_name = genus_name)
      breaks = c(-20, -1.96, -1, 1, 1.96, 20)
      palette=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080")
      par(oma = c(4,1,1,1), mfrow = c(2, 2), mar = c(4, 4, 1, 1))
      grid_size=2
      col = palette[cut(getisgrid_2$HOTSPOT, breaks)]
      plot(getisgrid_2, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_2),"regions"))
      grid_size=3
      col = palette[cut(getisgrid_3$HOTSPOT, breaks)]
      plot(getisgrid_3, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_3),"regions"))
      grid_size=6
      col = palette[cut(getisgrid_6$HOTSPOT, breaks)]
      plot(getisgrid_6, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_6),"regions"))
      grid_size=12
      col = palette[cut(getisgrid_12$HOTSPOT, breaks)]
      plot(getisgrid_12, col=col,main=paste(genus_name,grid_size,"x",grid_size,"-",nrow(getisgrid_12),"regions"))
      par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
      legend("bottom", inset=.02, title="Z score (local G analysis (Getis-Ord GI*))",
             legend =c("<-1.96","-1","1","1.96",">1.96"), 
             fill=c("#0000FF80", "#8080FF80", "#FFFFFF80", "#FF808080", "#FF000080"),
             xpd = TRUE, horiz  = TRUE, cex = 1, seg.len=1, bty = 'n')
    }
  }
  dev.off()
}


library(Seurat)
library(GEOquery)
library(tidyverse)
library(SeuratDisk)
library(ggplot2)
library(sctransform)
library(future)
library(MAST)
library(qs)
library(EnhancedVolcano)
library(xlsx)
library(gridExtra)
library(STRINGdb)
library(igraph)

#helper functions
source("covidutils.R")

##multiprocess settings (requires a ton of RAM)
plan("multiprocess")
options(future.globals.maxSize = 15000 *1024^2)
options(mc.cores = 8)


setwd("/home/kai/Documents/COVID_Data_Mine")

PBMC <- readRDS("/home/kai/Documents/COVID_Data_Mine/Blish PBMCs/blish_covid.seu.rds")

## Reapply all analysis steps and more stringent filtering

PBMC <- PercentageFeatureSet(PBMC, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
# VlnPlot(PBMC, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# # # plot1 <- FeatureScatter(PBMC, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # plot2 <- FeatureScatter(PBMC, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))

PBMC <- subset(PBMC, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 1000)
PBMC <- SCTransform(PBMC, vars.to.regress = "percent.mt", verbose = FALSE)

#Loading GSE data info for BAL dataset
# setwd("/media/kai/860 SSD/Single Cell/COVID")
covidGSE145926<-getGEO(filename = "GSE145926_series_matrix.txt.gz", GSEMatrix = TRUE)
covidGSE145926<-covidGSE145926@phenoData
covidGSE145926@data<-covidGSE145926@data[,c("title","geo_accession","cell subsets:ch1","patient group:ch1","supplementary_file_1")]

#Additional metadata from github code release
meta <- read.delim("meta.txt")
meta$sample_new_old[c(5,6,7)]<-meta$sample_new[c(5,6,7)]
meta$newgroup<-meta$group
meta$newgroup[9:13]<-"C"
covidGSE145926@data$title <- sub("BALF, ","",covidGSE145926@data$title)
covidGSE145926@data$title <- sub(" .*","",covidGSE145926@data$title)
covidGSE145926@data$patientseverity <- meta$sample_new_old[match(covidGSE145926@data$title,meta$sample)]
covidGSE145926@data$group <- meta$newgroup[match(covidGSE145926@data$title,meta$sample)]


files <- paste0(list.files(pattern=".*.h5", full.names=TRUE, recursive=FALSE))
name <- files

covid.list <- vector(mode = "list", length = length(files))
for(ii in 1:length(files)) {
  
  data.path <- files[ii]
  name[ii] <- sub("_C.*","",data.path)
  name[ii] <- sub("./","", name[ii])
  raw.data <-Read10X_h5(data.path)
  #datamatrix <- Matrix(data, sparse = TRUE)
  
  # BAL <- AddMetaData(BAL, metadata = predictions)
  #updated loading into list to facilitate seurat v3 cca integration
  covid.list[[ii]] <- CreateSeuratObject(counts = raw.data, min.cells = 0, min.features = 0, 
                                         project = name[ii])
  names(covid.list)[ii] <- name[ii]
  covid.list[[ii]] <- AddMetaData(covid.list[[ii]],metadata = covidGSE145926@data[name[ii],"title"],col.name = "patientcode")
  covid.list[[ii]] <- AddMetaData(covid.list[[ii]],metadata = covidGSE145926@data[name[ii],"cell subsets:ch1"],col.name = "accession")
  covid.list[[ii]] <- AddMetaData(covid.list[[ii]],metadata = covidGSE145926@data[name[ii],"patient group:ch1"],col.name = "condition")
  covid.list[[ii]] <- AddMetaData(covid.list[[ii]],metadata = covidGSE145926@data[name[ii],"patientseverity"],col.name = "patientseverity")
  covid.list[[ii]] <- AddMetaData(covid.list[[ii]],metadata = covidGSE145926@data[name[ii],"group"],col.name = "group")
  print(paste0("Loaded ", name[ii]))
}


#remove individual seurat objects
rm(raw.data)


for(ii in 1:length(covid.list)){
  ##QC loop on each individual run before integration
  covid.list[[ii]][["percent.mt"]] <- PercentageFeatureSet(covid.list[[ii]], pattern = "^MT-")
  seuratdata<-covid.list[[ii]]
  
  # Visualize QC metrics as a violin plot
  VlnPlot(seuratdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
  # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
  
  plot1 <- FeatureScatter(seuratdata, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(seuratdata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  CombinePlots(plots = list(plot1, plot2))
  
  
  BAL <- subset(covid.list[[ii]], subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 1000)
  
  BAL <- SCTransform(BAL, vars.to.regress = "percent.mt", verbose = FALSE)
  ##Normalization
  # BAL <- NormalizeData(BAL, normalization.method = "LogNormalize", scale.factor = 10000) #default settings
  
  ##Feature Selection
  # BAL <- FindVariableFeatures(BAL, selection.method = "vst", nfeatures = 2000)
  
  # Identify the 10 most highly variable genes
  # top10 <- head(VariableFeatures(BAL), 10)
  
  # plot variable features with and without labels
  # # plot1 <- VariableFeaturePlot(BAL)
  # plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
  # CombinePlots(plots = list(plot1, plot2))
  
  covid.list[[ii]] <- BAL
}

covid.list[[13]]<-PBMC
names(covid.list)[13] <- c("blishPBMC")


##SCT preintegration steps
covid.features <- SelectIntegrationFeatures(object.list = covid.list, nfeatures = 3000)
covid.list <- PrepSCTIntegration(object.list = covid.list, anchor.features = covid.features, 
                                    verbose = FALSE)

### Seurat integtration (removes batch effects and generates global variable features as well as CCA reduction)
covid.anchors <- FindIntegrationAnchors(object.list = covid.list, dims = 1:50, normalization.method = "SCT", anchor.features = covid.features)

save.image(file="covidpreint.RData")
gc()
options(mc.cores = 2)
covid.integrated <- IntegrateData(anchorset = covid.anchors, dims = 1:50, normalization.method = "SCT")

save.image(file="covid.RData")


##adjust condition naming conventions. note that S1 and M1 are the  same patient
covid.integrated$Donor.strat<-recode(covid.integrated$Donor.full,"C1 A" = "M1", "C1 B" = "S1","C2" = "M2", "C3" = "S2","C4" = "S3",
                                     "C5" = "M3","C6" = "S4","C7" = "M4")
covid.integrated$Status.strat<-gsub("[0-9]","",covid.integrated$Donor.strat)


##aggregate patient IDs #Critical cases renamed S2-S6
aggpatient<-gsub("HC","H",covid.integrated$patientseverity)
aggpatient<-recode(aggpatient, "C1" = "S2","C2" = "S3","C3" = "S4","C4"="S5","C5"="S6")
aggpatient[is.na(aggpatient)]<-paste0("PBMC_",covid.integrated$Donor.strat[is.na(aggpatient)])
covid.integrated$aggpatient<-aggpatient

aggcondition<-gsub("HC","H",covid.integrated$group)
aggcondition<-gsub("C","S",aggcondition)
aggcondition[is.na(aggcondition)]<-paste0("PBMC_",covid.integrated$Status.strat[is.na(aggcondition)])
covid.integrated$aggcondition<-aggcondition

## dim reduction
options(mc.cores = 4)

covid.integrated <- RunPCA(covid.integrated, verbose = FALSE)
covid.integrated <- FindNeighbors(covid.integrated, dims = 1:30)
covid.integrated <- FindClusters(covid.integrated, resolution = 0.5)

ElbowPlot(covid.integrated)
covid.integrated <- RunUMAP(covid.integrated, dims = 1:30)
DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("group","Status"),pt.size = .5)
DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("cell.type.fine"),pt.size = 1,label.size = 8,repel = TRUE)
DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("seurat_clusters_.5"),pt.size = .5)

###### Marker identification #####

##multiprocess settings Only for differential expression. these settigns will overload ram on preprocessing
plan("multiprocess")
options(future.globals.maxSize = 3000 *1024^2)
options(mc.cores = 16)

## first we need to either normalize the RNA assay or rerun SCT normalization on the combined dataset
DefaultAssay(covid.integrated)<-"RNA"
covid.integrated <- SCTransform(covid.integrated, vars.to.regress = "percent.mt", verbose = FALSE)
DefaultAssay(covid.integrated)<-"SCT"
covid.integrated.markers <- FindAllMarkers(covid.integrated, assay = "SCT", 
                                           slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

# ##cross reference with labeled PBMC dataset, prep dataset via SCT
# pbmctbcyto <- readRDS("pbmctbcytocleaned")
# pbmctbcyto <- SCTransform(pbmctbcyto, vars.to.regress = "percent.mt", verbose = FALSE)
# 
# 
# # #SCT preintegration
# # covid.transfer.features <- SelectIntegrationFeatures(object.list = covid.list, nfeatures = 3000)
# # covid.list <- PrepSCTIntegration(object.list = covid.list, anchor.features = covid.features,
# #                                  verbose = FALSE)
# options(mc.cores = 4)
# ##Seurat label transfer
# pbmctbcyto.anchors <- FindTransferAnchors(reference = pbmctbcyto, query = covid.integrated,
#                                           dims = 1:30, project.query = TRUE)
# predictions <- TransferData(anchorset = pbmctbcyto.anchors, refdata = pbmctbcyto$orig.ident,
#                             dims = 1:30)
# covid.integrated <- AddMetaData(covid.integrated, metadata = predictions)
# DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("predicted.id","cell.type.fine"),pt.size = 1,label.size = 8,repel = TRUE)
# DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("seurat_clusters","cell.type.fine"),pt.size = 1,label.size = 8,repel = TRUE)
# DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("cell.type.fine"),pt.size = .5,label.size = 8,repel = TRUE)

#### Per cluster identification ####
options(mc.cores = 16)
marker.list<-dlply(covid.integrated.markers, .(cluster))

##Cluster 0 - CD14 Monocytes, FCN1+
FeaturePlot(covid.integrated,features=marker.list$"0"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = "CD14",group.by = "seurat_clusters", pt.size = 0)

##Cluster 1 - CD14 Monocytes, FCN1-
FeaturePlot(covid.integrated,features=marker.list$"1"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = c("CD14","FCN1","CTSS"),group.by = "seurat_clusters", pt.size = 0)

 ##Cluster 2 - TCF7+ Memory T cells CD4+
FeaturePlot(covid.integrated,features=marker.list$"2"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = c("IL7R","LTB","ETS1","TCF7"),group.by = "seurat_clusters", pt.size = 0)
VlnPlot(covid.integrated,features = c("CD4"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 3 - Gamma Delta T cells CD8+
FeaturePlot(covid.integrated,features=marker.list$"3"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = c("CCL5","CD8A","IL32","GZMA"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 4 - CD14 Monocytes, FCN1+, CCL4+
FeaturePlot(covid.integrated,features=marker.list$"4"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"4"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("CCL4L2","CCL4","CXCL8"),group.by = "seurat_clusters", pt.size = 0)
VlnPlot(covid.integrated,features = c("CD14","FCN1","CTSS"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 5 - CD14 Monocytes, FCN1+, CCL8+,CCL2+,CCL7+
FeaturePlot(covid.integrated,features=marker.list$"5"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = c("CCL8","CCL2","CCL7"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 6 - CD14++ Monocytes, FCN1+, VCAN++, implicated in innate immune respose to lung infection - Chang 2017
FeaturePlot(covid.integrated,features=marker.list$"6"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = c("S100A8","S100A9","VCAN","CD14"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 7 - NK Cells -- updated after subclustering
FeaturePlot(covid.integrated,features=marker.list$"7"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = c("GNLY","GZMB","PRF1","NKG7","SYNE1"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 8 - CD14 Monocytes FCN1+ -- updated after subclustering FKBP4+, IER5+
FeaturePlot(covid.integrated,features=marker.list$"8"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"8"$gene[10:18],pt.size = .5)


##Cluster 9 - B cells
FeaturePlot(covid.integrated,features=marker.list$"9"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"9"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("MS4A1","BANK1","CD79A","PAX5"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 10 - Gamma Delta T cell CD8+-
FeaturePlot(covid.integrated,features=marker.list$"10"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = c("PRF1","FGFBP2","CD8A"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 11 - T cell - TCF7+, CD4- CD8-
markers11 <- FindMarkers(covid.integrated, ident.1 = 11, assay = "SCT", 
                                           slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")
FeaturePlot(covid.integrated,features=marker.list$"11"$gene[1:9],pt.size = .5)
VlnPlot(covid.integrated,features = c("IL7R","LTB","ETS1","TCF7"),group.by = "seurat_clusters", pt.size = 0)
VlnPlot(covid.integrated,features = c("CD4","CD8A","CD3E","TRBV1"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 12 - CCL2+ CD14+ FCN1+  Monocytes
FeaturePlot(covid.integrated,features=marker.list$"12"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=c("CD14","FCN1","CCL2"),pt.size = .1)
FeaturePlot(covid.integrated,features=c("CCL2"),pt.size = .1)
VlnPlot(covid.integrated,features = c("CCL2"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 13 - Lung epithelium and Basal Cells 
FeaturePlot(covid.integrated,features=marker.list$"13"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"13"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("SLPI","SCGB1A1","WFDC2","KRT19"),group.by = "seurat_clusters", pt.size = 0)

pneumomarkers<-FindMarkers(covid.integrated, ident.1 = "Epithelium/Basal",ident.2 = "Epithelium/Pneumocytes/Ciliary",assay = "SCT",slot="data",test.use="MAST",latent.vars="nCount_RNA")


##Cluster 14 - Neutrophils
FeaturePlot(covid.integrated,features=marker.list$"14"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"14"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("CXCL8","NAMPT","FCGR3B","S100A8","IFITM2","IL1RN"),group.by = "seurat_clusters", pt.size = 0)
neutrophilmarkers<-FindMarkers(covid.integrated, ident.1 = "Neutrophils",assay = "SCT",slot="data",test.use="MAST",latent.vars="nCount_RNA")


##Cluster 15 - Plasmablasts
FeaturePlot(covid.integrated,features=marker.list$"15"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"15"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("IGJ","IGHG1","IGHA1","MZB1"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 16 - T cell CD8+ CD3E+ IL7R-, possible treg - LAG3+, --listig as CD8 effector MKI67+, 
FeaturePlot(covid.integrated,features=marker.list$"16"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"16"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("CD8A","CD3E", "CD4","IL32","GZMA","LAG3","MKI67"),group.by = "seurat_clusters", pt.size = 0)
VlnPlot(covid.integrated,features = c("GNLY","GZMB","PRF1","NKG7","SYNE1"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 17 - CD14+/CD4+ Monocytes -- updated after subclustering
FeaturePlot(covid.integrated,features=marker.list$"17"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"17"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("CD4","CD8A"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 18 - Tregs CD4+, IL2RA+
FeaturePlot(covid.integrated,features=marker.list$"18"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"18"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("CD4","IL2RA","LAG3"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 19 - Monocytes FCN1+ CD14+ no distinct markers --updated, labeling with cluster0
FeaturePlot(covid.integrated,features=marker.list$"19"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"19"$gene[10:18],pt.size = .5)

##Cluster 20 - Monocytes CCL18+ FCN1-
FeaturePlot(covid.integrated,features=marker.list$"20"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"20"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("CD14","FCN1","CCL18","CCL2"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 21 - Lung Epithelium - Ependymal, Type II Pneumocytes, Cilited cells
FeaturePlot(covid.integrated,features=marker.list$"21"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"21"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("SLPI","CAPS","WFDC2","CD24","KRT19"),group.by = "seurat_clusters", pt.size = 0)
VlnPlot(covid.integrated,features = c("RSPH1","GSTA1","FAM183A","TPPP3"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 22 - Plasmacytoid Dendritic Cells
FeaturePlot(covid.integrated,features=marker.list$"22"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"22"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("IRF8","UGCG","TSPAN13","SERPINF1","LILRA"),group.by = "seurat_clusters", pt.size = 0)
VlnPlot(covid.integrated,features = c("PLD4","MAP1A","PTPRS","CLEC4C"),group.by = "seurat_clusters", pt.size = 0)

##Cluster 23 - Gamma Delta T Cells (memory?) TCF7+ 
FeaturePlot(covid.integrated,features=marker.list$"23"$gene[1:9],pt.size = .5)
FeaturePlot(covid.integrated,features=marker.list$"23"$gene[10:18],pt.size = .5)
VlnPlot(covid.integrated,features = c("IL7R","LTB","ETS1","TCF7","CD4"),group.by = "seurat_clusters", pt.size = 0)
VlnPlot(covid.integrated,features = c("CCL5","IL32","GZMA","TCF7"),group.by = "seurat_clusters", pt.size = 0)

qsave(covid.integrated,"covid.integrated.qs",nthreads = 16)
cat( "Compression Ratio: ", as.numeric(object.size(covid.integrated)) / file.info("covid.integrated.qs")$size, "\n" )

covid.integrated <- RenameIdents(covid.integrated, `0` = "FCN1+/VCAN+ Monocytes", `1` = "APOE+/CCL18+/FCN1- Monocytes", `2` = "CD4+/TCF7+ Memory T Cells", `3` = "CD8+ Gamma Delta T Cells",
                          `4` = "FCN1+/CCL4+ Monocytes", `5` = "FCN1+/CCL2/7/8+ Monocytes", `6` = "CD14+/FCN1+/VCAN++/S100A8++/FCGR3A- Monocytes", `7` = "NK", `8` = "CD14-/FCN1+/FKBP4+/IER5+ Monocytes", `9` = "B cells", `10` = "CD8+- Gamma Delta T Cells",
                          `11` = "CD4-/CD8-/TCF7+ T Cells", `12` = "CD14-/FCN1+/CCL2+ Monocytes", `13` = "Epithelium/Basal Cells", `14` = "Neutrophils",
                          `15` = "Plasmablasts", `16` = "CD8+/MKI67+ T Effector Cells", `17` = "MT+ Monocytes", `18` = "T Regulatory Cells", `19` = "FN1+/MRC1+/FCN1- Monocytes", `20` = "FN1+/CD14-/CCL18++/FCN1- Monocytes",
                          `21` = "Epithelium/Pneumocytes/Ciliary Cells", `22` = "Plasmacytoid Dendritic Cells", `23` = "TCF7+ Gamma Delta T Cells")
covid.integrated$cell.type.0.5 <- Idents(covid.integrated)

monocyteclusters<-c(0,1,4,5,6,8,12,17,19,20)

Idents(covid.integrated)<-covid.integrated$seurat_clusters_.5
covid.integrated <- RenameIdents(covid.integrated, `0` = "Monocytes", `1` = "Monocytes", `2` = "Memory T Cells", `3` = "CD8+ T Cells",
                                 `4` = "Monocytes", `5` = "Monocytes", `6` = "Monocytes", `7` = "NK", `8` = "Monocytes", `9` = "B cells", `10` = "CD8+ T Cells",
                                 `11` = "DN T Cells", `12` = "Monocytes", `13` = "Epithelium/Basal Cells", `14` = "Neutrophils",
                                 `15` = "Plasmablasts", `16` = "T Effector Cells", `17` = "Monocytes", `18` = "T Regulatory Cells", `19` = "Monocytes", `20` = "Monocytes",
                                 `21` = "Epithelium/Pneumocytes/Ciliary Cells", `22` = "Plasmacytoid Dendritic Cells", `23` = "Gamma Delta T Cells")
covid.integrated$cell.type.coarse.0.5<-Idents(covid.integrated)

uplot<-DimPlot(covid.integrated, reduction = "umap",label=FALSE,group.by = c("cell.type.merge"),
        pt.size = .5,label.size = 5,repel = TRUE)+NoLegend()+scale_color_viridis_d(begin = 0,end =1)
LabelPlotBackground(plot = uplot, id = "cell.type.merge", color = "black")

VlnPlot(covid.integrated,features = c("CD14","FCN1","CCL18","CD4","CD8A","MS4A1"),group.by = "seurat_clusters_.5", pt.size = 0)

DimPlot(covid.integrated, reduction = "umap",label=FALSE,group.by = c("aggcondition"),
        pt.size = .5,label.size = 5,repel = TRUE)+scale_color_viridis_d()


##### Second round of clustering and identification #####
covid.integrated$seurat_clusters_.5 <- covid.integrated$seurat_clusters
DefaultAssay(covid.integrated)<-"integrated"
covid.integrated <- FindClusters(covid.integrated, resolution = 1)
DefaultAssay(covid.integrated)<-"SCT"
DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("seurat_clusters","cell.type.0.5"),pt.size = 1)
DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("cell.type.0.5"),pt.size = 1)
clustertable<-table(as.data.frame(cbind(covid.integrated$seurat_clusters,covid.integrated$seurat_clusters_.5)))
colnames(clustertable)<-levels(covid.integrated$seurat_clusters_.5)
rownames(clustertable)<-levels(covid.integrated$seurat_clusters)

DimPlot(covid.integrated, reduction = "umap",label=TRUE,group.by = c("seurat_clusters","seurat_clusters_.5"),pt.size = .5,label.size = 8,repel = TRUE)+NoLegend()

#Parsing Cluster 7 - NK cells
Idents(covid.integrated)<-covid.integrated$seurat_clusters_.5
cluster7<-subset(covid.integrated, idents = "7")
DimPlot(cluster7, reduction = "umap", group.by = "cell.type.fine",cells.highlight = list(colnames(cluster7[,cluster7$cell.type.fine == "NK"])))


#### initial monocyte exploration ####
#Parsing Cluster 8
Idents(covid.integrated)<-covid.integrated$seurat_clusters_.5
cluster8<-subset(covid.integrated, idents = "8")
DimPlot(cluster8, reduction = "umap", group.by = "cell.type.fine")
DimPlot(cluster8, reduction = "umap", group.by = "seurat_clusters")

cluster8markers<-FindMarkers(covid.integrated,ident.1 = "0",ident.2 = "8",assay = "SCT", 
            slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
VlnPlot(covid.integrated,features = c("IER5","FKBP4","HSPA1A"),group.by = "seurat_clusters_.5", pt.size = 0)

#Parsing Cluster 16
Idents(covid.integrated)<-covid.integrated$seurat_clusters_.5
cluster16<-subset(covid.integrated, idents = "16")
DimPlot(cluster16, reduction = "umap", group.by = "cell.type.fine")
DimPlot(cluster16, reduction = "umap", group.by = "seurat_clusters")

#Parsing Cluster 17
Idents(covid.integrated)<-covid.integrated$seurat_clusters_.5
cluster17<-subset(covid.integrated, idents = "17")
DimPlot(cluster17, reduction = "umap", group.by = "cell.type.fine")
DimPlot(cluster17, reduction = "umap", group.by = "seurat_clusters")

#Parsing Cluster 19
Idents(covid.integrated)<-covid.integrated$seurat_clusters_.5
cluster19<-subset(covid.integrated, idents = "19")
DimPlot(cluster19, reduction = "umap", group.by = "cell.type.fine")
DimPlot(cluster19, reduction = "umap", group.by = "seurat_clusters")
DimPlot(covid.integrated, reduction = "umap", group.by = "seurat_clusters_.5", pt.size = .5, label = TRUE)


##### further macrophage investigation - FCGR3 is CD16, subunits a and b ####
FeaturePlot(covid.integrated,features=c("FN1","MRC1"),pt.size = .1)
VlnPlot(covid.integrated,features = c("FCGR3A","CD14","FN1","MRC1"),group.by = "seurat_clusters_.5", pt.size = 0)

monocytemarkers<-list()

for(ii in 1:length(monocyteclusters)){
monocytemarkers[[ii]]<-FindMarkers(covid.integrated, ident.1 = monocyteclusters[ii], ident.2 = monocyteclusters[-ii], 
            features = VariableFeatures(covid.integrated), assay = "SCT", 
            slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")
}

#checking cd14
names(monocytemarkers)<-monocyteclusters
cd14groups<-grep("CD14",sapply(monocytemarkers,rownames))

for(ii in cd14groups){
  print(names(monocytemarkers)[ii])
  print(monocytemarkers[[ii]][grep("CD14",rownames(monocytemarkers[[ii]])),])

  }

#checking cd18
grep("FCGR",sapply(monocytemarkers,rownames))

markers1.20<-FindMarkers(covid.integrated, ident.1 = c(1,20), ident.2 = monocyteclusters[!(monocyteclusters %in% c(1,20))], 
            features = VariableFeatures(covid.integrated), assay = "SCT", 
            slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

markers1vs20<-FindMarkers(covid.integrated, ident.1 = 1, ident.2 = 20, 
                         features = VariableFeatures(covid.integrated), assay = "SCT", 
                         slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

##S100high and FCGR3A low - highly inflammatory phenotype
markers6<-FindMarkers(covid.integrated, ident.1 = 6, ident.2 = monocyteclusters[-6], features = VariableFeatures(covid.integrated), assay = "SCT", 
            slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

##
markers046<-FindMarkers(covid.integrated, ident.1 = c(0,4,6), ident.2 = monocyteclusters[!(monocyteclusters %in% c(0,4,6))], features = VariableFeatures(covid.integrated), assay = "SCT", 
                      slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

##other groups
needreview<-c(0,1,17,19)
monocytemarkers$'1'[abs(monocytemarkers$'1'$avg_logFC)>=.4,]

monocytemarkers$'17'[abs(monocytemarkers$'17'$avg_logFC)>=.4,]
markers17<-FindMarkers(covid.integrated, ident.1 = 17, features = VariableFeatures(covid.integrated), assay = "SCT", 
            slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

#FN1 and MCR1 are highly specific canonical markers for M2 macrophages
monocytemarkers$'19'[abs(monocytemarkers$'19'$avg_logFC)>=.4,]

##### Merging better defined clusters into old clusters as letters #####
clusterstomerge<-c(12,3,7,30,26,22,5,14,27,15,21,28,25) ###from 1 resolution 30cluster set
merger<-covid.integrated$seurat_clusters %in% clusterstomerge
seurat_clusters.merge<-as.character(covid.integrated$seurat_clusters_.5)
seurat_clusters.merge[merger]<- paste0(covid.integrated$seurat_clusters[merger],"b")
covid.integrated$seurat_clusters.merge<-as.factor(seurat_clusters.merge)

DimPlot(covid.integrated, reduction = "umap", group.by = "seurat_clusters.merge",pt.size = .5,label = TRUE,label.size = 8, repel = TRUE)
mergetable<-table(covid.integrated$seurat_clusters.merge)
orphans<-mergetable[mergetable<300]

# #combine orphan clusters based on original assignments
# ##23
# DimPlot(covid.integrated, reduction = "umap", group.by = "seurat_clusters.merge",pt.size = .5,label = TRUE,label.size = 8, 
#         repel = TRUE, cells.highlight = list(colnames(covid.integrated[,covid.integrated$seurat_clusters.merge == "23"])))
# covid.integrated$seurat_clusters.merge<-recode(covid.integrated$seurat_clusters.merge, "23" = "20")
options(mc.cores = 1)

##combine using integration 
orphanselect<-vector(mode = "numeric",length = dim(covid.integrated)[2])
orphanselect[covid.integrated$seurat_clusters.merge %in% names(orphans)]<-1
covid.integrated$orphans<-as.factor(orphanselect)

qsaveworkspace()

rm(list = ls()[! (ls() %in% c('covid.integrated'))] )
covidorphans<-SplitObject(covid.integrated,split.by = "orphans")
names(covidorphans)<-c("main","orphan")
DimPlot(covidorphans$main, reduction = "umap", group.by = "seurat_clusters.merge",pt.size = .5,label = TRUE,label.size = 8, repel = TRUE)
DimPlot(covidorphans$orphan, reduction = "umap", group.by = "seurat_clusters.merge",pt.size = .5,label = TRUE,label.size = 8, repel = TRUE)

rm(covid.integrated)

DefaultAssay(covidorphans$main)<-"RNA"
DefaultAssay(covidorphans$orphan)<-"RNA"

covidorphans$main <- SCTransform(covidorphans$main, vars.to.regress = "percent.mt", verbose = TRUE)
covidorphans$orphan <- SCTransform(covidorphans$orphan, vars.to.regress = "percent.mt", verbose = FALSE)

#run on commandline
options(mc.cores = 4)
# #SCT preintegration
orphan.transfer.features <- SelectIntegrationFeatures(object.list = covidorphans, nfeatures = 3000)
covidorphans <- PrepSCTIntegration(object.list = covidorphans, anchor.features = orphan.transfer.features, 
                                  verbose = FALSE)

##Seurat label transfer
orphan.anchors <- FindTransferAnchors(reference = covidorphans$main, query = covidorphans$orphan, 
                                          dims = 1:30, project.query = TRUE)
orphanpredic <- TransferData(anchorset = orphan.anchors, refdata = covidorphans$main$seurat_clusters.merge, 
                            dims = 1:30)

##reload and apply predictions
orphanpredic<-qread("orphanpredic")
qloadworkspace()
covid.integrated$seurat_clusters.merge[rownames(orphanpredic)]<-orphanpredic$predicted.id

DimPlot(covid.integrated, reduction = "umap", group.by = c("seurat_clusters.merge"),pt.size = .5,label = TRUE,label.size = 3, 
                repel = TRUE)+NoLegend()

##### markers and identification for remaining clusters #####

##signature matrix from Monaco 2019 normalized deconvolution
sigmatrixRNAseq <- read.delim("~/Documents/COVID_Data_Mine/sigmatrixRNAseq.txt")

Idents(covid.integrated)<-covid.integrated$seurat_clusters.merge
options(mc.cores = 8)
merge.markers <- FindAllMarkers(covid.integrated, assay = "SCT", 
                                           slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")
merge.markers<-dlply(merge.markers, .(cluster))

tclusters<-c("3b","5b","7b","12b","14b","16","21b","22b","25b","26b")

tplots<-list()
for(ii in 1:length(tclusters)){
tplots[[ii]]<-DimPlot(covid.integrated, reduction = "umap", group.by = "cell.type.merge",cells.highlight = list(colnames(covid.integrated[,covid.integrated$seurat_clusters.merge %in% tclusters[ii]])))
}

clusterplots<-list()
for(ii in 1:length(levels(covid.integrated$cell.type.merge))){
  clusterplots[[ii]]<-DimPlot(covid.integrated, reduction = "umap", group.by = "cell.type.merge",cells.highlight = list(colnames(covid.integrated[,covid.integrated$cell.type.merge %in% levels(covid.integrated$cell.type.merge)[ii]])))
}

#check 25b and 3b and 5b

##22b - Tregs - canonical markers IL2RA and CTLA4, 
markers22b<-FindMarkers(covid.integrated, ident.1 = "22b", features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

table(subset(covid.integrated,idents = "22b")$predicted.id)
VlnPlot(covid.integrated,features = c("LAG3","IL2RA","CTLA4"),group.by = "seurat_clusters.merge", pt.size = 0)


##3b - γδ T cell - non-vd2 by CD8A-, GZMH+ SPON2+, renamed as NK cells, not distinct from other NK cluster
markers3b<-FindMarkers(covid.integrated, ident.1 = c("3b"), ident.2 = tclusters[!(tclusters %in% c("3b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

table(subset(covid.integrated,idents = "3b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("GNLY","SPON2","CD8A","TRDC","KLRF1","NCAM1"),group.by = "seurat_clusters.merge", pt.size = 0)
#subsequent marker analysis with NK is too close to gamma delta cells here. naming is preserved for now, but will be set to NK and regenerated later


##5b - CD4+ T cells TYROBP+, CST3+
markers5b<-FindMarkers(covid.integrated, ident.1 = c("5b"), ident.2 = tclusters[!(tclusters %in% c("5b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                       slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers5b14<-FindMarkers(covid.integrated, ident.1 = c("5b"), ident.2 = "14b",features = VariableFeatures(covid.integrated), assay = "SCT", 
                       slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")
markers5b14$diff<-markers5b14$pct.1-markers5b14$pct.2
arrange(markers5b14, desc(diff))

table(subset(covid.integrated,idents = "5b")$cell.type.fine)
table(subset(covid.integrated,idents = "14b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("CD4","MAL"),group.by = "seurat_clusters.merge", pt.size = 0)
FeaturePlot(covid.integrated,features=c("MAL"),pt.size = .1)

##7b CD8+ memory t cells CCL5+/GZMH+
markers7b<-FindMarkers(covid.integrated, ident.1 = c("7b"), ident.2 = tclusters[!(tclusters %in% c("7b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

table(subset(covid.integrated,idents = "7b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("CCL5","GZMH","CD8A","CD8B"),group.by = "seurat_clusters.merge", pt.size = 0)



##12b - NK cells -GNLY+,NKG7+
markers12b<-FindMarkers(covid.integrated, ident.1 = c("12b"), ident.2 = tclusters[!(tclusters %in% c("12b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                         slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

table(subset(covid.integrated,idents = "12b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("GNLY","NKG7","CD8A","TRDC","KLRF1","NCAM1"),group.by = "seurat_clusters.merge", pt.size = 0)

##14b - CD4+ T cells, LTB+ TRAT1+
markers14b<-FindMarkers(covid.integrated, ident.1 = c("14b"), ident.2 = tclusters[!(tclusters %in% c("14b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers14b5<-FindMarkers(covid.integrated, ident.1 = c("14b"), ident.2 = "5b",features = VariableFeatures(covid.integrated), assay = "SCT", 
                         slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

table(subset(covid.integrated,idents = "14b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("LTB","CD7","TRAT1","LEF1"),group.by = "seurat_clusters.merge", pt.size = 0)


##16 - CD8+ T cells HIST1H4C+ MKI67+
markers16<-FindMarkers(covid.integrated, ident.1 = c("16"), ident.2 = tclusters[!(tclusters %in% c("16"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers16<-FindMarkers(covid.integrated, ident.1 = c("16"),features = VariableFeatures(covid.integrated), assay = "SCT", 
                       slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

sig16<-sigmatrixRNAseq[rownames(sigmatrixRNAseq) %in% rownames(markers16),]
colSums(sig16*markers16$avg_logFC[rownames(markers16) %in% rownames(sig16)])

table(subset(covid.integrated,idents = "16")$cell.type.fine)
VlnPlot(covid.integrated,features = c("HIST1H4C","STMN1","HMGB2","MKI67"),group.by = "seurat_clusters.merge", pt.size = 0)


##21b - Plasmablasts 
markers21b<-FindMarkers(covid.integrated, ident.1 = c("21b"), ident.2 = tclusters[!(tclusters %in% c("21b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers21b15<-FindMarkers(covid.integrated, ident.1 = c("21b"), ident.2 = "15",features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")


sig21b<-sigmatrixRNAseq[rownames(sigmatrixRNAseq) %in% rownames(markers21b),]

table(subset(covid.integrated,idents = "21b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("FTH1","IGHM"),group.by = "seurat_clusters.merge", pt.size = 0)


##25b Monocytes - MALAT1+ - SLE related inflammatory RNA according to Yang et al 2017 in oncotarget
markers25b<-FindMarkers(covid.integrated, ident.1 = c("25b"), ident.2 = monocyteclusters,features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers25b<-FindMarkers(covid.integrated, ident.1 = c("25b"),features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")


sig25b<-sigmatrixRNAseq[rownames(sigmatrixRNAseq) %in% rownames(markers25b),]
colSums(sig25b*markers25b$avg_logFC[rownames(markers25b) %in% rownames(sig25b)])

table(subset(covid.integrated,idents = "25b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("MALAT1","NKTR","TTN","TCF7"),group.by = "seurat_clusters.merge", pt.size = 0)


##26b - CD8+ Memory T cells LTB+
markers26b<-FindMarkers(covid.integrated, ident.1 = c("26b"), ident.2 = tclusters[!(tclusters %in% c("26b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers26b<-FindMarkers(covid.integrated, ident.1 = c("26b"), ident.2 = "7b",features = VariableFeatures(covid.integrated), assay = "SCT", 
                          slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")


sig26b<-sigmatrixRNAseq[rownames(sigmatrixRNAseq) %in% rownames(markers26b),]
colSums(sig26b*markers26b$avg_logFC[rownames(markers26b) %in% rownames(sig26b)])

table(subset(covid.integrated,idents = "26b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("CCL5","CD8A"),group.by = "seurat_clusters.merge", pt.size = 0)

###B and other cells checking
bclusters<-c("15","15b","21b","27b")

##15 - Plasmablasts MZB1+,XBP1+
markers15<-FindMarkers(covid.integrated, ident.1 = c("15"), ident.2 = bclusters[!(bclusters %in% c("15"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers15<-FindMarkers(covid.integrated, ident.1 = c("15"),features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")


sig15<-sigmatrixRNAseq[rownames(sigmatrixRNAseq) %in% rownames(markers15),]
colSums(sig15*markers15$avg_logFC[rownames(markers15) %in% rownames(sig15)])

table(subset(covid.integrated,idents = "15")$cell.type.fine)
VlnPlot(covid.integrated,features = c("IGHG4","IGHG1","MZB1","XBP1"),group.by = "seurat_clusters.merge", pt.size = 0)

##15b - MS4A1+/CD79A+ Naive B Cells
markers15b<-FindMarkers(covid.integrated, ident.1 = c("15b"), ident.2 = bclusters[!(bclusters %in% c("15b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                       slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers15b<-FindMarkers(covid.integrated, ident.1 = c("15b"),features = VariableFeatures(covid.integrated), assay = "SCT", 
                       slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")


sig15b<-sigmatrixRNAseq[rownames(sigmatrixRNAseq) %in% rownames(markers15b),]
colSums(sig15b*markers15b$avg_logFC[rownames(markers15b) %in% rownames(sig15b)])

table(subset(covid.integrated,idents = "15b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("MS4A1","CD79A"),group.by = "seurat_clusters.merge", pt.size = 0)

##21b - IGHM+ IGJ+ Plasmablasts
markers21b<-FindMarkers(covid.integrated, ident.1 = c("21b"), ident.2 = bclusters[!(bclusters %in% c("21b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")

markers21b<-FindMarkers(covid.integrated, ident.1 = c("21b"),features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")


sig21b<-sigmatrixRNAseq[rownames(sigmatrixRNAseq) %in% rownames(markers21b),]
colSums(sig21b*markers21b$avg_logFC[rownames(markers21b) %in% rownames(sig21b)])

table(subset(covid.integrated,idents = "21b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("IGHM","IGHG2","IGHG3"),group.by = "seurat_clusters.merge", pt.size = 0)

##27b - myeloid dendritic cells CD1C+/LGALS2+
markers27b<-FindMarkers(covid.integrated, ident.1 = c("27b"), ident.2 = bclusters[!(bclusters %in% c("27b"))],features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = FALSE, latent.vars = "nCount_RNA")

markers27b<-FindMarkers(covid.integrated, ident.1 = c("27b"),features = VariableFeatures(covid.integrated), assay = "SCT", 
                        slot = "data",test.use = "MAST",only.pos = TRUE, latent.vars = "nCount_RNA")


sig27b<-sigmatrixRNAseq[rownames(sigmatrixRNAseq) %in% rownames(markers27b),]
colSums(sig27b*markers27b$avg_logFC[rownames(markers27b) %in% rownames(sig27b)])

table(subset(covid.integrated,idents = "27b")$cell.type.fine)
VlnPlot(covid.integrated,features = c("CD1C","LGALS2","CD1E"),group.by = "seurat_clusters.merge", pt.size = 0)

##2 NK cell merged in final merge IDs

Idents(covid.integrated)<-covid.integrated$seurat_clusters.merge

covid.integrated <- RenameIdents(covid.integrated, `0` = "FCN1+/VCAN+ Monocytes", `1` = "APOE+/CCL18+/FCN1- Monocytes", `3b` = "GNLY+/NKG7+ NK Cells",
                                 `4` = "FCN1+/CCL4+ Monocytes", `5` = "FCN1+/CCL2/7/8+ Monocytes", `5b` = "CD4+/TYROBP+/CST3+ T Cells", `6` = "CD14+/FCN1+/VCAN++/S100A8++/FCGR3A- Monocytes", `7b` = "CD8+/CCL5+/GZHM+ Memory T Cells", `8` = "CD14-/FCN1+/FKBP4+/IER5+ Monocytes",
                                 `12` = "CD14-/FCN1+/CCL2+ Monocytes",`12b` = "GNLY+/NKG7+ NK Cells", `13` = "Epithelium/Basal Cells", `14` = "Neutrophils", `14b` ="CD4+/LTB+/TRAT1+ T Cells",
                                 `15` = "MZB1+/XBP1+ Plasmablasts", `15b` ="MS4A1+/CD79A+ Naive B cells", `16` = "CD8+/MKI67+/HIST1H4C+ T Cells", `17` = "MT+ Monocytes", `19` = "FN1+/MRC1+/FCN1- Monocytes", `20` = "FN1+/CD14-/CCL18++/FCN1- Monocytes",
                                 `21b` = "IGJ+/IGHM+ Plasmablasts", `22` = "Plasmacytoid Dendritic Cells", `22b` = "IL32+/LAG3+ T Regulatory Cells", `25b` = "MALAT1+ Monocytes",`26b` = "CD8+/LTB+ Memory T Cells",`27b` = "Myeloid Dendritic Cells",
                                 `28b` = "Epithelium/Pneumocytes/Ciliary Cells")

covid.integrated$cell.type.merge<- Idents(covid.integrated)

DimPlot(covid.integrated, reduction = "umap", group.by = "cell.type.merge",pt.size = .5,label = TRUE,label.size = 5, 
        repel = TRUE)+NoLegend()

##### Checking Labels #####
DimPlot(covid.integrated, reduction = "umap", label = TRUE, pt.size = .5)+NoLegend()
#Plot Monocyte markers
FeaturePlot(covid.integrated,features=c("VCAN","FCN1","CCL2","CCL4","CCL7","FKBP4+","IER5","CD4","CCL18"),pt.size = .5)
VlnPlot(covid.integrated,features = c("VCAN","FCN1","CCL2","CCL4","CCL7","FKBP4","IER5","CD4","CCL18"),group.by = "seurat_clusters_.5", pt.size = 0)

##### cell proportions #####
#Note: removing M3 due to low cell count - full seurat object saved
qsave(covid.integrated,"covid.integrated.full",nthreads = 16)
Idents(covid.integrated)<-covid.integrated$aggpatient
covid.integrated<-subset(covid.integrated, idents = c("M3"), invert = TRUE) #remove M3

barplot(table(covid.integrated.full$aggpatient))

finetable<-table(as.data.frame(cbind(covid.integrated$cell.type.merge.coarse,covid.integrated$aggpatient)))
rownames(finetable)<-levels(covid.integrated$cell.type.merge)[as.integer(rownames(finetable))]
fineprop<-sweep(finetable,2,colSums(finetable),`/`)

coarsetable<-table(covid.integrated$cell.type.merge.coarse,covid.integrated$aggcondition.rn)
coarsetable<-as.data.frame.matrix(coarsetable)
coarsetable["Total",]<-colSums(coarsetable)
coarseprop<-as.data.frame(sweep(coarsetable,2,as.matrix(coarsetable["Total",]),"/")[1:15,])
coarsediff<-data.frame("M"=coarseprop$M-coarseprop$H,"S"=coarseprop$S-coarseprop$H,
                       "PBMC_M"=coarseprop$PBMC_M-coarseprop$PBMC_H,"PBMC_S"=coarseprop$PBMC_S-coarseprop$PBMC_H,
                       row.names = rownames(coarseprop))
coarsediffpct<-data.frame("M"=coarsediff$M/coarseprop$H,"S"=coarsediff$S/coarseprop$H,
                          "PBMC_M"=coarsediff$PBMC_M/coarseprop$PBMC_H,"PBMC_S"=coarsediff$PBMC_S/coarseprop$PBMC_H,rownames(coarsediff))

coarsematrix<-as.matrix(coarsetable)

coarseproptest<-list()
coarsepropsig<-list()
totals<-coarsetable["Total",]
for (ii in levels(covid.integrated$cell.type.merge.coarse)){
  props<-coarsetable[ii,]
  coarseproptest[[ii]]<-pairwise.prop.test(t(rbind(props,totals)), p.adjust.method = "BH", alternative = "two.sided",correct = FALSE)$p.value
  coarseproplist<-as.data.frame(coarseproptest[[ii]])
  coarseproplist$condition1<-rownames(coarseproplist)
  coarseproplist<-pivot_longer(coarseproplist,cols = c("BAL_H","BAL_M","BAL_S","PBMC_H","PBMC_M") )
  colnames(coarseproplist)<-c("condition1","condition2","adj_pval")
  coarseproplist<-coarseproplist[!(is.na(coarseproplist$adj_pval)),]
  print(ii)
  # print(coarseproplist[coarseproplist$adj_pval<.05,])
  # coarsepropsig[[ii]]<-coarseproplist[coarseproplist$adj_pval<.05,]
  coarsepropsig[[ii]]<-coarseproplist
  }



fineprop<-as.data.frame(fineprop)
fineprop$condition<-gsub("[0-9]","",fineprop$V2)

fineprop$condition <- factor(fineprop$condition,levels = c("H","M","S","PBMC_H","PBMC_M","PBMC_S"))

ggplot(fineprop, aes(fill=V1, y=Freq, x=condition))+ scale_fill_viridis_d() + 
  geom_bar(position="fill", stat="identity")

ggplot(subset(fineprop, fineprop$V1 %in% c("FCN1+/CCL4+ Monocytes","CD14+/FCN1+/VCAN++/S100A8++/FCGR3A- Monocytes",
                                           "Plasmacytoid Dendritic","FN1+/CD14-/CCL18++/FCN1- Monocytes")), 
       aes(fill=V1, y=Freq, x=condition))+ scale_fill_viridis_d() + geom_boxplot()+scale_y_sqrt()


#####  Differential Expression ##### initial code attempt #depreciated 
library(EnhancedVolcano)
options(mc.cores = 8)

celltypes<-levels(covid.integrated$cell.type.merge)
DEmetadata<-data.frame(BAL.sample = c("M-H","S-H","S-M"),PBMC.sample = c("PBMC_M-PBMC_H","PBMC_S-PBMC_H","PBMC_S-PBMC_M"),
                       p_val_adj = c("g2_1p_val_adj","g3_1p_val_adj","g3_2p_val_adj"), avg_logFC = c("avg_logFC_g2_1","avg_logFC_g3_1","avg_logFC_g3_2"))

# DE.bal<-group_split(DE.bal,cluster)
covid.DE.results<-list()
for (ii in celltypes){
modules<-as.data.frame(covid.module.results[[ii]]$modules)
rownames(modules)<-modules$id
subcells<-subset(covid.integrated, subset = cell.type.merge == ii)
Idents(subcells)<-subcells$aggcondition

DE.bal.combined<-FindDEseurat(subset(subcells,subset = aggcondition %in% c("H","M","S")),ident.1 = "H", ident.2 = "M",ident.3 = "S", assay = "SCT", 
                    slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
DE.bal.combined$module<-modules$module[rownames(DE.bal.combined)]

DE.pbmc.combined<-FindDEseurat(subset(subcells,subset = aggcondition %in% c("PBMC_H","PBMC_M","PBMC_S")),ident.1 = "PBMC_H", ident.2 = "PBMC_M",ident.3 = "PBMC_S", assay = "SCT", 
                               slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
DE.pbmc.combined$module<-modules$module[rownames(DE.pbmc.combined)]
  
volcano.plots<-list()
volcano.plots.modules<-list()
  for (jj in 1:3){
    bal.sample<-DEmetadata$BAL.sample[jj]
    pbmc.sample<-DEmetadata$PBMC.sample[jj]
    pval<-DEmetadata$p_val_adj[jj]
    logfc<-DEmetadata$avg_logFC[jj]
    
    volcano.plots[[bal.sample]]<-EnhancedVolcano(as.data.frame(DE.bal.combined),
                  lab = rownames(DE.bal.combined),
                  x = logfc, xlab = "ln fold change", 
                  legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
                  y = pval,
                  title = paste0("Differential Expression for ",ii),
                  subtitle = bal.sample,
                  FCcutoff = log(2),
                  pCutoff = 10e-6,
                  selectLab = rownames(DE.bal.combined)[rownames(DE.bal.combined) %in% bal.DE.cell$balDEs])
    volcano.plots[[pbmc.sample]]<-EnhancedVolcano(as.data.frame(DE.pbmc.combined),
                  lab = rownames(DE.pbmc.combined),
                  x = logfc, xlab = "ln fold change",
                  legendLabels = c("NS", "ln fold change > .4","p-val < .05","lnFC > .4 & p-val < .05"),
                  y = pval,
                  title = paste0("Differential Expression for ",ii),
                  subtitle = pbmc.sample,
                  FCcutoff = log(3)-log(2),
                  pCutoff = .05,
                  selectLab = rownames(DE.pbmc.combined)[rownames(DE.pbmc.combined) %in% pbmc.DE.cell$pbmcDEs])

    volcano.plots.modules[[bal.sample]]<-EnhancedVolcano(as.data.frame(DE.bal.combined),
                                                 lab = DE.bal.combined$module,
                                                 x = logfc, xlab = "ln fold change", 
                                                 legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
                                                 y = pval,
                                                 title = paste0("Differential Expression for ",ii),
                                                 subtitle = bal.sample,
                                                 FCcutoff = log(2),
                                                 pCutoff = 10e-6)
    volcano.plots.modules[[pbmc.sample]]<-EnhancedVolcano(as.data.frame(DE.pbmc.combined),
                                                 lab = DE.pbmc.combined$module,
                                                 x = logfc, xlab = "ln fold change",
                                                 legendLabels = c("NS", "ln fold change > .4","p-val < .05","lnFC > .4 & p-val < .05"),
                                                 y = pval,
                                                 title = paste0("Differential Expression for ",ii),
                                                 subtitle = pbmc.sample,
                                                 FCcutoff = log(3)-log(2),
                                                 pCutoff = .05)

    
    
  }

covid.DE.results[[ii]]<-list(DE.bal=DE.bal.combined,DE.pbmc=DE.pbmc.combined,volcano.plots=volcano.plots,volcano.plots.modules=volcano.plots.modules)
}

##### Differential Expression for coarse set ##### #actual code used for final dataset
library(EnhancedVolcano)
options(mc.cores = 8)

celltypes<-levels(covid.integrated$cell.type.merge.coarse)
DEmetadata<-data.frame(BAL.sample = c("M-H","S-H","S-M"),PBMC.sample = c("PBMC_M-PBMC_H","PBMC_S-PBMC_H","PBMC_S-PBMC_M"),
                       p_val_adj = c("g2_1p_val_adj","g3_1p_val_adj","g3_2p_val_adj"), avg_logFC = c("avg_logFC_g2_1","avg_logFC_g3_1","avg_logFC_g3_2"))

# DE.bal<-group_split(DE.bal,cluster)
covid.DE.results.coarse<-list()
for (ii in celltypes){
  # modules<-as.data.frame(covid.module.results[[ii]]$modules)
  # rownames(modules)<-modules$id
  subcells<-subset(covid.integrated, subset = cell.type.merge.coarse == ii)
  Idents(subcells)<-subcells$aggcondition
  
  DE.bal.combined<-FindDEseurat(subset(subcells,subset = aggcondition %in% c("H","M","S")),ident.1 = "H", ident.2 = "M",ident.3 = "S", assay = "SCT", 
                                slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
  # DE.bal.combined$module<-modules$module[rownames(DE.bal.combined)]
  
  DE.pbmc.combined<-FindDEseurat(subset(subcells,subset = aggcondition %in% c("PBMC_H","PBMC_M","PBMC_S")),ident.1 = "PBMC_H", ident.2 = "PBMC_M",ident.3 = "PBMC_S", assay = "SCT", 
                                 slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
  # DE.pbmc.combined$module<-modules$module[rownames(DE.pbmc.combined)]
  
  volcano.plots<-list()
  # volcano.plots.modules<-list()
  for (jj in 1:3){
    bal.sample<-DEmetadata$BAL.sample[jj]
    pbmc.sample<-DEmetadata$PBMC.sample[jj]
    pval<-DEmetadata$p_val_adj[jj]
    logfc<-DEmetadata$avg_logFC[jj]
    
    volcano.plots[[bal.sample]]<-EnhancedVolcano(as.data.frame(DE.bal.combined),
                                                 lab = rownames(DE.bal.combined),
                                                 x = logfc, xlab = "ln fold change", 
                                                 legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
                                                 y = pval,
                                                 title = paste0("Differential Expression for ",ii),
                                                 subtitle = bal.sample,
                                                 FCcutoff = log(2),
                                                 pCutoff = 10e-6,
                                                 selectLab = rownames(DE.bal.combined)[rownames(DE.bal.combined) %in% bal.DE.cell$balDEs])
    volcano.plots[[pbmc.sample]]<-EnhancedVolcano(as.data.frame(DE.pbmc.combined),
                                                  lab = rownames(DE.pbmc.combined),
                                                  x = logfc, xlab = "ln fold change",
                                                  legendLabels = c("NS", "ln fold change > .4","p-val < .05","lnFC > .4 & p-val < .05"),
                                                  y = pval,
                                                  title = paste0("Differential Expression for ",ii),
                                                  subtitle = pbmc.sample,
                                                  FCcutoff = log(3)-log(2),
                                                  pCutoff = .05,
                                                  selectLab = rownames(DE.pbmc.combined)[rownames(DE.pbmc.combined) %in% pbmc.DE.cell$pbmcDEs])
    
    # volcano.plots.modules[[bal.sample]]<-EnhancedVolcano(as.data.frame(DE.bal.combined),
    #                                                      lab = DE.bal.combined$module,
    #                                                      x = logfc, xlab = "ln fold change", 
    #                                                      legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
    #                                                      y = pval,
    #                                                      title = paste0("Differential Expression for ",ii),
    #                                                      subtitle = bal.sample,
    #                                                      FCcutoff = log(2),
    #                                                      pCutoff = 10e-6)
    # volcano.plots.modules[[pbmc.sample]]<-EnhancedVolcano(as.data.frame(DE.pbmc.combined),
    #                                                       lab = DE.pbmc.combined$module,
    #                                                       x = logfc, xlab = "ln fold change",
    #                                                       legendLabels = c("NS", "ln fold change > .4","p-val < .05","lnFC > .4 & p-val < .05"),
    #                                                       y = pval,
    #                                                       title = paste0("Differential Expression for ",ii),
    #                                                       subtitle = pbmc.sample,
    #                                                       FCcutoff = log(3)-log(2),
    #                                                       pCutoff = .05)
    
    
    
  }
  
  covid.DE.results.coarse[[ii]]<-list(DE.bal=DE.bal.combined,DE.pbmc=DE.pbmc.combined,volcano.plots=volcano.plots)
}
##### Differential Expression for subgroup set #####
library(EnhancedVolcano)
options(mc.cores = 8)

celltypes<-levels(covid.integrated$cell.type.merge.subgroup)
DEmetadata<-data.frame(BAL.sample = c("M-H","S-H","S-M"),PBMC.sample = c("PBMC_M-PBMC_H","PBMC_S-PBMC_H","PBMC_S-PBMC_M"),
                       p_val_adj = c("g2_1p_val_adj","g3_1p_val_adj","g3_2p_val_adj"), avg_logFC = c("avg_logFC_g2_1","avg_logFC_g3_1","avg_logFC_g3_2"))

# DE.bal<-group_split(DE.bal,cluster)
covid.DE.results.subgroup<-list()
for (ii in celltypes){
  # modules<-as.data.frame(covid.module.results[[ii]]$modules)
  # rownames(modules)<-modules$id
  subcells<-subset(covid.integrated, subset = cell.type.merge.subgroup == ii)
  Idents(subcells)<-subcells$aggcondition
  
  DE.bal.combined<-FindDEseurat(subset(subcells,subset = aggcondition %in% c("H","M","S")),ident.1 = "H", ident.2 = "M",ident.3 = "S", assay = "SCT", 
                                slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
  # DE.bal.combined$module<-modules$module[rownames(DE.bal.combined)]
  
  DE.pbmc.combined<-FindDEseurat(subset(subcells,subset = aggcondition %in% c("PBMC_H","PBMC_M","PBMC_S")),ident.1 = "PBMC_H", ident.2 = "PBMC_M",ident.3 = "PBMC_S", assay = "SCT", 
                                 slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
  # DE.pbmc.combined$module<-modules$module[rownames(DE.pbmc.combined)]
  
  volcano.plots<-list()
  # volcano.plots.modules<-list()
  for (jj in 1:3){
    bal.sample<-DEmetadata$BAL.sample[jj]
    pbmc.sample<-DEmetadata$PBMC.sample[jj]
    pval<-DEmetadata$p_val_adj[jj]
    logfc<-DEmetadata$avg_logFC[jj]
    
    volcano.plots[[bal.sample]]<-EnhancedVolcano(as.data.frame(DE.bal.combined),
                                                 lab = rownames(DE.bal.combined),
                                                 x = logfc, xlab = "ln fold change", 
                                                 legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
                                                 y = pval,
                                                 title = paste0("Differential Expression for ",ii),
                                                 subtitle = bal.sample,
                                                 FCcutoff = log(2),
                                                 pCutoff = 10e-6,
                                                 selectLab = rownames(DE.bal.combined)[rownames(DE.bal.combined) %in% bal.DE.cell$balDEs])
    volcano.plots[[pbmc.sample]]<-EnhancedVolcano(as.data.frame(DE.pbmc.combined),
                                                  lab = rownames(DE.pbmc.combined),
                                                  x = logfc, xlab = "ln fold change",
                                                  legendLabels = c("NS", "ln fold change > .4","p-val < .05","lnFC > .4 & p-val < .05"),
                                                  y = pval,
                                                  title = paste0("Differential Expression for ",ii),
                                                  subtitle = pbmc.sample,
                                                  FCcutoff = log(3)-log(2),
                                                  pCutoff = .05,
                                                  selectLab = rownames(DE.pbmc.combined)[rownames(DE.pbmc.combined) %in% pbmc.DE.cell$pbmcDEs])
    
    # volcano.plots.modules[[bal.sample]]<-EnhancedVolcano(as.data.frame(DE.bal.combined),
    #                                                      lab = DE.bal.combined$module,
    #                                                      x = logfc, xlab = "ln fold change", 
    #                                                      legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
    #                                                      y = pval,
    #                                                      title = paste0("Differential Expression for ",ii),
    #                                                      subtitle = bal.sample,
    #                                                      FCcutoff = log(2),
    #                                                      pCutoff = 10e-6)
    # volcano.plots.modules[[pbmc.sample]]<-EnhancedVolcano(as.data.frame(DE.pbmc.combined),
    #                                                       lab = DE.pbmc.combined$module,
    #                                                       x = logfc, xlab = "ln fold change",
    #                                                       legendLabels = c("NS", "ln fold change > .4","p-val < .05","lnFC > .4 & p-val < .05"),
    #                                                       y = pval,
    #                                                       title = paste0("Differential Expression for ",ii),
    #                                                       subtitle = pbmc.sample,
    #                                                       FCcutoff = log(3)-log(2),
    #                                                       pCutoff = .05)
    
    
    
  }
  
  covid.DE.results.subgroup[[ii]]<-list(DE.bal=DE.bal.combined,DE.pbmc=DE.pbmc.combined,volcano.plots=volcano.plots)
}




##### count DE occurrences in S-M cases #####
covid.DE.results<-qread("covid.DE.results.BH")



baltable<-table(balgenes)
bal.sig.genes<-arrange(as.data.frame(baltable),desc(Freq))[1:50,]


pbmctable<-table(pbmcgenes)
pbmc.sig.genes<-arrange(as.data.frame(pbmctable),desc(Freq))[1:50,]

DoHeatmap(subset(covid.integrated, downsample = 1000), group.by= "aggcondition",features =DE.rec.genes.var, size = 1)

##### counts in coarse cases #####


balgenescoarse<-NULL
balallgenes<-NULL
pbmcgenescoarse<-NULL
pbmcallgenes<-NULL
for (ii in names(covid.DE.results.coarse)){
  balpval<-summary(rowSums(covid.DE.results.coarse[[ii]]$DE.bal[,c("g2_1p_val_adj","g3_1p_val_adj","g3_2p_val_adj")]))
  pbmcpval<-10^-27.361
  
  balgenescoarse<-c(balgenescoarse,rownames(filter(covid.DE.results.coarse[[ii]]$DE.bal,g3_2p_val_adj <1e-7,g2_1p_val_adj <1e-7,g3_1p_val_adj <1e-7)))
  pbmcgenescoarse<-c(pbmcgenescoarse,rownames(filter(covid.DE.results.coarse[[ii]]$DE.pbmc,g3_2p_val_adj <.05,g2_1p_val_adj <.05,g3_1p_val_adj <.05)))
  balallgenes<-c(balallgenes,rownames(covid.DE.results.coarse[[ii]]$DE.bal))
  pbmcallgenes<-c(pbmcallgenes,rownames(covid.DE.results.coarse[[ii]]$DE.pbmc))
  
  if (sum(rownames(filter(covid.DE.results.coarse[[ii]]$DE.bal,g3_2p_val_adj <1e-6,g2_1p_val_adj <1e-6,g3_1p_val_adj <1e-6)) %in% "NEAT1")){
    print(ii)
  }
}

length(balgenescoarse)/length(balallgenes)
length(pbmcgenescoarse)/length(pbmcallgenes)

length(pbmcgenescoarse)/length(balgenescoarse)

baltablecoarse<-table(balgenescoarse)
bal.sig.genes.coarse<-arrange(as.data.frame(baltablecoarse),desc(Freq))[1:91,]


pbmctablecoarse<-table(pbmcgenescoarse)
pbmc.sig.genes.coarse<-arrange(as.data.frame(pbmctablecoarse),desc(Freq))[1:17,]

##### testing ##### NEED TO REMOVE genes with extremely high residual variance and rerun analyses #####
for (ii in levels(covid.integrated$cell.type.merge.subgroup)){
  print(ii)
  print(covid.DE.results.subgroup[[ii]]$DE.bal[c("IGLC3","IGHG1","MALAT1","NEAT1"),c("pct.1","pct.2","pct.3")])
  print(covid.DE.results.subgroup[[ii]]$DE.pbmc[c("IGLC3","IGHG1","MALAT1","NEAT1"),c("pct.1","pct.2","pct.3")])
}

varplot<-VariableFeaturePlot(covid.integrated)
varplot<-LabelPoints(varplot,points = DE.rec.genes.var)
varplot + scale_y_log10()
vargenedata<-as.data.frame(varplot$data)
vargenedata<-arrange(vargenedata,desc(residual_variance))
plot(vargenedata$residual_variance[1:100])

RidgePlot(covid.integrated,features = rownames(vargenedata)[1:9],group.by = "aggcondition",ncol = 3)
RidgePlot(covid.integrated,features = c("ISG20","IL1RN","HLA-DQB1","CD68","APOE","NUPR1","NEXN","IFITM1","OASL"),group.by = "aggcondition",ncol = 3)

genestoremove<-rownames(vargenedata)[1:21]

filteredvargenes<- covid.integrated@assays$SCT@var.features[!(covid.integrated@assays$SCT@var.features %in% genestoremove)]

##### Final list of genes of interest for module building #####
# DE.rec.genes<-c(as.character(pbmc.sig.genes$pbmcgenes),as.character(bal.sig.genes$balgenes[!(bal.sig.genes$balgenes %in% pbmc.sig.genes$pbmcgenes)]))
DE.rec.genes<-c(as.character(pbmc.sig.genes.coarse$pbmcgenes),as.character(bal.sig.genes.coarse$balgenes[!(bal.sig.genes.coarse$balgenes %in% pbmc.sig.genes.coarse$pbmcgenes)]))
DE.rec.genes.var<-DE.rec.genes[DE.rec.genes %in% filteredvargenes]
bal.sig.genes.coarse[bal.sig.genes.coarse$balgenescoarse %in% filteredvargenes,]
pbmc.sig.genes.coarse[pbmc.sig.genes.coarse$pbmcgenescoarse %in% filteredvargenes,]

DE.rec.cds<-covid.SCT.cds[rowData(covid.SCT.cds)$gene_short_name %in% DE.rec.genes.var,]
# DE.rec.modules<-moduletest(DE.rec.cds)

gene_module_df <- find_gene_modules(DE.rec.cds,cores =16,max_components =30, resolution = .8,random_seed = 1,leiden_iter = 100,k=14,verbose = TRUE)
agg_mat_cell <- as.data.frame(aggregate_gene_expression(DE.rec.cds, gene_module_df,scale_agg_values = FALSE,norm_method = "log"), "dgCMatrix")
row.names(agg_mat_cell) <- stringr::str_c("Module ", row.names(agg_mat_cell))
# agg_mat_cell <- t(scale(t(agg_mat_cell),center = FALSE))
# pheatmap::pheatmap(agg_mat_cell,
#                                  scale="column", clustering_method="ward.D2")

# aggcondition[aggcondition == "C"] <- "S"

metadata<-list()
##how many are there?
metadata$patientcellcount<-table(DE.rec.cds@colData$aggpatient)
metadata$conditioncellcount<-table(DE.rec.cds@colData$aggcondition)
metadata$colData<-colData(DE.rec.cds)

#module expression heatmap
cell_group_df <- tibble::tibble(cell=row.names(colData(DE.rec.cds)), 
                                cell_group=DE.rec.cds@colData$aggcondition)
agg_mat_gene <- aggregate_gene_expression(DE.rec.cds,gene_group_df = NULL,cell_group_df = cell_group_df)
conditionmap<-pheatmap::pheatmap(agg_mat_gene,
                                 scale="column", clustering_method="ward.D2")

cell_group_df <- tibble::tibble(cell=row.names(colData(DE.rec.cds)), 
                                cell_group=DE.rec.cds@colData$aggcondition)
agg_mat <- aggregate_gene_expression(DE.rec.cds,gene_group_df = gene_module_df,cell_group_df = cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
conditionmap<-pheatmap::pheatmap(agg_mat,
                                 scale="column", clustering_method="ward.D2")


cell_group_df <- tibble::tibble(cell=row.names(colData(DE.rec.cds)), 
                                cell_group=DE.rec.cds@colData$aggpatient)
agg_mat_patient <- aggregate_gene_expression(DE.rec.cds, gene_module_df, cell_group_df)
row.names(agg_mat_patient) <- stringr::str_c("Module ", row.names(agg_mat_patient))
patientmap<-pheatmap::pheatmap(agg_mat_patient,
                               scale="column", clustering_method="ward.D2")

agg_mat_patient_noscale <- aggregate_gene_expression(DE.rec.cds, gene_module_df, cell_group_df,scale_agg_values = FALSE)

#ANOVA
patient2condition<- gsub("[0-9]","",colnames(agg_mat_patient))


moduleanova<-list()
tukey<-list()

for (ii in 1:length(levels(gene_module_df$module))){
  
  module<-as.data.frame(cbind(agg_mat_patient[ii,],patient2condition))
  colnames(module)<-c("exp","condition")
  module$exp<-as.numeric(module$exp)
  modulePBMC<-module[grepl("PBMC",module$condition),]
  module<-module[!grepl("PBMC",module$condition),]
  
  #logFC
  balmean<-aggregate(module$exp, by = list(module$condition), FUN = mean)
  pbmcmean<-aggregate(modulePBMC$exp, by = list(modulePBMC$condition), FUN = mean)
  rownames(balmean)<-balmean$Group.1
  rownames(pbmcmean)<-pbmcmean$Group.1
  
  
  # Compute the analysis of variance
  BAL.aov <- aov(exp ~ condition, data = module)
  PBMC.aov <- aov(exp ~ condition, data = modulePBMC)
  # PostHoc via Tukey
  BAL.ph<-TukeyHSD(BAL.aov)
  PBMC.ph<-TukeyHSD(PBMC.aov)
  tukeyres<-as.data.frame(rbind(BAL.ph$condition,PBMC.ph$condition))
  
  ##logFC
  module<-as.data.frame(cbind(agg_mat_patient[ii,],patient2condition))
  colnames(module)<-c("exp","condition")
  module$exp<-as.numeric(module$exp)
  
  #logFC
  logfc<-data.frame(rep(0,6),row.names = rownames(tukeyres))
  colnames(logfc)<-"fc"
  means<-aggregate(module$exp, by = list(module$condition), FUN = mean)
  rownames(means)<-means$Group.1
  means$log<-log(means$x+1)
  
  logfc["M-H","fc"]<-means["M","log"]-means["H","log"]
  logfc["S-H","fc"]<-means["S","log"]-means["H","log"]
  logfc["S-M","fc"]<-means["S","log"]-means["M","log"]
  
  logfc["PBMC_M-PBMC_H","fc"]<-means["PBMC_M","log"]-means["PBMC_H","log"]
  logfc["PBMC_S-PBMC_H","fc"]<-means["PBMC_S","log"]-means["PBMC_H","log"]
  logfc["PBMC_S-PBMC_M","fc"]<-means["PBMC_S","log"]-means["PBMC_M","log"]
  tukeyres$avg_logFC<-logfc[rownames(tukeyres),"fc"]
  
  moduleanova[[ii]]<-list(BAL.aov,PBMC.aov)
  
  
  tukey[[ii]]<-tukeyres
}


#topGO module analysis


module_enrichment<-list()

for (ii in 1:length(levels(gene_module_df$module))){
  
  # genelist <- numeric(length = dim(gene_module_df)[1]) 
  # names(genelist) <- gene_module_df$id
  # 
  # genelist[gene_module_df[gene_module_df$module == ii,1]$id]<-1
  
  
  genelist <- numeric(length = length(covid.integrated@assays$SCT@var.features))
  names(genelist) <- covid.integrated@assays$SCT@var.features

  genelist[gene_module_df[gene_module_df$module == ii,1]$id]<-1
  
  
  
  GOdata <- new("topGOdata", description = paste0("Module ",ii), ontology = "BP",
                allGenes = factor(genelist), 
                nodeSize = 10,
                annot = annFUN.org, mapping = "org.Hs.eg.db", ID="symbol")  
  
  # Test for enrichment using Fisher's Exact Test
  resultFisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")  
  module_enrichment[[ii]] <- GenTable(GOdata, Fisher = resultFisher, topNodes = 20, numChar = 60)
  names(module_enrichment)[ii] <- paste0("Module ",ii)
}

results<-list()
results[["module_enrichment"]]<-module_enrichment
results[["module_anova"]]<-moduleanova
results[["modules"]]<-gene_module_df
results[["condition_heatmap"]]<-conditionmap
results[["patient_heatmap"]]<-patientmap
names(results$module_anova)<-names(results$module_enrichment)
results[["tukey"]]<-tukey
results[["sigtukey"]]<-map(tukey, ~ .x %>% filter(`p adj` < .05))
results[["metadata"]]<-metadata
results[["aggregates"]]<-list(cell = agg_mat_cell, condition = agg_mat, patient = agg_mat_patient)
DEG.module.results<-results

bal.sig.genes.coarse[bal.sig.genes.coarse$balgenescoarse %in% DEG.module.results$modules$id,]
pbmc.sig.genes.coarse[pbmc.sig.genes.coarse$pbmcgenescoarse %in% DEG.module.results$modules$id,]

##### Automated Module testing with 4 top gene modules for all cell types #####
celltypes<-gsub("\\+","\\\\+",levels(covid.SCT.cds@colData$cell.type.merge))

covid.moduleDE.results<-list()  

for (ii in celltypes){
  subset.cds<- covid.SCT.cds[,grepl(ii, colData(covid.SCT.cds)$cell.type.merge, ignore.case=TRUE)]

  covid.moduleDE.results[[gsub("\\\\","\\",ii)]]<-moduleDE(subset.cds,DEG.module.results$modules)
  
  print(ii)
  print(covid.moduleDE.results[[ii]]$sigtukey)
}

##### rDEGs extracted per cell type and consolidated #####

rDEG.all<-list(bal = NULL, pbmc = NULL)
for (ii in names(covid.DE.results.coarse)){

    celltype<-ii
    
    balgenes<-rownames(filter(covid.DE.results.coarse[[ii]]$DE.bal,g3_2p_val_adj <1e-7,g2_1p_val_adj <1e-7,g3_1p_val_adj <1e-7))
    balgenes<-balgenes[balgenes %in% DEG.module.results$modules$id]
    balDEs<-covid.DE.results.coarse[[ii]]$DE.bal[balgenes,]
    balDEs$gene<-rownames(balDEs)
    rownames(balDEs)<-NULL
    
    if (dim(balDEs)[1]>0){
    balDEs<-cbind(balDEs,celltype)
    rDEG.all$bal<-rbind(rDEG.all$bal,balDEs)
    }
    
    pbmcgenes<-rownames(filter(covid.DE.results.coarse[[ii]]$DE.pbmc,g3_2p_val_adj <.05,g2_1p_val_adj <.05,g3_1p_val_adj <.05))
    pbmcgenes<-pbmcgenes[pbmcgenes %in% DEG.module.results$modules$id]
    pbmcDEs<-covid.DE.results.coarse[[ii]]$DE.pbmc[pbmcgenes,]
    pbmcDEs$gene<-rownames(pbmcDEs)
    rownames(pbmcDEs)<-NULL
    
    if (dim(pbmcDEs)[1]>0){
    pbmcDEs<-cbind(pbmcDEs,celltype)
    rDEG.all$pbmc<-rbind(rDEG.all$pbmc,pbmcDEs)
    }
  
}

rDEG.all$combined<-cbind(rDEG.all$bal,sample="BAL")
rDEG.all$combined<-rbind(rDEG.all$combined,cbind(rDEG.all$pbmc,sample="PBMC"))

##### number of DEG in the simplified set in each cell group #####

bal.DE.cell<-data.frame()
pbmc.DE.cell<-data.frame()
for (ii in names(covid.DE.results)){
  balDEs<-DE.rec.genes.var[DE.rec.genes.var %in% rownames(filter(covid.DE.results[[ii]]$DE.bal,g3_2p_val_adj <1e-6,g2_1p_val_adj <1e-6,g3_1p_val_adj <1e-6))]
  pbmcDEs<-DE.rec.genes.var[DE.rec.genes.var %in% rownames(filter(covid.DE.results[[ii]]$DE.pbmc,g3_2p_val_adj <.05,g2_1p_val_adj <.05,g3_1p_val_adj <.05))]
  
  if(!(is_empty(balDEs))){
  balDEs<-cbind(balDEs,celltype = ii)
  bal.DE.cell<-rbind(bal.DE.cell,balDEs)
  }
  
  if(!(is_empty(pbmcDEs))){
  pbmcDEs<-cbind(pbmcDEs,celltype = ii)
  pbmc.DE.cell<-rbind(pbmc.DE.cell,pbmcDEs)
  }
}

arrange(as.data.frame(table(bal.DE.cell$celltype)),desc(Freq))
arrange(as.data.frame(table(pbmc.DE.cell$celltype)),desc(Freq))

arrange(as.data.frame(table(bal.DE.cell$balDEs)),desc(Freq))
arrange(as.data.frame(table(pbmc.DE.cell$pbmcDEs)),desc(Freq))


bal.DE.cell[bal.DE.cell$balDEs == "NEAT1",]
pbmc.DE.cell[pbmc.DE.cell$pbmcDEs == "NEAT1",]
pbmc.DE.cell[pbmc.DE.cell$pbmcDEs == "IGLC3",]


##### Cross-Compartment Differential Expression #####
compartment<-gsub("_[A-Z]","",covid.integrated$aggcondition.rn)
names(compartment)<-names(covid.integrated$aggcondition.rn)
covid.integrated$compartment<-compartment

library(EnhancedVolcano)
options(mc.cores = 8)
 
celltypes<-levels(covid.integrated$cell.type.merge.coarse)
DEmetadata<-data.frame(cross.sample = c("BAL_H-PBMC_H","BAL_M-PBMC_M","BAL_S-PBMC_S","BAL-PBMC"),
                       p_val_adj = c("H_p_val_adj","M_p_val_adj","S_p_val_adj","all_p_val_adj"), avg_logFC = c("H_avg_logFC","M_avg_logFC","S_avg_logFC","all_avg_logFC"))


covid.DE.cross.coarse<-list()
for (ii in celltypes){

  subcells<-subset(covid.integrated, subset = cell.type.merge.coarse == ii)
  Idents(subcells)<-subcells$aggcondition.rn
  DE.cross.list<-list()
  
  DE.cross.list[["H"]]<-FindMarkers(subcells,ident.1 = "BAL_H", ident.2 = "PBMC_H", assay = "SCT", 
                                slot = "data",test.use = "MAST", latent.vars = "nCount_RNA", features = DE.rec.genes.var)
  
  
  DE.cross.list[["M"]]<-FindMarkers(subcells,ident.1 = "BAL_M", ident.2 = "PBMC_M", assay = "SCT", 
                         slot = "data",test.use = "MAST", latent.vars = "nCount_RNA", features = DE.rec.genes.var)
  
  DE.cross.list[["S"]]<-FindMarkers(subcells,ident.1 = "BAL_S", ident.2 = "PBMC_S", assay = "SCT", 
                         slot = "data",test.use = "MAST", latent.vars = "nCount_RNA", features = DE.rec.genes.var)
  
  DE.cross.list[["All"]]<-FindMarkers(subcells, ident.1 = "BAL", group.by = "compartment", assay = "SCT", 
                         slot = "data",test.use = "MAST", latent.vars = "nCount_RNA", features = DE.rec.genes.var)
  
  volcano.plots<-list()

  for (jj in 1:4){
    cross.sample<-DEmetadata$cross.sample[jj]
    pval<-DEmetadata$p_val_adj[jj]
    logfc<-DEmetadata$avg_logFC[jj]
    colnames(DE.cross.list[[jj]])[c(2,5)]<-c(logfc,pval)
    
    volcano.plots[[cross.sample]]<-EnhancedVolcano(as.data.frame(DE.cross.list[[jj]]),
                                                 lab = rownames(DE.cross.list[[jj]]),
                                                 x = logfc, xlab = "ln fold change", 
                                                 legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
                                                 y = pval,
                                                 title = paste0("Differential Expression for ",ii),
                                                 subtitle = cross.sample,
                                                 FCcutoff = log(2),
                                                 pCutoff = 10e-6
                                                 )
    
    
  }
  covid.DE.cross.coarse[[ii]]<-list(DE.cross=DE.cross.list,volcano.plots=volcano.plots)
}

##### Checking cross-DEGs #####

for (ii in names(covid.DE.cross.coarse)){
  gridExtra::grid.arrange(covid.DE.cross.coarse[[ii]]$volcano.plots$`BAL_H-PBMC_H`,
            covid.DE.cross.coarse[[ii]]$volcano.plots$`BAL_M-PBMC_M`,
            covid.DE.cross.coarse[[ii]]$volcano.plots$`BAL_S-PBMC_S`,
            covid.DE.cross.coarse[[ii]]$volcano.plots$`BAL-PBMC`
            ,nrow=2,ncol=2)
  
  }

##### DE processings######
DE<-covid.DE.results$`CD8+/CCL5+/GZHM+ Memory T Cells`
grid.arrange(DE$volcano.plots$`M-H`,DE$volcano.plots$`S-H`,DE$volcano.plots$`S-M`,ncol=3)
grid.arrange(DE$volcano.plots$`PBMC_M-PBMC_H`,DE$volcano.plots$`PBMC_S-PBMC_H`,DE$volcano.plots$`PBMC_S-PBMC_M`,ncol=3)
grid.arrange(DE$volcano.plots.modules$`M-H`,DE$volcano.plots.modules$`S-H`,DE$volcano.plots.modules$`S-M`,ncol=3)
grid.arrange(DE$volcano.plots.modules$`PBMC_M-PBMC_H`,DE$volcano.plots.modules$`PBMC_S-PBMC_H`,DE$volcano.plots.modules$`PBMC_S-PBMC_M`,ncol=3)

DEsum<-summarize.DE(DE)
saveDExlsx(DEsum,"/mnt/g/My Drive/Finn-Perkins Lab/COVID/Monocle Modules/GZMH+SPON2+ GD T cells DE.xlsx")


VlnPlot(covid.integrated,features = c("NEAT1","C1QC","IGLC3","IGHG1"),group.by = "aggcondition", pt.size = 0)
VlnPlot(covid.integrated,features = c("MALAT1","IGLC3","FABP4"),group.by = "aggcondition", pt.size = 0)
VlnPlot(subset(covid.integrated,idents = "APOE+/CCL18+/FCN1- Monocytes"),features = c("MALAT1","IGLC3"),group.by = "aggpatient", pt.size = 0)
VlnPlot(covid.integrated,features = c("IGHG2","IGHA1","IGLC3","IGHG1"),group.by = "aggcondition", pt.size = 0)
VlnPlot(covid.integrated,features = c("SAT1","FOS","MX1"),group.by = "aggcondition", pt.size = 0)
FeaturePlot(covid.integrated,features=c("NEAT1","MS4A1"),pt.size = .5)
FeaturePlot(covid.integrated,features=c("IGLC3","IGHG1"),pt.size = .5)

DotPlot(subset(covid.integrated,idents = "APOE+/CCL18+/FCN1- Monocytes"),features = c("MALAT1","IGLC3"),group.by = "aggcondition")


DotPlot(subset(covid.integrated,idents = c("MZB1+/XBP1+ Plasmablasts")),features = c("MS4A1","IGLC3","MZB1"),group.by = "aggcondition")
DotPlot(subset(covid.integrated,idents = "APOE+/CCL18+/FCN1- Monocytes"),features = c("MS4A1","IGLC3","MZB1"),group.by = "aggcondition")
DotPlot(subset(covid.integrated,idents = "IGJ+/IGHM+ Plasmablasts"),features = c("MS4A1","IGLC3","MZB1"),group.by = "aggcondition")


VlnPlot(covid.integrated,features = c("NEAT1"),group.by = "cell.type.merge", pt.size = 0)

##resident macrophage marker cd68 from tissue resident macrophages, davies et al nature immunology
VlnPlot(covid.integrated, features = "CD68",group.by = "cell.type.merge",pt.size = 0)+NoLegend()
dot<-DotPlot(subset(covid.integrated,idents = "CD14+/FCN1+/VCAN++/S100A8++/FCGR3A- Monocytes"), features = "CD68",group.by = "aggcondition")

#how many DEs per cell type
balDE<-NULL
pbmcDE<-NULL
for (ii in names(covid.DE.results)){
  balDE[ii]<-dim(covid.DE.results[[ii]]$bal)[1]
  pbmcDE[ii]<-dim(covid.DE.results[[ii]]$pbmc)[1]

  if(is.null(pbmcDE[ii])){
    pbmcDE[ii]<-0
  }
  
}
DE.by.type<-data.frame(row.names = names(covid.DE.results))
DE.by.type$bal<-balDE
DE.by.type$pbmc<-pbmcDE
colnames(DE.by.type)<-c("bal","pbmc")


##### inverted matrix Gene network test #####
tyro<-subset(covid.integrated, idents = "CD4+/TYROBP+/CST3+ T Cells")
tyro<-RunTSNE(tyro)
test<-FindNeighbors(tyro,reduction = "pca",dims = 1:10,do.plot = TRUE)

##list of unique DEGs
genes<-NULL
for(ii in "covid.integrated$cell.type.merge"){
  #bal
  names<-rownames(covid.DE.results[[ii]]$bal)
  genes<-c(genes,names[!(names %in% genes)])
  
  #pbmc
  names<-rownames(covid.DE.results[[ii]]$pbmc)
  genes<-c(genes,names[!(names %in% genes)])
}

#removing misc non-coding genes and MT genes
genes<-genes[!(grepl("A(C|L)[0-9]{6}",genes))]
genes<-genes[!(grepl("MT-",genes))]


covid.invert<-CreateSeuratObject(t(covid.integrated@assays$SCT@counts[genes,]),assay= "RNA", min.cells = 0, min.features = 0)
covid.invert <- NormalizeData(covid.invert, normalization.method = "LogNormalize", scale.factor = 10000)
covid.invert <- FindVariableFeatures(covid.invert, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(covid.invert)
covid.invert <- ScaleData(covid.invert, features = all.genes)
covid.invert <- RunPCA(covid.invert, features = VariableFeatures(object = covid.invert))
ElbowPlot(covid.invert)
covid.invert <- RunUMAP(covid.invert, dims = 1:20)

covid.invert <- FindNeighbors(covid.invert, dims = 1:20,prune.SNN = .25)
covid.invert <- FindClusters(covid.invert, resolution = 5)

DimPlot(covid.invert,reduction = "umap")


plot(sort(colSums(covid.invert@graphs$RNA_snn),decreasing=TRUE))

nodescores<-colSums(covid.invert@graphs$RNA_snn)
covid.invert$nodescores<-nodescores

hubs<-data.frame(gene = NULL, cluster = NULL,score = NULL)
for(ii in levels(covid.invert$seurat_clusters)){
  subset<-subset(covid.invert,idents = ii)
  scores<-sort(nodescores[colnames(subset)],decreasing = TRUE)[1:2]
  scores<-data.frame(gene = names(scores),cluster = as.character(ii),score = scores)
  hubs<-rbind(hubs,scores)
  
}

nodelabel<-names(nodescores)
# for(ii in 1:length(nodescores)){
#   if(!(names(nodescores[ii])%in%"NEAT1")){
#     nodelabel[ii]<-""
#   }
# }
for(ii in 1:length(nodescores)){
  if(!(names(nodescores[ii])%in%hubs$gene)){
  nodelabel[ii]<-""
  }
}


net <- graph.adjacency(adjmatrix = as.matrix(x = covid.invert@graphs$RNA_snn),
                       mode = "undirected", weighted = TRUE, diag = FALSE)
plot.igraph(x = net, layout = as.matrix(x = Embeddings(object = covid.invert[["umap"]])),
            edge.width = E(graph = net)$weight, vertex.label = nodelabel,
            vertex.size = 0)

###string


for_gen <- unique(sort(genes))
new_genes <- as.data.frame(for_gen, stringsAsFactors = FALSE)
names(new_genes) <- "gene"

string_db<-STRINGdb$new(version="11",species=9606,score_threshold=400,input_directory = "")
genes.mapped<-string_db$map(new_genes,"gene",removeUnmappedRows = TRUE)

final_graph<-ppi.string(new_genes)
betweenness(final_graph)
subgraph<-subgraph(graph,names(sort(betweenness(graph),decreasing = TRUE))[1:20])
size<-sort(betweenness(graph),decreasing = TRUE)[1:20]/mean(sort(betweenness(graph),decreasing = TRUE))
size<-scale(size,center = FALSE)+1
plot(subgraph,layout = layout_with_lgl,vertex.size = size*10)


##### aggregation and bulk DEGs from DEseq2 #####
cell_group_df <- tibble::tibble(cell=row.names(colData(covid.SCT.cds)),
                                cell_group=covid.SCT.cds@colData$aggpatient)
agg_mat_bulk_patient <- aggregate_gene_expression(covid.SCT.cds,cell_group_df =  cell_group_df)

coldatabulk<-data.frame(colnames(agg_mat_bulk_patient),row.names = colnames(agg_mat_bulk_patient))
colnames(coldatabulk)<-"condition"
coldatabulk$condition<-gsub("[0-9]","",coldatabulk$condition)
coldatabulk$condition<-recode(coldatabulk$condition, "H" = "BAL_H","M" = "BAL_M","S" = "BAL_S")

agg_mat_bulk_patient<-agg_mat_bulk_patient[c("NEAT1","MALAT1","MX1","ISG15"),]

#ANOVA
patient2condition<- gsub("[0-9]","",colnames(agg_mat_bulk_patient))


moduleanova<-list()
tukey<-list()

for (ii in rownames(agg_mat_bulk_patient)){
  
  module<-as.data.frame(cbind(agg_mat_bulk_patient[ii,],patient2condition))
  colnames(module)<-c("exp","condition")
  module$exp<-as.numeric(module$exp)
  modulePBMC<-module[grepl("PBMC",module$condition),]
  module<-module[!grepl("PBMC",module$condition),]
  
  #logFC
  balmean<-aggregate(module$exp, by = list(module$condition), FUN = mean)
  pbmcmean<-aggregate(modulePBMC$exp, by = list(modulePBMC$condition), FUN = mean)
  rownames(balmean)<-balmean$Group.1
  rownames(pbmcmean)<-pbmcmean$Group.1
  
  
  # Compute the analysis of variance
  BAL.aov <- aov(exp ~ condition, data = module)
  PBMC.aov <- aov(exp ~ condition, data = modulePBMC)
  # PostHoc via Tukey
  BAL.ph<-TukeyHSD(BAL.aov)
  PBMC.ph<-TukeyHSD(PBMC.aov)
  tukeyres<-as.data.frame(rbind(BAL.ph$condition,PBMC.ph$condition))
  
  ##logFC
  module<-as.data.frame(cbind(agg_mat_bulk_patient[ii,],patient2condition))
  colnames(module)<-c("exp","condition")
  module$exp<-as.numeric(module$exp)
  
  #logFC
  logfc<-data.frame(rep(0,6),row.names = rownames(tukeyres))
  colnames(logfc)<-"fc"
  means<-aggregate(module$exp, by = list(module$condition), FUN = mean)
  rownames(means)<-means$Group.1
  means$log<-log(means$x+1)
  
  logfc["M-H","fc"]<-means["M","log"]-means["H","log"]
  logfc["S-H","fc"]<-means["S","log"]-means["H","log"]
  logfc["S-M","fc"]<-means["S","log"]-means["M","log"]
  
  logfc["PBMC_M-PBMC_H","fc"]<-means["PBMC_M","log"]-means["PBMC_H","log"]
  logfc["PBMC_S-PBMC_H","fc"]<-means["PBMC_S","log"]-means["PBMC_H","log"]
  logfc["PBMC_S-PBMC_M","fc"]<-means["PBMC_S","log"]-means["PBMC_M","log"]
  tukeyres$avg_logFC<-logfc[rownames(tukeyres),"fc"]
  
  moduleanova[[ii]]<-list(BAL.aov,PBMC.aov)
  
  tukeyres$p_FDR<-c(p.adjust(tukeyres$`p adj`[1:3],method = "fdr",n = dim(agg_mat_bulk_patient)[1]),
                    p.adjust(tukeyres$`p adj`[1:3],method = "fdr",n = dim(agg_mat_bulk_patient)[1]))
  
  tukey[[ii]]<-tukeyres
}

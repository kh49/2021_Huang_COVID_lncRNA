library(Seurat)
library(tidyverse)
library(ggplot2)
library(sctransform)
library(future)
library(MAST)
library(qs)
library(EnhancedVolcano)
library(xlsx)
library(gridExtra)




##import Chua et al. nbt samples with bal, all severe
covid.main.loc<-readRDS("covid_nbt_main.rds")

##multiprocess settings (requires a ton of RAM)
plan("multiprocess")
options(future.globals.maxSize = 15000 *1024^2)
options(mc.cores = 4)

## Reapply all analysis steps and more stringent filtering

covid.main.loc <- PercentageFeatureSet(covid.main.loc, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
# VlnPlot(covid.main.loc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# # # plot1 <- FeatureScatter(covid.main.loc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # plot2 <- FeatureScatter(covid.main.loc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))

covid.main.loc <- subset(covid.main.loc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 1000)
covid.main.loc <- SCTransform(covid.main.loc, vars.to.regress = "percent.mt", verbose = FALSE)


options(mc.cores = 4)
covid.main.loc <- RunPCA(covid.main.loc, verbose = FALSE)
covid.main.loc <- RunUMAP(covid.main.loc, dims =1:30)
saveRDS(covid.main.loc,"covidmainloc.rds")

##Seurat label transfer ##note that there is a bug with seurat SCtransform and anchor finding -
#if there are row.name errors, SCtransform was ran on one or both datasets in a previous version with the bug. rerun SCtransform on covid.integrated if error occurs

covid.main.loc<-readRDS("covidmainloc.rds")
covid.integrated <- readRDS("covidintegratedS4.rds")
DE.rec.genes.var <- readRDS("~/Documents/COVID_Data_Mine/derecgenesvar.rds")

# covid.integrated <- subset(covid.integrated,subset = compartment =="BAL")

#following steps are replaced by SCtransform
# DefaultAssay(covid.main.loc)<-"RNA"
# DefaultAssay(covid.integrated)<-"RNA"
# covid.integrated<-NormalizeData(covid.integrated)
# covid.integrated<-ScaleData(covid.integrated)
# covid.main.loc<-NormalizeData(covid.main.loc)
# covid.main.loc<-ScaleData(covid.main.loc)
# covid.integrated<-FindVariableFeatures(covid.integrated)
covid.integrated <- RunPCA(covid.integrated, verbose = FALSE)
covid.main.loc <- FindVariableFeatures(covid.main.loc)
covid.main.loc <- RunPCA(covid.main.loc, verbose = FALSE)
covid.main.loc <- RunUMAP(covid.main.loc, dims =1:30)

covid.integrated.anchors <- FindTransferAnchors(reference = covid.integrated, query = covid.main.loc, 
                                                dims = 1:50, project.query = TRUE, normalization.method = "LogNormalize",verbose = TRUE)
predictions <- TransferData(anchorset = covid.integrated.anchors, refdata = covid.integrated$cell.type.merge.coarse, 
                            dims = 1:50)
covid.main.loc <- AddMetaData(covid.main.loc, metadata = predictions)

Idents(covid.main.loc)<-covid.main.loc$celltype
p1<-DimPlot(covid.main.loc,label = TRUE)
Idents(covid.main.loc)<-covid.main.loc$predicted.id
p2<-DimPlot(covid.main.loc,label=TRUE)
p1+p2

saveRDS(covid.main.loc,"covidmainlocpred.rds")

covid.main.loc <- readRDS("covidmainlocpred.rds")


##### cleanup ######
Idents(covid.main.loc)<-covid.main.loc$severity
covid.main.loc<-RenameIdents(covid.main.loc, `control` = "NP_H",`moderate` = "NP_M", `critical` = "NP_S")
covid.main.loc$aggcondition.rn<-Idents(covid.main.loc)
Idents(covid.main.loc)<-covid.main.loc$aggcondition.rn

##### Differential Expression of nbt set #####
library(EnhancedVolcano)
options(mc.cores = 8)

DefaultAssay(covid.main.loc) <- "RNA"

# celltypes<-levels(as.factor(covid.main.loc$predicted.id))
celltypes<-c("M1 MoMa","M2 MoMa","NK","CD4 T", "CD8 Memory T")

DEmetadata<-data.frame(NBT.sample = c("M-H","S-H","S-M"),
                       p_val_adj = c("g2_1p_val_adj","g3_1p_val_adj","g3_2p_val_adj"), avg_logFC = c("avg_logFC_g2_1","avg_logFC_g3_1","avg_logFC_g3_2"))


covid.DE.results.nbt<-list()
for (ii in celltypes){
  subcells<-subset(covid.main.loc, subset = predicted.id == ii)
  Idents(subcells)<-subcells$severity
  
  # DE.nbt.combined<-FindDEseurat(subcells,ident.1 = "control", ident.2 = "moderate",ident.3 = "critical", assay = "SCT", 
                                # slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
  DE.nbt.combined<-FindDEseurat(subcells,ident.1 = "control", ident.2 = "moderate",ident.3 = "critical", assay = "RNA",
                                slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
  # DE.bal.combined$module<-modules$module[rownames(DE.bal.combined)]
  
  
  volcano.plots<-list()
  # volcano.plots.modules<-list()
  for (jj in 1:3){
    nbt.sample<-DEmetadata$NBT.sample[jj]
    pval<-DEmetadata$p_val_adj[jj]
    logfc<-DEmetadata$avg_logFC[jj]
    
    volcano.plots[[nbt.sample]]<-EnhancedVolcano(as.data.frame(DE.nbt.combined),
                                                 lab = rownames(DE.nbt.combined),
                                                 x = logfc, xlab = "ln fold change", 
                                                 legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
                                                 y = pval,
                                                 title = paste0("Differential Expression for ",ii),
                                                 subtitle = nbt.sample,
                                                 FCcutoff = log(2),
                                                 pCutoff = 10e-6)
  
  
    
  }
  
  covid.DE.results.nbt[[ii]]<-list(DE.nbt=DE.nbt.combined,volcano.plots=volcano.plots)
}


nbtgenescoarse<-NULL
nbtallgenes<-NULL
nbt.val.genes <- list()
manuscript.genes <- c("NEAT1","MALAT1","MX1","IFIT1","IFIT2","IFIT3","NFKBIA","NUPR1","BCL2A1","MTRNR2L12","CTSL","CTSB")

for (ii in names(covid.DE.results.nbt)){
  
  nbtgenescoarse<-c(nbtgenescoarse,rownames(filter(covid.DE.results.nbt[[ii]]$DE.nbt,g3_2p_val_adj <1e-7,g2_1p_val_adj <1e-7,g3_1p_val_adj <1e-7)))
  nbtallgenes<-c(nbtallgenes,rownames(covid.DE.results.nbt[[ii]]$DE.nbt))
  
  for (jj in manuscript.genes){
    if (jj %in% nbtallgenes){
      nbt.val.genes[[jj]]<-rbind(nbt.val.genes[[jj]],covid.DE.results.nbt[[ii]]$DE.nbt[jj,])
      rownames(nbt.val.genes[[jj]])[dim(nbt.val.genes[[jj]])[1]]<-ii
    }
    
  }
}



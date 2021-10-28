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
library(ComplexHeatmap)
library(ggpubr)
library(circlize)
library(RColorBrewer)



##import Chua et al. pbmc samples with bal, all severe
covid.pbmc<-readRDS("seurat_COVID19_PBMC_cohort1_10x_jonas_FG_2020-08-15.rds")

##multiprocess settings (requires a ton of RAM)
plan("multiprocess")
options(future.globals.maxSize = 15000 *1024^2)
options(mc.cores = 4)

## Reapply all analysis steps and more stringent filtering

covid.pbmc <- PercentageFeatureSet(covid.pbmc, pattern = "^MT-", col.name = "percent.mt")
# Visualize QC metrics as a violin plot
# VlnPlot(covid.pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 1)

# FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

# # # plot1 <- FeatureScatter(covid.pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
# # plot2 <- FeatureScatter(covid.pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# CombinePlots(plots = list(plot1, plot2))

covid.pbmc <- subset(covid.pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA > 1000)
covid.pbmc <- SCTransform(covid.pbmc, vars.to.regress = "percent.mt", verbose = FALSE)


options(mc.cores = 4)
covid.pbmc <- RunPCA(covid.pbmc, verbose = FALSE)
covid.pbmc <- RunUMAP(covid.pbmc, dims =1:30)
saveRDS(covid.pbmc,"covidpbmc.rds")

##Seurat label transfer ##note that there is a bug with seurat SCtransform and anchor finding -
#if there are row.name errors, SCtransform was ran on one or both datasets in a previous version with the bug. rerun SCtransform on covid.integrated if error occurs


# covid.integrated <- subset(covid.integrated,subset = compartment =="PBMC")

#following steps are replaced by SCtransform
# DefaultAssay(covid.pbmc)<-"RNA"
# DefaultAssay(covid.integrated)<-"RNA"
# covid.integrated<-NormalizeData(covid.integrated)
# covid.integrated<-ScaleData(covid.integrated)
# covid.pbmc<-NormalizeData(covid.pbmc)
# covid.pbmc<-ScaleData(covid.pbmc)
# covid.integrated<-FindVariableFeatures(covid.integrated)
covid.integrated <- RunPCA(covid.integrated, verbose = FALSE)
covid.pbmc <- RunPCA(covid.pbmc, verbose = FALSE)
covid.pbmc <- RunUMAP(covid.pbmc, dims =1:15)



covid.integrated.anchors <- FindTransferAnchors(reference = covid.integrated, query = covid.pbmc, 
                                                dims = 1:50, project.query = TRUE, normalization.method = "LogNormalize",verbose = TRUE)
predictions <- TransferData(anchorset = covid.integrated.anchors, refdata = covid.integrated$cell.type.merge.coarse, 
                            dims = 1:50)
covid.pbmc <- AddMetaData(covid.pbmc, metadata = predictions)

Idents(covid.pbmc)<-covid.pbmc$id.celltype
p1<-DimPlot(covid.pbmc,label = TRUE)
Idents(covid.pbmc)<-covid.pbmc$predicted.id
p2<-DimPlot(covid.pbmc,label=TRUE)
p1+p2

saveRDS(covid.integrated,"covidintegratedS4.rds")
saveRDS(covid.pbmc,"covidpbmc.rds")

covid.pbmc <- readRDS("covidpbmc.rds")
# covid.integrated <- readRDS("covidintegratedS4.rds")
DE.rec.genes.var <- readRDS("~/Documents/COVID_Data_Mine/derecgenesvar.rds")

##### cleanup ######
Idents(covid.pbmc)<-covid.pbmc$group_per_sample
covid.pbmc<-RenameIdents(covid.pbmc, `control` = "PBMC_H",`mild` = "PBMC_M", `severe` = "PBMC_S")
covid.pbmc$aggcondition.rn<-Idents(covid.pbmc)
Idents(covid.pbmc)<-covid.pbmc$aggcondition.rn

##### Differential Expression of pbmc set #####
library(EnhancedVolcano)
options(mc.cores = 8)

DefaultAssay(covid.pbmc) <- "RNA"
# covid.pbmc <- NormalizeData(covid.pbmc)

# celltypes<-levels(as.factor(covid.pbmc$predicted.id))
celltypes<-c("M1 MoMa","M2 MoMa","NK","CD4 T", "CD8 Memory T")

DEmetadata<-data.frame(PBMC.sample = c("M-H","S-H","S-M"),
                       p_val_adj = c("g2_1p_val_adj","g3_1p_val_adj","g3_2p_val_adj"), avg_logFC = c("avg_logFC_g2_1","avg_logFC_g3_1","avg_logFC_g3_2"))


covid.DE.results.pbmc<-list()
for (ii in celltypes){
  subcells<-subset(covid.pbmc, subset = predicted.id == ii)
  Idents(subcells)<-subcells$group_per_sample
  
  # DE.pbmc.combined<-FindDEseurat(subcells,ident.1 = "control", ident.2 = "moderate",ident.3 = "critical", assay = "SCT", 
  # slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
  DE.pbmc.combined<-FindDEseurat(subcells,ident.1 = "control", ident.2 = "mild",ident.3 = "severe", assay = "RNA",
                                slot = "data",test.use = "MAST", latent.vars = "nCount_RNA")
  # DE.bal.combined$module<-modules$module[rownames(DE.bal.combined)]
  
  
  volcano.plots<-list()
  # volcano.plots.modules<-list()
  for (jj in 1:3){
    pbmc.sample<-DEmetadata$PBMC.sample[jj]
    pval<-DEmetadata$p_val_adj[jj]
    logfc<-DEmetadata$avg_logFC[jj]
    
    volcano.plots[[pbmc.sample]]<-EnhancedVolcano(as.data.frame(DE.pbmc.combined),
                                                 lab = rownames(DE.pbmc.combined),
                                                 x = logfc, xlab = "ln fold change", 
                                                 legendLabels = c("NS", "ln fold change > ln(2)","p-val < 10e-6","lnFC > ln(2) & p-val < 10e-6"),
                                                 y = pval,
                                                 title = paste0("Differential Expression for ",ii),
                                                 subtitle = pbmc.sample,
                                                 FCcutoff = log(2),
                                                 pCutoff = 10e-6)
    
    
    
  }
  
  covid.DE.results.pbmc[[ii]]<-list(DE.pbmc=DE.pbmc.combined,volcano.plots=volcano.plots)
}
save.image()

pbmcgenescoarse<-NULL
pbmcallgenes<-NULL
pbmc.val.genes <- list()
# MALAT1.pbmc <- NULL
manuscript.genes <- c("NEAT1","MALAT1","MX1","IFIT1","IFIT2","IFIT3","NFKBIA","NUPR1","BCL2A1","MTRNR2L12","CTSL","CTSB")

for (ii in names(covid.DE.results.pbmc)){
  
  pbmcgenescoarse<-c(pbmcgenescoarse,rownames(filter(covid.DE.results.pbmc[[ii]]$DE.pbmc,g3_2p_val_adj <1e-7,g2_1p_val_adj <1e-7,g3_1p_val_adj <1e-7)))
  pbmcallgenes<-c(pbmcallgenes,rownames(covid.DE.results.pbmc[[ii]]$DE.pbmc))
  
  for (jj in manuscript.genes){
    if (jj %in% pbmcallgenes){
      pbmc.val.genes[[jj]]<-rbind(pbmc.val.genes[[jj]],covid.DE.results.pbmc[[ii]]$DE.pbmc[jj,])
      rownames(pbmc.val.genes[[jj]])[dim(pbmc.val.genes[[jj]])[1]]<-ii
    }

  }
  
  # if ("MALAT1" %in% pbmcallgenes){
  #   MALAT1.pbmc<-rbind(MALAT1.pbmc,covid.DE.results.pbmc[[ii]]$DE.pbmc["MALAT1",])
  #   rownames(MALAT1.pbmc)[dim(MALAT1.pbmc)[1]]<-ii
  # }
}


##### plotting #####

celltypes<-c("M1 MoMa","M2 MoMa","NK","CD4 T", "CD8 Memory T")
Idents(covid.pbmc)<-covid.pbmc$predicted.id
Idents(covid.main.loc)<-covid.main.loc$predicted.id

val.heatmaps<-list()
heatmapDEs<-list()
for (ii in celltypes){
  balsubset<-subset(covid.main.loc,idents = ii)
  pbmcsubset<-subset(covid.pbmc,idents = ii)
  Idents(balsubset)<-balsubset$aggcondition.rn
  Idents(pbmcsubset)<-pbmcsubset$aggcondition.rn
  
  subset<-merge(balsubset, y = pbmcsubset, add.cell.ids = c("nbt", "pbmc"), project = "validation")
  
  
  
  balsubset<-AverageExpression(balsubset,return.seurat = TRUE,assays = "RNA")
  
  pbmcsubset<-AverageExpression(pbmcsubset,return.seurat = TRUE,assays = "RNA")
  
  subset<-AverageExpression(subset,return.seurat = TRUE,assays = "RNA")
  
  # using gene results from full group testing to filter differential genes at the cell type level
  balgenes<-c(rownames(filter(covid.DE.results.nbt[[ii]]$DE.nbt,g3_2p_val_adj <.05,g2_1p_val_adj <.05,g3_1p_val_adj <.05)))
  pbmcgenes<-c(rownames(filter(covid.DE.results.pbmc[[ii]]$DE.pbmc,g3_2p_val_adj <.05,g2_1p_val_adj <.05,g3_1p_val_adj <.05)))
  
  DE.genes<-c(balgenes,pbmcgenes[!(pbmcgenes %in% balgenes)])
  DE.genes<- DE.genes[DE.genes %in% manuscript.genes]
  
  # DE.genes<-c(DE.genes)
  DE.genes <- manuscript.genes
  

  # heatmapDEs[[ii]]<-cbind(balDEs,pbmcDEs)
  
  # write.xlsx(heatmapDEs[[ii]], file="heatmapDEs.xlsx",
  # col.names=TRUE, row.names=TRUE, append = T,sheetName = ii)  
  
  # sweep(typetotal,2,as.matrix(typetotal["Total",]),"/")[1:15,]
  
  #send to coarseheatmap to get clustering, then back to seurat for nice looking figure
  zbal<-DoHeatmap(balsubset,features =c(DE.genes), size = 3, draw.lines = FALSE)
  zpbmc<-DoHeatmap(pbmcsubset,features =c(DE.genes), size = 3, draw.lines = FALSE)
  
  # this will use separately scaled BAL and PBMC d
  z<-DoHeatmap(subset,features =c(DE.genes), size = 3, draw.lines = FALSE)
  
  z$data<-rbind(zbal$data,zpbmc$data)
  
  y <- z$data %>% drop_na()
  x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
    tidyr::spread(key = Feature, value = Expression)
  w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
    tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()
  
  rdeganno<-data.frame(rDEG=matrix(NA,nrow=nrow(w),ncol = 1))
  rownames(rdeganno)<-rownames(w)
  rdeganno[rownames(rdeganno) %in% balgenes,1]<-"nbt"
  rdeganno[rownames(rdeganno) %in% pbmcgenes,1]<-"pbmc"
  rdeganno[(rownames(rdeganno) %in% pbmcgenes)&(rownames(rdeganno) %in% balgenes),1]<-"both"
  avgexp<-subset@assays$RNA@counts[rownames(rdeganno),]
  avgexp<-data.frame(nbt=rowMeans(avgexp[,c("NP_H","NP_M","NP_S")]),pbmc=rowMeans(avgexp[,c("PBMC_H","PBMC_M","PBMC_S")]))
  avgexp$avg<-rowMeans(avgexp)
  
  ##set annotation after clustering order has been set in w
  ordered.DE.genes<-rownames(w)
  blank<-matrix(NA,nrow=length(ordered.DE.genes))
  manu.gene.anno<-data.frame(nbt.HvsM=blank,nbt.HvsS=blank,nbt.SvsM=blank,pbmc.HvsM=blank,pbmc.HvsS=blank,pbmc.SvsM=blank)
  rownames(manu.gene.anno)<-ordered.DE.genes
  
  balDEs<-covid.DE.results.nbt[[ii]]$DE.nbt[ordered.DE.genes,]
  balDEs<-balDEs[,-(1:3)]
  colnames(balDEs)<-paste0("NBT_",colnames(balDEs))
  rownames(balDEs)<-ordered.DE.genes
  manu.gene.anno[,1:3][balDEs[,7:9]<.001]<-"*"
  manu.gene.anno[,1:3][.001<=balDEs[,7:9] & balDEs[,7:9]<.01]<-""
  manu.gene.anno[,1:3][.01<=balDEs[,7:9] & balDEs[,7:9]<.05]<-""
  manu.gene.anno[,1:3][.05<=balDEs[,7:9]]<-""
  # manu.gene.anno[,1:3][is.na(balDEs[,9])]<-""
  
  
  pbmcDEs<-covid.DE.results.pbmc[[ii]]$DE.pbmc[ordered.DE.genes,]
  pbmcDEs<-pbmcDEs[,-(1:3)]
  colnames(pbmcDEs)<-paste0("PBMC_",colnames(pbmcDEs))
  rownames(pbmcDEs)<-ordered.DE.genes
  manu.gene.anno[,4:6][pbmcDEs[,7:9]<.001]<-"*"
  manu.gene.anno[,4:6][.001<=pbmcDEs[,7:9] & pbmcDEs[,7:9]<.01]<-""
  manu.gene.anno[,4:6][.01<=pbmcDEs[,7:9] & pbmcDEs[,7:9]<.05]<-""
  manu.gene.anno[,4:6][.05<=pbmcDEs[,7:9]]<-""
  # manu.gene.anno[,4:6][is.na(pbmcDEs[,9])]<-""
  
  manu.gene.anno[,7:9]<-balDEs[,7:9]
  manu.gene.anno[,10:12]<-pbmcDEs[,7:9]
  colnames(manu.gene.anno)[7:12]<-c("nbt.HvsM.p","nbt.HvsS.p","nbt.SvsM.p","pbmc.HvsM.p","pbmc.HvsS.p","pbmc.SvsM.p")
  
  #Z transform
  avgexp<-sweep(avgexp[,1:2],1,as.matrix(avgexp$avg),"/")
  
  
  rdeganno$"Exp Z Diff"<-avgexp$nbt-avgexp$pbmc
  
  col_fun<-colorRamp2(c(-2,0,2),c("purple","black","gold"))
  col_fun2<-colorRamp2(c(0,.05,.1),RColorBrewer::brewer.pal(9,"PiYG")[c(1,5,9)])
  
  ha.l<-rowAnnotation(NP_HvsM = anno_simple(manu.gene.anno$nbt.HvsM.p,col =col_fun2, pch = manu.gene.anno$nbt.HvsM),
                    NP_HvsS = anno_simple(manu.gene.anno$nbt.HvsS.p,col =col_fun2, pch = manu.gene.anno$nbt.HvsS),
                    NP_SvsM = anno_simple(manu.gene.anno$nbt.SvsM.p,col =col_fun2, pch = manu.gene.anno$nbt.SvsM)
                    )
  ha.r<-rowAnnotation(PBMC_HvsM = anno_simple(manu.gene.anno$pbmc.HvsM.p,col =col_fun2, pch = manu.gene.anno$pbmc.HvsM),
                      PBMC_HvsS = anno_simple(manu.gene.anno$pbmc.HvsS.p,col =col_fun2, pch = manu.gene.anno$pbmc.HvsS),
                      PBMC_SvsM = anno_simple(manu.gene.anno$pbmc.SvsM.p,col =col_fun2, pch = manu.gene.anno$pbmc.SvsM)
  )
  
  # ha<-rowAnnotation(df=manu.gene.anno,col=list(nbt.p=c("nbt"="blue","pbmc"="red","both"="purple"),`Exp Z Diff` = col_fun2))
  H<-Heatmap(w, left_annotation = ha.l,right_annotation = ha.r,col=col_fun,cluster_columns = FALSE,
             column_split = c("NBT","NBT","NBT","PBMC","PBMC","PBMC"),name = "Exp",
             column_title = ii,width = unit(5, "cm"))
  
  
  # val.heatmaps[[ii]]<-DoHeatmap(subset,features = rev(DE.genes)[row_order(H)],size = 3, draw.lines = FALSE)
  lgd_pvalue = Legend(title = "p-value", col_fun = col_fun2, at = c(.001, .05, .1), 
                      labels = c("0.001", "0.05", ".1"))
  # and one for the significant p-values
  lgd_sig = Legend(pch = "*", type = "points", labels = "< 0.001")
  # these two self-defined legends are added to the plot by `annotation_legend_list`
  Hanno<-grid.grabExpr(draw(H, annotation_legend_list = list(lgd_pvalue, lgd_sig)))
  
  
  val.heatmaps[[ii]]<-Hanno
  

}

for (plot in val.heatmaps){
  grid.newpage()
  grid.draw(plot)
  
}

##### UMAP plotting #####
p1<-DimPlot(covid.main.loc, reduction = "umap",label=FALSE,group.by = c("celltype"),
        pt.size = .5,label.size = 5,repel = TRUE)+scale_color_viridis_d(begin = 0,end =1)+NoLegend()
p1<-LabelPlotBackground(plot = p1, id = "celltype", color = "black",size = 5)
p1<-p1+ggtitle("NP Original")

p2<-DimPlot(covid.main.loc, reduction = "umap",label=FALSE,group.by = c("predicted.id"),
            pt.size = .5,label.size = 5,repel = TRUE)+scale_color_viridis_d(begin = 0,end =1)+NoLegend()
p2<-LabelPlotBackground(plot = p2, id = "predicted.id", color = "black",size = 5)
p2<-p2+ggtitle("NP Predicted")
p1+p2

p3<-DimPlot(covid.pbmc, reduction = "umap",label=FALSE,group.by = c("id.celltype"),
            pt.size = .5,label.size = 5,repel = TRUE)+scale_color_viridis_d(begin = 0,end =1)+NoLegend()
p3<-LabelPlotBackground(plot = p3, id = "id.celltype", color = "black",size = 5)
p3<-p3+ggtitle("PBMC Original")

p4<-DimPlot(covid.pbmc, reduction = "umap",label=FALSE,group.by = c("predicted.id"),
            pt.size = .5,label.size = 5,repel = TRUE)+scale_color_viridis_d(begin = 0,end =1)+NoLegend()
p4<-LabelPlotBackground(plot = p4, id = "predicted.id", color = "black",size = 5)
p4<-p4+ggtitle("PBMC Predicted")
p3+p4

#####logfc comparison #####
covid.DE.results.coarse<-readRDS("covidDEresultscoarse.rds")
covid.DE.results.coarse<-covid.DE.results.coarse[1:5]
names(covid.DE.results.coarse)[1:2]<-c("M1 MoMa","M2 MoMa")

DE.fc<-list()
DE.fch<-list()
DE.agree<-list()


DE.genes<-manuscript.genes
# DE.genes<-DE.rec.genes.var

for (ii in celltypes){
  defc<-data.frame(cbind(covid.DE.results.coarse[[ii]]$DE.bal[DE.genes,4:6],
                         covid.DE.results.coarse[[ii]]$DE.pbmc[DE.genes,4:6],
                         covid.DE.results.nbt[[ii]]$DE.nbt[DE.genes,4:6],
                         covid.DE.results.pbmc[[ii]]$DE.pbmc[DE.genes,4:6]))
  p.high.select<-data.frame(cbind(covid.DE.results.coarse[[ii]]$DE.bal[DE.genes,10:12],
                             covid.DE.results.coarse[[ii]]$DE.pbmc[DE.genes,10:12],
                             covid.DE.results.nbt[[ii]]$DE.nbt[DE.genes,10:12],
                             covid.DE.results.pbmc[[ii]]$DE.pbmc[DE.genes,10:12]))
  p.high.select<-p.high.select > .05
  defc[p.high.select]<-NA
  
  colnames(defc)<-c("balMH","balSH","balSM","pbmcMH","pbmcSH","pbmcSM","npMH","npSH","npSM","pbMH","pbSH","pbSM")
  defc.bal<-defc[!(is.na(defc[,1])),c(1,2,3,7,8,9)]
  defc.pbmc<-defc[!(is.na(defc[,4])),c(4,5,6,10,11,12)]
  
  DE.fc[[ii]]<-defc
  bal.sign<-sign(defc.bal)[,1:3]
  bal.sign[is.na(bal.sign)]<- -3
  pbmc.sign<-sign(defc.pbmc)[,1:3]
  pbmc.sign[is.na(pbmc.sign)]<- -3
  nbt.sign<-sign(defc.bal)[,4:6]
  nbt.sign[is.na(nbt.sign)]<- -2
  pb.sign<-sign(defc.pbmc)[,4:6]
  pb.sign[is.na(pb.sign)]<- -2
  DE.agree$bal[[ii]]<-bal.sign*nbt.sign
  DE.agree$pbmc[[ii]]<-pbmc.sign*pb.sign

  # DE.fch[[ii]]<-Heatmap(defc,col=col_fun,cluster_columns = FALSE,
             # column_split = c("BAL","BAL","BAL","PBMC","PBMC","PBMC","NP","NP","NP","PB","PB","PB"),name = "logFC",
              # column_title = ii)
  # defc<-rownames_to_column(defc,var = "gene") %>% pivot_longer(!gene,names_to = "condition",values_to = "log2FC")%>%
  #   mutate(sign=sign(log2FC),compartment=ifelse(grepl("bal/|np",condition),"Local","Peripheral"))
  
}

##summarizing DE direction trends
table(unlist(DE.agree))
table(unlist(DE.agree$bal))

bal.dir<-data.frame(matrix(0,ncol = 8, nrow = 5), row.names = celltypes)
colnames(bal.dir)<-c("-3","3","-2","2","-1","1","-6","6")

pbmc.dir<-bal.dir
bal.dir.control<-bal.dir
pbmc.dir.control<-bal.dir
bal.dir.cases<-bal.dir
pbmc.dir.cases<-bal.dir

bal.genes.dir<-DE.agree$bal
pbmc.genes.dir<-DE.agree$pbmc

for (ii in celltypes){
bal.dir[ii,names(table(unlist(DE.agree$bal[[ii]])))]<-table(unlist(DE.agree$bal[[ii]]))
pbmc.dir[ii,names(table(unlist(DE.agree$pbmc[[ii]])))]<-table(unlist(DE.agree$pbmc[[ii]]))  

bal.dir.control[ii,names(table(unlist(DE.agree$bal[[ii]][,-3])))]<-table(unlist(DE.agree$bal[[ii]][,-3]))
pbmc.dir.control[ii,names(table(unlist(DE.agree$pbmc[[ii]][,-3])))]<-table(unlist(DE.agree$pbmc[[ii]][,-3]))  

bal.dir.cases[ii,names(table(unlist(DE.agree$bal[[ii]][,3])))]<-table(unlist(DE.agree$bal[[ii]][,3]))
pbmc.dir.cases[ii,names(table(unlist(DE.agree$pbmc[[ii]][,3])))]<-table(unlist(DE.agree$pbmc[[ii]][,3])) 

bal.genes.dir[[ii]]<-bal.genes.dir[[ii]]%>% mutate_all(function(x){recode(x,`3` = "na.orig",`-3` = "na.orig",
                                                                          `2`="na.val",`-2`="na.val",
                                                                          `6`="na.all",`-6`="na.all",
                                                                          `1`="agree",`-1`="disagree")})
pbmc.genes.dir[[ii]] <- pbmc.genes.dir[[ii]] %>% mutate_all(function(x){recode(x,`3` = "na.orig",`-3` = "na.orig",
                                                                               `2`="na.val",`-2`="na.val",
                                                                               `6`="na.all",`-6`="na.all",
                                                                               `1`="agree",`-1`="disagree")})
write.xlsx(bal.genes.dir[[ii]], file="balvalgene.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = ii)

write.xlsx(pbmc.genes.dir[[ii]], file="pbmcvalgene.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = ii)

}

bal.dir <- bal.dir %>% mutate(na.orig = `-3`+`3`,na.val = `-2`+`2`,na.all = `-6`+`6`,disagree = `-1`,agree = `1`) %>%
  select(last_col(4):last_col())
pbmc.dir <- pbmc.dir %>% mutate(na.orig = `-3`+`3`,na.val = `-2`+`2`,na.all = `-6`+`6`,disagree = `-1`,agree = `1`) %>%
  select(last_col(4):last_col())

bal.dir.control <- bal.dir.control %>% mutate(na.orig = `-3`+`3`,na.val = `-2`+`2`,na.all = `-6`+`6`,disagree = `-1`,agree = `1`) %>%
  select(last_col(4):last_col())
pbmc.dir.control <- pbmc.dir.control %>% mutate(na.orig = `-3`+`3`,na.val = `-2`+`2`,na.all = `-6`+`6`,disagree = `-1`,agree = `1`) %>%
  select(last_col(4):last_col())

bal.dir.cases <- bal.dir.cases %>% mutate(na.orig = `-3`+`3`,na.val = `-2`+`2`,na.all = `-6`+`6`,disagree = `-1`,agree = `1`) %>%
  select(last_col(4):last_col())
pbmc.dir.cases <- pbmc.dir.cases %>% mutate(na.orig = `-3`+`3`,na.val = `-2`+`2`,na.all = `-6`+`6`,disagree = `-1`,agree = `1`) %>%
  select(last_col(4):last_col())


write.xlsx(bal.dir, file="DEvalidation.xlsx",
col.names=TRUE, row.names=TRUE, append = T,sheetName = "BAL Validation")

write.xlsx(bal.dir.control, file="DEvalidation.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = "BAL Val - Controls Only")

write.xlsx(bal.dir.cases, file="DEvalidation.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = "BAL Val - Cases Only")

write.xlsx(pbmc.dir, file="DEvalidation.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = "PBMC Validation")

write.xlsx(pbmc.dir.control, file="DEvalidation.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = "PBMC Val - Controls Only")

write.xlsx(pbmc.dir.cases, file="DEvalidation.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = "PBMC Val - Cases Only")
        

##### dotplot ####
orderednamescoarse<-c("M1 MoMa","M2 MoMa","Int MoMa","CD4 T",
                      "Treg","CD8 T","CD8 Memory T","NK","Neutrophils","Naive B",
                      "Plasmablasts","Plasmacytoid Dendritic","Myeloid Dendritic",
                      "Epithelium/Granulocytes","Epithelium/Pneumocytes/Ciliary")

orderedpbmc <- orderednamescoarse[orderednamescoarse %in% covid.pbmc$predicted.id]
covid.pbmc$predicted.id<-factor(covid.pbmc$predicted.id,levels = orderedpbmc)

orderednbt <- orderednamescoarse[orderednamescoarse %in% covid.main.loc$predicted.id]
covid.main.loc$predicted.id<-factor(covid.main.loc$predicted.id,levels = orderednbt)

#reorder names in validation sets for ease of graph interpretation
orderednamesnbt<-c("nrMa","MoD-Ma","moDC","rMa","CTL","Treg","NKT","NKT-p","NK","Neu","B cell",
                   "pDC","Basal","Ciliated","Ciliated-diff","FOXN4","Ionocyte","IRC","MC",
                   "outliers_epithelial","Secretory","Secretory-diff","Squamous","unknown_epithelial")
covid.main.loc$celltype<-factor(covid.main.loc$celltype,levels = orderednamesnbt)

orderednamespbmc<-c("0: Classical Monocytes","1: HLA-DR+ CD83+ Monocytes ","2: CD163+ Monocytes","3: HLA-DR- S100A+ monocytes",
                    "4: Non-classical Monocytes","9: CD4+ T cells_1","10: CD4+ T cells_2","11: CD4+ T cells_3","12: CD8+ T cells_1",
                    "13: CD8+ T cells_2","14: CD8+ T cells_3","15: NK cells","5: Neutrophils","6: Immature Neutrophils","16: B cells_1",
                    "17: B cells_2","18: B cells_3","19: Plasmablasts","7: mDCs","8: pDCs","20: Megakaryocyte","21: mixed","22: undefined")
covid.pbmc$id.celltype<-factor(covid.pbmc$id.celltype,levels = orderednamespbmc)

my_breaks<-c(10,100,1000,5000,10000)

covid.main.cells<-as.data.frame(as.data.frame(table(covid.main.loc$celltype,covid.main.loc$predicted.id)))
nbtcelltotals<-table(covid.main.loc$predicted.id)
covid.main.cells$Total<-as.matrix(nbtcelltotals[covid.main.cells$Var2])
covid.main.cells %>%  mutate(Pct = (Freq/Total) * 100) %>%
  ggplot(aes(x=Var1, y = Var2, color = Freq, size = Pct)) + 
  geom_point()+theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1.05,hjust = 1) )+
  scale_color_gradient(name = "Freq", trans = "log",
                     breaks = my_breaks, labels = my_breaks)+xlab("Nasopharyngeal")+ylab("Bronchoalveolar")

covid.pbmc.cells<-as.data.frame(as.data.frame(table(covid.pbmc$id.celltype,covid.pbmc$predicted.id)))
pbmccelltotals<-table(covid.pbmc$predicted.id)
covid.pbmc.cells$Total<-as.matrix(pbmccelltotals[covid.pbmc.cells$Var2])
covid.pbmc.cells %>%  mutate(Pct = (Freq/Total) * 100) %>%
  ggplot(aes(x=Var1, y = Var2,color = Freq, size = Pct)) + 
  geom_point()+theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1.05,hjust = 1) )+
  scale_color_gradient(name = "Freq", trans = "log",
                      breaks = my_breaks, labels = my_breaks)+xlab("PBMC (Validation)")+ylab("PBMC")

covid.main.cells %>%  mutate(Pct = (Freq/Total) * 100) %>%
  ggplot(aes(x=Var1, y = Var2, color = log10(Freq), size = Pct)) + 
  geom_point()+theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1.05,hjust = 1) )+
  scale_color_viridis_c()+xlab("Nasopharyngeal")+ylab("Bronchoalveolar")

covid.pbmc.cells %>%  mutate(Pct = (Freq/Total) * 100) %>%
  ggplot(aes(x=Var1, y = Var2,color = log10(Freq), size = Pct)) + 
  geom_point()+theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1.05,hjust = 1) )+
  scale_color_viridis_c()+xlab("PBMC (Validation)")+ylab("PBMC")

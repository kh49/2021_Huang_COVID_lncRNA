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
library(ComplexHeatmap)
library(ggpubr)
library(circlize)
library(RColorBrewer)


#preprocessing
balpatients<-!(grepl("PBMC",covid.integrated$aggpatient))

covid.integrated$aggpatient.rn<-covid.integrated$aggpatient
covid.integrated$aggpatient.rn[balpatients]<-paste0("BAL_",covid.integrated$aggpatient.rn[balpatients])

#simplified naming scheme #includes correction for 3b gamma delta T to NK
Idents(covid.integrated)<-covid.integrated$seurat_clusters.merge

covid.integrated <- RenameIdents(covid.integrated, `0` = "M1 MoMa", `1` = "M2 MoMa", `3b` = "NK",
                                 `4` = "M1 MoMa", `5` = "M1 MoMa", `5b` = "CD4 T", `6` = "M1 MoMa", `7b` = "CD8 Memory T", `8` = "M1 MoMa",
                                 `12` = "M1 MoMa",`12b` = "NK", `13` = "Epithelium/Granulocytes", `14` = "Neutrophils", `14b` ="CD4 T",
                                 `15` = "Plasmablasts", `15b` ="Naive B", `16` = "CD8 T", `17` = "Int MoMa", `19` = "M2 MoMa", `20` = "M2 MoMa",
                                 `21b` = "Plasmablasts", `22` = "Plasmacytoid Dendritic", `22b` = "Treg", `25b` = "Int MoMa",`26b` = "CD8 Memory T",`27b` = "Myeloid Dendritic",
                                 `28b` = "Epithelium/Pneumocytes/Ciliary")

covid.integrated$cell.type.merge.coarse<-Idents(covid.integrated)
orderednamescoarse<-c("M1 MoMa","M2 MoMa","Int MoMa","CD4 T",
                      "Treg","CD8 T","CD8 Memory T","NK","Neutrophils","Naive B",
                      "Plasmablasts","Plasmacytoid Dendritic","Myeloid Dendritic",
                      "Epithelium/Granulocytes","Epithelium/Pneumocytes/Ciliary")

covid.integrated$cell.type.merge.coarse<-factor(covid.integrated$cell.type.merge.coarse,levels = orderednamescoarse)


Idents(covid.integrated)<-covid.integrated$seurat_clusters.merge

covid.integrated <- RenameIdents(covid.integrated, `0` = "M1-a", `1` = "M2-a", `3b` = "NK",
                                 `4` = "M1-b", `5` = "M1-c", `5b` = "CD4 T-a", `6` = "M1-d", `7b` = "CD8 Memory T-a", `8` = "M1-e",
                                 `12` = "M1-f",`12b` = "NK", `13` = "E/G", `14` = "Neutrophils", `14b` ="CD4 T-b",
                                 `15` = "Plasmablasts-a", `15b` ="Naive B", `16` = "CD8 T", `17` = "Int-a", `19` = "M2-b", `20` = "M2-c",
                                 `21b` = "Plasmablasts-b", `22` = "Plasmacytoid Dendritic", `22b` = "Treg", `25b` = "Int-b",`26b` = "CD8 Memory T-b",`27b` = "Myeloid Dendritic",
                                 `28b` = "E/P/C")

covid.integrated$cell.type.merge.subgroup<-Idents(covid.integrated)
covid.integrated$cell.type.merge.subgroup<-factor(covid.integrated$cell.type.merge.subgroup,levels = c("M1-a","M1-b","M1-c","M1-d",
                                                                                                       "M1-e","M1-f","M2-a","M2-b","M2-c","Int-a","Int-b",
                                                                                                       "CD4 T-a","CD4 T-b","Treg","CD8 T","CD8 Memory T-a","CD8 Memory T-b",
                                                                                                       "NK","Neutrophils","Naive B", "Plasmablasts-a",
                                                                                                       "Plasmablasts-b","Plasmacytoid Dendritic","Myeloid Dendritic","E/G",
                                                                                                       "E/P/C"))
Idents(covid.integrated)<-covid.integrated$cell.type.merge.coarse

VlnPlot(subset(covid.integrated, idents = c("M1 MoMa","M2 MoMa") ),features = c("NEAT1","MALAT1","MX1","CTSL"),
        group.by = "aggcondition.rn",split.by = "cell.type.merge.coarse", pt.size = 0)
RidgePlot(subset(covid.integrated, idents = c("M1 MoMa","M2 MoMa") ),features = c("NEAT1","MALAT1"),
        group.by = "aggcondition.rn")


#append BAL to naming scheme
Idents(covid.integrated)<-covid.integrated$aggcondition

covid.integrated <- RenameIdents(covid.integrated, `H`="BAL_H",`M`="BAL_M",`S`="BAL_S")
covid.integrated$aggcondition.rn<-Idents(covid.integrated)

#change macrophage naming in DE tables
names(covid.DE.results.coarse)<-gsub("Macrophage","MoMa",names(covid.DE.results.coarse))
names(covid.DE.results.subgroup)<-gsub(" Macrophage","",names(covid.DE.results.subgroup))


##### Figure 1 #####
# covid.integrated$aggcondition<-factor(covid.integrated$aggcondition, levels = c("H","M","S","PBMC_H","PBMC_M","PBMC_S"))
Idents(covid.integrated)<-covid.integrated$cell.type.merge.coarse
DimPlot(covid.integrated, reduction = "umap",label=FALSE,
        group.by = c("aggcondition.rn"),pt.size = .5)+scale_color_viridis_d(begin = 0,end =1)

DimPlot(covid.integrated, reduction = "umap",label=FALSE,
        group.by = c("cell.type.merge.coarse"),pt.size = .5,
        split.by = "aggcondition.rn",ncol = 3)+scale_color_viridis_d(begin = 0,end =1)


uplot<-DimPlot(covid.integrated, reduction = "umap",label=FALSE,group.by = c("cell.type.merge.coarse"),
               pt.size = .5,label.size = 5,repel = TRUE)+scale_color_viridis_d(begin = 0,end =1)+NoLegend()

LabelPlotBackground(plot = uplot, id = "cell.type.merge.coarse", color = "black",size = 5)

uplot<-DimPlot(covid.integrated, reduction = "umap",label=FALSE,group.by = c("cell.type.merge.subgroup"),
               pt.size = .5,label.size = 5,repel = TRUE)+scale_color_viridis_d(begin = 0,end =1)+NoLegend()

LabelPlotBackground(plot = uplot, id = "cell.type.merge.subgroup", color = "black",size = 5)

#Dotplot
cellmarkers<- c("FCGR3A","CD14","FCN1","VCAN","FN1","IL7R","CD4","IL2RA","LAG3","CD8A","CCL5","GZMH","SPON2",
                "NCAM1","NAMPT","CXCR2","MS4A1","IGJ","MZB1","IRF8","PLD4","CD1C","LGALS2","KRT19","SLPI","LCN2","S100A2","CFAP300","PPIL6")

VlnPlot(covid.integrated,features = c("IGHG1","IGHL3"), pt.size = 0)

# Idents(covid.integrated)<-covid.integrated$cell.type.merge.subgroup
DotPlot(covid.integrated,features = cellmarkers, group.by = "cell.type.merge.subgroup", cols = c("darkorchid3","chartreuse3"))+
  theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1.05,hjust = 1))

FeaturePlot(covid.integrated,features=c("FCN1","FN1","IL7R","LAG3","CD8A","SPON2","NAMPT","MS4A1","MZB1","IRF8","CD1C","SLPI"),ncol = 2)


##### Figure 1 Supplemental Identification #####
Idents(covid.integrated)<-covid.integrated$aggcondition.rn
balsubset_ident<-subset(covid.integrated,idents = c("BAL_H","BAL_M","BAL_S"))
pbmcsubset_ident<-subset(covid.integrated,idents = c("PBMC_H","PBMC_M","PBMC_S"))

#loading in the celltypes from BAL data
bal.cell.idents<-read.delim("all.cell.annotation.meta.txt")
balnum2patient<-data.frame(num = c(9,1:3,5,6,10,11,12,7,8),patientcode = names(table(balsubset_ident$patientcode)) )
rownames(balnum2patient)<-balnum2patient$patientcode
bal.cell.idents$updateID<-base::apply(bal.cell.idents,1,function(X) gsub("_[0-9]*",paste0("-1_",balnum2patient[X[2],1]),X[1]))
balnamedidents<-bal.cell.idents$celltype
names(balnamedidents)<-bal.cell.idents$updateID
balnamedidents<-balnamedidents[names(balnamedidents) %in% colnames(balsubset_ident)]
balsubset_ident$balcelltype<-balnamedidents


balcoarsetable<-table(balsubset_ident$cell.type.merge.coarse,balsubset_ident$balcelltype)
pbmccoarsetable<-table(pbmcsubset_ident$cell.type.merge.coarse,pbmcsubset_ident$cell.type)

balsubtable<-table(balsubset_ident$cell.type.merge.subgroup,balsubset_ident$balcelltype)
pbmcsubtable<-table(pbmcsubset_ident$cell.type.merge.subgroup,pbmcsubset_ident$cell.type)

write.xlsx(balcoarsetable, file="celltable.xlsx",
           col.names=TRUE, row.names=TRUE, append = F,sheetName = "BAL Coarse vs. Original")
write.xlsx(balsubtable, file="celltable.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = "BAL Subgroup vs. Original")

write.xlsx(pbmccoarsetable, file="celltable.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = "PBMC Coarse vs. Original")
write.xlsx(pbmcsubtable, file="celltable.xlsx",
           col.names=TRUE, row.names=TRUE, append = T,sheetName = "PBMC Subgroup vs. Original")


##### Figure 2 proportions#####
# coarseprop$celltype<-rownames(coarseprop)
rownames(coarseprop)<-gsub("Macrophage","MoMa",rownames(coarseprop))
Heatmap(t(scale(t(coarseprop))))
coarseproptotal<-coarseprop
coarseproptotal$celltype<-rownames(coarseproptotal)
coarseproptotal<- pivot_longer(coarseproptotal,cols = c("BAL_H","BAL_M","BAL_S","PBMC_H","PBMC_M","PBMC_S"))
colnames(coarseproptotal)<-c("celltype","condition","prop")
coarseproptotal$sample<-"BAL"
coarseproptotal$condition<-factor(coarseproptotal$condition, levels = c("BAL_H","BAL_M","BAL_S","PBMC_H","PBMC_M","PBMC_S"))

coarseproptotal$sample<-apply(coarseproptotal,1, function(x){ifelse(grepl("PBMC",x[2]),"PBMC","BAL")})
coarseproptotal$celltype<-factor(coarseproptotal$celltype, levels = orderednamescoarse)

theme_set(theme_minimal(base_size=18))
ggplot(coarseproptotal, aes(fill=condition, y=prop, x=sample))+ scale_fill_viridis_d() + 
  geom_bar(position="dodge", stat="identity")+facet_wrap(. ~ celltype,ncol = 5, scales = "free")+labs(x="Sample",y="Abundance")


#ordered abundance of cell types
typetotal<-as.data.frame(cbind(rowSums(coarsetable[,1:3]),rowSums(coarsetable[,4:6])))
rownames(typetotal)<-gsub("Macrophage","MoMa",rownames(typetotal))
colnames(typetotal)<-c("BAL","PBMC")
typeprop<-sweep(typetotal,2,as.matrix(typetotal["Total",]),"/")[1:15,]
typeprop$celltype<-rownames(typeprop)
typeprop<- pivot_longer(typeprop,cols = c("BAL","PBMC"))
colnames(typeprop)<-c("celltype","condition","prop")
typeprop$celltype<-factor(typeprop$celltype,levels = orderednamescoarse)

ggplot(typeprop,aes(fill=condition,y=prop,x=celltype))+scale_fill_brewer(palette = "Set1")+geom_col(position="dodge")+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+labs(x="Cell Type",y="Abundance")


#per patient cell types
patientcoarsetable<-table(covid.integrated$cell.type.merge.coarse,covid.integrated$aggpatient.rn)
patientcoarsetable<-as.data.frame.matrix(patientcoarsetable)
patientcoarsetable["Total",]<-colSums(patientcoarsetable)
patientcoarseprop<-as.data.frame(sweep(patientcoarsetable,2,as.matrix(patientcoarsetable["Total",]),"/")[1:15,])
patientcoarseprop$celltype<-rownames(patientcoarseprop)
patientcoarseprop<-pivot_longer(patientcoarseprop,cols = 1:25)
colnames(patientcoarseprop)<-c("Cell Type","Patient","Abundance")
patientcoarseprop$Patient<-factor(patientcoarseprop$Patient,levels=colnames(patientcoarsetable))

ggplot(patientcoarseprop,aes(fill=`Cell Type`,y=Abundance,x=Patient))+scale_fill_viridis_d()+
  geom_col(position = "fill")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


##### Figure 3 #####
Idents(covid.integrated)<-covid.integrated$cell.type.merge.coarse

balde<-c(0)
pbmcde<-c(0)
lencoarse<-length(levels(covid.integrated$cell.type.merge.coarse))
for (ii in levels(covid.integrated$cell.type.merge.coarse)){
  balde<-sum(balde,dim(covid.DE.results.coarse[[ii]]$DE.bal))
  pbmcde<-sum(pbmcde,dim(covid.DE.results.coarse[[ii]]$DE.pbmc))
  
}
balde/lencoarse
pbmcde/lencoarse


agg_mat<-DEG.module.results$aggregates$condition
colnames(agg_mat) <- recode(colnames(agg_mat),H="BAL_H",M="BAL_M",S="BAL_S")
conditionmap<-pheatmap::pheatmap(agg_mat,
                                 scale="column", clustering_method="ward.D2")
modulegoplots<-list()
for (ii in names(DEG.module.results$module_enrichment)){
plotdata<-DEG.module.results$module_enrichment[[ii]][,c("GO.ID","Term","Fisher")]
plotdata$nlogfisher<--log10(as.numeric(plotdata$Fisher))
plotdata$fullterm<-paste0(plotdata$GO.ID,": ",plotdata$Term)
modulegoplots[[ii]]<-ggplot(plotdata,aes(reorder(fullterm,nlogfisher),nlogfisher))+geom_col(fill="#56B4E9")+
          geom_hline(yintercept = -log10(.05),linetype="dashed",color="red")+theme(axis.text = element_text(size = 10))+coord_flip()+
          labs(title = paste0(ii," GO Enrichment"),y= "-log10 P-val", x = "Enriched GO Term")
}

do.call("ggarrange",c(modulegoplots,ncol=2,nrow=2))


## module membership circle
modules<-DEG.module.results$modules[,1:2]
color = brewer.pal(4,"Set3")
circos.par(gap.degree=0,start.degree=90)
circos.initialize(modules$module,xlim=matrix(c(1,1,1,1,15,15,14,6),ncol=2))
circos.track(ylim=c(0,1),bg.border=NA,track.height=.25)
circos.track(ylim = c(0, 1),bg.col=color,panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ycenter, paste0("Module ",CELL_META$sector.index), 
              facing = "downward")
  circos.axis(h = 1.2, major.at = seq(0, round(CELL_META$xlim[2])), minor.ticks = 0,
              labels.cex = 1,labels = modules$id[modules$module==CELL_META$sector.index],labels.facing="clockwise")
}, track.height = .25)
circos.clear()

##### FIgure 4 DEGs by cell type #####

celltypes<-c("M1 MoMa","M2 MoMa","NK","CD4 T", "CD8 Memory T")
DE.rec.genes.module <- dplyr::select(DEG.module.results$modules,id,module) %>% arrange(module)

coarseheatmaps<-list()
heatmapDEs<-list()
for (ii in celltypes){
subset<-subset(covid.integrated,idents = ii)
Idents(subset)<-subset$aggcondition.rn
balsubset<-subset(subset,idents = c("BAL_H","BAL_M","BAL_S"))
pbmcsubset<-subset(subset,idents = c("PBMC_H","PBMC_M","PBMC_S"))



balsubset<-AverageExpression(balsubset,return.seurat = TRUE,assays = "SCT")

pbmcsubset<-AverageExpression(pbmcsubset,return.seurat = TRUE,assays = "SCT")

subset<-AverageExpression(subset,return.seurat = TRUE,assays = "SCT")

# using gene results from full group testing to filter differential genes at the cell type level
balgenes<-c(rownames(filter(covid.DE.results.coarse[[ii]]$DE.bal,g3_2p_val_adj <1e-7,g2_1p_val_adj <1e-7,g3_1p_val_adj <1e-7)))
pbmcgenes<-c(rownames(filter(covid.DE.results.coarse[[ii]]$DE.pbmc,g3_2p_val_adj <.05,g2_1p_val_adj <.05,g3_1p_val_adj <.05)))

DE.genes<-c(balgenes,pbmcgenes[!(pbmcgenes %in% balgenes)])
DE.genes<- DE.genes[DE.genes %in% DEG.module.results$modules$id]

DE.genes<-c(DE.genes)

##store DE info from heatmaps for a seperate table

balDEs<-covid.DE.results.coarse[[ii]]$DE.bal[DE.genes,]
  balDEs<-balDEs[,-(1:3)]
  colnames(balDEs)<-paste0("BAL_",colnames(balDEs))


pbmcDEs<-covid.DE.results.coarse[[ii]]$DE.pbmc[DE.genes,]
  pbmcDEs<-pbmcDEs[,-(1:3)]
  colnames(pbmcDEs)<-paste0("PBMC_",colnames(pbmcDEs))

heatmapDEs[[ii]]<-cbind(balDEs,pbmcDEs)

# write.xlsx(heatmapDEs[[ii]], file="heatmapDEs.xlsx",
                 # col.names=TRUE, row.names=TRUE, append = T,sheetName = ii)  

# sweep(typetotal,2,as.matrix(typetotal["Total",]),"/")[1:15,]

#send to coarseheatmap to get clustering, then back to seurat for nice looking figure
zbal<-DoHeatmap(balsubset,features =c(DE.genes), size = 3, draw.lines = FALSE)
zpbmc<-DoHeatmap(pbmcsubset,features =c(DE.genes), size = 3, draw.lines = FALSE)

z<-DoHeatmap(subset,features =c(DE.genes), size = 3, draw.lines = FALSE)

# this will use separately scaled BAL and PBMC data
z$data<-rbind(zbal$data,zpbmc$data)

y <- z$data %>% drop_na()
x <- y %>% group_by(Identity) %>% dplyr::select(Feature, Cell, Identity, Expression) %>%
  tidyr::spread(key = Feature, value = Expression)
w <- y %>% dplyr::select(Feature, Cell, Expression) %>%
  tidyr::spread(key = Cell, value = Expression) %>% column_to_rownames("Feature") %>% as.matrix()

rdeganno<-data.frame(rDEG=matrix(NA,nrow=nrow(w),ncol = 1))
rownames(rdeganno)<-rownames(w)
rdeganno[rownames(rdeganno) %in% balgenes,1]<-"bal"
rdeganno[rownames(rdeganno) %in% pbmcgenes,1]<-"pbmc"
rdeganno[(rownames(rdeganno) %in% pbmcgenes)&(rownames(rdeganno) %in% balgenes),1]<-"both"
avgexp<-subset@assays$SCT@counts[rownames(rdeganno),]
avgexp<-data.frame(bal=rowMeans(avgexp[,c("BAL_H","BAL_M","BAL_S")]),pbmc=rowMeans(avgexp[,c("PBMC_H","PBMC_M","PBMC_S")]))
avgexp$avg<-rowMeans(avgexp)
#Z transform
avgexp<-sweep(avgexp[,1:2],1,as.matrix(avgexp$avg),"/")


rdeganno$"Exp Z Diff"<-avgexp$bal-avgexp$pbmc

col_fun<-colorRamp2(c(-2,0,2),c("purple","black","gold"))
col_fun2<-colorRamp2(c(-2,0,2),RColorBrewer::brewer.pal(9,"PiYG")[c(1,5,9)])

ha<-rowAnnotation(df=rdeganno,col=list(rDEG=c("bal"="blue","pbmc"="red","both"="purple"),`Exp Z Diff` = col_fun2))
H<-Heatmap(w, right_annotation = ha,col=col_fun,cluster_columns = FALSE,
           column_split = c("BAL","BAL","BAL","PBMC","PBMC","PBMC"),name = "Exp",column_title = ii)


coarseheatmaps[[ii]]<-DoHeatmap(subset,features = rev(DE.genes)[row_order(H)],size = 3, draw.lines = FALSE)

coarseheatmaps[[ii]]<-H
}

## featureplots

genestoplot<-c("NFKBIA","NUPR1","BCL2A1","MTRNR2L12","CTSL","CTSB")
FeaturePlot(covid.integrated,features = genestoplot,split.by = "aggcondition.rn")&NoAxes()


#deg exploration
celltypes<-c("M1 MoMa","M2 MoMa","NK","CD4 T", "CD8 Memory T")
geneoi<-"NUPR1"
geneofinterest<-data.frame(row.names = colnames(covid.DE.results.coarse$`M1 MoMa`$DE.bal)[-c(1:3)])

for (ii in celltypes){
  ballabel<-paste0("BAL ",ii)
  pbmclabel<-paste0("PBMC ",ii)
  geneofinterest[,ballabel]<-t(covid.DE.results.coarse[[ii]]$DE.bal[geneoi,-c(1:3)])
  geneofinterest[,pbmclabel]<-t(covid.DE.results.coarse[[ii]]$DE.pbmc[geneoi,-c(1:3)])
  
}  

t(geneofinterest)


##### Figure 4 #####

#ridgeplots
VlnPlot(covid.integrated,features=c("NEAT1","MALAT1"),ncol = 1,group.by = "aggcondition.rn",pt.size = 0)
RidgePlot(covid.integrated,features=c("NEAT1","MALAT1"),ncol = 1,group.by = "aggcondition.rn")

highlightplots<-list()
#NEAT1's high logFC
neatfc<-rDEG.all$bal[,c("gene","avg_logFC_g3_2")]
neatfc<-aggregate(neatfc[,2],list(neatfc$gene),mean) %>% arrange(desc(abs(x)))
neatfc$Group.1<-factor(neatfc$Group.1, levels = neatfc$Group.1)
neatfc$fill<-"nofill"
neatfc$fill[neatfc$Group.1=="NEAT1"]<-"fill"

#covert to log2
neatfc$x<-log2(exp(neatfc$x))

highlightplots$neatfc<-ggplot(neatfc[1:10,],aes(x=Group.1,y=x,fill=fill))+geom_col()+scale_fill_manual(values = c("nofill"="gray","fill"="tomato"),guide = FALSE)+
  labs(x="",y="Severe - Mild logFC",title="BAL rDEG Top LogFC")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

#MALAT1 in pbmc
malatfc<-rDEG.all$pbmc[,c("gene","avg_logFC_g3_2")]
malatfc<-aggregate(malatfc[,2],list(malatfc$gene),mean) %>% arrange(desc(abs(x)))
malatfc$Group.1<-factor(malatfc$Group.1,levels = malatfc$Group.1)
malatfc$fill<-"nofill"
malatfc$fill[malatfc$Group.1=="MALAT1"]<-"fill"

#convert natural log to log2 for graphing
malatfc$x<-log2(exp(malatfc$x))

highlightplots$malatfc<-ggplot(malatfc[1:10,],aes(x=Group.1,y=x,fill=fill))+geom_col()+scale_fill_manual(values = c("nofill"="gray","fill"="tomato"),guide = FALSE)+
  labs(x="",y="Severe - Mild logFC",title="PBMC rDEG Top LogFC")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


###rDEG counts plotted
balcount<-bal.sig.genes.coarse[bal.sig.genes.coarse$balgenescoarse %in% DEG.module.results$modules$id,]
pbmccount<-pbmc.sig.genes.coarse[pbmc.sig.genes.coarse$pbmcgenescoarse %in% DEG.module.results$modules$id,]

balcount$balgenescoarse<-factor(balcount$balgenescoarse,levels = balcount$balgenescoarse)
pbmccount$pbmcgenescoarse<-factor(pbmccount$pbmcgenescoarse,levels = pbmccount$pbmcgenescoarse)

balcount$fill<-"nofill"
balcount$fill[balcount$balgenescoarse=="NEAT1"]<-"fill"

pbmccount$fill<-"nofill"
pbmccount$fill[pbmccount$pbmcgenescoarse=="MALAT1"]<-"fill"


highlightplots$neatcount<-ggplot(balcount[1:5,],aes(x=balgenescoarse,y=Freq,fill=fill))+geom_col()+scale_y_continuous(breaks = c(1:9))+labs(x="",y="Frequency",title="BAL rDEG Top Frequency")+scale_fill_manual(values = c("nofill"="gray","fill"="tomato"),guide = FALSE)

highlightplots$malatcount<-ggplot(pbmccount[1:5,],aes(x=pbmcgenescoarse,y=Freq,fill=fill))+geom_col()+scale_y_continuous(breaks = c(1:9))+labs(x="",y="Frequency",title="PBMC rDEG Top Frequency")+scale_fill_manual(values = c("nofill"="gray","fill"="tomato"),guide = FALSE)

do.call("ggarrange",c(highlightplots[c(3,4,1,2)],ncol=2,nrow=2))

##
lncdotplots<-list()
for (ii in 1:6){
  subset<-covid.integrated[,covid.integrated$aggcondition.rn == levels(covid.integrated$aggcondition.rn)[ii]]
  if (ii==1){
    lncdotplots[[ii]]<-DotPlot(subset,features = c("NEAT1","MALAT1"),group.by = "cell.type.merge.subgroup")+NoLegend()+labs(title=levels(covid.integrated$aggcondition.rn)[ii])
  }else if(ii==6){
    lncdotplots[[ii]]<-DotPlot(subset,features = c("NEAT1","MALAT1"),group.by = "cell.type.merge.subgroup")+labs(title=levels(covid.integrated$aggcondition.rn)[ii])+
      theme(axis.line.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())
  }else{
    lncdotplots[[ii]]<-DotPlot(subset,features = c("NEAT1","MALAT1"),group.by = "cell.type.merge.subgroup")+NoLegend()+labs(title=levels(covid.integrated$aggcondition.rn)[ii])+
    theme(axis.line.y = element_blank(),axis.text.y = element_blank(),axis.title.y = element_blank(),axis.ticks.y = element_blank())
  }
}

do.call("ggarrange",c(lncdotplots,ncol=6))

FeaturePlot(covid.integrated,features = c("NEAT1","MALAT1"),split.by = "aggcondition.rn")&NoAxes()

##
Idents(covid.integrated)<-covid.integrated$cell.type.merge.subgroup
celltypes<-c("M1-a","M1-b","M1-c","M1-d",
             "M1-e","M1-f","M2-a","M2-b","M2-c","Int-a","Int-b",
             "CD4 T-a","CD4 T-b","Treg","CD8 T","CD8 Memory T-a","CD8 Memory T-b",
             "NK")
detailgenes<-c("NEAT1","MALAT1","NFKBIA","FOS","CTSL","BCL2A1","MTRNR2L12","NUPR1")

subgroupheatmaps<-list()
subgroupridgeplots<-list()
subgroupviolinplots<-list()
for (ii in celltypes){
  subset<-subset(covid.integrated,idents = ii)
  Idents(subset)<-subset$aggcondition.rn
  subsetname<-unique(covid.integrated$cell.type.merge[covid.integrated$cell.type.merge.subgroup == ii])
  
  # using gene results from full group testing to filter differential genes at the cell type level
  balgenes<-c(rownames(filter(covid.DE.results[[subsetname]]$DE.bal,g3_2p_val_adj <1e-6,g2_1p_val_adj <1e-6,g3_1p_val_adj <1e-6)))
  pbmcgenes<-c(rownames(filter(covid.DE.results[[subsetname]]$DE.pbmc,g3_2p_val_adj <.05,g2_1p_val_adj <.05,g3_1p_val_adj <.05)))
  
  DE.genes<-c(balgenes,pbmcgenes[!(pbmcgenes %in% balgenes)])
  DE.genes<- DE.genes[DE.genes %in% detailgenes]
  
  #send to subgroupheatmap to get clustering, then back to seurat for nice looking figure
  if (length(DE.genes > 0)){
    subgroupridgeplots[[ii]]<-RidgePlot(subset,features =DE.genes, ncol = 2)+ylab(ii)
    subgroupviolinplots[[ii]]<-VlnPlot(subset,features = DE.genes, ncol = 2,pt.size = 0)+ylab(ii)
    subset<-AverageExpression(subset,return.seurat = TRUE,assays = "SCT")
    subgroupheatmaps[[ii]]<-DoHeatmap(subset,features = DE.genes,  size = 3, draw.lines = FALSE )+ylab(ii)
    }
}

grid.arrange(subgroupheatmaps[[1:9]],ncol = 3)


#cell group dotplot
subset<-subset(covid.integrated,idents = c("M1 MoMa","M2 MoMa","NK","CD4 T", "CD8 Memory T"))
DotPlot(subset,features = detailgenes, group.by = "cell.type.merge.subgroup", cols = c("darkorchid3","chartreuse3"))+
  theme(axis.text.x = element_text(size = 10,angle = 45,vjust = 1.05,hjust = 1))


##### Supplementals #####

### file 1 markers
write.xlsx(covid.integrated.markers, file="celltypemarkers.xlsx",
col.names=TRUE, row.names=TRUE, append = F,sheetName = "Cell Type Markers")

### file 2 
names(coarsepropsig)<-gsub("Macrophage","MoMa",names(coarsepropsig))
for (ii in names(coarsepropsig)){
  write.xlsx(coarsepropsig[[ii]], file="cellproportionpvalues.xlsx",
             col.names=TRUE, row.names=TRUE, append = T,sheetName = gsub("/","_",ii))
  
}

### coarse DEs printed
for (ii in names(covid.DE.results.coarse)){
  write.xlsx(covid.DE.results.coarse[[ii]]$DE.bal, file="coarseDEsbal.xlsx",
             col.names=TRUE, row.names=TRUE, append = T,sheetName = paste0("bal ",gsub("/","_",ii)))
  write.xlsx(covid.DE.results.coarse[[ii]]$DE.pbmc, file="coarseDEspbmc.xlsx",
             col.names=TRUE, row.names=TRUE, append = T,sheetName = paste0("pbmc ",gsub("/","_",ii)))
  
}



##### Cross DEG Supplemental #####
names(covid.DE.cross.coarse)<-gsub("Macrophage","MoMa",names(covid.DE.cross.coarse))

for (ii in names(covid.DE.cross.coarse)){
  for (jj in names(covid.DE.cross.coarse$`M1 MoMa`$DE.cross)){
  write.xlsx(covid.DE.cross.coarse[[ii]]$DE.cross[[jj]], file="coarseDEcross.xlsx",
             col.names=TRUE, row.names=TRUE, append = T,sheetName = paste0(jj," ",gsub("/","_",ii)))
  }
}


## Table 6
allrDEG<-rDEG.all$combined
allrDEG$celltype<-gsub("Macrophage","MoMa",allrDEG$celltype)
write.xlsx(allrDEG, file="rDEGcombined.xlsx",
           col.names=TRUE, row.names=TRUE, append = F, sheetName = "All rDEGs by Cell Type")

##### Validation DEG Heatmaps #####
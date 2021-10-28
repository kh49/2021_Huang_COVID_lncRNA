library(monocle3)
library(Seurat)
library(tidyverse)
library(topGO)

set.seed(0)
#Extract count data, phenotype data, and feature data from the Seurat Object.
counts.data <- covid.integrated@assays$SCT@data
pheno.data <- covid.integrated@meta.data
feature.data <- data.frame(gene_short_name = row.names(counts.data), row.names = row.names(counts.data))


#Construct a CellDataSet.
#construction function now includes size factor estimation
covid.SCT.cds <- new_cell_data_set(counts.data, cell_metadata = pheno.data, gene_metadata = feature.data)

#this includes log normalization
covid.SCT.cds <- preprocess_cds(covid.SCT.cds,method = "PCA",
                                num_dim =50, norm_method = "none")

covid.SCT.cds@preprocess_aux@listData$gene_loadings <- covid.integrated@reductions$pca@feature.loadings
covid.SCT.cds@preprocess_aux@listData$prop_var_expl<-covid.integrated@reductions$pca@stdev

reducedDims(covid.SCT.cds)[["Aligned"]] <- reducedDims(covid.SCT.cds)[["PCA"]]

covid.SCT.cds <- reduce_dimension(covid.SCT.cds)

#Import UMAP coordinates from Seurat
s.umap <- covid.integrated@reductions$umap@cell.embeddings
reducedDims(covid.SCT.cds)$UMAP<- s.umap

plot_cells(covid.SCT.cds, color_cells_by="cell.type.merge")

# #generate colData with all conditions from both studies this section has been added to the main seurat analysis and is depreciated
# aggcondition<-colData(covid.SCT.cds)$group
# # aggcondition[aggcondition == "C"] <- "S"
# aggcondition[is.na(aggcondition)]<-colData(covid.SCT.cds)$Status[is.na(aggcondition)]
# colData(covid.SCT.cds)$aggcondition<-aggcondition
# 
# aggpatient<-colData(covid.SCT.cds)$patientseverity
# aggpatient[is.na(aggpatient)]<-paste0("PBMC_",colData(covid.SCT.cds)$Donor[is.na(aggpatient)])
# colData(covid.SCT.cds)$aggpatient<-aggpatient

##removing m3 because it only has 364 cells - already completed in seurat
# covid.SCT.cds<-covid.SCT.cds[,covid.SCT.cds@colData$aggpatient!="M3"]

### exporting aggregates for WGCNA
### aggregation takes raw counts, defaults to log transformed and scaled data
aggexpfine<-list()

for(ii in levels(covid.SCT.cds@colData$cell.type.0.5)){
 
  export.cds<-covid.SCT.cds[,colData(covid.SCT.cds)$cell.type.0.5 %in% ii]

  cell_group_df <- tibble::tibble(cell=row.names(colData(export.cds)), 
                                  cell_group=export.cds@colData$aggpatient)
  aggexpfine[[ii]]<-aggregate_gene_expression(export.cds,cell_group_df = cell_group_df)
  
}




##### modules for T regs #####

treg.cds<- covid.SCT.cds[,grepl("T Regulatory", colData(covid.SCT.cds)$cell.type.coarse.merge, ignore.case=TRUE)]
tregresults<-moduletest(treg.cds)
map(tregresults$tukey, ~ .x %>% filter(`p adj` < .05))
#modules significant between M and S none, SM and H 17,16,15,9,6,5,2   , SH 11,8,1

##### modules for CD8+ T cells#####

cd8pos.cds<- covid.SCT.cds[,grepl("CD8+", colData(covid.SCT.cds)$cell.type.merge, ignore.case=TRUE)]

cd8posresults<-moduletest(cd8pos.cds)
map(cd8posresults$tukey, ~ .x %>% filter(`p adj` < .05))
#modules significant between M and S ,SM and H  18,16,15,10,9,6,2, SH 12,11

##### modules for double negative t cells#####

dnt.cds<- covid.SCT.cds[,grepl("TCF7+", colData(covid.SCT.cds)$cell.type.merge, ignore.case=TRUE)]

dntresults<-moduletest(dnt.cds)
map(dntresults$tukey, ~ .x %>% filter(`p adj` < .05))
#modules significant between M and S , SM and H 14,12,10,7,2  , SH  11 , MH 4

##### modules for FCN1- Monocytes#####

fcn1neg.cds<- covid.SCT.cds[,grepl("FCN1-", colData(covid.SCT.cds)$cell.type.merge, ignore.case=TRUE)]

fcn1negresults<-moduletest(fcn1neg.cds)

map(fcn1negresults$tukey, ~ .x %>% filter(`p adj` < .05))
#modules significant between M and S 14,4 SM and H 18,5  , SH 10,14,41
#11 is sig between pbmc H nad MS



##### modules for FCN1+ Monocytes#####

fcn1pos.cds<- covid.SCT.cds[,grepl("FCN1+", colData(covid.SCT.cds)$cell.type.merge, ignore.case=TRUE)]

fcn1posresults<-moduletest(fcn1pos.cds)
map(fcn1posresults$tukey, ~ .x %>% filter(`p adj` < .05))
#modules significant between M and S 16,14,9,7 SM and H  18,12,5,1 , SH 16,14,10,8,4,2
#pbmc s and h - 6

##### modules for MT+ Monocytes#####

mtpos.cds<- covid.SCT.cds[,grepl("MT+", colData(covid.SCT.cds)$cell.type.merge, ignore.case=TRUE)]

mtposresults<-moduletest(mtpos.cds)
map(mtposresults$tukey, ~ .x %>% filter(`p adj` < .05))

#modules significant between M and S 14,8, SM and H  7, SH 14,8


##### modules for FCN1- Monocytes#####

ccl2macro.cds<- covid.SCT.cds[,grepl("CCL2", colData(covid.SCT.cds)$cell.type.merge, ignore.case=TRUE)]

ccl2macroresults<-moduletest(ccl2macro.cds)

map(ccl2macroresults$tukey, ~ .x %>% filter(`p adj` < .05))

##### modules for CD14-/FCN1\\+/FKBP4\\+/IER5\\+ Monocytes#####

FKBP4IER5mono.cds<- covid.SCT.cds[,grepl("CD14-/FCN1\\+/FKBP4\\+/IER5\\+ Monocytes", colData(covid.SCT.cds)$cell.type.merge, ignore.case=TRUE)]

FKBP4IER5monoresults<-moduletest(FKBP4IER5mono.cds)

FKBP4IER5monoresults$sigtukey<-map(FKBP4IER5monoresults$tukey, ~ .x %>% filter(`p adj` < .05))

##### Automated Module Generation for all cell types #####
celltypes<-gsub("\\+","\\\\+",levels(covid.SCT.cds@colData$cell.type.merge.coarse))

covid.module.results<-list()  

for (ii in celltypes){
  subset.cds<- covid.SCT.cds[,grepl(ii, colData(covid.SCT.cds)$cell.type.merge.coarse, ignore.case=TRUE)]
  
  covid.module.results[[gsub("\\\\","\\",ii)]]<-moduletest(subset.cds)
  
}
qsave(covid.module.results,"covid.module.results")



##### cell processing #####
tukey<-bind_rows(covid.module.results$`MALAT1+ Monocytes`$tukey, .id = "column_label")
tukey<-tukey[tukey$`p adj`<=.075,]
write.xlsx(tukey,"/mnt/g/My Drive/Finn-Perkins Lab/COVID/Monocle Modules/MALAT1+ Monocytes/MALAT1+ Monocytes.xlsx")


##### module correlation analysis #####
library(Hmisc)
library(corrplot)
library(circlize)
library(igraph)

#generate modules based on UMAP graph of genes
# gene_module_df <- DEG.module.results$modules
gene_module_df <- find_gene_modules(covid.SCT.cds[,colData(covid.SCT.cds)$cell.type.merge.coarse %in%
                                                    c("M1 Macrophage","M2 Macrophage", "Int Macrophage")],core =16,max_components =50, resolution = .001,random_seed = 1)

module_enrichment<-list()

for (ii in 1:length(levels(gene_module_df$module))){
  
  genelist <- numeric(length = dim(gene_module_df)[1]) 
  names(genelist) <- gene_module_df$id
  
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

#first just genes into modules to build correlation matrices, one for each macrophage group
agg_mod <- aggregate_gene_expression(covid.SCT.cds, gene_module_df)
row.names(agg_mod) <- stringr::str_c("Module ", row.names(agg_mod))
modulenames<- c("Cell Motility","Neutrophil/Interferon/Inflammation",
                "Redox Metabolism","Cell Division","Small Molecule Metabolism",
                "Protein Processing","Viral Transcription")

cor <- list(NA,NA,NA)
heatmaps<-list(NA,NA,NA)
anova<-list(NA,NA,NA)
for (ii in c("M1 Macrophage","M2 Macrophage","Int Macrophage")){
    
    par(mfrow=c(1,1))
    cor[[ii]]<-list(NA,NA,NA)

    sub.cds<-covid.SCT.cds[,(colData(covid.SCT.cds)$cell.type.merge.coarse == ii)&
                             (colData(covid.SCT.cds)$aggcondition %in% c("H","M","S"))]
    
    cell_group_df <- tibble::tibble(cell=row.names(colData(sub.cds)), 
                                    cell_group=sub.cds@colData$aggcondition)
    
    agg_heat<-aggregate_gene_expression(sub.cds,
                                        gene_group_df = gene_module_df,cell_group_df = cell_group_df)
    
    row.names(agg_heat) <- modulenames
    heatmaps[[ii]]<-pheatmap::pheatmap(agg_heat,
                       scale="column", clustering_method="ward.D2", cluster_rows = FALSE, cluster_cols = FALSE)
    
    #ANOVA
    
    cell_group_df <- tibble::tibble(cell=row.names(colData(sub.cds)), 
                                    cell_group=sub.cds@colData$aggpatient)
    agg_mat_patient <- aggregate_gene_expression(sub.cds, gene_module_df, cell_group_df)
    row.names(agg_mat_patient) <- stringr::str_c("Module ", row.names(agg_mat_patient))
    
    patient2condition<- gsub("[0-9]","",colnames(agg_mat_patient))
    
    agg_mat_patient_noscale <- aggregate_gene_expression(sub.cds, gene_module_df, cell_group_df,scale_agg_values = FALSE)
    
    
    moduleanova<-list()
    tukey<-list()
    
    for (ll in 1:length(levels(gene_module_df$module))){
      
      module<-as.data.frame(cbind(agg_mat_patient[ll,],patient2condition))
      colnames(module)<-c("exp","condition")
      module$exp<-as.numeric(module$exp)
      
      #logFC
      balmean<-aggregate(module$exp, by = list(module$condition), FUN = mean)
      rownames(balmean)<-balmean$Group.1
      
      
      
      # Compute the analysis of variance
      BAL.aov <- aov(exp ~ condition, data = module)
      
      # PostHoc via Tukey
      BAL.ph<-TukeyHSD(BAL.aov)
 
      tukeyres<-as.data.frame(BAL.ph$condition)
      
      ##switch to un-scaled modules for logFC
      module<-as.data.frame(cbind(agg_mat_patient_noscale[ll,],patient2condition))
      colnames(module)<-c("exp","condition")
      module$exp<-as.numeric(module$exp)
      
      #logFC
      logfc<-data.frame(rep(0,3),row.names = rownames(tukeyres))
      colnames(logfc)<-"fc"
      means<-aggregate(module$exp, by = list(module$condition), FUN = mean)
      rownames(means)<-means$Group.1
      means$log<-log(means$x+1)
      
      logfc["M-H","fc"]<-means["M","log"]-means["H","log"]
      logfc["S-H","fc"]<-means["S","log"]-means["H","log"]
      logfc["S-M","fc"]<-means["S","log"]-means["M","log"]
      
      tukeyres$avg_logFC<-logfc[rownames(tukeyres),"fc"]
      
      moduleanova[[ll]]<-BAL.aov
      
      
      tukey[[ll]]<-tukeyres
    }
    anova[[ii]]<-list(moduleanova,tukey)
    
    par(mfrow=c(1,3))
    for (jj in c("H","M","S")){
  cor[[ii]][[jj]] <- rcorr(t(agg_mod[,(colData(covid.SCT.cds)$cell.type.merge.coarse == ii)&
                                 (colData(covid.SCT.cds)$aggcondition == jj)]),type = "pearson")
  
  crit_agg<-aggregate_gene_expression(covid.SCT.cds[,(colData(covid.SCT.cds)$cell.type.merge.coarse == ii)&
                                                      (colData(covid.SCT.cds)$aggcondition == jj)],
                                                      gene_module_df)
  
  row.names(crit_agg) <- modulenames

  
  M <- cor[[ii]][[jj]]$r
  colnames(M)<-rownames(M)<-modulenames ##generated later in module analysis
  p_mat <- cor[[ii]][[jj]]$P
  p_adj <- matrix(p.adjust(p_mat, method = "fdr"), nrow = nrow(p_mat))
  diag(p_adj)<-1
  logp<--log10(p_adj)
  logp[logp==Inf]<-10
  logp[logp<=1.30103]<-0
  diag(logp)<-1
  
  # corrplot(M, type = "upper", order = "hclust",
  #          p.mat = p_adj, sig.level = 0.05,title = names(macrophagecor)[ii], diag = FALSE)
  
  # col_fun = colorRamp2(range(M), c("#4B0082", "#FFD700"), transparency = .025)
  # # col_ramp = colorRampPalette(colors = c("#4B008220", "#FFD70080"),alpha = TRUE)(9)
  # 
  # 
  # col_M <- col_fun(M)
  # 
  # # col_M[p_adj>.05] <-  sub("..$","10",col_M[p_adj>.05])
  # col_M[abs(M)<.25] <-  sub("..$","00",col_M[abs(M)<.25])
  # 
  # # p_col<-p_mat
  # # diag(p_col) <- 0
  # # 
  # # for (ii in seq(0,1,.125)){
  # #   p_col[p_col<=ii]<-col_ramp[1+ii*8]
  # # }
  # # p_col[p_col>=1]<-col_ramp[1]
  # 
  # chordDiagram(M, col = col_M, symmetric = TRUE)
  # Create igraph
  M_graph<-graph_from_adjacency_matrix(
    M, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  p_graph<-graph_from_adjacency_matrix(
    logp/4+1, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  
  # Function for mapping correlations to colors
  make_colors <- function(wgt) {  
    cr <- colorRamp(c("#AF8DC3", "white", "#7FBF7B"), space = "Lab")
    v <- cr((wgt + 1)/2) # map linearly to [0, 1]
    return(rgb(v[,1], v[,2], v[, 3], maxColorValue = 255))
  }
  
  # Tweaks for some pretty plotting 
  E(M_graph)$width <- 2*(E(p_graph)$weight-1)^2
  E(M_graph)$color <- make_colors(E(M_graph)$weight)
  V(M_graph)$label.color <- "black"
  V(M_graph)$label.cex <- 1.2
  V(M_graph)$size <- 10
  V(M_graph)$color <- "white"
  plot(M_graph, layout = layout_in_circle(M_graph))
  title(main = jj,sub = ii)
  
  }
}

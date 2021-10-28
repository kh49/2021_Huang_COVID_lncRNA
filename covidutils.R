library(dplyr)
library(ggplot2)
library(ggrepel)
library(qs)
library(topGO)
library(monocle3)
library(Seurat)
library(tidyverse)
library(xlsx)
library(gridExtra)
library(STRINGdb)
library(igraph)


saveDExlsx<-function(DE,file){
  wb <- createWorkbook()
  sheet1 <- createSheet(wb, sheetName = "BAL.DE.sig")
  sheet2 <- createSheet(wb, sheetName = "PBMC.DE.sig")
  addDataFrame(DE$bal, sheet = sheet1)
  addDataFrame(DE$pbmc, sheet = sheet2)
  saveWorkbook(wb, file = file)
  
}

qsaveworkspace<-function(){
  files<-sapply(setdiff(ls(1), lsf.str(1)),as.name)
  for (ii in 1:length(files)){
      qsave(eval(files[[ii]]), file = paste0("./qData/",as.character(files[[ii]])), nthreads = 20)
  }
  
  qsave(files,file = "./qData/files")
}

qloadworkspace<-function(){
  files<-as.character(qread("./qData/files"))
  for(ii in 1:length(files))
  assign(files[ii],qread(paste0("./qData/",gsub("`","",files[ii])),nthreads = 20),envir = globalenv())
  
}
##### per cell type gene SNN graph exploration #####
gene.graph<-function(covid.integrated,covid.DE.results,idents){
  ##### inverted matrix Gene network test #####
  subset<-subset(covid.integrated, idents = idents)
  
  ##list of unique DEGs
  genes<-NULL
  for(ii in idents){
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
  
  cells<-WhichCells(covid.integrated,idents = idents)
  
  covid.invert<-CreateSeuratObject(t(subset@assays$SCT@counts[genes,cells]),assay= "RNA", min.cells = 0, min.features = 0)
  covid.invert <- NormalizeData(covid.invert, normalization.method = "LogNormalize", scale.factor = 10000)
  covid.invert <- FindVariableFeatures(covid.invert, selection.method = "vst", nfeatures = 5000)
  all.genes <- rownames(covid.invert)
  covid.invert <- ScaleData(covid.invert, features = all.genes)
  covid.invert <- RunPCA(covid.invert, features = VariableFeatures(object = covid.invert))
  ElbowPlot(covid.invert)
  covid.invert <- RunUMAP(covid.invert, dims = 1:10)
  
  covid.invert <- FindNeighbors(covid.invert, dims = 1:10,prune.SNN = .35)
  covid.invert <- FindClusters(covid.invert, resolution = .5)
  
  # DimPlot(covid.invert,reduction = "umap")
  # 
  # 
  # plot(sort(colSums(covid.invert@graphs$RNA_snn),decreasing=TRUE))
  
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
  plot<-plot.igraph(x = net, layout = as.matrix(x = Embeddings(object = covid.invert[["umap"]])),
              edge.width = E(graph = net)$weight, vertex.label = nodelabel,
              vertex.size = 0)
  
  return(list(nodescores=nodescores,graph=net,plot=plot))
}

ppi.string<-function (covid.DE.results,idents, score_threshold = 400, speciesID = 9606) 
{
    ##list of unique DEGs
    genes<-NULL
    for(ii in idents){
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
    
    
    for_gen <- unique(sort(genes))
    new_genes <- as.data.frame(for_gen, stringsAsFactors = FALSE)
    names(new_genes) <- "gene"
    database <- STRINGdb$new(version = "11", species = speciesID, 
                             score_threshold = score_threshold, input_directory = "./")
    mapped <- database$map(new_genes, "gene", removeUnmappedRows = TRUE)
    mapped <- mapped[!duplicated(mapped$STRING_id), ]
    interactions <- database$get_interactions(mapped$STRING_id)
    graph_ppi <- data.frame(interactions$from, interactions$to, 
                                  stringsAsFactors = FALSE)

    for (n in seq_len(nrow(graph_ppi))) {
      graph_ppi$interactions.from[n] <- mapped[graph_ppi$interactions.from[n] == 
                                                 mapped$STRING_id, ][1]
    }
    for (n in seq_len(nrow(graph_ppi))) {
      graph_ppi$interactions.to[n] <- mapped[graph_ppi$interactions.to[n] == 
                                               mapped$STRING_id, ][1]
    }
    edge_list <- as.matrix(as.vector(graph_ppi[, 1], mode = "character"))
    edge_list <- cbind(edge_list, as.vector(graph_ppi[, 2], 
                                            mode = "character"))
    final_graph <- graph.edgelist(edge_list, directed = FALSE)

  return(final_graph)
}

##### 3factor MAST DE functions #####
threefactorMASTDETest <- function(
  data.use,
  cells.1,
  cells.2,
  cells.3,
  latent.vars = NULL,
  verbose = TRUE,
  ...
) {
  
  if (length(x = latent.vars) > 0) {
    latent.vars <- scale(x = latent.vars)
  }
  group.info <- data.frame(row.names = c(cells.1, cells.2, cells.3))
  latent.vars <- latent.vars %||% group.info
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[cells.3, "group"] <- "Group3"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  latent.vars.names <- c("condition", colnames(x = latent.vars))
  latent.vars <- cbind(latent.vars, group.info)
  latent.vars$wellKey <- rownames(x = latent.vars)
  fdat <- data.frame(rownames(x = data.use))
  colnames(x = fdat)[1] <- "primerid"
  rownames(x = fdat) <- fdat[, 1]
  sca <- MAST::FromMatrix(
    exprsArray = as.matrix(x = data.use),
    cData = latent.vars,
    fData = fdat
  )
  to.return<-list()
  
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group1")
  SummarizedExperiment::colData(sca)$condition <- cond
  fmla <- as.formula(
    object = paste0(" ~", paste(latent.vars.names, collapse = "+"))
  )
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  
  summaryCond <- summary(object = zlmCond, doLRT = 'conditionGroup2')
  summaryDt <- as.data.frame(summaryCond$datatable)
  p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
  genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  to.return$g2_1 <- data.frame(p_val, row.names = genes.return)
  
  summaryCond <- summary(object = zlmCond, doLRT = 'conditionGroup3')
  summaryDt <- as.data.frame(summaryCond$datatable)
  p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
  genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  to.return$g3_1 <- data.frame(p_val, row.names = genes.return)
  
  ##reset ZLM with alternate reference to do group3-group2
  cond <- factor(x = SummarizedExperiment::colData(sca)$group)
  cond <- relevel(x = cond, ref = "Group2")
  SummarizedExperiment::colData(sca)$condition <- cond
  zlmCond <- MAST::zlm(formula = fmla, sca = sca, ...)
  
  summaryCond <- summary(object = zlmCond, doLRT = 'conditionGroup3')
  summaryDt <- as.data.frame(summaryCond$datatable)
  p_val <- summaryDt[summaryDt[, "component"] == "H", 4]
  genes.return <- summaryDt[summaryDt[, "component"] == "H", 1]
  to.return$g3_2 <- data.frame(p_val, row.names = genes.return)
  return(to.return)
}

FindDEseurat <- function(
  object,
  ident.1 = NULL,
  ident.2 = NULL,
  ident.3 = NULL,
  group.by = NULL,
  subset.ident = NULL,
  assay = NULL,
  slot = 'data',
  reduction = NULL,
  features = NULL,
  logfc.threshold = 0.25,
  test.use = "MAST",
  min.pct = 0.1,
  min.diff.pct = -Inf,
  verbose = TRUE,
  max.cells.per.ident = Inf,
  random.seed = 1,
  latent.vars = NULL,
  min.cells.feature = 3,
  min.cells.group = 3,
  pseudocount.use = 1,
  ...
) {
  if (!is.null(x = group.by)) {
    if (!is.null(x = subset.ident)) {
      object <- subset(x = object, idents = subset.ident)
    }
    Idents(object = object) <- group.by
  }
  if (!is.null(x = assay) && !is.null(x = reduction)) {
    stop("Please only specify either assay or reduction.")
  }
  data.slot <- ifelse(
    test = test.use %in% c("negbinom", "poisson", "DESeq2"),
    yes = 'counts',
    no = slot
  )
  if (is.null(x = reduction)) {
    assay <- assay %||% DefaultAssay(object = object)
    data.use <-  GetAssayData(object = object[[assay]], slot = data.slot)
  } else {
    if (data.slot == "counts") {
      stop("The following tests cannot be used when specifying a reduction as they assume a count model: negbinom, poisson, DESeq2")
    }
    data.use <- t(x = Embeddings(object = object, reduction = reduction))
  }
  if (is.null(x = ident.1)) {
    stop("Please provide ident.1")
  } 
  if (length(x = as.vector(x = ident.1)) > 1 &&
      any(as.character(x = ident.1) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(x = !as.character(x = ident.1) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.1 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    ident.1 <- WhichCells(object = object, idents = ident.1)
  }
  
  if (length(x = as.vector(x = ident.2)) > 1 &&
      any(as.character(x = ident.2) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.2) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.2 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.2)) {
      stop("Please provide ident.2")
    } else {
      ident.2 <- WhichCells(object = object, idents = ident.2)
    }
  }
  
  if (length(x = as.vector(x = ident.3)) > 1 &&
      any(as.character(x = ident.3) %in% colnames(x = data.use))) {
    bad.cells <- colnames(x = data.use)[which(!as.character(x = ident.3) %in% colnames(x = data.use))]
    if (length(x = bad.cells) > 0) {
      stop(paste0("The following cell names provided to ident.3 are not present in the object: ", paste(bad.cells, collapse = ", ")))
    }
  } else {
    if (is.null(x = ident.3)) {
      stop("Please provide ident.3")
    } else {
      ident.3 <- WhichCells(object = object, idents = ident.3)
    }
  }
  
  if (!is.null(x = latent.vars)) {
    latent.vars <- FetchData(
      object = object,
      vars = latent.vars,
      cells = c(ident.1, ident.2, ident.3)
    )
  }
  counts <- switch(
    EXPR = data.slot,
    'scale.data' = GetAssayData(object = object[[assay]], slot = "counts"),
    numeric()
  )
  
  cells.1<-ident.1
  cells.2<-ident.2
  cells.3<-ident.3
  object<-data.use
  
  features <- features %||% rownames(x = object)
  
  # error checking
  if (length(x = cells.1) == 0) {
    stop("Cell group 1 is empty - no cells with identity class ", cells.1)
  } else if (length(x = cells.2) == 0) {
    stop("Cell group 2 is empty - no cells with identity class ", cells.2)
    return(NULL)
  } else if (length(x = cells.3) == 0) {
    stop("Cell group 3 is empty - no cells with identity class ", cells.3)
    return(NULL)
  }	else if (length(x = cells.1) < min.cells.group) {
    stop("Cell group 1 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.2) < min.cells.group) {
    stop("Cell group 2 has fewer than ", min.cells.group, " cells")
  } else if (length(x = cells.3) < min.cells.group) {
    stop("Cell group 3 has fewer than ", min.cells.group, " cells")
  } else if (any(!cells.1 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.1) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.1 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.2 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.2) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.2 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  } else if (any(!cells.3 %in% colnames(x = object))) {
    bad.cells <- colnames(x = object)[which(x = !as.character(x = cells.3) %in% colnames(x = object))]
    stop(
      "The following cell names provided to cells.3 are not present: ",
      paste(bad.cells, collapse = ", ")
    )
  }
  # feature selection (based on percentages)
  data <- switch(
    EXPR = slot,
    'scale.data' = counts,
    object
  )
  if (is.null(x = reduction)) {
    thresh.min <- 0
    pct.1 <- round(
      x = rowSums(x = data[features, cells.1, drop = FALSE] > thresh.min) /
        length(x = cells.1),
      digits = 3
    )
    pct.2 <- round(
      x = rowSums(x = data[features, cells.2, drop = FALSE] > thresh.min) /
        length(x = cells.2),
      digits = 3
    )
    pct.3 <- round(
      x = rowSums(x = data[features, cells.3, drop = FALSE] > thresh.min) /
        length(x = cells.3),
      digits = 3
    )
    data.alpha <- cbind(pct.1, pct.2, pct.3)
    colnames(x = data.alpha) <- c("pct.1", "pct.2", "pct.3")
    alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
    names(x = alpha.min) <- rownames(x = data.alpha)
    features <- names(x = which(x = alpha.min > min.pct))
    if (length(x = features) == 0) {
      stop("No features pass min.pct threshold")
    }
    alpha.diff <- alpha.min - apply(X = data.alpha, MARGIN = 1, FUN = min)
    features <- names(
      x = which(x = alpha.min > min.pct & alpha.diff > min.diff.pct)
    )
    if (length(x = features) == 0) {
      stop("No features pass min.diff.pct threshold")
    }
  } else {
    data.alpha <- data.frame(
      pct.1 = rep(x = NA, times = length(x = features)),
      pct.2 = rep(x = NA, times = length(x = features)),
      pct.3 = rep(x = NA, times = length(x = features))
    )
  }
  # feature selection (based on average difference)
  mean.fxn <- if (is.null(x = reduction) && slot != "scale.data") {
    switch(
      EXPR = slot,
      'data' = function(x) {
        return(log(x = rowMeans(x = expm1(x = x)) + pseudocount.use))
      },
      function(x) {
        return(log(x = rowMeans(x = x) + pseudocount.use))
      }
    )
  } else {
    rowMeans
  }
  data.1 <- mean.fxn(data[features, cells.1, drop = FALSE])
  data.2 <- mean.fxn(data[features, cells.2, drop = FALSE])
  data.3 <- mean.fxn(data[features, cells.3, drop = FALSE])
  
  total.diff <- apply(cbind((data.2 - data.1),(data.3-data.1),(data.3-data.2)),1,function(x){x[which.max(abs(x))]})
  
  if (is.null(x = reduction) && slot != "scale.data") {
    features.diff <- names(x = which(x = abs(x = total.diff) > logfc.threshold))
    
    features <- intersect(x = features, y = features.diff)
    if (length(x = features) == 0) {
      stop("No features pass logfc.threshold threshold")
    }
  }
  if (max.cells.per.ident < Inf) {
    set.seed(seed = random.seed)
    # Should be cells.1 and cells.2?
    if (length(x = cells.1) > max.cells.per.ident) {
      cells.1 <- sample(x = cells.1, size = max.cells.per.ident)
    }
    if (length(x = cells.2) > max.cells.per.ident) {
      cells.2 <- sample(x = cells.2, size = max.cells.per.ident)
    }
    if (length(x = cells.3) > max.cells.per.ident) {
      cells.3 <- sample(x = cells.3, size = max.cells.per.ident)
    }
    if (!is.null(x = latent.vars)) {
      latent.vars <- latent.vars[c(cells.1, cells.2, cells.3), , drop = FALSE]
    }
  }
  # perform DE
  if (!(test.use %in% c('negbinom', 'poisson', 'MAST', "LR")) && !is.null(x = latent.vars)) {
    warning(
      "'latent.vars' is only used for 'negbinom', 'poisson', 'LR', and 'MAST' tests",
      call. = FALSE,
      immediate. = TRUE
    )
  }
  if (!test.use %in% c('wilcox', 'MAST', 'DESeq2')) {
    CheckDots(...)
  }
  de.results <- switch(
    EXPR = test.use,
    'wilcox' = stop("Only MAST is supported", test.use),
    'bimod' = stop("Only MAST is supported", test.use),
    'roc' = stop("Only MAST is supported", test.use),
    't' = stop("Only MAST is supported", test.use),
    'negbinom' = stop("Only MAST is supported", test.use),
    'poisson' = stop("Only MAST is supported", test.use),
    'MAST' = threefactorMASTDETest(
      data.use = object[features, c(cells.1, cells.2, cells.3), drop = FALSE],
      cells.1 = cells.1,
      cells.2 = cells.2,
      cells.3 = cells.3,
      latent.vars = latent.vars,
      verbose = verbose,
      ...
    ),
    "DESeq2" = stop("Only MAST is supported", test.use),
    "LR" = stop("Only MAST is supported", test.use),
    stop("Unknown test: ", test.use)
  )
  if (is.null(x = reduction)) {
    diff.col <- ifelse(
      test = slot == "scale.data" || test.use == 'roc',
      yes = "avg_diff",
      no = "avg_logFC"
    )
    de.results.tmp<-arrange(de.results$g2_1,p_val)
    colnames(de.results.tmp)<-"g2_1_p_val"
    de.results.tmp$g3_1_p_val<-de.results$g3_1[rownames(de.results.tmp),]
    de.results.tmp$g3_2_p_val<-de.results$g3_2[rownames(de.results.tmp),]
    de.results <- de.results.tmp
    de.results <- cbind(de.results,avg_logFC_g2_1=as.data.frame(data.2-data.1)[rownames(de.results),],avg_logFC_g3_1=as.data.frame(data.3-data.1)[rownames(de.results),]
                        ,avg_logFC_g3_2=as.data.frame(data.3-data.2)[rownames(de.results),])
    de.results <- cbind(de.results, data.alpha[rownames(x = de.results), , drop = FALSE])
  } else {
    diff.col <- "avg_diff"
    de.results[, diff.col] <- total.diff[rownames(x = de.results)]
  }
  if (test.use == "roc") {
    de.results <- de.results[order(-de.results$power, -de.results[, diff.col]), ]
  } else {
    de.results <- de.results[order(de.results$g3_2_p_val, -de.results[, "avg_logFC_g3_2"]), ]
    de.results$g2_1p_val_adj = p.adjust(
      p = de.results$g2_1_p_val,
      method = "BH",
      n = nrow(x = object)
      
    )
    de.results$g3_1p_val_adj = p.adjust(
      p = de.results$g3_1_p_val,
      method = "BH",
      n = nrow(x = object)
      
    )
    de.results$g3_2p_val_adj = p.adjust(
      p = de.results$g3_2_p_val,
      method = "BH",
      n = nrow(x = object)
      
    )
  }
  
  return(de.results)
}

summarize.DE<-function(DE){
  sigbal<-DE$DE.bal
  sigbal$'M-H'<-DE$volcano.plots$`M-H`$data$Sig
  sigbal$'S-H'<-DE$volcano.plots$`S-H`$data$Sig
  sigbal$'S-M'<-DE$volcano.plots$`S-M`$data$Sig
  sigbal$importance<-0
  sigbal$importance<-apply(sigbal,1,DE.importance)
  sigbal<-sigbal[sigbal$importance>0,]
  
  sigpbmc<-DE$DE.pbmc
  sigpbmc$'PBMC_M-PBMC_H'<-DE$volcano.plots$`PBMC_M-PBMC_H`$data$Sig
  sigpbmc$'PBMC_S-PBMC_H'<-DE$volcano.plots$`PBMC_S-PBMC_H`$data$Sig
  sigpbmc$'PBMC_S-PBMC_M'<-DE$volcano.plots$`PBMC_S-PBMC_M`$data$Sig
  sigpbmc$importance<-0
  sigpbmc$importance<-apply(sigpbmc,1,DE.importance)
  sigpbmc<-sigpbmc[sigpbmc$importance>0,]  
  
  results<-list()
  results$bal<-arrange(sigbal,desc(importance),g2_1p_val_adj+g3_1p_val_adj)
  results$pbmc<-arrange(sigpbmc,desc(importance),g2_1p_val_adj+g3_1p_val_adj)
  return(results)
}

DE.importance<-function(gene){
  
  if(!is.na(gene[14])|!is.na(gene[15])){
    if((gene[14] == "FC_P")&&(gene[15] == "FC_P")){
      gene[17] <- as.numeric(gene[17]) + 1
    
    }
  }
  if((gene[16]=="FC_P")&&!(is.na(gene[16]))){
    gene[17] <- as.numeric(gene[17]) + 2
  }
  
  return(gene[17])
}
  

moduletest<-function(cds){
  gene_module_df <- find_gene_modules(cds,core =16,max_components = 50, resolution = .001)
  agg_mat_cell <- as(aggregate_gene_expression(cds, gene_module_df,scale_agg_values = FALSE,norm_method = "log"), "dgCMatrix")
  row.names(agg_mat_cell) <- stringr::str_c("Module ", row.names(agg_mat_cell))
  # agg_mat_cell <- t(scale(t(agg_mat_cell),center = FALSE))
  # pheatmap::pheatmap(agg_mat_cell,
  #                                  scale="column", clustering_method="ward.D2")
  
  # aggcondition[aggcondition == "C"] <- "S"
  
  metadata<-list()
  ##how many are there?
  metadata$patientcellcount<-table(cds@colData$aggpatient)
  metadata$conditioncellcount<-table(cds@colData$aggcondition)
  metadata$colData<-colData(cds)
  
  #module expression heatmap
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=cds@colData$aggcondition)
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  conditionmap<-pheatmap::pheatmap(agg_mat,
                     scale="column", clustering_method="ward.D2")
  
  
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=cds@colData$aggpatient)
  agg_mat_patient <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  row.names(agg_mat_patient) <- stringr::str_c("Module ", row.names(agg_mat_patient))
  patientmap<-pheatmap::pheatmap(agg_mat_patient,
                     scale="column", clustering_method="ward.D2")
  
  agg_mat_patient_noscale <- aggregate_gene_expression(cds, gene_module_df, cell_group_df,scale_agg_values = FALSE)
  
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
    
    ##switch to un-scaled modules for logFC
    module<-as.data.frame(cbind(agg_mat_patient_noscale[ii,],patient2condition))
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
  
  return(results)
}

moduleDE<-function(cds,gene_module_df){
  agg_mat_cell <- as(aggregate_gene_expression(cds, gene_module_df,scale_agg_values = FALSE,norm_method = "log"), "dgCMatrix")
  row.names(agg_mat_cell) <- stringr::str_c("Module ", row.names(agg_mat_cell))
  # agg_mat_cell <- t(scale(t(agg_mat_cell),center = FALSE))
  # pheatmap::pheatmap(agg_mat_cell,
  #                                  scale="column", clustering_method="ward.D2")
  
  # aggcondition[aggcondition == "C"] <- "S"
  
  metadata<-list()
  ##how many are there?
  metadata$patientcellcount<-table(cds@colData$aggpatient)
  metadata$conditioncellcount<-table(cds@colData$aggcondition)
  metadata$colData<-colData(cds)
  
  #module expression heatmap
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=cds@colData$aggcondition)
  agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
  conditionmap<-pheatmap::pheatmap(agg_mat,
                                   scale="column", clustering_method="ward.D2")
  
  
  cell_group_df <- tibble::tibble(cell=row.names(colData(cds)), 
                                  cell_group=cds@colData$aggpatient)
  agg_mat_patient <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
  row.names(agg_mat_patient) <- stringr::str_c("Module ", row.names(agg_mat_patient))
  patientmap<-pheatmap::pheatmap(agg_mat_patient,
                                 scale="column", clustering_method="ward.D2")
  
  agg_mat_patient_noscale <- aggregate_gene_expression(cds, gene_module_df, cell_group_df,scale_agg_values = FALSE)
  
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
    
    ##switch to un-scaled modules for logFC
    module<-as.data.frame(cbind(agg_mat_patient_noscale[ii,],patient2condition))
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
  
  
  results<-list()
  results[["module_anova"]]<-moduleanova
  results[["condition_heatmap"]]<-conditionmap
  results[["patient_heatmap"]]<-patientmap
  names(results$module_anova)<-names(results$module_enrichment)
  results[["tukey"]]<-tukey
  results[["sigtukey"]]<-map(tukey, ~ .x %>% filter(`p adj` < .05))
  results[["metadata"]]<-metadata
  results[["aggregates"]]<-list(cell = agg_mat_cell, condition = agg_mat, patient = agg_mat_patient)
  
  return(results)
}

#from seurat LabelClusters
LabelPlotBackground<-function (plot, id, clusters = NULL, labels = NULL, split.by = NULL, 
                               repel = TRUE, ...) 
{
  xynames <- unlist(x = GetXYAesthetics(plot = plot), use.names = TRUE)
  if (!id %in% colnames(x = plot$data)) {
    stop("Cannot find variable ", id, " in plotting data")
  }
  if (!is.null(x = split.by) && !split.by %in% colnames(x = plot$data)) {
    warning("Cannot find splitting variable ", id, " in plotting data")
    split.by <- NULL
  }
  data <- plot$data[, c(xynames, id, split.by)]
  possible.clusters <- as.character(x = na.omit(object = unique(x = data[, 
                                                                         id])))
  groups <- clusters %||% as.character(x = na.omit(object = unique(x = data[, 
                                                                            id])))
  if (any(!groups %in% possible.clusters)) {
    stop("The following clusters were not found: ", paste(groups[!groups %in% 
                                                                   possible.clusters], collapse = ","))
  }
  labels.loc <- lapply(X = groups, FUN = function(group) {
    data.use <- data[data[, id] == group, , drop = FALSE]
    data.medians <- if (!is.null(x = split.by)) {
      do.call(what = "rbind", args = lapply(X = unique(x = data.use[, 
                                                                    split.by]), FUN = function(split) {
                                                                      medians <- apply(X = data.use[data.use[, split.by] == 
                                                                                                      split, xynames, drop = FALSE], MARGIN = 2, 
                                                                                       FUN = median, na.rm = TRUE)
                                                                      medians <- as.data.frame(x = t(x = medians))
                                                                      medians[, split.by] <- split
                                                                      return(medians)
                                                                    }))
    }
    else {
      as.data.frame(x = t(x = apply(X = data.use[, xynames, 
                                                 drop = FALSE], MARGIN = 2, FUN = median, na.rm = TRUE)))
    }
    data.medians[, id] <- group
    return(data.medians)
  })
  labels.loc <- do.call(what = "rbind", args = labels.loc)
  labels <- labels %||% groups
  if (length(x = unique(x = labels.loc[, id])) != length(x = labels)) {
    stop("Length of labels (", length(x = labels), ") must be equal to the number of clusters being labeled (", 
         length(x = labels.loc), ").")
  }
  names(x = labels) <- groups
  for (group in groups) {
    labels.loc[labels.loc[, id] == group, id] <- labels[group]
  }
  geom.use <- ifelse(test = repel, yes = geom_label_repel, no = geom_label)
  plot <- plot + geom.use(data = labels.loc, mapping = aes_string(x = xynames["x"], 
                                                                  y = xynames["y"], label = id), ...)
  return(plot)
}

"%||%" <- devtools:::`%||%`

###pulled vrom seurat visualization.R
GetXYAesthetics <- function(plot, geom = 'GeomPoint', plot.first = TRUE) {
  geoms <- sapply(
    X = plot$layers,
    FUN = function(layer) {
      return(class(x = layer$geom)[1])
    }
  )
  geoms <- which(x = geoms == geom)
  if (length(x = geoms) == 0) {
    stop("Cannot find a geom of class ", geom)
  }
  geoms <- min(geoms)
  if (plot.first) {
    x <- as.character(x = plot$mapping$x %||% plot$layers[[geoms]]$mapping$x)[2]
    y <- as.character(x = plot$mapping$y %||% plot$layers[[geoms]]$mapping$y)[2]
  } else {
    x <- as.character(x = plot$layers[[geoms]]$mapping$x %||% plot$mapping$x)[2]
    y <- as.character(x = plot$layers[[geoms]]$mapping$y %||% plot$mapping$y)[2]
  }
  return(list('x' = x, 'y' = y))
}


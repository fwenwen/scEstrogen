# Figure3
options(stringsAsFactors = FALSE)
library(plyr)
library(permute)
library(data.table)
library(hdf5r)
library(SCopeLoomR)
library(SCENIC)
library(data.table)
library(pbapply)
library(plyr)
library(philentropy)
library(ggplot2)
library(ggrepel)
library(latex2exp)

Organ.group.list <- SplitObject(seuratObj,split.by = "organ.group")
for(i in names(Organ.group.list)){
  Organ.group.list[[i]] <- droplevels(Organ.group.list[[i]])
  refGene <- seuratObj@misc$geneName
  if(!(any(is.na(rownames(Organ.group.list[[i]]@assays$SCT@counts))))){
    rownames(Organ.group.list[[i]]@assays$SCT@counts) <- refGene[rownames(Organ.group.list[[i]]@assays$SCT@counts)]
  }
}

for(i in names(Organ.group.list)){
  setwd(sprintf("./SCENIC/%s",i))
  dir.create("data")
  dir.create("output")
  subobj <- Organ.group.list[[i]]
  saveRDS(subobj,sprintf("./%s.rds",i))
  mat <- subobj@assays$SCT@counts
  mat[1:3,1:3]
  cellinfo <- subobj@meta.data
  head(cellinfo)
  
  length(colnames(mat))
  length(rownames(cellinfo))
  table(colnames(mat) %in% rownames(cellinfo))
  
  mat <- as.matrix(mat)
  dim(mat)
  saveRDS(mat,sprintf("./data/%s_mat.rds",i))
  
  colnames(cellinfo)
  cellinfo <- cellinfo[,c("organ.group","seurat_clusters","final","minor2")]
  head(cellinfo)
  saveRDS(cellinfo,sprintf("./data/cellinfo_%s.rds",i))
    loom <- build_loom("./data/exprMat.loom", dgem=mat)
  saveLoom(mat, "./data/s1_exprMat.loom")
}

myupset <- function (data, nsets = 5, nintersects = 40, sets = NULL, keep.order = F, 
                     set.metadata = NULL, intersections = NULL, matrix.color = "gray23", mat_col=NULL,
                     main.bar.color = "gray23", mainbar.y.label = "Intersection Size", 
                     mainbar.y.max = NULL, sets.bar.color = "#1170AA", sets.x.label = "Set Size", 
                     point.size = 2.2, line.size = 0.7, mb.ratio = c(0.7, 0.3), 
                     expression = NULL, att.pos = NULL, att.color = main.bar.color, 
                     order.by = c("freq", "degree"), decreasing = c(T, F), show.numbers = "yes", 
                     number.angles = 0, group.by = "degree", cutoff = NULL, queries = NULL, 
                     query.legend = "none", shade.color = "gray88", shade.alpha = 0.25, 
                     matrix.dot.alpha = 0.5, empty.intersections = NULL, color.pal = 1, 
                     boxplot.summary = NULL, attribute.plots = NULL, scale.intersections = "identity", 
                     scale.sets = "identity", text.scale = 1, set_size.angles = 0, 
                     set_size.show = FALSE, set_size.numbers_size = NULL, set_size.scale_max = NULL)  {
  
  startend <- UpSetR:::FindStartEnd(data)
  first.col <- startend[1]
  last.col <- startend[2]
  if (color.pal == 1) {
    palette <- c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", 
                 "#9467BD", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22", 
                 "#17BECF")
  }
  else {
    palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                 "#0072B2", "#D55E00", "#CC79A7")
  }
  if (is.null(intersections) == F) {
    Set_names <- unique((unlist(intersections)))
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <-UpSetR:::Number_of_sets(Set_names)
    if (keep.order == F) {
      Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::specific_intersections(data, first.col, 
                                                 last.col, intersections, order.by, group.by, decreasing, 
                                                 cutoff, main.bar.color, Set_names)
  }
  else if (is.null(intersections) == T) {
    Set_names <- sets
    if (is.null(Set_names) == T || length(Set_names) == 0) {
      Set_names <- UpSetR:::FindMostFreq(data, first.col, last.col, 
                                         nsets)
    }
    Sets_to_remove <- UpSetR:::Remove(data, first.col, last.col, Set_names)
    New_data <- UpSetR:::Wanted(data, Sets_to_remove)
    Num_of_set <- UpSetR:::Number_of_sets(Set_names)
    if (keep.order == F) {
      Set_names <- UpSetR:::order_sets(New_data, Set_names)
    }
    All_Freqs <- UpSetR:::Counter(New_data, Num_of_set, first.col, 
                                  Set_names, nintersects, main.bar.color, order.by, 
                                  group.by, cutoff, empty.intersections, decreasing)
  }
  Matrix_setup <- UpSetR:::Create_matrix(All_Freqs)
  labels <- UpSetR:::Make_labels(Matrix_setup)
  att.x <- c()
  att.y <- c()
  if (is.null(attribute.plots) == F) {
    for (i in seq_along(attribute.plots$plots)) {
      if (length(attribute.plots$plots[[i]]$x) != 0) {
        att.x[i] <- attribute.plots$plots[[i]]$x
      }
      else if (length(attribute.plots$plots[[i]]$x) == 
               0) {
        att.x[i] <- NA
      }
      if (length(attribute.plots$plots[[i]]$y) != 0) {
        att.y[i] <- attribute.plots$plots[[i]]$y
      }
      else if (length(attribute.plots$plots[[i]]$y) == 
               0) {
        att.y[i] <- NA
      }
    }
  }
  BoxPlots <- NULL
  if (is.null(boxplot.summary) == F) {
    BoxData <- UpSetR:::IntersectionBoxPlot(All_Freqs, New_data, first.col, 
                                            Set_names)
    BoxPlots <- list()
    for (i in seq_along(boxplot.summary)) {
      BoxPlots[[i]] <- UpSetR:::BoxPlotsPlot(BoxData, boxplot.summary[i], 
                                             att.color)
    }
  }
  customAttDat <- NULL
  customQBar <- NULL
  Intersection <- NULL
  Element <- NULL
  legend <- NULL
  EBar_data <- NULL
  if (is.null(queries) == F) {
    custom.queries <- UpSetR:::SeperateQueries(queries, 2, palette)
    customDat <- UpSetR:::customQueries(New_data, custom.queries, 
                                        Set_names)
    legend <- UpSetR:::GuideGenerator(queries, palette)
    legend <- UpSetR:::Make_legend(legend)
    if (is.null(att.x) == F && is.null(customDat) == F) {
      customAttDat <- UpSetR:::CustomAttData(customDat, Set_names)
    }
    customQBar <- UpSetR:::customQueriesBar(customDat, Set_names, 
                                            All_Freqs, custom.queries)
  }
  if (is.null(queries) == F) {
    Intersection <- UpSetR:::SeperateQueries(queries, 1, palette)
    Matrix_col <- intersects(UpSetR:::QuerieInterData, Intersection, 
                             New_data, first.col, Num_of_set, All_Freqs, expression, 
                             Set_names, palette)
    Element <- UpSetR:::SeperateQueries(queries, 1, palette)
    EBar_data <- UpSetR:::ElemBarDat(Element, New_data, first.col, 
                                     expression, Set_names, palette, All_Freqs)
  }
  else {
    Matrix_col <- NULL
  }
  if (!is.null(mat_col)) {
    Matrix_col <- mat_col
  }
  Matrix_layout <- UpSetR:::Create_layout(Matrix_setup, matrix.color, 
                                          Matrix_col, matrix.dot.alpha)
  Set_sizes <- UpSetR:::FindSetFreqs(New_data, first.col, Num_of_set, 
                                     Set_names, keep.order)
  Bar_Q <- NULL
  if (is.null(queries) == F) {
    Bar_Q <- intersects(UpSetR:::QuerieInterBar, Intersection, New_data, 
                        first.col, Num_of_set, All_Freqs, expression, Set_names, 
                        palette)
  }
  QInter_att_data <- NULL
  QElem_att_data <- NULL
  if ((is.null(queries) == F) & (is.null(att.x) == F)) {
    QInter_att_data <- intersects(UpSetR:::QuerieInterAtt, Intersection, 
                                  New_data, first.col, Num_of_set, att.x, att.y, expression, 
                                  Set_names, palette)
    QElem_att_data <- elements(UpSetR:::QuerieElemAtt, Element, New_data, 
                               first.col, expression, Set_names, att.x, att.y, palette)
  }
  AllQueryData <- UpSetR:::combineQueriesData(QInter_att_data, QElem_att_data, 
                                              customAttDat, att.x, att.y)
  ShadingData <- NULL
  if (is.null(set.metadata) == F) {
    ShadingData <- UpSetR:::get_shade_groups(set.metadata, Set_names, 
                                             Matrix_layout, shade.alpha)
    output <- UpSetR:::Make_set_metadata_plot(set.metadata, Set_names)
    set.metadata.plots <- output[[1]]
    set.metadata <- output[[2]]
    if (is.null(ShadingData) == FALSE) {
      shade.alpha <- unique(ShadingData$alpha)
    }
  }
  else {
    set.metadata.plots <- NULL
  }
  if (is.null(ShadingData) == TRUE) {
    ShadingData <- UpSetR:::MakeShading(Matrix_layout, shade.color)
  }
  Main_bar <- suppressMessages(UpSetR:::Make_main_bar(All_Freqs, Bar_Q, 
                                                      show.numbers, mb.ratio, customQBar, number.angles, EBar_data, 
                                                      mainbar.y.label, mainbar.y.max, scale.intersections, 
                                                      text.scale, attribute.plots))
  Matrix <- UpSetR:::Make_matrix_plot(Matrix_layout, Set_sizes, All_Freqs, 
                                      point.size, line.size, text.scale, labels, ShadingData, 
                                      shade.alpha)
  Sizes <- UpSetR:::Make_size_plot(Set_sizes, sets.bar.color, mb.ratio, 
                                   sets.x.label, scale.sets, text.scale, set_size.angles, 
                                   set_size.show, set_size.scale_max, set_size.numbers_size)
  structure(class = "upset", .Data = list(Main_bar = Main_bar, 
                                          Matrix = Matrix, Sizes = Sizes, labels = labels, mb.ratio = mb.ratio, 
                                          att.x = att.x, att.y = att.y, New_data = New_data, expression = expression, 
                                          att.pos = att.pos, first.col = first.col, att.color = att.color, 
                                          AllQueryData = AllQueryData, attribute.plots = attribute.plots, 
                                          legend = legend, query.legend = query.legend, BoxPlots = BoxPlots, 
                                          Set_names = Set_names, set.metadata = set.metadata, set.metadata.plots = set.metadata.plots))
}
library(ggthemes)
library(UpSetR)
# head(fromList(UnciliatedAll.list))
library(ggpubr)
pair.list <- list(Normal=Normal.pair,EN=EN.pair,EH=EH.pair)

u1 <- myupset(fromList(link.list), order.by = "freq",
              sets = names(link.list),
              point.size = 3, line.size = 1,
              mainbar.y.label = "Intersectionh size", 
              sets.x.label = "Gene set size",
              # sets.bar.color=subgroupcols[c(4,2,1,7,3,5,6)],
              main.bar.color= "#9C755F",
              matrix.color ="#FF9DA7")
u1
pdf("./upsetplot.pdf", height=4.27, width=6.27)
print(u1)
dev.off()

Idents(seuratObj) <- seuratObj$organ.group
ras_exp_scatter_sample <- function(seuratObj, rasMat, gene, reduction = 'umap') {
  Idents(seuratObj) <- seuratObj$organ.group
  plist1 <- list()
  for(i in levels(seuratObj$organ.group)){
    plist1[[i]] <- FeaturePlot(seuratObj, reduction = reduction, order=TRUE,
                               cols=c('lightgrey', pals::brewer.blues(10)),
                               coord.fixed = F,cells = WhichCells(seuratObj,idents = i),
                               pt.size=0.2, combine=T, features = gene) + NoAxes() + NoLegend() +
      theme(plot.title = element_blank())
      ggtitle(gene) + guides(color=guide_colorbar(title="Expression"))
  }
  p1 <- FeaturePlot(seuratObj, reduction = reduction, order=TRUE,
                    cols=c('lightgrey', pals::brewer.blues(10)),
                    coord.fixed = F,cells = WhichCells(seuratObj,idents = "Normal"),
                    pt.size=0.2, combine=T, features = gene) + NoAxes() + #NoLegend() +
    theme(plot.title = element_blank())
    ggtitle(gene) + guides(color=guide_colorbar(title="Expression"))

  p2 <- FeaturePlot(seuratObj, reduction = reduction, order=TRUE,
                    cols=c('lightgrey', pals::brewer.blues(10)),
                    coord.fixed = F,cells = WhichCells(seuratObj,idents = "Estrogen"),
                    pt.size=0.2, combine=T, features = gene) + NoAxes() + #NoLegend() +
    theme(plot.title = element_blank())
    ggtitle(gene) + guides(color=guide_colorbar(title="Expression"))
  
  seuratObj$Regulon <- rasMat[rownames(seuratObj@meta.data), gene]
  plist2 <- list()
  for(j in levels(seuratObj$organ.group)){
    plist2[[j]] <- FeaturePlot(seuratObj, reduction = reduction, order=TRUE, coord.fixed = F,
                               cols=c("lightgrey", pals::brewer.orrd(10)),
                               pt.size=0.2, combine=T,cells = WhichCells(seuratObj,idents = j),
                               features="Regulon") + NoAxes() + NoLegend() +
      theme(plot.title = element_blank())
      ggtitle(gene) + guides(color=guide_colorbar(title="Regulon Activity"))
    
  }
  p3 <- FeaturePlot(seuratObj, reduction = reduction, order=TRUE, coord.fixed = F,
                    cols=c("lightgrey", pals::brewer.orrd(10)),
                    pt.size=0.2, combine=T,cells = WhichCells(seuratObj,idents = "Normal"),
                    features="Regulon") + NoAxes() + #NoLegend() +
    theme(plot.title = element_blank())
    ggtitle(gene) + guides(color=guide_colorbar(title="Regulon Activity"))
  p4 <- FeaturePlot(seuratObj, reduction = reduction, order=TRUE, coord.fixed = F,
                    cols=c("lightgrey", pals::brewer.orrd(10)),
                    pt.size=0.2, combine=T,cells = WhichCells(seuratObj,idents = "Estrogen"),
                    features="Regulon") + NoAxes() + #NoLegend() +
    theme(plot.title = element_blank())
    ggtitle(gene) + guides(color=guide_colorbar(title="Regulon Activity"))
  pp <- ggarrange(plotlist = c(plist2,plist1),ncol = 10,nrow = 2)
  return(pp)
}
ras_exp_scatter(seuratObj = seuratObj, rasMat = rasMat_merge, gene = "SOX13", reduction = "umap")

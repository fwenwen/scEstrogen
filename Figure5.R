library(Seurat)
library(dplyr)
library(future)
library(glmGamPoi)
library(ggforce)
library(pals)
library(RColorBrewer)
library(ggplot2)
library(scatterpie)
library(ggpubr)
library(readxl)
library(stringr)
library(ggthemes)
library(pals)

# geneID transform
head(rownames(Lung@assays$SCT@data))
refGenes <- Lung@misc$geneName
if(!any(is.na(refGenes[rownames(Lung@assays$SCT@data)]))){
  rownames(Lung@assays$SCT@data) <- refGenes[rownames(Lung@assays$SCT@data)]
}
head(rownames(Lung@assays$SCT@data))

head(rownames(Lung@assays$SCT@counts))
refGenes <- Lung@misc$geneName
if(!any(is.na(refGenes[rownames(Lung@assays$SCT@counts)]))){
  rownames(Lung@assays$SCT@counts) <- refGenes[rownames(Lung@assays$SCT@counts)]
}
head(rownames(Lung@assays$SCT@counts))
###-------------------------5A
Lung <- readRDS("/public/workspace/wanggenjun/work/scRNA/estrogen/zhao/LungAnalysis/Lung_integrated.rds")
DefaultAssay(Lung) <- 'SCT'
Idents(Lung) <- 'final'
Lung <- RunUMAP(Lung, dims=1:50, min.dist=0.4, n.neighbors=80, verbose=FALSE) ##  dims=1:20, min.dist=0.5

cell_eb <- data.frame(Lung@reductions$umap@cell.embeddings[,1:2],cluster=as.character(Lung$final),
                      group=as.character(Lung$group),stringsAsFactors = F) 
colnames(cell_eb) <- c('x','y','cluster','group') 

pp2 <- ggplot(data=cell_eb, aes(x,y)) + geom_point(aes(colour = factor(cluster)), size=0.25) +
  scale_color_manual(values =Colors$Cluster.cols) + theme_pubr() +
  NoAxes() + NoLegend() 
pp2

#medium pie
df <- Lung@meta.data %>% dplyr::count(final, group) %>% dplyr::rename(cluster = final) %>% dplyr::rename(group = group)
df$n <- df$n/as.numeric(table(Lung@meta.data$group)[df$group]) # normal 组的细胞数量
coords <- cell_eb %>% group_by(cluster) %>% dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% as.data.frame
piedf <- reshape2::acast(df, cluster ~ group, value.var = "n")
piedf[is.na(piedf)] <- 0
piedf <- cbind(piedf, coords[match(rownames(piedf), coords$cluster), -1])
cellnumbers <- dplyr::count(cell_eb, cluster)
sizee <- cellnumbers[match(rownames(piedf), cellnumbers[,1]), 2]
piedf$radius <- scales::rescale(sizee, c(0.4,0.9))
DimPieM <- pp2 + ggnewscale::new_scale_fill() + 
  geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=piedf,
                  cols=colnames(piedf)[1:(ncol(piedf)-3)],color=NA) +
  coord_equal()+ theme_pubr() + NoAxes()  + scale_fill_manual(values = Colors$group.cols)
DimPieM
###-------------------------5B
library(pals)
orange9 <- brewer.pal(9,"YlOrRd")
show_col(orange9)
library(RColorBrewer)
library(igraph)
library(reshape2)
library(scales)
library(tidyverse)
library(xlsx)
library(cpplot)
library(tidyverse)
pvalues_Nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_Normal/pvalues.txt", check.names = FALSE)
means_Nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_Normal/means.txt", check.names = FALSE)
sig.means_Nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_Normal/significant_means.txt", check.names = FALSE)
deconvoluted_Nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_Normal/deconvoluted.txt", check.names = FALSE)

pvalues_EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_EH/pvalues.txt", check.names = FALSE)
means_EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_EH/means.txt", check.names = FALSE)
sig.means_EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_EH/significant_means.txt", check.names = FALSE)
deconvoluted_EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_EH/deconvoluted.txt", check.names = FALSE)

celltype <- c('Epithelial cells','Ciliated epithelial cells','Endothelial cells','Fibroblasts','Smooth muscle cells','Naive T cells','Proliferating T cells',
              'CD4+ T cells','CD8+ T cells','NK cells','B cells','Plasma cells','Neutrophils','pDCs','Mac/Mono')


cellinter_stat <- NULL
library(dplyr)
#Normal
op <- par(mfrow=c(1,1))
for(db in c("Nor")){
  sig_means <- read.delim(sprintf("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_Normal/significant_means.txt",db), check.names = F)
  sum(!is.na(sig_means[,-c(1:12)]))
  intersum <- sort(apply(sig_means[,-c(1:12)], 2, function(x){
    sum(!is.na(x))
  }))
  matx <- matrix(0, nrow=length(Colors$Cluster.cols), ncol=length(Colors$Cluster.cols))
  rownames(matx) <- colnames(matx) <- names(Colors$Cluster.cols)
  intersum <- intersum / sum(intersum)
  lapply(names(intersum), function(i){
    x <- unlist(strsplit(i, split = "\\|"))
    matx[x[1],x[2]] <<- intersum[i]
  })
  ## matx[upper.tri(matx)] <- -matx[lower.tri(matx)]
  corrplot(matx, is.corr = FALSE, 
           method = "square", type = "full", 
           tl.col = "black", order = "original", 
           col = orange9)
}
par(op)

#estrogen
cellinter_stat <- NULL
op <- par(mfrow=c(1,1))
for(db in c("Estrogen")){
  sig_means <- read.delim(sprintf("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Lung_EH/significant_means.txt",db), check.names = F)
  sub_data <- sig_means[rowSums(sig_means[,-c(1:12)], na.rm = T) > 0, c(2,13:ncol(sig_means))] %>% 
    droplevels() %>% dplyr::mutate(group=db)
  cellinter_stat <- rbind(cellinter_stat, sub_data)
  intersum <- sort(apply(sig_means[,-c(1:12)], 2, function(x){
    sum(!is.na(x))
  }))
  matx <- matrix(0, nrow=length(Colors$Cluster.cols), ncol=length(Colors$Cluster.cols))
  rownames(matx) <- colnames(matx) <- names(Colors$Cluster.cols)
  intersum <- intersum / sum(intersum)
  lapply(names(intersum), function(i){
    x <- unlist(strsplit(i, split = "\\|"))
    matx[x[1],x[2]] <<- intersum[i]
  })
corrplot(matx, is.corr = FALSE, 
         method = "square", type = "full", 
         tl.col = "black", order = "original", 
         col = orange9)
}

###-------------------------5C
library(dplyr)
group <- levels(Lung$group) ## Revise
AKData <- list()
for(x in group){
  for(y in group){
    i <- paste(x,y,sep="-")
    if(x != y){
      AKData[[i]] <- list()
      subcells <- Cells(Lung)[Lung$group %in% c(x,y)]  ## Revise
      subObj <- subset(Lung, cells = subcells)  ## Revise
      subObj@meta.data <- droplevels(subObj@meta.data)
      Idents(subObj) <- subObj$group  ## Revise
      subMarkers <- FindAllMarkers(subObj, min.pct=0.1, 
                                   logfc.threshold=0.1, 
                                   return.thresh=1, only.pos=TRUE) %>%
        dplyr::mutate(group=ifelse(cluster==y, "Up", "Down"))
      subMarkers <- subMarkers[-grep("ENSMFAG",subMarkers$gene),]
      subMarkers <- subMarkers[-grep("MT",subMarkers$gene),]
      rownames(subMarkers) <- subMarkers$gene
      DEgenes <- subMarkers %>% dplyr::filter(p_val < 0.05 & avg_log2FC > 0.5) 
      topDE <- DEgenes %>% dplyr::filter(!grepl("MT-",gene)) %>%
        dplyr::group_by(cluster) %>% dplyr::top_n(n=10, wt=avg_log2FC)
      avg <- AverageExpression(subObj, assays="SCT")$SCT[, c(x,y)] %>% log1p()
      colnames(avg) <- c("x","y")
      avg <- avg[rowSums(avg) > 0, ] %>% as.data.frame() %>% rownames_to_column(var='gene') %>%
        dplyr::mutate(group=ifelse(gene %in% rownames(DEgenes), DEgenes[gene, 'group'],"NS")) %>%
        dplyr::mutate(size=ifelse(group=="NS", 0.4, DEgenes[gene,'avg_log2FC'])) %>%
        dplyr::mutate(shape=ifelse(group=="NS", 18, 16))
      rownames(avg) <- avg$gene
      labs <- avg[topDE$gene,]
      AKData[[i]][['avg']] <- avg
      AKData[[i]][['DEgenes']] <- DEgenes
      AKData[[i]][['topDE']] <- topDE
      AKData[[i]][['subMarkers']] <- subMarkers
    }
  }
}
Idents(Lung) <- 'final'
SubExp <- AverageExpression(Lung, assays = 'SCT', slot = 'counts')$SCT
AKTop <- AKTop[AKTop %in% rownames(SubExp)]
SMeanExp <- SubExp[AKTop, ]
table(AKTop %in% rownames(SubExp))
library(ggplot2)
apply(SMeanExp, 1, function(x){(x - min(x)) / max(x)})
SMeanExp <- t(apply(SMeanExp, 1, function(x){(x - min(x)) / max(x)}))
SMeanExp <- reshape2::melt(t(SMeanExp)) %>% dplyr::rename(Stage="Var1", Gene="Var2")
SMeanExp
#visualization
library(ggthemes)
library(ggplot2)
library(ggpubr)
library(ggrepel)
AKScaterPlots <- list()
y <- group[1]
AKTop <- AKDEgenes <- c()
for(x in group[-1]){
  i <- paste(x,y,sep="-")
  avg <- AKData[[i]][['avg']] 
  topDE <- AKData[[i]][['topDE']]
  DEgenes <- AKData[[i]][['DEgenes']]
  AKDEgenes <- rbind(AKDEgenes, DEgenes[,c('gene','group')])
  AKTop <- unique(c(AKTop, topDE$gene))
  labs <- avg[topDE$gene,]
  set.seed(20210106)
  ids <- sort(c(sample(which(avg$group=="NS"),1000),which(avg$group!="NS")))
  AKScaterPlots[[i]] <- avg[ids, ] %>% ggplot(aes(x=x, y=y, color=group, size=size)) + 
    geom_point() + xlab(x) + ylab(y) + theme_few() +
    geom_text_repel(aes(x, y, label=gene, color=group), size=3, data=labs) +
    scale_color_manual(values = c(NS='lightgrey',Down="#EFE118",Up="#4E889B")) #Normal'#4E889B'Estrogen'#EFE118'
}
AKScaterPlots


###-------------------------5D   
data <- read.table('/public/workspace/wanggenjun/work/scRNA/estrogen/zhao/LungAnalysis/data_genes3_expression.txt',header = T,sep='\t')
head(data)
p <- ggplot() +
  geom_half_violin(data = data2[data2$group == 'Normal',],
                   aes(x = gene, y = Exp, fill = group),
                   color = "black",
                   linewidth = 0.4, 
                   draw_quantiles = c(0.5), 
                   scale = 'width') + 
  facet_grid(rows = vars(final), scales = 'free_y') 
p

p2 <- p +
  geom_half_violin(data = data2[data2$group == 'Estrogen',],
                   aes(x = gene, y = Exp, fill = group),
                   color = "black",
                   linewidth = 0.4, 
                   draw_quantiles = c(0.5), 
                   scale = 'width',
                   side='r') + 
  facet_grid(rows = vars(final), scales = 'free_y') 
p2

mytheme <- theme_bw() +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 12, angle = 60,hjust = 1),
        axis.title.y = element_text(size = 12), 
        axis.title.x = element_blank(), 
        strip.background = element_blank(), 
        strip.text.y = element_text(size = 12, angle = 0), 
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 12)) 
p3 <- p2 +
  mytheme +
  scale_fill_manual(values = Colors$group.cols) +
  scale_y_continuous(breaks = seq(0, 11, by = 5)) + 
  labs(y = 'Log Normalized Expression') 
p3

###-------------------------5E/5F
library(foreach)
library(enrichplot)
library(clusterProfiler)
library(org.Mfascicularis.eg.db)
library(data.table)
library(dplyr)
library(stringr)
library(DOSE)
library(ComplexHeatmap)
top.markers_Lung <- FindMarkers(Lung,ident.1='Estrogen',ident.2='Normal',group.by = 'group', min.pct=0.1, logfc.threshold=0.1, return.thresh=0.1)#, only.pos=TRUE) ## return.thresh=0.01, test.use="wilcox" 
write.table(top.markers_Lung,'top.markers_Lung_EvsN.txt',row.names =TRUE,, quote = F, sep = "\t" )
submarkers <- top.markers_Lung %>% filter( pct.1 >= 0.1)
geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$gene), decreasing = T)
ego <- gseGO(geneList     = geneList,
             OrgDb        = org.Mfascicularis.eg.db,
             ont          = "BP",
             keyType      = "GID",
             nPerm        = 10000,
             minGSSize    = 10,
             maxGSSize    = 500,
             pvalueCutoff = 1,
             verbose      = FALSE)

gsea <- gseaplot2(ego, geneSetID = rownames(ego@result)[head(order(ego@result$enrichmentScore))])+
  gseaplot2(ego, geneSetID = rownames(ego@result)[tail(order(ego@result$enrichmentScore))]) #尾
gsea

gseaplot2(ego,
          title = "immune response",  
          "GO:0006955", 
          color="red",
          base_size = 20, 
          subplots = 1:2,
          pvalue_table = T) 
gseaplot2(ego,
          title = "inflammatory response",  
          "GO:0006954",
          color="red", 
          base_size = 20, 
          subplots = 1:2, 
          pvalue_table = T) 

###-------------------------5G
selected_rows <- cellinter_stat2 <- NULL
 
for(db in c("Lung_Normal","Lung_EH")){
  sig_means <- read.delim(sprintf("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/%s/significant_means.txt",db), check.names = F)
  sub_data <- sig_means[rowSums(sig_means[,-c(1:12)], na.rm = T) > 0, c(2,13:ncol(sig_means))] %>%
    droplevels() %>% dplyr::mutate(group=db)
  cellinter_stat2 <- cellinter_stat2 %>% dplyr::bind_rows(sub_data)
  selected_rows <- c(selected_rows, as.character(cellinter_stat2[which(rowSums(apply(cellinter_stat2[,colnames(cellinter_stat2)[-c(1,6)]],1,as.numeric), na.rm=TRUE) > 0),'interacting_pair']))
}

x <- table(selected_rows)
selected_rows <- names(x[x>1])
# 
length(selected_rows)
write(selected_rows, '/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/Lung_rows.1.txt')

selected_columns <- colnames(cellinter_stat2)[grep("Mac/Mono",colnames(cellinter_stat2))]
selected_columns
selected_columns <- c(selected_columns[-grep("Mac/Mono",gsub("\\|.*$","",selected_columns))],
                      selected_columns[grep("Mac/Mono",gsub("\\|.*$","",selected_columns))])


library(ggplot2)
dot_plot <- function(selected_rows = NULL,
                     selected_columns = NULL,
                     means_path = './means.txt',
                     pvalues_path = './pvalues.txt',
                     means_separator = '\t',
                     pvalues_separator = '\t', 
                     silence=F){
  all_pval = read.table(pvalues_path, header=T, stringsAsFactors = F, sep=means_separator, comment.char = '', check.names=F)
  all_means = read.table(means_path, header=T, stringsAsFactors = F, sep=pvalues_separator, comment.char = '', check.names=F)
  
  intr_pairs = all_pval$interacting_pair
  all_pval = all_pval[,-c(1:11)]
  all_means = all_means[,-c(1:11)]
  
  if(is.null(selected_rows)){
    selected_rows = intr_pairs
  }
  
  if(is.null(selected_columns)){
    selected_columns = colnames(all_pval)
  }
  selected_columns <- intersect(selected_columns, colnames(all_pval))
  sel_pval = all_pval[match(selected_rows, intr_pairs), selected_columns]
  sel_means = all_means[match(selected_rows, intr_pairs), selected_columns]
  
  df_names = expand.grid(selected_rows, selected_columns)
  pval = unlist(sel_pval)
  pval[pval==0] = 0.0009
  plot.data = cbind(df_names,pval)
  pr = unlist(as.data.frame(sel_means))
  ## pr[pr==0] = 1
  plot.data = cbind(plot.data,pr) ## log2(pr)
  colnames(plot.data) = c('pair', 'clusters', 'pvalue', 'mean')
  plot.data$clusters <- factor(plot.data$clusters,levels = selected_columns)
  
  if(!silence){
    my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399)
    ## my_palette <- colorRampPalette(rev(brewer.rdylbu(20)), alpha=TRUE)(n=500)
    p <- ggplot(plot.data,aes(x=clusters,y=pair)) +
      geom_point(aes(size=-log10(pvalue),color=mean)) +
      scale_color_gradientn('Mean (Molecule 1, Molecule 2)', colors=my_palette) +
      theme_bw() +
      theme(panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            axis.text=element_text(size=14, colour = "black"),
            axis.text.x = element_text(angle = 90, hjust = 1),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title=element_blank(),
            panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
    p
  }else{
    plot.data
  }
}

cellinter_stat1 <- NULL 
for(db in c("Lung_Normal","Lung_EH")){
  cellinter_stat1 <- rbind(cellinter_stat1, 
                           dot_plot(selected_rows, selected_columns, silence = T, 
                                    means_path = sprintf("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/%s/means.txt",db), 
                                    pvalues_path = sprintf("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/%s/pvalues.txt",db)) %>%
                             dplyr::mutate(clusters=sprintf("%s(%s)",clusters,db))
  )
}

column_order <- unlist(lapply(selected_columns, function(x){
  sprintf("%s(%s)",x,c("Lung_Normal","Lung_EH"))
}))
column_order

cellinter_stat1$clusters <- factor(cellinter_stat1$clusters, levels=column_order)

pdat <- cellinter_stat1 %>% dplyr::select(-mean) %>% 
  tidyr::pivot_wider(names_from=c(clusters), values_from=pvalue) %>% 
  as.data.frame()
dim(pdat)

rownames(pdat) <- pdat$pair

pdat <- pdat %>% dplyr::select(-pair) %>% log10()
pdat <- pdat[rowSums(pdat)!=0,]
pdat[is.na(pdat)]=0
phtm <- pheatmap::pheatmap(-pdat, silent = T, cluster_cols = F)
pair_levels <- rev(phtm$tree_row$labels[phtm$tree_row$order])
cellinter_stat1 <- cellinter_stat1 %>% dplyr::filter(pair %in% pair_levels) %>% 
  dplyr::mutate(pair=factor(pair, levels=pair_levels))
my_palette <- colorRampPalette(c("black", "blue", "yellow", "red"), alpha=TRUE)(n=399) 
my_palette <- colorRampPalette(c("grey","black","blue",rep(kovesi.rainbow(5),time=1:5)), alpha=TRUE)(n=399)
emptyRows <- levels(cellinter_stat1$clusters)[!levels(cellinter_stat1$clusters) %in% cellinter_stat1$clusters]
addItems <- expand.grid(levels(cellinter_stat1$pair), emptyRows) %>% dplyr::rename(pair=Var1, clusters=Var2) %>% 
  dplyr::mutate(pvalue=1, mean=min(cellinter_stat1$mean))
cellinter_stat1 <- rbind(cellinter_stat1, addItems)

p <- ggplot(cellinter_stat1,aes(x=clusters,y=pair)) +
  geom_point(aes(size=-log10(pvalue),color=log1p(mean))) +
  scale_color_gradientn('Log mean (Molecule 1, Molecule 2)', colors=my_palette) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text=element_text(size=14, colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_blank(),
        panel.border = element_rect(size = 0.7, linetype = "solid", colour = "black"))
p



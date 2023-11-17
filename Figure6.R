# pkgs
library(dplyr)
library(ggpie)
library(ggplot2)
library(ggpubr)
library(ggthemes)
library(pals)
library(RColorBrewer)
library(Seurat)
library(stringr)

## GSEGO
library(clusterProfiler)
library(enrichplot)
library(foreach)
library(ggstatsplot)
library(org.Mfascicularis.eg.db)
library(tidyverse)

meta <- read.delim('/public/workspace/wanggenjun/work/scRNA/estrogen/cao/meta.txt',head=F,stringsAsFactors=F)
Breast <- readRDS("/public/workspace/wanggenjun/work/scRNA/estrogen/data/Organ/Breast.rds")

groupCols <- readRDS("/public/workspace/wanggenjun/work/scRNA/estrogen/data/color.list.rds")
minor2.cols <- readRDS("/public/workspace/wanggenjun/work/scRNA/estrogen/data/minor_id.col.rds")

Breast <- RunUMAP(Breast, dims=1:30, min.dist=0.5, 
                  n.neighbors=50, verbose=FALSE, umap.method = "umap-learn",metric = "correlation")

saveRDS(object = Breast, file = "/public/workspace/wanggenjun/work/scRNA/estrogen/cao/Breast/Breast.rds")

setwd('/public/workspace/wanggenjun/work/scRNA/estrogen/cao/Breast')

Breast.mm <- subset(Breast, idents=c('Mye_0','Mye_1','Mye_2','Mye_3','Mye_4'))
Breast.mm@meta.data <- droplevels(Breast.mm@meta.data)
Breast.mm <- RunUMAP(Breast.mm, dims=1:30, min.dist=0.5, n.neighbors=50, verbose=FALSE, umap.method = "umap-learn",metric = "correlation")

# celltype
clusters <- Breast$final
celltype = data.frame(ClusterID=Breast$seurat_clusters, minorID=Breast$minorID, celltype.final=Breast$final, celltype.minor=Breast$minor2, stringsAsFactors = F)
write.csv(celltype,"celltype.csv",row.names = T)
head(celltype)

final.marker <- read.csv("/public/workspace/wanggenjun/work/scRNA/estrogen/data/final_marker.csv")

options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 100)
library(pals)
library(RColorBrewer)

# final
unique(Breast$final)
final.cols <- c('#B6992D','#6B468F','#9E67A0','#5D4EA2','#DDA0A7','#BB4922','#C3784D','#9D0142',
                '#E15759','#4E79A7','#4ea7b0','#245728','#377F32','#8FBA2A')
final <- DimPlot(Breast,label = F,cols = final.cols, 
                       group.by = 'final', pt.size = 0.5, seed = 1)
final

Breast.mm <- RunUMAP(Breast.mm, dims=1:30, min.dist=0.1, n.neighbors=50, verbose=FALSE, umap.method = "umap-learn",metric = "correlation")

options(repr.plot.width = 6, repr.plot.height = 5, repr.plot.res = 100)

minorcols.mm <- c('#FBCF35FF','#ED4C1CFF','#9C7E70FF','#5AC2F1FF','#11776CFF')
umap.mm <- DimPlot(Breast.mm,label = F,cols = minorcols.mm, 
                        group.by = 'minor2', pt.size = 3) + ggtitle("mm")

samplegroup <- data.frame(sample=meta$V1, group=gsub("^.*-","",meta$V1))
samplegroup$group <- str_sub(samplegroup$group,1,6)
samplegroup$group[samplegroup$group == "EH"] <- "Estrogen"
rownames(samplegroup) <- samplegroup$sample
head(samplegroup)

options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 100)

Breast$group <- factor(samplegroup[as.character(Breast$orig.ident),"group"],
                       levels = c("Normal","Estrogen"))
group.cols <- setNames(c("#4E889B","#EFE118"),c("Normal","Estrogen"))

bygroup <- DimPlot(Breast, group.by = "group", cols = groupCols$group.cols, pt.size = 0.5,, seed = 1) +
                    NoAxes() + ggtitle("Breast-bygroup")

options(repr.plot.width = 14, repr.plot.height = 5, repr.plot.res = 100)

mm.gene <- c(
            'S100A9','IL1B',#mye0
             'PSAP','VMO1','FCER1G','GALK2',#mye1
             'S100A10','CD68','FABP4',#mye2
             'C1QC','C1QA','CD74',#mye3
             'CST3',#mye4
             'MYL6','COX1','KRT16','ND5','NKG7',#mye5
             'IRF8','RAMP1','GADD45B')#Mye6
p.mm <- DotPlot(Breast.mm, group.by = "minor2", 
                            features = mm.gene, dot.scale = 12)+
  scale_color_gradientn(colours=rev(brewer.rdylbu(20)), guide = "colourbar") +
  theme(axis.text.x=element_text(angle=60, hjust=1),panel.grid.major = element_line(colour = "grey70", size = 0.2))
p.mm
pdf("markers_mm.pdf", height=5, width=14)
print(p.mm)
dev.off()

plot <- dplyr::count(Breast@meta.data, final, orig.ident)

temp_labels.final <- Breast@meta.data %>%
     group_by(final) %>%
     tally() %>%
     dplyr::rename('cluster' = final)
 temp_labels.final

# by group
p.group.scalediagram <- ggplot(data = plot, aes(fill = orig.ident, n, final))+ 
                     geom_bar(stat = "identity", position="fill", width=0.8,size=0.25)+
                    scale_fill_manual(values = c("#4E889B", "#EFE118"))+ theme_classic()

# by final
p.final.scalediagram <- ggplot(data = plot, aes(fill = final, n, orig.ident))+ 
                geom_bar(stat = "identity", position="fill", width=0.5,size=0.25)+
                scale_fill_manual(values = groupCols$Cluster.cols) + theme_classic()

minor2.mm.pie <- read.delim('/public/workspace/wanggenjun/work/scRNA/estrogen/cao/Breast/minor2.mm.pie.txt',head=T,stringsAsFactors=F)
head(minor2.mm.pie)
unique(minor2.mm.pie)

minor2.mm.pie <- minor2.mm.pie[1:319,]

cols.mm <- c('#FBCF35FF','#ED4C1CFF','#9C7E70FF','#5AC2F1FF','#11776CFF')

rosepie.mm <- ggrosepie(data = minor2.mm.pie, group_key = c("minor2","orig.ident"), 
          count_type = "full",label_info = "all",show_tick = F, 
                        donut_label = F,
                        fill_color = c('#EFE118','#4E889B'),
                       donut_frac=0.2) + theme_pubr() + NoLegend() + NoAxes()
    
rosepie.mm

library(cpplot)
library(Cairo)
Cairo.capabilities() # 检查当前电脑所支持的格式
library(qgraph)
library(igraph)
library(pals)
library(psych)
library(RColorBrewer)
library(reshape2)
library(scales)
library(tidyverse)
library(xlsx)

setwd = '/public/workspace/wanggenjun/work/scRNA/estrogen/cao/Breast'

orange9 <- brewer.pal(9,"YlOrRd")
show_col(orange9)

##Normal
# counts_Nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/Liver_EH_counts.txt", check.names = FALSE)

pvalues.nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_Normal/pvalues.txt", check.names = FALSE)
means.nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_Normal/means.txt", check.names = FALSE)
sig.means.nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_Normal/significant_means.txt", check.names = FALSE)
deconvoluted.nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_Normal/deconvoluted.txt", check.names = FALSE)

##EH
pvalues.EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_EH/pvalues.txt", check.names = FALSE)
means.EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_EH/means.txt", check.names = FALSE)
sig.means.EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_EH/significant_means.txt", check.names = FALSE)
deconvoluted.EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_EH/deconvoluted.txt", check.names = FALSE)

colnames(means.EH)

levels(Breast@meta.data$final)

celltype <- c('Epithelial cells','Endothelial cells','Fibroblasts','Smooth muscle cells','Naive T cells','Proliferating T cells',
              'CD4+ T cells','CD8+ T cells','NK cells','B cells','Plasma cells','Neutrophils','pDCs','Mac/Mono')
head(pvalues.nor)

# Normal
op <- par(mfrow=c(1,1))
for(db in c("Nor")){
  sig_means <- read.delim(sprintf("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_Normal/significant_means.txt",db), check.names = F)
  sum(!is.na(sig_means[,-c(1:12)]))
  intersum <- sort(apply(sig_means[,-c(1:12)], 2, function(x){
    sum(!is.na(x))
  }))
  matx <- matrix(0, nrow=length(Cluster.cols_NoCiliated), ncol=length(Cluster.cols_NoCiliated))
  rownames(matx) <- colnames(matx) <- names(Cluster.cols_NoCiliated)
  intersum <- intersum / sum(intersum)
  lapply(names(intersum), function(i){
    x <- unlist(strsplit(i, split = "\\|"))
    matx[x[1],x[2]] <<- intersum[i]
  })
  ## matx[upper.tri(matx)] <- -matx[lower.tri(matx)]
  corrplot(matx, is.corr = FALSE, 
           method = "square", type = "full", 
           tl.col = "black", order = "original", 
           col = c('#FEE090','#F46D43','#A50026'))
}
par(op)


op <- par(mfrow=c(1,1))
for(db in c("EH")){
  sig_means <- read.delim(sprintf("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_EH/significant_means.txt",db), check.names = F)
  sum(!is.na(sig_means[,-c(1:12)]))
  intersum <- sort(apply(sig_means[,-c(1:12)], 2, function(x){
    sum(!is.na(x))
  }))
  matx <- matrix(0, nrow=length(Cluster.cols_NoCiliated), ncol=length(Cluster.cols_NoCiliated))
  rownames(matx) <- colnames(matx) <- names(Cluster.cols_NoCiliated)
  intersum <- intersum / sum(intersum)
  lapply(names(intersum), function(i){
    x <- unlist(strsplit(i, split = "\\|"))
    matx[x[1],x[2]] <<- intersum[i]
  })
  ## matx[upper.tri(matx)] <- -matx[lower.tri(matx)]
  corrplot(matx, is.corr = FALSE, 
           method = "square", type = "full", 
           tl.col = "black", order = "original", 
           col = c('#FEE090','#F46D43','#A50026'))
}
par(op)


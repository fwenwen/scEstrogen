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
library(cowplot)
library(tidyverse)

## GSEGO
library(clusterProfiler)
library(enrichplot)
library(foreach)
library(ggstatsplot)
library(org.Mfascicularis.eg.db)  #library(org.Hs.eg.db)
library(tidyverse)

meta <- read.delim('/public/workspace/wanggenjun/work/scRNA/estrogen/cao/meta.txt',head=F,stringsAsFactors=F)
Colon <- readRDS("/public/workspace/wanggenjun/work/scRNA/estrogen/data/Organ/Colon.rds")
groupCols <- readRDS("/public/workspace/wanggenjun/work/scRNA/estrogen/data/color.list.rds")
minor2.cols <- readRDS("/public/workspace/wanggenjun/work/scRNA/estrogen/data/minor_id.col.rds")
# final marker
final.marker <- read.csv("/public/workspace/wanggenjun/work/scRNA/estrogen/data/final_marker.csv")

Colon <- RunUMAP(Colon, dims=1:30, min.dist=0.5, 
                  n.neighbors=50, verbose=FALSE, umap.method = "umap-learn",metric = "correlation")

# subset fib
Colon.fib <- subset(Colon, idents=c('Fib_0','Fib_1'))
Colon.fib <- RunUMAP(Colon.fib, dims=1:30, min.dist=0.5, n.neighbors=50, verbose=FALSE, umap.method = "umap-learn",metric = "correlation")
Colon.fib@meta.data <- droplevels(Colon.fib@meta.data)

# subset
nColon <- subset(Colon, idents=c('Epi_0','Epi_1','Epi_2','Epi_3','Epi_4','Epi_5','Epi_6','Epi_7','Epi_8','Epi_10',
                                 'Endo_0','Endo_1','Endo_2','Endo_3',
                                 'Fib_0','Fib_1',
                                 'SMC_0','SMC_1','SMC_2',
                                 'Tn','ProliferatingT',
                                 'CD4T_0','CD4T_1','CD4T_2','CD8T_0','CD8T_1','CD8T_2','CD8T_3','CD8T_4',
                                 'NK_0','NK_1','NK_2','NK_3','NK_4',
                                 'B_0','B_1','PlasmaB_0','PlasmaB_1',
                                 'Neu_0','Neu_2','pDC',
                                 'Mye_0','Mye_1','Mye_2','Mye_3','Mye_4','Mye_5','Mye_6'))
nColon <- RunUMAP(nColon, dims=1:30, min.dist=0.5, n.neighbors=50, verbose=FALSE, umap.method = "umap-learn",metric = "correlation")
nColon@meta.data <- droplevels(nColon@meta.data)

saveRDS(object = Colon.fib, file = "Colon.fib.rds") 

groupCols$minor2.cols <- minor2.cols

final.cols <- c('#B6992D','#6B468F','#9E67A0','#5D4EA2','#DDA0A7','#BB4922','#C3784D','#9D0142',
                '#E15759','#4E79A7','#4ea7b0','#245728','#377F32','#8FBA2A')

minor.cols <- c('#FED439','#709AE1','#8A9197','#D2AF81','#FD7446',
                '#D5E4A2','#197EC0','#F05C3B','#46732E','#1A9993', # Epi
  '#2aa198','#b58900','#6c71c4','#d33682', #Endo
  '#5FB233FF','#6A7F93FF','#F57206FF',#Fib
  '#EB0F13FF','#8F2F8BFF','#1396DBFF', #SMC
  '#BA6222FF',# Tn
  '#FB82BEFF',# ProliferatingT
  '#026CCBFF','#F51E02FF','#05B102FF',# CD4T
  '#2F86FFFF','#EBAB16FF','#22C408FF','#FECDAAFF','#F14809FF', #CD8T
  '#0C5BB0FF','#15983DFF','#EC579AFF','#FA6B09FF','#FEC10BFF',# NK
  '#3D79F3FF','#E6352FFF', # B
  '#EB5291FF','#1794CEFF', # PlasmaB
  '#1d457f','#c36377',# Neu
  '#f2af4a',# pDC
  '#FBCF35FF','#ED4C1CFF','#9C7E70FF','#5AC2F1FF','#11776CFF','#2366C0FF','#A3DA4BFF' # Mye
               )

# final
Idents(Colon) <- 'final'

unique(Colon$final)
options(repr.plot.width = 8, repr.plot.height = 6, repr.plot.res = 100)
umap.final <- DimPlot(Colon,label = F,cols = final.cols, 
                       group.by = 'final', pt.size = 0.5, seed = 1)#+ NoLegend() + NoAxes()
umap.final

## save as pdf
pdf("umap_final.pdf", height = 6,width = 6)
print(umap.final + NoLegend() + NoAxes() + ggtitle(NULL, subtitle = NULL))
dev.off()

# fib
options(repr.plot.width = 6, repr.plot.height = 6, repr.plot.res = 100)
cols.fib <- c('#5FB233FF','#6A7F93FF')
umap.fib <- DimPlot(Colon.fib,label = F,cols = cols.fib, 
                        group.by = 'minor2', pt.size = 3) + ggtitle("fib")
umap.fib
## save as pdf
pdf("umap_fib.pdf", height = 6,width = 6)
print(umap.fib + NoLegend() + NoAxes() + ggtitle(NULL, subtitle = NULL))
dev.off()

samplegroup <- data.frame(sample=meta$V1, group=gsub("^.*-","",meta$V1))
samplegroup$group <- str_sub(samplegroup$group,1,6)
samplegroup$group[samplegroup$group == "EH"] <- "Estrogen"
rownames(samplegroup) <- samplegroup$sample
head(samplegroup)

cell_eb <- data.frame(Colon@reductions[["umap"]]@cell.embeddings[,1:2], 
                      cluster = as.character(Colon$final),
                      group = as.character(Colon$group), 
                      stringsAsFactors = F)
colnames(cell_eb) <- c('x','y','cluster', 'group')

# add pie
df <- Colon@meta.data %>% 
    dplyr::count(final, group) %>% 
    dplyr::rename(cluster = final) %>% 
    dplyr::rename(group = group)
df$n <- df$n/as.numeric(table(Colon@meta.data$group)[df$group])
coords <- cell_eb %>% 
            group_by(cluster) %>% 
            dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% 
            as.data.frame

piedf <- reshape2::acast(df, cluster ~ group, value.var = "n")
piedf[is.na(piedf)] <- 0
piedf <- cbind(piedf, coords[match(rownames(piedf), coords$cluster), -1])

cellnumbers <- dplyr::count(cell_eb, cluster)
sizee <- cellnumbers[match(rownames(piedf), cellnumbers[,1]), 2]
piedf$radius <- scales::rescale(sizee, c(0.6,0.8))

pie <- umap.final + ggnewscale::new_scale_fill() + 
    geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=piedf,
                  cols=colnames(piedf)[1:(ncol(piedf)-3)],color=NA) +
    scale_fill_manual(values = groupCols$group.cols) + 
    coord_equal() + theme_classic() + NoAxes()+ NoLegend()
pie

# save as pdf
pdf("umap+pie.pdf", height = 6,width = 6)
print(pie + NoAxes() + NoLegend())
dev.off()

cell_eb.fib <- data.frame(Colon.fib@reductions[["umap"]]@cell.embeddings[,1:2], 
                      cluster = as.character(Colon.fib$minor2),
                      group = as.character(Colon.fib$group), 
                      stringsAsFactors = F)
colnames(cell_eb.fib) <- c('x','y','cluster', 'group')

# add pie
df <- Colon.fib@meta.data %>% 
    dplyr::count(minor2, group) %>% 
    dplyr::rename(cluster = minor2) %>% 
    dplyr::rename(group = group)
df$n <- df$n/as.numeric(table(Colon.fib@meta.data$group)[df$group])
coords <- cell_eb.fib %>% 
            group_by(cluster) %>% 
            dplyr::summarise(x=median(x), y=median(y), .groups = "keep") %>% 
            as.data.frame
options(repr.plot.width = 6, repr.plot.height = 6, repr.plot.res = 100)
piedf <- reshape2::acast(df, cluster ~ group, value.var = "n")
piedf[is.na(piedf)] <- 0
piedf <- cbind(piedf, coords[match(rownames(piedf), coords$cluster), -1])

cellnumbers <- dplyr::count(cell_eb.fib, cluster)
sizee <- cellnumbers[match(rownames(piedf), cellnumbers[,1]), 2]
piedf$radius <- scales::rescale(sizee, c(0.6,0.8))

pie.fib <- umap.fib + ggnewscale::new_scale_fill() + 
    geom_scatterpie(aes_(x=~x,y=~y,r=~radius), data=piedf,
                  cols=colnames(piedf)[1:(ncol(piedf)-3)],color=NA) +
    scale_fill_manual(values = groupCols$group.cols) + 
    coord_equal() + theme_pubr() + NoAxes()+ NoLegend()
pie.fib

# save as pdf
pdf("umap+pie_fib.pdf", height = 6,width = 6)
print(pie.fib + NoAxes() + NoLegend())
dev.off()

Colon@assays$SCT@data[1:3,1:3]
DefaultAssay(Colon) <- "SCT"

# Ensembl to Symbol -----------
refGenes <- Colon@misc$geneName
if(!any(is.na(refGenes[rownames(Colon@assays$SCT@data)]))){
  rownames(Colon@assays$SCT@data) <- refGenes[rownames(Colon@assays$SCT@data)]
}

head(rownames(Colon@assays$SCT@data))

Idents(Colon.fib) <- 'minor2'

# fib
options(repr.plot.width = 6, repr.plot.height = 6, repr.plot.res = 100)
feat.fib.FIBIN <- FeaturePlot(Colon.fib, reduction = 'umap', order=TRUE,
            cols=c("lightgrey", brewer.orrd(10)),slot = "data",
            min.cutoff='q1', pt.size=3, combine=T,
            features=c("FIBIN")) + NoLegend() + NoAxes()
feat.fib.ADH1B <- FeaturePlot(Colon.fib, reduction = 'umap', order=TRUE,
            cols=c("lightgrey", brewer.orrd(10)),slot = "data",
            min.cutoff='q1', pt.size=3, combine=T,
            features=c("ADH1B")) + NoLegend() + NoAxes()

feat.fib.ADH1B <- FeaturePlot(Colon.fib, reduction = 'umap', order=TRUE,
            cols=c("lightgrey", brewer.orrd(10)),slot = "data",
            min.cutoff='q1', pt.size=3, combine=T,
            features=c("TNC")) + NoLegend() + NoAxes()

# Ensembl to Symbol -----------
refGenes <- Colon@misc$geneName
if(!any(is.na(refGenes[rownames(Colon@assays$SCT@data)]))){
  rownames(Colon@assays$SCT@data) <- refGenes[rownames(Colon@assays$SCT@data)]
}

head(rownames(Colon@assays$SCT@data))

DefaultAssay(Colon) <- "SCT"
# final
Idents(Colon) <- Colon$final

genes_to_check = c('WFDC2','EPCAM','MMP7',# Epithelial cells 
                   #'DNAH7','C9orf116','C20orf85',# Ciliated epithelial cells
                   'PLVAP','PECAM1','VWF',# Endothelial cells 
                   'DCN','COL1A1','LUM', # Fibroblasts 
                   'MYH11','MYL9','TAGLN',  # Smooth muscle cells
                   'CCR7' ,#Naive T cells
                   'MKI67','TOP2A','CD3D',  #Proliferating T cells
                   'CD4',# CD4+ T cells
                   'CD3E','CD8A',# CD8+ T cells  #CD3D
                   'NKG7', # NK cells
                   'BANK1','MS4A1','CD79B', # B cells
                   'JCHAIN',#Plasma cells
                   'CSF3R','S100A8','S100A9', ## Neutrophils
                   'POLD1','SCT','IRF7', # pDCs
                   'C1QA','PSAP','C1QC','CD68') #Mac/Mono
options(repr.plot.width = 12, repr.plot.height = 4, repr.plot.res = 100)
p.final.marker <- DotPlot(Colon, features = genes_to_check) + 
    scale_color_gradientn(colours=rev(brewer.rdylbu(20)), guide = "colourbar") + 
    theme(axis.text.x=element_text(angle=90, hjust=1),panel.grid.major = element_line(colour = "grey70", size = 0.2))
p.final.marker

ggsave("marker_final.pdf", p.final.marker, limitsize = FALSE, width=11 ,height=4)

a <- rep(c('Fib_0','Fib_1'), times=c(670,691))
b <- rep(c('Colon-EH','Colon-Normal','Colon-EH','Colon-Normal'), times=c(289,381,595,96))

fib.pie <- data.frame(minor2=a,orig.ident=b)
head(fib.pie)
dim(fib.pie)

library(plotrix)

unique(plot$final)

labels<-c('Epithelial cells','Endothelial cells',
             'Fibroblasts','Smooth muscle cells',
             'Naive T cells','Proliferating T cells','CD4+ T cells','CD8+ T cells',
             'NK cells','B cells','Plasma cells','Neutrophils','pDCs','Mac/Mono')

xy.pop <- c(8.61,0.52,5.30,0.24,8.24,8.00,16.86,13.17,3.04,30.91,4.03,0.02,0.03,1.02)
xx.pop <- c(20.86,5.75,16.04,2.43,1.78,7.69,6.66,13.75,4.14,6.24,13.06,0.29,0.00,1.31)
xx.pop

head(final.cols)

options(repr.plot.width = 10, repr.plot.height = 8, repr.plot.res = 100)
pdf("pyramid_final.pdf", height=8, width=10)
par(mar=pyramid.plot(xy.pop,xx.pop,labels=F,                    
                     lxcol=final.cols,rxcol=final.cols,
                                      laxlab=c(0,10,20,30,40),
                                      raxlab=c(0,10,20,30,40),
                     gap=2,space=0.1,show.values=F))
dev.off()

Colon.fib$samplefinal <- factor(sprintf("%s-%s", Colon.fib$minor2, Colon.fib$group), 
                                levels=unlist(lapply(levels(Colon.fib$minor2), 
                                                     function(x){sprintf("%s-%s",x, levels(Colon.fib$group))})))

Colon.fib.01$samplefinal <- factor(sprintf("%s-%s", Colon.fib.01$minor2, Colon.fib.01$group), 
                                levels=unlist(lapply(levels(Colon.fib.01$minor2), 
                                                     function(x){sprintf("%s-%s",x, levels(Colon.fib.01$group))})))
head(Colon.fib.01$samplefinal)
levels(Colon.fib.01$samplefinal)

Idents(Colon.fib) <- Colon.fib$minor2

unique(Colon.fib$minor2)

# top markers for fib
# Symbol to Ensembl ------------
refGenes <- Colon.fib@misc$geneName
rownames(Colon.fib@assays$SCT@data) <- names(refGenes)[match(rownames(Colon.fib@assays$SCT@data),refGenes)]

top.markers.fib <- FindAllMarkers(Colon.fib, assay="SCT", 
                              min.pct=0.1, logfc.threshold=0.1, return.thresh=1, only.pos=T) 
top.markers.fib$pct.diff <- top.markers.fib$pct.1 - top.markers.fib$pct.2

write.table(top.markers.fib, "./Fib.top.markers.txt", row.names = F, quote = F, sep = "\t")

library(foreach)
library(enrichplot)
library(clusterProfiler)

minpct <- .1
gseOut.fib <- foreach(cls=levels(top.markers.fib$cluster)) %do% {
  submarkers <- top.markers.fib %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
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
  ego@result <- cbind(ego@result, cluster=cls)
  ego
}
gseResult.fib <- do.call("rbind", lapply(gseOut.fib, function(x){x@result}))
head(gseResult.fib)

test <- gseResult.fib %>% filter(pvalue <= 0.05)

topPath.fib <- gseResult.fib %>% filter(pvalue <= 0.05) %>% dplyr::group_by(cluster) %>% 
  dplyr::top_n(n=10, wt=enrichmentScore)
head(topPath.fib)
table(topPath.fib$cluster)

write.table(topPath.fib, "top10Path.fib.txt", row.names = F, quote = F, sep = "\t")

topPath.fib$enrichmentScore1 <- ifelse(topPath.fib$cluster=="Fib_0",-topPath.fib$enrichmentScore,topPath.fib$enrichmentScore)


library(ggplot2)
library(forcats)
options(repr.plot.width = 9, repr.plot.height = 7, repr.plot.res = 100)

pdf("go_BP.pdf", height = 6.5,width = 9)
ggplot(topPath.fib, showCategory = 10, 
       aes(enrichmentScore1, fct_reorder(Description, enrichmentScore1))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=-log10(pvalue), size = enrichmentScore)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=F)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  ylab(NULL) 

dev.off()

library(SCORPIUS)

expression <- t(as.matrix(Colon.fib@assays$SCT@scale.data))

group_name =  as.factor(as.character(Colon.fib$minor2))

table(group_name)
dim(expression)

expression[1:4,1:4]

space <- reduce_dimensionality(expression, "spearman")
draw_trajectory_plot(space, progression_group=group_name, contour = T,
                     progression_group_palette = setNames(c('#5FB233FF','#6A7F93FF'),
                                      levels(group_name)))

space <- reduce_dimensionality(expression, "spearman")
sp1 <- draw_trajectory_plot(space, group_name, contour = TRUE,progression_group_palette = setNames(c('#5FB233FF','#6A7F93FF'),
                                      levels(group_name)))
traj <- infer_trajectory(space)
sp2 <- draw_trajectory_plot(space, group_name, traj$path, contour = TRUE,path_size = 1,
                     progression_group_palette = setNames(c('#5FB233FF','#6A7F93FF'),
                                      levels(group_name)))

gimp <- gene_importances(
  expression, 
  traj$time, 
  num_permutations = 10, 
  num_threads = 8, 
  ntree = 10000,
  ntree_perm = 1000
) 
gimp$qvalue <- p.adjust(gimp$pvalue, "BH", length(gimp$pvalue))

gene_sel_n <- gimp$gene[gimp$qvalue < .05 & gimp$importance>0.9]

expr_sel <- scale_quantile(expression[,gene_sel_n])

# Ensembl to Symbol -----------
refGenes <- Colon.fib@misc$geneName
if(!any(is.na(refGenes[colnames(expr_sel)]))){
  colnames(expr_sel) <- refGenes[colnames(expr_sel)]
}

head(colnames(expr_sel))

expr_sel_n <- expr_sel[, -grep("^ENSMFAG", colnames(expr_sel))]

# Draw a time series heatmap
time <- traj$time
draw_trajectory_heatmap(expr_sel, time)

## Also show the progression groupings
draw_trajectory_heatmap(expr_sel, time, 
                        progression_group=group_name)

hmcols = colorRampPalette(rev(c("#CF384D","#ED6345","#FA9A58","#F3FAAD","#D1EC9C","#96D5A4","#5BB6A9","#3682BA")))(62)

sp4 <- draw_trajectory_heatmap(
  expr_sel_n, time, progression_group=group_name,show_labels_row =T,color = hmcols,
  progression_group_palette = setNames(c('#5FB233FF','#6A7F93FF'),
                                      levels(group_name))
)

library(cpplot)
library(Cairo)
Cairo.capabilities() 
library(qgraph)
library(igraph)
library(pals)
library(psych)
library(RColorBrewer)
library(reshape2)
library(scales)
library(tidyverse)
library(xlsx)
setwd = '/public/workspace/wanggenjun/work/scRNA/estrogen/cao/Colon'

orange9 <- brewer.pal(9,"YlOrRd")
show_col(orange9)

pvalues.nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_Normal/pvalues.txt", check.names = FALSE)
means.nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_Normal/means.txt", check.names = FALSE)
sig.means.nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_Normal/significant_means.txt", check.names = FALSE)
deconvoluted.nor <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_Normal/deconvoluted.txt", check.names = FALSE)

pvalues.EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_EH/pvalues.txt", check.names = FALSE)
means.EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_EH/means.txt", check.names = FALSE)
sig.means.EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_EH/significant_means.txt", check.names = FALSE)
deconvoluted.EH <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_EH/deconvoluted.txt", check.names = FALSE)

celltype <- c('Epithelial cells','Endothelial cells','Fibroblasts','Smooth muscle cells','Naive T cells','Proliferating T cells',
              'CD4+ T cells','CD8+ T cells','NK cells','B cells','Plasma cells','Neutrophils','pDCs','Mac/Mono')

ccc_compare2(group1.name = "Normal",group2.name = "Estrogen",
             group1.pfile = "/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_Normal/pvalues.txt",
             group1.mfile="/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_Normal/means.txt",
             group2.pfile="/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_EH/pvalues.txt",
             group2.mfile="/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Breast_EH/means.txt",
             p.threshold = 0.05,thre=0.5,
             cell.pair="Neutrophils|B cells", #指定ligand产生的细胞|receptor产生的细胞
             plot.width=15,plot.height=30,filename = "0914_ccinter"
)

library(circlize)
library(corrplot)


Cluster.cols_Ciliated <- Cluster.cols_NoCiliated[-13]

options(repr.plot.width = 9, repr.plot.height = 9, repr.plot.res = 100)
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


# EH
op <- par(mfrow=c(1,1))
for(db in c("EH")){
  sig_means <- read.delim(sprintf("/public/workspace/wanggenjun/work/scRNA/estrogen/data/cellphonedb/out/Colon_EH/significant_means.txt",db), check.names = F)
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


refGenes <- Colon.B@misc$geneName
if(!any(is.na(refGenes[rownames(Colon.B@assays$SCT@data)]))){
  rownames(Colon.B@assays$SCT@data) <- refGenes[rownames(Colon.B@assays$SCT@data)]
}
head(rownames(Colon.B@assays$SCT@data))

options(repr.plot.width = 15, repr.plot.height = 5, repr.plot.res = 100)
vln.fib <- VlnPlot(Colon.fib, features = unique(top5.fib$gene), split.by = "group", group.by = 'minor2',
                 cols=c("#4E889B","#EFE118"),pt.size = 0, combine = True, stack = TRUE)
vln.fib

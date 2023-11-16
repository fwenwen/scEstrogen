# Figure1 
options(stringsAsFactors = F)
library(cowplot)
library(ggpubr)
library(dplyr)
library(readr)
library(tidyr)
library(ggforce)
library(pals)
library(pheatmap)
library(scales)
library(ggthemes)
library(Seurat)
# meta
p1 <- DimPlot(seuratObj, group.by = "final",reduction = "umap", 
              cols = Cluster.cols, pt.size = 0.2,order = T)

p2 <- DimPlot(seuratObj, group.by = "organ", cols = organ.cols,#order = T, 
              pt.size = 0.2,reduction = "umap")
p3 <- DimPlot(seuratObj, group.by = "group", cols = group.cols,#order = T, 
              pt.size = 0.2,reduction = "umap")

p4 <- DimPlot(seuratObj, group.by = "organ.group", cols = sample.cols,#order = T, 
              pt.size = 0.2,reduction = "umap")

ggarrange(p1+NoAxes()+NoLegend(),p2+NoAxes()+NoLegend(),
          p3+NoAxes()+NoLegend(),p4+NoAxes()+NoLegend(),
          ncol = 2,nrow = 2)

library(ggpubr)
ggarrange(as_ggplot(get_legend(p1)),as_ggplot(get_legend(p2)),
          as_ggplot(get_legend(p3)),as_ggplot(get_legend(p4)),
          ncol = 2,nrow = 2)

pdf("meta.pdf", height=8.27, width=8.27)
print(ggarrange(p1+NoAxes()+NoLegend(),p2+NoAxes()+NoLegend(),
                p3+NoAxes()+NoLegend(),p4+NoAxes()+NoLegend(),
                ncol = 2,nrow = 2))
print(ggarrange(as_ggplot(get_legend(p1)),as_ggplot(get_legend(p2)),
                as_ggplot(get_legend(p3)),as_ggplot(get_legend(p4)),
                ncol = 2,nrow = 2))
dev.off()

pdf("meta_legend.pdf", height=8.27, width=8.27)
print(ggarrange(as_ggplot(get_legend(p1)),as_ggplot(get_legend(p2)),
                as_ggplot(get_legend(p3)),as_ggplot(get_legend(p4)),
                ncol = 2,nrow = 2))
dev.off()

# percentage
cellstat <- table(seuratObj$group, seuratObj$final)
cellstat <- (cellstat * 100 / rowSums(cellstat))
cellstat <- reshape2::melt(t(cellstat) * 100 / colSums(cellstat))
cellstat <- as.data.frame(cellstat)
head(cellstat)
colnames(cellstat) <- c("cellType", "group", "percentage")

unique(seuratObj$organ)
seuratObj.list <- SplitObject(seuratObj,split.by = "organ")

cellstat_split <- data.frame()
for(i in unique(seuratObj$organ)){
  cellstat <- table(seuratObj.list[[i]]$group, seuratObj.list[[i]]$final)
  cellstat <- (cellstat * 100 / rowSums(cellstat))
  cellstat <- reshape2::melt(t(cellstat) * 100 / colSums(cellstat))
  cellstat <- as.data.frame(cellstat)
  colnames(cellstat) <- c("cellType", "group", "percentage")
  cellstat$organ <- i
  cellstat_split <- rbind(cellstat_split,cellstat)
}
organ.number <- data.frame(organ=names(organ.cols),num=1:5)
rownames(organ.number) <- organ.number$organ
cellstat_split$y <- organ.number[as.character(cellstat_split$organ),"num"]
cellstat_split$y <- as.numeric(cellstat_split$y)

final.number <- data.frame(final=names(Cluster.cols),num=1:15)
rownames(final.number) <- final.number$final
cellstat_split$x <- final.number[as.character(cellstat_split$cellType),"num"]
cellstat_split$x <- as.numeric(cellstat_split$x)

cellstat_split_Normal <- cellstat_split[cellstat_split$group=="Normal",] 
cellstat_split_Estrogen <- cellstat_split[cellstat_split$group=="Estrogen",]

head(cellstat_split_Normal)
cellstat1 <- merge(cellstat_split_Normal,cellstat_split_Estrogen,by=c("cellType","organ")) %>% 
  dplyr::select(.,c(1,2,4,8:10))
head(cellstat1)
colnames(cellstat1) <- c("cellType", "organ","Normal", "Estrogen", "x","y"     )
library(scatterpie)
p5 <- ggplot() +
  geom_scatterpie(aes(x,y,r=0.4), cellstat1,
                  cols = c("Normal","Estrogen")) 


p6 <- p5+coord_equal()+
  theme_void()+
  ggtitle("Percentage of cells in each group")+
  theme(legend.position = "",
        plot.title = element_text(vjust =0,hjust =0.5,size = 14))+
  scale_fill_manual(values = group.cols)

pdf("pie.pdf",height = 5.27,width = 10.27)
print(p6)
dev.off()


# *organ_marker
organ_marker <- data.frame()
for(i in unique(seuratObj$organ)){
  marker <- read.delim(sprintf("./%s_marker1.txt",i))
  organ_marker <- rbind(organ_marker,marker)
}
# 
head(organ_marker)
dim(organ_marker)
organ_marker <- organ_marker %>% dplyr::filter(p_val<0.05 & avg_log2FC > 0.25 & pct.1 > pct.2)
dim(organ_marker)

head(organ_marker)
colnames(organ_marker)[6] <- "group"
head(organ_marker %>% dplyr::count(organ,celltype,group))
organ_marker_number <- organ_marker %>% dplyr::count(organ,celltype,group)
head(organ_marker_number)

organ_marker_number$group[organ_marker_number$group=="EH"] <- "Estrogen"
organ_marker_number$celltype <- factor(organ_marker_number$celltype,levels = rev(names(Cluster.cols)))
organ_marker_number$organ <- factor(organ_marker_number$organ,levels = names(organ.cols))
# plot 
p7 <- organ_marker_number %>% dplyr::arrange(celltype) %>% 
  ggplot(aes(y=celltype,x = ifelse(group=="Normal",-n,n),fill = group)) + 
  geom_bar(stat = "identity") +
  scale_fill_manual(values = group.cols)+
  theme_classic()+
  theme(panel.spacing.x = unit(0.1, "cm"),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        axis.title = element_blank(),
        
        strip.text.x = element_text(size=10,color="black",face="bold",vjust=0.5),
        strip.background.x = element_blank(),
        axis.text = element_text(color="black"),
        plot.margin=unit(c(0.2,0.5,0.2,0.2),units=,"cm"),
        legend.text=element_text(color="black",size=9),
        legend.title = element_blank(),
        legend.key=element_blank(),  
        legend.background=element_blank(),
        legend.box.margin = margin(1,1,1,1),
        legend.position = c(0.16,0.4), legend.justification = c(0.1,0.6))


pdf("bar_gene.pdf", height=2.57, width=8.27)
print(p7)
dev.off()

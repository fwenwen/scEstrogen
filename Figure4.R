# Figure4
library(dplyr)
library(Seurat)
library(ggforce)
DefaultAssay(seuratObj) <- 'SCT'

unique(seuratObj$organ)
Idents(seuratObj) <- seuratObj$organ
Breast <- subset(seuratObj,ident="Breast")
Breast@meta.data <- droplevels(Breast@meta.data)
unique(Breast$final)
Idents(Breast) <- Breast$group
unique(Breast$group)

top.markers <- FindAllMarkers(Breast, assay="SCT", min.pct=0.1, logfc.threshold=0.1, return.thresh=1, only.pos=TRUE) 

top.markers$pct.diff <- top.markers$pct.1 - top.markers$pct.2
write.table(top.markers, "./group.markers.txt", row.names = F, quote = F, sep = "\t")
head(top.markers)
sig.markers <- top.markers[!duplicated(top.markers$gene),] %>% dplyr::filter(p_val_adj < 0.01)
head(sig.markers)

BC <- read.delim("./Breast/MONDO_0007254_associations_export.tsv")
head(BC)
unique(BC$mappedGenes)
BC <- BC %>% dplyr::mutate(as.data.frame(str_split_fixed(BC$mappedGenes,",",2))) %>% 
  dplyr::select(,c("pValue","mappedGenes","V1","V2")) %>% 
  magrittr::set_colnames(c("pValue","mappedGenes","gene1","gene2"))
head(BC)

BC_dup <- BC %>% dplyr::group_by(gene) %>% 
  dplyr::summarise(pValue = pValue[which.max(pValue)])
head(BC_dup)
BC_dup <- BC_dup[-1,]
saveRDS(BC_dup,"./Breast/BC_dup_max.rds")
BC_dup <- as.data.frame(BC_dup)
rownames(BC_dup) <- BC_dup$gene

top.markers[top.markers$gene %in% intersect(top.markers$gene,BC_dup$gene),]

BC_vol <- top.markers[top.markers$gene %in% intersect(top.markers$gene,BC_dup$gene),]

BC_vol <- BC_vol[,c("avg_log2FC","cluster","gene")] %>% 
  dplyr::mutate(pval=-log10(BC_dup[as.character(BC_vol$gene),"pValue"]))
head(BC_vol)
BC_vol$avg_log2FC <- ifelse(BC_vol$cluster=="Normal",-BC_vol$avg_log2FC,BC_vol$avg_log2FC)
# 
genes <- BC_vol$gene[abs(BC_vol$avg_log2FC)>log2(1.1)]
pltd <- BC_vol %>%
  mutate(label=ifelse(gene%in%genes,gene,""),
         category=ifelse(abs(avg_log2FC)< log2(1.1), "NS", ifelse(avg_log2FC>0,"Estrogen","Normal"))) 

# 
c2 <- pltd %>%
  ggplot(aes(x=pval,y=avg_log2FC)) + 
  geom_point(aes(fill=category,color=category),
             size=ifelse(pltd$label=="",1,3),shape=21) +
  geom_text(aes(label=label), hjust=0, vjust=-0.5, size=3) +
  ## ggrepel::geom_text_repel(aes(label=label)) +
  scale_fill_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
  scale_color_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
  # theme_minimal()
  theme_pubr()
c2
pdf("./Breast/group.vol.pdf",height = 4.27,width = 4.27)
print(c2)
dev.off()
# 

Idents(Breast) <- Breast$final
for(i in levels(Breast$final)[-13]){
  subobj <- subset(Breast,idents=i)
  Idents(subobj) <- subobj$group
  markers <- FindAllMarkers(subobj, assay="SCT", min.pct=0.1, logfc.threshold=0, return.thresh=1, only.pos=TRUE) 
  markers$pct.diff <- markers$pct.1 - markers$pct.2
  markers$final <- i
  write.table(markers, sprintf("./Breast/%s.group.markers.txt",i), row.names = F, quote = F, sep = "\t")
}


GWAS_vol <- function(BC,organ,final){
  BC <- read.delim("./Breast/MONDO_0007254_associations_export.tsv")
  BC <- BC[-grep(";",BC$mappedGenes),]
  BC <- BC %>% dplyr::mutate(as.data.frame(str_split_fixed(BC$mappedGenes,",",2))) %>% 
    dplyr::select(,c("pValue","mappedGenes","V1","V2")) %>% 
    magrittr::set_colnames(c("pValue","mappedGenes","gene1","gene2"))
  # 
  head(BC)
  BC_dup <- BC %>% dplyr::group_by(gene) %>% 
    dplyr::summarise(pValue = pValue[which.max(pValue)])
  head(BC_dup)
  BC_dup <- BC_dup[-1,]
  BC_dup <- as.data.frame(BC_dup)
  rownames(BC_dup) <- BC_dup$gene
  saveRDS(BC_dup,sprintf("./%s/BC_dup_max.rds",organ))
  
  top.markers <- read.delim(sprintf("./%s/%s.group.markers.txt",organ,final))
  BC_vol <- top.markers[top.markers$gene %in% intersect(top.markers$gene,BC_dup$gene),]
  BC_vol <- BC_vol[,c("avg_log2FC","cluster","gene")] %>% 
    dplyr::mutate(pval=-log10(BC_dup[as.character(BC_vol$gene),"pValue"]))
  BC_vol$avg_log2FC <- ifelse(BC_vol$cluster=="Normal",-BC_vol$avg_log2FC,BC_vol$avg_log2FC)
  # 
  genes <- BC_vol$gene[abs(BC_vol$avg_log2FC)>log2(1.1)]
  pltd <- BC_vol %>%
    mutate(label=ifelse(gene%in%genes,gene,""),
           category=ifelse(abs(avg_log2FC)< log2(1.1), "NS", ifelse(avg_log2FC>0,"Estrogen","Normal"))) 
  
  # 
  c2 <- pltd %>%
    ggplot(aes(x=pval,y=avg_log2FC)) + 
    geom_point(aes(fill=category,color=category),
               size=ifelse(pltd$label=="",1,3),shape=21) +
    geom_text(aes(label=label), hjust=0, vjust=-0.5, size=3) +
    ## ggrepel::geom_text_repel(aes(label=label)) +
    scale_fill_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
    scale_color_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
    # theme_minimal()
    theme_pubr()+ggtitle(final)
  pdf(sprintf("./%s/fig/%s.vol.pdf",organ,final),
      height = 4.27,width = 4.27)
  print(c2)
  dev.off()
  return(c2)
}

GWAS_vol(BC = BC,organ = "Breast",final = "B cells")

c.list <- list()
for(i in c(levels(Breast$final)[-(13:14)],"Mac.Mono")){
  c.list[[i]] <- GWAS_vol(BC = BC,organ = "Breast",final = i)
}
i

length(c.list)

ggarrange(plotlist = c.list, ncol = 5, nrow = 3)

pdf("./Breast/fig/final.vol.pdf",
    height = 14.27,width = 15.27)
print(ggarrange(plotlist = c.list, ncol = 5, nrow = 3))
dev.off()

unique(seuratObj$majorcluster)

DimPlot(seuratObj,group.by = "finalmajor",cols = major.cols,pt.size = 0.2)+NoAxes()
# 
unique(seuratObj$organ)
Idents(seuratObj) <- seuratObj$organ
Organ.list <- SplitObject(seuratObj,split.by = "organ")
# seuratObj@meta.data
for(i in names(Organ.list)){
  Organ.list[[i]]@meta.data <- droplevels(Organ.list[[i]]@meta.data)
  # saveRDS(Organ.list[[i]],paste0(i,".rds"))
}
# 
Breast <- Organ.list[["Breast"]]
# 
Idents(Breast) <- Breast$finalmajor
for(i in levels(Breast$finalmajor)){
  subobj <- subset(Breast,idents=i)
  Idents(subobj) <- subobj$group
  markers <- FindAllMarkers(subobj, assay="SCT", min.pct=0.1, logfc.threshold=0, return.thresh=1, only.pos=TRUE) 
  markers$pct.diff <- markers$pct.1 - markers$pct.2
  markers$finalmajor <- i
  write.table(markers, sprintf("./%s.group.markers.txt",i), row.names = F, quote = F, sep = "\t")
}
# 
GWAS_major_vol <- function(BC,finalmajor){
  BC <- read.delim("/public/workspace202011/singlecell/Monkey_estrogen/RNA/analysis/New/Figure9/Breast/breast_cancer.tsv")
  BC <- BC[-grep(";",BC$mappedGenes),]
  BC <- BC %>% dplyr::mutate(as.data.frame(str_split_fixed(BC$mappedGenes,",",2))) %>% 
    dplyr::select(,c("pValue","mappedGenes","V1","V2")) %>% 
    magrittr::set_colnames(c("pValue","mappedGenes","gene1","gene2"))
  # 
  m1 <- BC[,c("pValue","mappedGenes","gene1")] %>% magrittr::set_colnames(c("pValue","mappedGenes","gene"))
  m2 <- BC[,c("pValue","mappedGenes","gene1")] %>% magrittr::set_colnames(c("pValue","mappedGenes","gene"))
  BC <- rbind(m1,m2)
  head(BC)
  BC_dup <- BC %>% dplyr::group_by(gene) %>% 
    dplyr::summarise(pValue = pValue[which.max(pValue)])
  head(BC_dup)
  BC_dup <- BC_dup[-1,]
  BC_dup <- as.data.frame(BC_dup)
  rownames(BC_dup) <- BC_dup$gene
  # saveRDS(BC_dup,sprintf("./%s/BC_dup_max.rds",organ))
  
  top.markers <- read.delim(sprintf("./%s.group.markers.txt",finalmajor))
  BC_vol <- top.markers[top.markers$gene %in% intersect(top.markers$gene,BC_dup$gene),]
  BC_vol <- BC_vol[,c("avg_log2FC","cluster","gene")] %>% 
    dplyr::mutate(pval=-log10(BC_dup[as.character(BC_vol$gene),"pValue"]))
  BC_vol$avg_log2FC <- ifelse(BC_vol$cluster=="Normal",-BC_vol$avg_log2FC,BC_vol$avg_log2FC)
  # 
  genes <- BC_vol$gene[abs(BC_vol$avg_log2FC)>log2(1.1)]
  pltd <- BC_vol %>%
    mutate(label=ifelse(gene%in%genes,gene,""),
           category=ifelse(abs(avg_log2FC)< log2(1.1), "NS", ifelse(avg_log2FC>0,"Estrogen","Normal"))) 
  
  # 
  c2 <- pltd %>%
    ggplot(aes(x=pval,y=avg_log2FC)) + 
    geom_point(aes(fill=category,color=category),
               size=ifelse(pltd$label=="",1,3),shape=21) +
    geom_text(aes(label=label), hjust=0, vjust=-0.5, size=3) +
    ## ggrepel::geom_text_repel(aes(label=label)) +
    scale_fill_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
    scale_color_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
    # theme_minimal()
    theme_bw()+ggtitle(finalmajor)
  pdf(sprintf("./%s.vol.pdf",finalmajor),
      height = 4.27,width = 4.27)
  print(c2)
  dev.off()
  return(c2)
}

GWAS_major_vol(BC = BC,final = "B cells")

c.list <- list()
for(i in levels(Breast$finalmajor)){
  c.list[[i]] <- GWAS_major_vol(BC = BC, final = i)
}

ggarrange(plotlist = c.list, ncol = 6, nrow = 1)

pdf("./final.vol.pdf",
    height = 8.27,width = 20.27)
print(ggarrange(plotlist = c.list, ncol = 4, nrow = 2))
dev.off()

# Colorectal_Cancer
Colon <- Organ.list[["Colon"]]
# 
Idents(Colon) <- Colon$finalmajor
for(i in levels(Colon$finalmajor)){
  subobj <- subset(Colon,idents=i)
  Idents(subobj) <- subobj$group
  markers <- FindAllMarkers(subobj, assay="SCT", min.pct=0.1, logfc.threshold=0, return.thresh=1, only.pos=TRUE) 
  markers$pct.diff <- markers$pct.1 - markers$pct.2
  markers$finalmajor <- i
  write.table(markers, sprintf("./%s.group.markers.txt",i), row.names = F, quote = F, sep = "\t")
}
# 
GWAS_major_vol <- function(BC,finalmajor){
  BC <- read.delim("/public/workspace202011/singlecell/Monkey_estrogen/RNA/analysis/New/Figure9/Colon/Colorectal_Cancer.tsv")
  BC <- BC[-grep(";",BC$mappedGenes),]
  BC <- BC %>% dplyr::mutate(as.data.frame(str_split_fixed(BC$mappedGenes,",",2))) %>% 
    dplyr::select(,c("pValue","mappedGenes","V1","V2")) %>% 
    magrittr::set_colnames(c("pValue","mappedGenes","gene1","gene2"))
  # 
  head(BC)
  BC_dup <- BC %>% dplyr::group_by(gene) %>% 
    dplyr::summarise(pValue = pValue[which.max(pValue)])
  head(BC_dup)
  BC_dup <- BC_dup[-1,]
  BC_dup <- as.data.frame(BC_dup)
  rownames(BC_dup) <- BC_dup$gene
  # saveRDS(BC_dup,sprintf("./%s/BC_dup_max.rds",organ))
  
  top.markers <- read.delim(sprintf("./%s.group.markers.txt",finalmajor))
  BC_vol <- top.markers[top.markers$gene %in% intersect(top.markers$gene,BC_dup$gene),]
  BC_vol <- BC_vol[,c("avg_log2FC","cluster","gene")] %>% 
    dplyr::mutate(pval=-log10(BC_dup[as.character(BC_vol$gene),"pValue"]))
  BC_vol$avg_log2FC <- ifelse(BC_vol$cluster=="Normal",-BC_vol$avg_log2FC,BC_vol$avg_log2FC)
  # 
  genes <- BC_vol$gene[abs(BC_vol$avg_log2FC)>log2(1.1)]
  pltd <- BC_vol %>%
    mutate(label=ifelse(gene%in%genes,gene,""),
           category=ifelse(abs(avg_log2FC)< log2(1.1), "NS", ifelse(avg_log2FC>0,"Estrogen","Normal"))) 
  
  # 
  c2 <- pltd %>%
    ggplot(aes(x=pval,y=avg_log2FC)) + 
    geom_point(aes(fill=category,color=category),
               size=ifelse(pltd$label=="",1,3),shape=21) +
    geom_text(aes(label=label), hjust=0, vjust=-0.5, size=3) +
    ## ggrepel::geom_text_repel(aes(label=label)) +
    scale_fill_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
    scale_color_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
    # theme_minimal()
    theme_bw()+ggtitle(finalmajor)
  pdf(sprintf("./%s.vol.pdf",finalmajor),
      height = 4.27,width = 4.27)
  print(c2)
  dev.off()
  return(c2)
}

# GWAS_major_vol(BC = BC,final = "B cells")

c.list <- list()
for(i in levels(Colon$finalmajor)){
  c.list[[i]] <- GWAS_major_vol(BC = BC, final = i)
}

ggarrange(plotlist = c.list, ncol = 6, nrow = 1)

pdf("./final.vol.pdf",
    height = 8.27,width = 20.27)
print(ggarrange(plotlist = c.list, ncol = 4, nrow = 2))
dev.off()

unique(seuratObj$majorcluster)

# **endometrial cancer ---------------------------------
setwd("/public/workspace202011/singlecell/Monkey_estrogen/RNA/analysis/New/Figure9/Uterus/major2/")
# 
Uterus <- Organ.list[["Uterus"]]
Idents(Uterus) <- Uterus$orig.ident
unique(Uterus$orig.ident)
Uterus <- subset(Uterus,idents=c("Uterus-EH","Uterus-Normal1"))
Uterus@meta.data <- droplevels(Uterus@meta.data)
# 

# 
Idents(Uterus) <- Uterus$finalmajor
for(i in levels(Uterus$finalmajor)){
  subobj <- subset(Uterus,idents=i)
  Idents(subobj) <- subobj$group
  markers <- FindAllMarkers(subobj, assay="SCT", min.pct=0.1, logfc.threshold=0, return.thresh=1, only.pos=TRUE) 
  markers$pct.diff <- markers$pct.1 - markers$pct.2
  markers$finalmajor <- i
  write.table(markers, sprintf("./%s.group.markers.txt",i), row.names = F, quote = F, sep = "\t")
}
# 
GWAS_major_vol <- function(BC,finalmajor){
  BC <- read.delim("./endometrial_cancer.tsv")
  # BC <- BC[-grep(";",BC$mappedGenes),]
  BC <- BC %>% dplyr::mutate(as.data.frame(str_split_fixed(BC$mappedGenes,",",2))) %>% 
    dplyr::select(,c("pValue","mappedGenes","V1","V2")) %>% 
    magrittr::set_colnames(c("pValue","mappedGenes","gene1","gene2"))
  # 
  head(BC)
  BC_dup <- BC %>% dplyr::group_by(gene) %>% 
    dplyr::summarise(pValue = pValue[which.max(pValue)])
  head(BC_dup)
  # BC_dup <- BC_dup[-1,]
  BC_dup <- as.data.frame(BC_dup)
  rownames(BC_dup) <- BC_dup$gene
  # saveRDS(BC_dup,sprintf("./%s/BC_dup_max.rds",organ))
  
  top.markers <- read.delim(sprintf("./%s.group.markers.txt",finalmajor))
  BC_vol <- top.markers[top.markers$gene %in% intersect(top.markers$gene,BC_dup$gene),]
  BC_vol <- BC_vol[,c("avg_log2FC","cluster","gene")] %>% 
    dplyr::mutate(pval=-log10(BC_dup[as.character(BC_vol$gene),"pValue"]))
  BC_vol$avg_log2FC <- ifelse(BC_vol$cluster=="Normal",-BC_vol$avg_log2FC,BC_vol$avg_log2FC)
  # 
  genes <- BC_vol$gene[abs(BC_vol$avg_log2FC)>log2(1.1)]
  pltd <- BC_vol %>%
    mutate(label=ifelse(gene%in%genes,gene,""),
           category=ifelse(abs(avg_log2FC)< log2(1.1), "NS", ifelse(avg_log2FC>0,"Estrogen","Normal"))) 
  
  # 
  c2 <- pltd %>%
    ggplot(aes(x=pval,y=avg_log2FC)) + 
    geom_point(aes(fill=category,color=category),
               size=ifelse(pltd$label=="",1,3),shape=21) +
    geom_text(aes(label=label), hjust=0, vjust=-0.5, size=3) +
    ## ggrepel::geom_text_repel(aes(label=label)) +
    scale_fill_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
    scale_color_manual(values=c(NS="#161919", Estrogen="#EFE118", Normal="#4E889B")) +
    # theme_minimal()
    theme_pubr()+ggtitle(finalmajor)
  pdf(sprintf("./%s.vol_endometrial_cancer.pdf",finalmajor),
      height = 4.27,width = 4.27)
  print(c2)
  dev.off()
  return(c2)
}

# GWAS_major_vol(BC = BC,final = "B cells")

c.list <- list()
for(i in levels(Uterus$finalmajor)){
  c.list[[i]] <- GWAS_major_vol(BC = BC, final = i)
}

ggarrange(plotlist = c.list, ncol = 6, nrow = 1)

pdf("./final.vol_endometrial_cancer.pdf",
    height = 4.27,width = 20.27)
print(ggarrange(plotlist = c.list, ncol = 6, nrow = 1))
dev.off()


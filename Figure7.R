###-------------------------7A
# the same as Fig. 5A

###-------------------------7B
plot <- dplyr::count(Liver@meta.data, final, orig.ident)
plot_N <- subset(plot,orig.ident=='Liver-Normal')
plot_E <- subset(plot,orig.ident=='Liver-EH')
scalediagram.bycluster_N <- ggplot(data = plot_N, aes(fill = final, n, orig.ident))+ 
  geom_bar(stat = "identity", position="fill", width=0.5,size=0.25)+
  scale_fill_manual(values = Colors$Cluster.cols) + coord_flip()+
  coord_polar(theta = 'y')

print(scalediagram.bycluster_N)
scalediagram.bycluster_E <- ggplot(data = plot_E, aes(fill = final, n, orig.ident))+ 
  geom_bar(stat = "identity", position="fill", width=0.5,size=0.25)+
  scale_fill_manual(values = Colors$Cluster.cols) + coord_flip()+
  coord_polar(theta = 'y')
print(scalediagram.bycluster_E)

###-------------------------7C
#the same as Fig. 5B

###-------------------------7D
par(mfrow = c(4,4), xpd=TRUE)
for (i in 1:nrow(matx)) {
  mat2 <- matrix(0, nrow = nrow(matx), ncol = ncol(matx), dimnames = dimnames(matx))
  mat2[i, ] <- matx[i, ]
  netVisual_circle(mat2, 
                   weight.scale = T, 
                   edge.weight.max = max(matx), 
                   title.name = rownames(matx)[i],
                   arrow.size=0.2,col=Cluster.cols_NoCiliated)
}

###-------------------------7E
#the same as Fig. 5G

###-------------------------7F
GO <- enrichGO(OrgDb=org.Mfascicularis.eg.db,
               gene = deg_neu$gene,
               pvalueCutoff = 0.1,
               qvalueCutoff = 1,
               keyType = 'GID',
               pAdjustMethod = 'fdr',
               ont = "All",
               readable=TRUE)
head(GO)
barplot(GO, x = "GeneRatio", color = "p.adjust", showCategory =10) 

###-------------------------7G
p1_neu <- DimPlot(scRNA.neu, reduction = "umap",
                  label = T,label.box = T)+xlim(-15,15) 
p1_neu
###-------------------------7H
library(monocle)
scRNA.neu@meta.data  <- droplevels(scRNA.neu@meta.data)
data <- as(as.matrix(scRNA.neu@assays$SCT@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA.neu_bycluster@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        lowerDetectionLimit = 0.5,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
saveRDS(mycds,"cds_neu_bycluster.rds")
mycds <- readRDS("/public/workspace/wanggenjun/work/scRNA/estrogen/zhao/LiverAnalysis/cds_neu_bycluster.rds")

Idents(scRNA.neu_bycluster) <- 'minor2'
top.markers_bycluster <- read.delim("/public/workspace/wanggenjun/work/scRNA/estrogen/zhao/LiverAnalysis/top.markers_neu.txt",header = T)
diff.genes <- top.markers_bycluster
head(diff.genes)
dim(diff.genes)
diff.genes <- subset(diff.genes,p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)$name
# diff.genes <- subset(diff.genes,p_val_adj < 0.01)$name
head(diff.genes)
length(diff.genes)
diff.genes <- unique(diff.genes)
length(unique(diff.genes))

mycds <- setOrderingFilter(mycds, diff.genes)
M1 <- plot_ordering_genes(mycds)
M1
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
mycds <- orderCells(mycds, root_state = 2 )
traj_neu1 <- plot_cell_trajectory(mycds, color_by = "State")+
  scale_colour_manual(values = cols_1)
traj_neu2 <- plot_cell_trajectory(mycds, color_by = "Pseudotime",)
traj_neu3 <- plot_cell_trajectory(mycds, color_by = "minor2")+
  scale_colour_manual(values = c('#666666','#FB8072','#BEBADA'))
traj_neu4 <- plot_cell_trajectory(mycds, color_by = "group") + scale_color_manual(values = Colors$group.cols)
traj_neu1
traj_neu2
traj_neu3
traj_neu4

###-------------------------7I
cs2 <- table(mycds$group, mycds$State)
cs1 <- as.data.frame(cs1 * 100 / rowSums(cs1)) 
cs2 <- as.data.frame(cs2 * 100 / rowSums(cs2))
n1 <- cs1 %>% mutate(State=factor(Var2)) %>% 
  ggplot(aes(x=Var1, y=Freq, fill=State)) + 
  geom_bar(stat="identity") + rotate() + xlab("") +
  scale_fill_manual(values = c('1'='#1B9E77','2'='#D95F02','3'='#7570B3')) + theme_pubr()
n1
n2 <- cs1 %>% mutate(State=Var2) %>% 
  ggplot(aes(x=Var1, y=Freq, fill=State)) + 
  geom_bar(stat="identity") + rotate() + xlab("") +
  scale_fill_manual(values=c('1'='#1B9E77','2'='#D95F02','3'='#7570B3')) + theme_pubr()
n2
n3 <- cs2 %>% mutate(State=Var2) %>% 
  ggplot(aes(x=Var1, y=Freq, fill=State)) + 
  geom_bar(stat="identity") + rotate() + xlab("") +
  scale_fill_manual(values=c('1'='#1B9E77','2'='#D95F02','3'='#7570B3')) + theme_pubr()
n3

###-------------------------7J
BEAM_res <- BEAM(mycds, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
heatmap_byPsedotime <- plot_genes_branched_heatmap(mycds[row.names(subset(BEAM_res,qval<0.001)),],
                                                   branch_point = 1,
                                                   num_clusters = 4, #这些基因被分成几个group
                                                   cores = 1,
                                                   branch_labels = c("Cell fate 1", "Cell fate 2"),
                                                   hmcols = colorRampPalette(rev(c("#CF384D","#ED6345","#FA9A58","#F3FAAD","#D1EC9C","#96D5A4","#5BB6A9","#3682BA")))(62),# colorRampPalette(rev(brewer.pal(9, "PRGn")))(62),
                                                   branch_colors = c('#D95F02', '#1B9E77', '#7570B3'), #pre-branch, Cell fate 1, Cell fate 2分别用什么颜色
                                                   use_gene_short_name = TRUE,
                                                   show_rownames = TRUE,
                                                   return_heatmap = T)
heatmap_byPsedotime$ph_res

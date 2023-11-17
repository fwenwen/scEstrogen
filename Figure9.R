###-------------------------9A
# the same as Fig. 5A

###-------------------------9B
gene_cell_exp <- AverageExpression(Lung,
                                   features = markergenes,
                                   group.by = 'final',
                                   slot = 'data') 
gene_cell_exp <- as.data.frame(gene_cell_exp$SCT)
head(gene_cell_exp)
marker_exp <- t(scale(t(gene_cell_exp),scale = T,center = T))
marker_exp
pheatmap(marker_exp,
         cluster_cols = F, cluster_rows = F, scale = "none",
         treeheight_col = 0, treeheight_row = 0,
         display_numbers = F,
         border_color = "black",
         color = col_heatmap)

###-------------------------9C
# the same as Fig. 5C

###-------------------------9D
FeaturePlot(Lung,features= topDEGs,
            cols=c("lightgrey", brewer.orrd(10)),min.cutoff='q1')

###-------------------------9E
library(simplifyEnrichment)
go_id <- read.table('D:/bioinformation/RData/UterusAnalysis/go_id_Uterus.txt',header = TRUE, sep = "\t")
head(go_id)
mat <- GO_similarity(unique(go_id$ID))
binary_cut(mat)
cluster_terms(mat, method = "binary_cut") #和上面一样：binary_cut(mat)
df = simplifyGO(mat)
df
head(df)
sort(table(df$cluster))
split(df, df$cluster)
compare_clustering_methods(mat) #比较
compare_clustering_methods(mat, plot_type = "heatmap")

###-------------------------9F
# the same as Fig. 7G,H

###-------------------------9G
library(CytoTRACE)
plot_complex_cell_trajectory(mycds, x = 1, y = 2,color_by = "minor3")+
  scale_colour_manual(values = cols_Epi)

###-------------------------9H
Epi_minor3 <- c('COX1+ Epi',
                'NME2+ Epi',
                'WFDC2+ Epi',
                'ARGLU1+ Epi',
                'KRT16+ Epi',
                'GADD45B+ Epi',
                'GALK2+ Epi',
                'S100A1+ Epi',
                'NDRG2+ Epi',
                'CSN3+ Epi')
#cs_Stateminor3$minor3<-factor(cs_Stateminor3$minor3,levels=Epi_minor3)
#cs_Stateminor3 <- cs_Stateminor3[order(cs_Stateminor3$minor3),] 
#head(cs_Stateminor3)
Epi$State <- mycds$State
head(Epi@meta.data)
Epi@meta.data$State[Epi$State == "1"] <- "S7"
Epi@meta.data$State[Epi$State == "2"] <- "S5"
Epi$State[Epi$State == "3"] <- "S6"
Epi$State[Epi$State == "4"] <- "S4"
Epi$State[Epi$State == "5"] <- "S3"
Epi$State[Epi$State == "6"] <- "S1"
Epi$State[Epi$State == "7"] <- "S2"
table(Epi$State)
cs_Stateminor3$State[cs_Stateminor3$State == "1"] <- "S7"
cs_Stateminor3$State[cs_Stateminor3$State == "2"] <- "S5"
cs_Stateminor3$State[cs_Stateminor3$State == "3"] <- "S6"
cs_Stateminor3$State[cs_Stateminor3$State == "4"] <- "S4"
cs_Stateminor3$State[cs_Stateminor3$State == "5"] <- "S3"
cs_Stateminor3$State[cs_Stateminor3$State == "6"] <- "S1"
cs_Stateminor3$State[cs_Stateminor3$State == "7"] <- "S2"
table(cs_Stateminor3$State)

cs_Stateminor3$State<-factor(cs_Stateminor3$State,levels=c('S1','S2','S3','S4','S5','S6','S7'))
cs_Stateminor3 <- cs_Stateminor3[order(cs_Stateminor3$State),] 
head(cs_Stateminor3)
bar <- cs_Stateminor3 %>% mutate(State) %>% 
  ggplot(aes(x=minor3, y=Freq, fill=State)) + 
  geom_bar(stat="identity") + rotate() + xlab("") +
  scale_fill_manual(values = cols_1A) + theme_pubr()
bar

###-------------------------9I
Idents(Epi) <- 'State'
dge.State <- FindAllMarkers(Episub, min.pct=0.1, logfc.threshold=0.1)
sig_dge.State <- subset(dge.State, p_val_adj<0.01&abs(avg_log2FC)>1)
table((sig_dge.State %>% dplyr::filter( pct.1 >= 0.1))$cluster)

minpct <- .1
gseOut <- foreach(cls=levels(sig_dge.State$cluster)[-5]) %do% {
  submarkers <- sig_dge.State %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
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
gseResult <- do.call("rbind", lapply(gseOut, function(x){x@result}))
go_id <- gseResult
go_id <- droplevels(go_id)
head(go_id)
mxid <- as.data.frame(go_id %>% group_by(ID) %>% 
                        dplyr::summarise(cellType = cluster[which.min(p.adjust)]))
rownames(mxid) <- mxid$ID
sgseResult <- go_id %>% arrange(p.adjust)

go_id$cluster<-factor(go_id$cluster,levels=c('S1','S2','S3','S4','S5','S6','S7'))
go_id <- go_id[order(go_id$cluster),] 
head(go_id)
topPath <- go_id %>% group_by(cluster) %>% filter(pvalue < 0.01)%>% 
  top_n(n=10, wt=abs(enrichmentScore))
head(topPath)
table(topPath$cluster)
gseMat <- reshape2::acast(go_id,
                          ID~cluster, value.var='enrichmentScore')
dim(gseMat)
gseMat[is.na(gseMat)] <- 0
head(gseMat)
topgo <- unique(topPath$ID)
topgo
ngseMat0 <- gseMat[topgo,]
ngseMat0[1:3,]
dim(ngseMat0)
library(GOfuncR)
labels_row <- get_names(rownames(ngseMat0))$go_name
ngseMat0[1:3,]
pheatmap::pheatmap(ngseMat0, cluster_cols = F, cluster_rows = T,
                   color=cols_btor, labels_row = labels_row,angle_col = 45,
                   fontsize_row = 6, border_color = NA)



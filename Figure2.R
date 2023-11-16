# Figure2
for (si in as.character(unique(epi$organ.group))) { 
  sma11.meta.data=epi@meta.data %>% filter(organ.group == si)
  sma11.count=as.data.frame(epi[["RNA"]]@counts[,rownames(sma11.meta.data)]) 
  sma11.count=sma11.count[rowSums(sma11.count) > 0,]
  sma11.count=data.frame(t(sma11.count)) 
  sma11.count=sma11.count[!str_detect(rownames(sma11.count),"^MT-"),]
  write.table(sma11.count,file = paste0(si,".count.txt"),quote = F,sep = "\t",row.names = T,col.names = T)
}

library(tidyverse)

mat.files=dir("./res2/",pattern = "dt_0_02.txt$")

all.mat=data.frame()
for (fi in mat.files) {
  tmp.mat=read.table(paste0("./res2/",fi),header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  tmp.mat=as.data.frame(t(tmp.mat))
  sampleid=str_replace(fi,"\\..*$","")
  colnames(tmp.mat)=paste(sampleid,colnames(tmp.mat),sep = "_")
  tmp.mat$gene=rownames(tmp.mat)
  
  if (sampleid == "Breast_EH") {
    all.mat=tmp.mat
  }else{
    all.mat=all.mat %>% full_join(tmp.mat,by="gene") #元素的并集进行合并
  }
}
head(all.mat)
# 
signature.programs=c("Colon_EH_1","Breast_EH_5","Uterus_EH_3","Lung_EH_16","Breast_EH_2","Lung_EH_13","Breast_EH_6","Lung_EH_9")
signature.loading=all.mat[,c("gene",signature.programs)]

used.gene=c()
for (pi in signature.programs) {
  tmp.df=signature.loading[,c("gene",pi)]
  tmp.loading=tmp.df[,2]
  names(tmp.loading)=tmp.df[,1]
  
  tmp.loading=tmp.loading[!is.na(tmp.loading)]
  used.gene=append(used.gene,names(tail(sort(tmp.loading),100)))
}
used.gene=unique(used.gene)

signature.loading=signature.loading[signature.loading$gene %in% used.gene,]
rownames(signature.loading)=signature.loading$gene
signature.loading$gene=NULL
signature.loading[is.na(signature.loading)]<-0
signature.loading$total_loading=rowSums(signature.loading)
signature.loading$average_loading=signature.loading$total_loading / length(signature.programs)

signature.loading=signature.loading%>%arrange(desc(average_loading))
head(rownames(signature.loading),30)
# 

signature.loading[1:30,]

SpecificityRank <- function(rssMat, cluster, topn = 30) {
  data <- tibble::enframe(sort(rssMat[, as.character(cluster)], decreasing = T)) %>%
    magrittr::set_colnames(c('regulon', 'RSS')) %>%
    tibble::rowid_to_column("index")
  data <- data.frame(gene=rownames(signature.loading),
                     average_score=signature.loading$average_loading) %>% 
    tibble::rowid_to_column("index")
  head(data)
  data$pt.col <- ifelse(data$index <= topn, "#007D9B", "#BECEE3")
  data <- head(data, n=100)
  data.label <- head(data, n=topn)
  
  pp <- ggplot(data, aes(index, average_score)) +
    geom_point(size=3, color=data$pt.col) +
    ggrepel::geom_text_repel(inherit.aes = FALSE, data = data.label, aes(index, average_score, label=gene), size=4,
                             max.overlaps = 30) +
    # ggtitle(as.character(cluster)) + 
    ylab("Average score") + xlab('Genes') +
    theme_bw(base_size = 12) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(color="black"),
          axis.ticks = element_line(color="black"),
          axis.text = element_text(color = "black"),
          plot.title = element_text(hjust = .5)
    )
  return(pp)
}

pdf("./gene_top30.pdf",height = 4.27,width = 4.87)
print(pp)
dev.off()

library(org.Hs.eg.db)
library(foreach)
library(enrichplot)
library(clusterProfiler)
library(dplyr)

library(org.Mfascicularis.eg.db)
keytypes(org.Mfascicularis.eg.db)
head(keys(org.Mfascicularis.eg.db,keytype = "GID"))

genes <- head(rownames(signature.loading),30)
refGenes <- seuratObj@misc$geneName
genes <- names(refGenes)[match(genes,refGenes)]

terms <- enrichGO(OrgDb="org.Mfascicularis.eg.db",
                  gene = genes,
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.1,
                  keyType = 'GID',
                  pAdjustMethod = 'fdr',
                  ont = "BP",
                  readable=TRUE)

head(terms)
dim(terms)
library(topGO)
library(enrichplot)
goplot(terms)

plotGOgraph(terms)

grep("T\\ cell",terms@result$Description)
grep("immu",terms@result$Description)

GOs
terms@result <- terms@result[terms@result$Description %in% GOs,]

Gresult <- terms@result
Gresult <- Gresult[-grep("T\\ cell",Gresult$Description)[-(1:2)],]
dim(Gresult)
Gresult <- Gresult[-grep("immu",Gresult$Description)[-(1:2)],]
Gresult$Count <- as.numeric(gsub("/.*$","",Gresult$BgRatio))
Gresult <- Gresult %>% mutate(logFDR=-log10(pvalue)) %>%  dplyr::top_n(n = 30,wt = logFDR)


get("Gresult")
#GO
GO_bar <- function(x){
  y <- get(x)
  ggplot(data = y,
         aes(x = Count,
             y = Description,
             fill = -log10(pvalue))) +
    scale_y_discrete(labels = function(y) str_wrap(y, width = 50) ) + #label换行，部分term描述太长
    geom_bar(stat = "identity",width = 0.8) +
    labs(x = "Gene Number",
         y = "Description",
         title = paste0(x," of GO enrichment barplot")) +
    theme_bw() 
}

#MF:
GO_bar("Gresult")+scale_fill_distiller(palette = "Reds",direction = 1)

pdf("./gene_top30terms.pdf",height = 5.27,width = 8.27)
print(GO_bar("Gresult")+scale_fill_distiller(palette = "Reds",direction = 1))
dev.off()


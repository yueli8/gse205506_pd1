library(monocle)
setwd("d:/R/data/GSE227744")
Neutrophil<-readRDS(file = "Neutrophil2_cluster_id.rds")
Neutrophil_matrix<- as(as.matrix(Neutrophil@assays$RNA@counts),'sparseMatrix')
Neutrophil_data<- Neutrophil@meta.data
Neutrophil_data$celltype<-Neutrophil@active.ident
Neutrophil_f_data <- data.frame(gene_short_name = row.names(Neutrophil),row.names = row.names(Neutrophil))
pd<-new('AnnotatedDataFrame',data=Neutrophil_data)
fd<- new('AnnotatedDataFrame',data=Neutrophil_f_data)
#你可以去看一下identical(rownames(fd),rownames(Neutrophil_matrix))返回的结果是否为TURE，若不是可以运行下面的代码
Neutrophil_matrix<- Neutrophil_matrix[rownames(fd), ]

cds<- newCellDataSet(Neutrophil_matrix,
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit =0.5,
                     expressionFamily =negbinomial.size())
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
cds <- detectGenes(cds,min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >=10))  #过滤掉在小于10个细胞中表达的基因，还剩11095个基因:

express_genes<-VariableFeatures(Neutrophil)
cds<-setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)

diff<-differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~celltype",cores = 1)
deg<-subset(diff,qval<0.01)
deg<-deg[order(deg$qval,decreasing = F),]
write.table(deg, file = "deg.xls",col.names = T,row.names = F,sep = "\t",quote = F)

ordergene<-rownames(deg)
cds<-setOrderingFilter(cds,ordergene)
pdf("ordergenes.pdf")
plot_ordering_genes(cds)

cds<-reduceDimension(cds,max_components = 2, method='DDRTree')
cds<- orderCells(cds) # In dfs(graph = graph, root = root, mode = mode, unreachable = unreachable,  :Argument `neimode' is deprecated; use `mode' instead报错没关系，可以继续跑
plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "celltype",size=1,show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone = TRUE)

time_diff<-differentialGeneTest(cds[ordergene,],cores = 1,fullModelFormulaStr = "~sm.ns(Pseudotime)")
library(magrittr)
library(tidyverse)
time_genes<-time_diff %>% pull(gene_short_name) %>% as.character()
plot_pseudotime_heatmap(cds[time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)

time_genes<-top_n(time_diff,n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character() #提取前100个基因
plot_pseudotime_heatmap(cds[time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)
#https://www.jianshu.com/p/5d6fd4561bc0
library(Seurat)
library(ggplot2)
library(cowplot)
library(scater)
library(scran)
library(BiocParallel)
library(BiocNeighbors)
library(data.table)
library(dplyr)
library(Matrix)
library(clustree)
library(monocle)
library(ggsci)
library(ggpubr)
setwd("~/gse205506/finalize")
Only_Epithelial<-readRDS("Epithelial_monocle2.rds")
DimPlot(Only_Epithelial, reduction = "umap", label = TRUE, pt.size = 0.01) 
DimPlot(Only_Epithelial, reduction = "umap", label = FALSE, pt.size = 0.01) 
DimPlot(Only_Epithelial, reduction = "umap")
Epithelial<-readRDS("Epithelial_monocle2.rds")

Epithelial_matrix<- as(as.matrix(Epithelial@assays$RNA@counts),'sparseMatrix')
Epithelial_data<- Epithelial@meta.data
Epithelial_data$celltype<-Epithelial@active.ident
Epithelial_f_data <- data.frame(gene_short_name = row.names(Epithelial),row.names = row.names(Epithelial))
#直接读取表达矩阵来创建
pd<-new('AnnotatedDataFrame',data=Epithelial_data)
fd<- new('AnnotatedDataFrame',data=Epithelial_f_data)
#你可以去看一下identical(rownames(fd),rownames(Epithelial_matrix))返回的结果是否为TURE，若不是可以运行下面的代码
Epithelial_matrix<- Epithelial_matrix[rownames(fd), ]
cds<- newCellDataSet(Epithelial_matrix,
                     phenoData = pd,
                     featureData = fd,
                     lowerDetectionLimit =0.5,
                     expressionFamily =negbinomial.size())
#估计size factor和离散度
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)
#过滤低质量的细胞
cds <- detectGenes(cds,min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed >=0))  #过滤掉在小于10个细胞中表达的基因，
#轨迹定义基因选择及可视化和构建轨迹
#使用seurat選擇的高變基因
express_genes<-VariableFeatures(Epithelial)
cds<-setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)
#使用cluster差異表達基因
#deg.cluster<-FindAllMarkers(Epithelial)
#express_genes<-subset(deg.cluster,p_val_adj<0.05)$gene
#cds<-setOrderingFilter(cds,express_genes)
#plot_ordering_genes(cds)
#使用monocle選擇的高變基因
#disp_table<-dispersionTable(cds)
#disp.genes<-subset(disp_table,mean_expression >=0.1 & dispersion_empirical >= 1*dispersion_fit)$gene_id
#cds<-setOrderingFilter(cds,disp.genes)
#plot_ordering_genes(cds)

diff<-differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~celltype",cores = 1)
head(diff)
deg<-subset(diff,qval<0.01)
deg<-deg[order(deg$qval,decreasing = F),]
head(deg)
write.table(deg, file = "deg.xls",col.names = T,row.names = F,sep = "\t",quote = F)

ordergene<-rownames(deg)
cds<-setOrderingFilter(cds,ordergene)
#pdf("ordergenes.pdf")
plot_ordering_genes(cds)
#使用DDRTree降維
cds<-reduceDimension(cds,max_components = 2, method='DDRTree')
cds<- orderCells(cds)
#使用root_state參數可以設置擬時間軸的根,如下面的擬時間着色圖中可以看出,左邊是根.
#根據state圖可以看出,根是state1,若要想把另一端設爲根,可以按如下操作
# In dfs(graph = graph, root = root, mode = mode, unreachable = unreachable,  :Argument `neimode' is deprecated; use `mode' instead报错没关系，可以继续跑
#pdf("train.monocle.pseudotime.pdf",width=7,height=7)
#cds<- orderCells(cds,root_state = 2)

plot_cell_trajectory(cds, color_by = "Pseudotime",size=1,show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "celltype",size=1,show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "State",size=1,show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "cell",size=1,show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "tech",size=1,show_backbone = TRUE)
plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap("~cell",nrow = 1)+ scale_color_npg()
plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap("~cell",nrow = 1)+scale_color_nejm()
colour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#14a0dc","#b414dc")
plot_cell_trajectory(cds, color_by = "celltype") + facet_wrap("~cell",nrow = 1)+scale_color_manual(values = colour)

#覺得這種規矩不太直觀,可以畫成樹形圖
p1<-plot_cell_trajectory(cds, x=1,y=2,color_by = "celltype") + 
  theme(legend.position = 'none',panel.border = element_blank())+#去掉第一個的legend
scale_color_manual(values = colour)
p2<-plot_complex_cell_trajectory(cds, x=1,y=2,color_by = "celltype") + 
  scale_color_manual(values = colour)+
  theme(legend.title = element_blank())
p1|p2 
 
p1<-plot_cell_trajectory(cds, x=1,y=2,color_by = "celltype") + facet_wrap("~cell",nrow = 1)+
   theme(legend.position = 'none',panel.border = element_blank())+#去掉第一個的legend
   scale_color_manual(values = colour)
p2<-plot_complex_cell_trajectory(cds, x=1,y=2,color_by = "celltype") +facet_wrap("~cell",nrow = 1)+
   scale_color_manual(values = colour)+
   theme(legend.title = element_blank())
p1|p2 
 
#畫沿時間軸的細胞密度
df<-pData(cds)
View(df)
ggplot(df,aes(Pseudotime,colour=celltype,fill=celltype))+facet_wrap("~cell",nrow = 1)+
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()

#手動設置顏色
ClusterName_color_panel<-c("goblet"="#DC143C","basal"="#0000FF","enterocyte"="#20B2AA",
"FGFR1_High_NME5"="#FFA500","cortical_thymus"="#9370DB","ionocytes"="#98FB98",
"FOXN4+"="#F08080","club"="#0000FF")

ggplot(df,aes(Pseudotime,colour=celltype,fill=celltype))+facet_wrap("~cell",nrow = 1)+
  geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic2()+
  scale_fill_manual(name="",values = ClusterName_color_panel)+
  scale_color_manual(name="",values = ClusterName_color_panel)

#指定基因的可視化
#選擇前4個基因並將其對象取出
keygenes<-head(ordergene,4)
cds_subset<-cds[keygenes,]
plot_genes_in_pseudotime(cds_subset,color_by = "State")
plot_genes_in_pseudotime(cds_subset,color_by = "celltype")
plot_genes_in_pseudotime(cds_subset,color_by = "Pseudotime")
#指定基因
s.genes<-c("ANPEP","CA2")

plot_genes_in_pseudotime(cds[s.genes,],color_by = "celltype")
plot_genes_jitter(cds[s.genes,],color_by = "celltype")
plot_genes_violin(cds[s.genes,],color_by = "celltype")

#擬時序展示單個基因表達量
colnames(pData(cds))
pData(cds)$CA2 =log2(exprs(cds)['CA2',]+1)
plot_cell_trajectory(cds,color_by = "CA2")+scale_color_gsea()

pData(cds)$ANPEP =log2(exprs(cds)['ANPEP',]+1)
plot_cell_trajectory(cds,color_by = "ANPEP")+scale_color_gsea()
#尋找擬時相關的基因
Time_diff<-differentialGeneTest(cds[ordergene,],cores = 1,
                                fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff<-Time_diff[,c(5,2,3,4,1,6,7)]
Time_genes<-Time_diff %>% pull(gene_short_name) %>%  as.character()
plot_pseudotime_heatmap(cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)
p=plot_pseudotime_heatmap(cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)
#前面通過設置num_clusters將熱圖聚成四個cluster,把每一個cluster的基因單獨提出來分析
p$tree_row
clusters<-cutree(p$tree_row,k=4)
clustering<-data.frame(clusters)
clustering[,1]<-as.character(clustering[,1])
colnames(clustering)<-"Gene_Clusters"
table(clustering)
#提取前100個
Time_genes<-top_n(Time_diff,n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)
#顯著差異基因按熱圖結果排序並保存
hp.genes<-p$tree_row$table[p$tree_row$order]
Time_diff_sig<-Time_diff[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"Time_diff_sig,csv",row.names = F)

#提取前100個
Time_genes<-top_n(Time_diff,n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)
#顯著差異基因按熱圖結果排序並保存
hp.genes<-p$tree_row$table[p$tree_row$order]
Time_diff_sig<-Time_diff[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"Time_diff_sig,csv",row.names = F)

#提取前50個
Time_genes<-top_n(Time_diff,n=25,desc(qval)) %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)
#顯著差異基因按熱圖結果排序並保存
hp.genes<-p$tree_row$table[p$tree_row$order]
Time_diff_sig<-Time_diff[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"Time_diff_sig,csv",row.names = F)

features<-c("CA2","MLEC","LGALS4")
VlnPlot(Epithelial,features=features,group.by = "cell",ncol=3)
plot_genes_in_pseudotime(cds[features,],color_by="celltype")
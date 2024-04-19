#install.packages('devtools')
#devtools::install_github("jinworks/CellChat")
library(CellChat)
library(tidyverse)
library(ggalluvial)
#install.packages("anndata")
library(anndata)
library(anndata)
library(Seurat)
library(patchwork)
library(uwot)
library(reticulate)
library(glmGamPoi)
library(umap)
library(ComplexHeatmap)
#install.packages('NMF')
library(NMF)
library(dplyr)
library(SeuratData)
library(ggplot2)
library(svglite)
library(wordcloud)
library(wordcloud2)
library(tm)
#devtools::install_github("jokergoo/circlize")
library(circlize)
#devtools::install_github("jokergoo/ComplexHeatmap")

#source('functional.R')#千萬不要加上
setwd("/home/r2309261027/integrated/GJ")
setwd("/home/r2309261027/integrated/CYTOTRACE_178318")
rm(list=ls()) #清空所有变量
options(stringsAsFactors = FALSE)#输入数据不自动转换成因子（防止数据格式错误）

#(A) Starting from a count data matrix
#load("data_humanSkin_CellChat.rda")
#data.input = data_humanSkin$data
#meta = data_humanSkin$meta
#cell.use = rownames(meta)[meta$condition == "LS"] # extract the c
#data.input = data.input[, cell.use]
#meta = meta[cell.use, ]

#(B) Starting from a Seurat object
#data.input <- GetAssayData(seurat_object, assay = "RNA", slot = "data") # normalized data matrix
#labels <- Idents(seurat_object)
#meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

#(C) Starting from a SingleCellExperiment object
#data.input <- SingleCellExperiment::logcounts(object) # normalized data matrix
#meta <- as.data.frame(SingleCellExperiment::colData(object)) # extract a dataframe of the cell labels
#meta$labels <- meta[["sce.clusters"]]

#(D) Starting from an Anndata object
# read the data into R using anndata R package
#ad <- read_h5ad("scanpy_object.h5ad")# access count data matrix
#counts <- t(as.matrix(ad$X)) # normalize the count data if the normalized data is not available in the .h5ad file
#library.size <- Matrix::colSums(counts)
#data.input <- as(log1p(Matrix::t(Matrix::t(counts)/library.size)* 10000), "dgCMatrix")
#meta <- ad$obs# access meta data
#meta$labels <- meta[["clusters"]]
#GJ的数据引入代码
#分组
sce <- readRDS("All_right_cluster_id_test.rds")
sce_LM <- subset(sce,subset = tech =="metastasis")
sce_Tumor <- subset(sce,subset = tech =="tumor")
sce_PBMC <- subset(sce,subset = tech =="pbmc")
saveRDS(sce_LM, file = "sce_LM.rds")
saveRDS(sce_Tumor, file = "sce_Tumor.rds")
saveRDS(sce_PBMC, file = "sce_PBMC.rds")
#2.Create a CellChat object by following option
#(A) Starting from the digital gene expression matrix and cell label information
#cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
# [1] "Create a CellChat object from a data matrix"
# Set cell identities for the new CellChat object 
# The cell groups used for CellChat analysis are  APOE+ FIB FBN1+ FIB COL11A1+ FIB Inflam. FIB cDC1 cDC2 LC Inflam. DC TC Inflam. TC CD40LG+ TC NKT 
#(B) Starting from a Seurat object
#cellChat <- createCellChat(object = seurat.obj, group.by = "ident", assay = "RNA")
#(C) Starting from a SingleCellExperiment object
#cellChat <- createCellChat(object = sce.obj, group.by = "sce.clusters")
#(D) Starting from an AnnData object
#sce <- zellkonverter::readH5AD(file = "adata.h5ad") # retrieve all the available assays within sce object assayNames(sce)
# add a new assay entry "logcounts" if not available
#counts <- assay(sce, "X") # make sure this is the original count data matrix
#library.size <- Matrix::colSums(counts)
#logcounts(sce) <- log1p(Matrix::t(Matrix::t(counts)/library.size) * 10000)
# extract a cell meta data
#meta <- as.data.frame(SingleCellExperiment::colData(sce)) #
#cellChat <- createCellChat(object = sce, group.by = "sce.clusters")
rm(list = ls())
sce_ac <- readRDS("Metastasis_cluster_id.rds")
#创建cellchat对象
sce_ac[['celltype']] <- Idents(sce_ac)#将celltype指定为细胞类型
cellchat <- createCellChat(sce_ac@assays$RNA@data, meta = sce_ac@meta.data, group.by = "celltype")
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents))
CellChatDB <- CellChatDB.human#导入配受体库
str(CellChatDB) #查看数据库信息
#以下是李老师代码
#CellChatDB <- CellChatDB.human#重复；导入配受体库
#showDatabaseCategory(CellChatDB)#查看数据库信息
# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)#显示数据集结构
# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")# use Secreted Signaling；选择特定的信号来进行分析，这里还可以选择ECM-receptor和Cell-Cell Contact。 
#CellChatDB.use <- CellChatDB # simply use the default CellChatDB to use all CellChatDB for cell-cell communication analysis
#这是使用所有 use all CellChatDB for cell-cell communication analysis
# set the used database in the object
cellchat@DB <- CellChatDB.use#这行代码将之前筛选得到的特定信号的数据库 CellChatDB.use 赋值给 cellchat 对象的 DB 属性。这样，cellchat 对象将使用这个特定信号的数据库进行细胞间通信分析。

#数据预处理
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
#这一步会首先需要选取在上一步选择的interaction database中的基因，这一步并不是一般意义的subset，因为它还会将对应基因的表达矩阵赋值给cellchat@data.Signaling，所以不是可选项，而是必须有的一步
future::plan("multisession", workers = 4) # do parallel一行用于设置并行计算的代码。
#接下来需要识别各个细胞类型过表达（差异）的基因和互作；
#各个细胞类型过表达的基因是基于表达该基因的细胞比例，差异倍数，和p值来判定的。默认参数是，细胞比例阈值为0（thresh.pc = 0），差异倍数为0（thresh.fc = 0，only.pos = TRUE），p值0.05（thresh.p = 0.05）（这里不是校正后的p值，虽然源码中有做bonferroni矫正，但是没有用这个值做判定）。最后过表达的基因会存于object@var.features，基因名存在形式为cellchat@var.features$features，差异基因计算结果的表格存在形式为cellchat@var.features$features.info然后这些过表达基因所在的互作即过表达互作。过表达互作会存于cellchat@LR$LRsig 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human) 
cellchat@data.project[1:4,1:4]
#为可选项，找到配体受体关系之后，projectData将配体受体对的表达值投射到PPI，来对@data.signaling表达值进行校正。结果保存在@data.project
#Inference of cell-cell communication network
#cellchat <- computeCommunProb(cellchat, type = "triMean")
#triMean is used for calculating the average gene expression per cell group. 
#[1] ">>> Run CellChat on sc/snRNA-seq data <<< [2023-12-15 12:49:11.48819]"
#|=============================================================================================================| 100%
#[1] ">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [2023-12-15 12:51:34.162272]"
#使用表达值推测细胞互作的概率，该步骤相对较耗时一些。
#注2：在假设细胞数较多的群 往往比 细胞数较少的群发送更强的信号的前提下，当population.size = TRUE时候，CellChat可以在概率计算中考虑每个细胞群中细胞比例的影响。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE)#如果不想用上一步PPI矫正的结果，raw.use = TRUE即可
#先保存所有的结果
#all the inferred cell-cell communications at the level of ligands/receptors
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "LM_cell-cell_communications.all.csv")

#access the the inferred communications at the level of signaling pathways
#df.net1 <- subsetCommunication(cellchat,slot.name = "netP")

#gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
#levels(cellchat@idents)
#df.net2 <- subsetCommunication(cellchat, sources.use = c("Epi"), targets.use = c("Fibroblast" ,"T")) 

#gives the inferred cell-cell communications mediated by signaling WNT and TGFb.
#df.net3 <- subsetCommunication(cellchat, signaling = c("CCL", "TGFb"))
###Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
#使用 computeCommunProbPathway 函数对 cellchat 对象进行通信通路水平的细胞间通信概率计算。
cellchat <- computeCommunProbPathway(cellchat)
#使用 aggregateNet 函数对 cellchat 对象进行网络聚合，将相同的配体-受体对进行合并，得到整体的细胞间通信网络。
cellchat <- aggregateNet(cellchat)
saveRDS(cellchat, file = "cellchat_LM.rds")
rm(list = ls())
cellchat <- readRDS("cellchat_LM.rds")
#Visualization of cell-cell communication network 
groupSize <- as.numeric(table(cellchat@idents))

netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")

netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength") 
#看不同细胞群相互作用的数量和强度
gg1 <- netVisual_heatmap(cellchat, object = 1)
gg2 <- netVisual_heatmap(cellchat, measure = "weight",comparison = c(1,2))
gg1 + gg2#热图展示
# 分细胞类型展示
# control the parameter edge.weight.max so that we can compare edge weights between differet networks.
weight_mat <- cellchat@net$weight
par(mfrow = c(4,4),mgp=c(0,0,0), xpd=TRUE)
par(mar = c(2, 2, 2, 2))  # 设置边距为2个字符
for (cel in unique(cellchat@idents)){
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle( cir_mat, vertex.weight= groupSize, weight.scale= T,edge.weight.max = max(weight_mat), vertex.label.cex=0.4,title.name=cel)
}
#细胞通路展示
library(CellChat)
library(ggalluvial)
library(ggplot2)
library(patchwork)
library(ggalluvial)
library(ComplexHeatmap)
library(igraph)
library(Seurat)
library(svglite)
# Signaling role analysis on the aggregated cell-cell communication network naling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",width = 6, height = 13)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",width = 6, height = 13)
draw(ht1 + ht2, ht_gap = unit(2, "cm"))
#通过整合特定通路的所有受配体对展示推断的信号网络 主要功能函数：netVisual_aggregate
pathways.show.all <- cellchat@netP$pathways
#select one pathway
#(A)Circle plot
pathways.show <- c("CXCL")
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", color.use = NULL, sources.use = NULL, targets.use = NULL)
#(B) Hierarchy plot
cellchat@meta$celltype %>% head()
vertex.receiver = seq(1,4)
netVisual_aggregate(cellchat, signaling = pathways.show, layout ="hierarchy", vertex.receiver = vertex.receiver)
#(C) Chord diagram#暂时画不出来
par(mar = c(5, 5, 4, 2))
netVisual_aggregate(cellchat, signaling = pathways.show, layout ="chord")
par(mfrow=c(1,1))
group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 8)) #grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group =group.cellType, title.name = paste0(pathways.show, " signaling network"))
#(D) Heatmap plot
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR.CXCL <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR.CXCL[3,]
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
#(A) Bubble plot
# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_bubble(cellchat, sources.use = 9, targets.use = c(1:16), remove.isolate = TRUE)
#以上是把所有的信号通路都指定了
# (2) show all the significant interactions (L-R pairs) associated with certain signaling pathways
netVisual_bubble(cellchat, sources.use = 9, targets.use = c(1:16), signaling = c("CCL","CXCL","MIF","ANNEXIN"), remove.isolate = TRUE)
#以上是指定信号通路
## (3) show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("CCL","CXCL","FGF","ANNEXIN"))
netVisual_bubble(cellchat, sources.use = 9, targets.use = c(1:16), pairLR.use = pairLR.use, remove.isolate = TRUE)
# set the order of interacting cell pairs on x-axis
# (4) Default: first sort cell pairs based on the appearance of sources in levels(object@idents), and then based on the appearance of targets in levels(object@idents)
# (5) sort cell pairs based on the targets.use defined by users
netVisual_bubble(cellchat, targets.use = c("Exhausted_CD8_T","Cytotoxic_CD8_T"), pairLR.use = pairLR.use, 
                 remove.isolate = TRUE, sort.by.target = T)
# (6) sort cell pairs based on the sources.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE,sort.by.source = T)
# (7) sort cell pairs based on the sources.use and then targets.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), pairLR.use = pairLR.use, remove.isolate = TRUE, 
                 sort.by.source = T, sort.by.target = T)
# (8) sort cell pairs based on the targets.use and then sources.use defined by users
netVisual_bubble(cellchat, sources.use = c("FBN1+ FIB","APOE+ FIB","Inflam. FIB"), targets.use = c("LC","Inflam. DC","cDC2","CD40LG+ TC"), 
                 pairLR.use = pairLR.use, remove.isolate = TRUE, sort.by.source = T, sort.by.target = T, sort.by.source.priority = FALSE)
#5\6\7\8待定

#(A)Compute and visualize the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")# the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groupsnetAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
#(B) Visualize dominant senders (sources) and receivers (targets) in a 2D space
netAnalysis_signalingRole_scatter(cellchat, signaling = NULL)
#(C) Identify signals contributing the most to outgoing or incoming signaling of certain cell groups
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2
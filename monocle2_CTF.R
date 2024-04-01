#monocle2

##系统报错改为英文
Sys.setenv(LANGUAGE = "en")
##禁止转化为因子
options(stringsAsFactors = FALSE)
##清空环境
rm(list=ls())

library(dplyr)
library(Seurat)
library(tidyverse)
library(patchwork)
#BiocManager::install("monocle")
library(monocle)

setwd("D:\\MM\\WBM\\HM28\\Plasma\\")
scRNA_harmony=readRDS("scRNA_harmony_HM28.Plasma.rds")

table(scRNA_harmony@meta.data$seurat_clusters)

#可以选择cluster，也可以选择细胞类型做轨迹分析
Idents(scRNA_harmony)="seurat_clusters"
scRNA.T=subset(scRNA_harmony,ident=c(2,6))

Idents(scRNA_harmony)="anno"
scRNA.T=subset(scRNA_harmony,ident="T cells")

T_matrix <- as(as.matrix(scRNA.T@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = scRNA.T@meta.data)
fData <- data.frame(gene_short_name = row.names(T_matrix), row.names = row.names(T_matrix))
fd <- new('AnnotatedDataFrame', data = fData)
cds <- newCellDataSet(T_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 0.1)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 10))

cds1=cds
##选择离散程度高的基因作为轨迹分析
disp_table <- dispersionTable(cds1)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds1 <- setOrderingFilter(cds1, disp.genes)
plot_ordering_genes(cds1)

cds1 <- reduceDimension(cds1, max_components = 2,
                        method = 'DDRTree')

cds1 <- orderCells(cds1,root_state = NULL)
plot_cell_trajectory(cds1, color_by = "seurat_clusters")
plot_cell_trajectory(cds1, color_by = "State")
plot_cell_trajectory(cds1, color_by = "Pseudotime")
plot_cell_trajectory(cds1, color_by = "group.a")##group.a为metadata中的行名group.a
plot_cell_trajectory(cds1, color_by = "Pseudotime") +facet_wrap(~group.a, nrow = 3)#可分组展示拟时间

cds1_expressed_genes <-  row.names(subset(fData(cds1),
                                          num_cells_expressed >= 10))
cds1_filtered <- cds1[cds1_expressed_genes,]
my_genes <- row.names(subset(fData(cds1_filtered),
                             gene_short_name %in% c("CSDE1","IGHG4", "IFITM1")))
cds_subset <- cds1_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters")
plot_genes_in_pseudotime(cds_subset, color_by =  "State")
plot_genes_in_pseudotime(cds_subset, color_by =  "anno")
plot_genes_in_pseudotime(cds_subset, color_by =  "group.a")


genes <-  c("CSDE1","IGHG4", "IFITM1")
p1 <- plot_genes_jitter(cds1[genes,], grouping = "State", color_by = "State")
p2 <- plot_genes_violin(cds1[genes,], grouping = "State", color_by = "State")
p3 <- plot_genes_in_pseudotime(cds1[genes,], color_by = "State")

marker_genes <- row.names(subset(fData(cds1),
                                 gene_short_name %in% c("CSDE1","IGHG4", "IFITM1")))
plot_genes_in_pseudotime(marker_genes, color_by ="State")
diff_test_res <- differentialGeneTest(cds1[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))

plot_pseudotime_heatmap(cds1[sig_gene_names,],
                        num_clusters = 5,
                        cores = 1,
                        show_rownames = T)

BEAM_cds1 <- BEAM(cds1, branch_point = 1, cores = 4)#运行时间较长
BEAM_cds1 <- BEAM_cds1[order(BEAM_cds1$qval),]
BEAM_cds1 <- BEAM_cds1[,c("gene_short_name", "pval", "qval")]
plot_genes_branched_heatmap(cds1[row.names(subset(BEAM_cds1,
                                                  qval < 1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)

genes <- row.names(subset(fData(cds1),
                          gene_short_name %in% c("CSDE1","IGHG4", "IFITM1")))

plot_genes_branched_pseudotime(cds1[genes,],
                               branch_point = 1,
                               color_by = "group.a",
                               ncol = 1)
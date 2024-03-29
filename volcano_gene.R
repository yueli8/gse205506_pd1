install.packages("devtools")
devtools::install_github("BioSenior/ggVolcano")
library(ggVolcano)


res <- read.csv("CXCL8_IL1B_Neutrophil.csv", header=TRUE)
res_1 <- add_regulate(res, log2FC_name = "log2FoldChange", fdr_name = "padj",log2FC = 2, fdr = 0.01)
ggvolcano(res_1, x = "log2FoldChange", y = "padj", label = "row", label_number = 10, output = FALSE)

#更改火山图颜色
ggvolcano(res_1, x = "log2FoldChange", y = "padj",
          fills = c("#002394","#ABABAB","#D00000"),
          colors = c("#002394","#ABABAB","#D00000"),
          pointsize = 2,  ##调节点大小
          label = "row", label_number = 10, output = FALSE)

Usage
ggvolcano(
  data,
  x = "log2FoldChange",
  y = "padj",
  pointSize = 1,
  pointShape = 21,
  fills = c("#00AFBB", "#999999", "#FC4E07"),
  colors = c("#00AFBB", "#999999", "#FC4E07"),
  x_lab = NULL,
  y_lab = NULL,
  legend_title = NULL,
  legend_position = "UL",
  log2FC_cut = 1,
  FDR_cut = 0.05,
  add_line = TRUE,
  add_label = TRUE,
  label = "row",
  label_number = 10,
  custom_label = c("目标基因1", "目标基因2", "目标基因3"),
  output = TRUE,
  filename = "volcano_plot"
)
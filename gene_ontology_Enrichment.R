library(ggplot2)
setwd("~/gse205506/gsea")
a<-read.table("epithelial_gsea", head=TRUE, sep="\t")
S1<- ggplot(a, aes(x= NES, y=reorder(GSEA, NES), size=SIZE,fill=FDR.q.val)) + geom_point(shape = 21) +theme_bw() +theme()
S1=S1+ scale_fill_continuous(low = '#d90424', high = '#374a89')+scale_x_continuous(
  labels = scales::number_format(accuracy = 0.1))+ theme(axis.text.y = element_text(size = 16,
                                                                                    family = "Arial", color = "black", face = "bold"))+xlab("NES")+scale_x_continuous(breaks = seq(0, 3, by = 0.2))
S1

a<-read.table("cd4_naive_gsea", head=TRUE, sep="\t")
S1<- ggplot(a, aes(x= NES, y=reorder(GSEA, NES), size=SIZE,fill=FDR.q.val)) + geom_point(shape = 21) +theme_bw() +theme()
S1=S1+ scale_fill_continuous(low = '#d90424', high = '#374a89')+scale_x_continuous(
  labels = scales::number_format(accuracy = 0.1))+ theme(axis.text.y = element_text(size = 16,
                                                                                    family = "Arial", color = "black", face = "bold"))+xlab("NES")+scale_x_continuous(breaks = seq(0, 3, by = 0.2))
S1

a<-read.table("cd4_chronic_activation_gsea", head=TRUE, sep="\t")
S1<- ggplot(a, aes(x= NES, y=reorder(GSEA, NES), size=SIZE,fill=FDR.q.val)) + geom_point(shape = 21) +theme_bw() +theme()
S1=S1+ scale_fill_continuous(low = '#d90424', high = '#374a89')+scale_x_continuous(
  labels = scales::number_format(accuracy = 0.1))+ theme(axis.text.y = element_text(size = 16,
                                                                                    family = "Arial", color = "black", face = "bold"))+xlab("NES")+scale_x_continuous(breaks = seq(0, 3, by = 0.2))
S1

a<-read.table("cd4_regulatory_gsea", head=TRUE, sep="\t")
S1<- ggplot(a, aes(x= NES, y=reorder(GSEA, NES), size=SIZE,fill=FDR.q.val)) + geom_point(shape = 21) +theme_bw() +theme()
S1=S1+ scale_fill_continuous(low = '#d90424', high = '#374a89')+scale_x_continuous(
  labels = scales::number_format(accuracy = 0.1))+ theme(axis.text.y = element_text(size = 16,
                                                                                    family = "Arial", color = "black", face = "bold"))+xlab("NES")+scale_x_continuous(breaks = seq(0, 3, by = 0.2))
S1

a<-read.table("cd8_cytotoxic_gsea", head=TRUE, sep="\t")
S1<- ggplot(a, aes(x= NES, y=reorder(GSEA, NES), size=SIZE,fill=FDR.q.val)) + geom_point(shape = 21) +theme_bw() +theme()
S1=S1+ scale_fill_continuous(low = '#d90424', high = '#374a89')+scale_x_continuous(
  labels = scales::number_format(accuracy = 0.1))+ theme(axis.text.y = element_text(size = 16,
                                                                                    family = "Arial", color = "black", face = "bold"))+xlab("NES")+scale_x_continuous(breaks = seq(0, 3, by = 0.2))
S1

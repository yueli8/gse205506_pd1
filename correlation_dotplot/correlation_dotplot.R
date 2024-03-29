library(ggplot2)
setwd("~/gse205506/gsea")
a<-read.table("ANPEP_CA2", head=TRUE, sep="\t")
ggplot(a,aes(x=features.plot,y=id))+geom_point()
s1<-ggplot(a,aes(x=features.plot,y=id))+geom_point(aes(size=`chisq.test`, color=`correlation`))  
s2<-s1+ 
  theme(axis.text.x=element_text(size = 16,angle=90,hjust = 1,vjust=0.5,family="Arial",color="black",face="bold"))+
  theme(axis.text.y.left =element_text(size = 16,family="Arial",color="black",face="bold"))+scale_color_gradient(low="blue",high="red")
labs(x=NULL,y=NULL)
s2

library(ggplot2)
setwd("~/gse205506/gsea")
a<-read.table("CDC42SE2_HLA-A_KLF6", head=TRUE, sep="\t")
ggplot(a,aes(x=features.plot,y=id))+geom_point()
s1<-ggplot(a,aes(x=features.plot,y=id))+geom_point(aes(size=`chisq.test`, color=`correlation`))  
s2<-s1+ 
  theme(axis.text.x=element_text(size = 16,angle=90,hjust = 1,vjust=0.5,family="Arial",color="black",face="bold"))+
  theme(axis.text.y.left =element_text(size = 16,family="Arial",color="black",face="bold"))+scale_color_gradient(low="blue",high="red")
labs(x=NULL,y=NULL)
s2

library(ggplot2)
setwd("~/gse205506/gsea")
a<-read.table("NINJ2_UTS2_MT2A", head=TRUE, sep="\t")
ggplot(a,aes(x=features.plot,y=id))+geom_point()
s1<-ggplot(a,aes(x=features.plot,y=id))+geom_point(aes(size=`chisq.test`, color=`correlation`))  
s2<-s1+ 
  theme(axis.text.x=element_text(size = 16,angle=90,hjust = 1,vjust=0.5,family="Arial",color="black",face="bold"))+
  theme(axis.text.y.left =element_text(size = 16,family="Arial",color="black",face="bold"))+scale_color_gradient(low="blue",high="red")
labs(x=NULL,y=NULL)
s2

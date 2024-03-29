library(ggplot2)
setwd("~/gse205506/finalize")
a<-read.table("ANPEP_CA2", head=TRUE, sep="\t")
ggplot(a,aes(x=features.plot,y=id))+geom_point()
ggplot(a,aes(x=features.plot,y=id))+geom_point(aes(size=`correlation`,
                 color=`chisq.test`))  

a<-read.table("CDC42SE2_HLA-A_KLF6", head=TRUE, sep="\t")
ggplot(a,aes(x=features.plot,y=id))+geom_point()
ggplot(a,aes(x=features.plot,y=id))+geom_point(aes(size=`correlation`,
                                                   color=`chisq.test`))  


a<-read.table("NINJ2_UTS2_MT2A", head=TRUE, sep="\t")
ggplot(a,aes(x=features.plot,y=id))+geom_point()
ggplot(a,aes(x=features.plot,y=id))+geom_point(aes(size=`correlation`,
                                                   color=`chisq.test`))  
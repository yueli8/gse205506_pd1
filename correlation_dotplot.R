library(ggplot2)
setwd("~/gse205506/finalize")
a<-read.table("ANPEP_CA2", head=TRUE, sep="\t")
ggplot(a,aes(x=features.plot,y=id))+geom_point()
ggplot(a,aes(x=features.plot,y=id))+geom_point(aes(size=`correlation`,
                 color=`chisq.test`))  

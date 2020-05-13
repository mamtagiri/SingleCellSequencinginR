# create histogram for your UMI count, RNA count, Percent mitochondrial genes for seurat object

library(Seurat)
library(ggplot)

g<-YOUROBJECT@metadata 

#Percent mitocondrial histo

ggplot(g, aes(x=percent.mt)) + geom_histogram( binwidth=2, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    ggtitle("Percentage Mitochondrial genes histo") 
    
# UMI count histo    
ggplot(g, aes(x=nCount_RNA)) + geom_histogram( binwidth=2, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    ggtitle("UMI count histo") 
    
# Gene count histo    
ggplot(g, aes(x=nCount_features)) + geom_histogram( binwidth=2, fill="#69b3a2", color="#e9ecef", alpha=0.9) +
    ggtitle("Total genes histo") 

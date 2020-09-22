## SEurat functions###

library(Seurat)
library(cowplot)


####average expression table

AvgExTable<-AverageExpression(YOUROBJECT,slot="data",assays = "RNA")

## Violin plot ##
pdf(file="Myviolin.pdf")
plots <- VlnPlot(object = a1, features = c('ADM', 'CPA3','EPCAM','TPSAB1','CAMP',"AGR2"), split.by = 'group', pt.size = 0, combine = FALSE, log = TRUE,adjust=TRUE)
plots <- lapply(X = plots,FUN = function(x) x +  theme(axis.text.x = element_text(size = 16, angle = 90, hjust = .5, vjust = .5),axis.text.y = element_text(size = 10, angle = 0, hjust = 1, vjust = 0),axis.title.x = element_text(size = 1, angle = 0, hjust = .5, vjust = 0),axis.title.y = element_text(size = 16, angle = 90, hjust = .5, vjust = .5)))

CombinePlots(plots = plots, legend = 'right',nrow=1)
dev.off()


###Heatmap

DoHeatmap(
  object = immune.WT, 
  genes.use =include, 
  slim.col.label = TRUE,disp.min = -Inf,
  disp.max = Inf
)


##UMAP
UMAPPlot(a,pt.size=0.7,label.size=13,repel=TRUE,label=TRUE)


## Scale all genes in scale data
all.set<-rownames(myfile)
myfile <- ScaleData(myfile, features=all.set)
DoHeatmap(myfile,features = YOURLIST) + theme(axis.text.y = element_text(size = 20)) + theme(text = element_text(size = 20))

##Add new meta data to the clustered seurat object

mygroup<-c("Untreated","Untreated","DSS only","DSS only","Untreated","DSS only","DSS+BZA","Untreated","DSS only","DSS+BZA") ##keep in same order as output of Idents(OBJECT)
b$cluster<-b@active.ident
Idents(b)<-"orig.ident"
names(mygroup)<-levels(b)
b<-RenameIdents(b,mygroup)
b$treatment<-b@active.ident
Idents(b)<-"cluster"

##create a histogram for UMI, Genes
meta<-YOUROBJECT@meta.data
ggplot(meta, aes(x=percent.mt)) + geom_histogram( binwidth=100, fill="#69b3a2", color="#e9ecef", alpha=0.9) + ggtitle("Histogram of Percent Mito")

##Create cluster based UMI, gene, mito percent plot
meta<-YOUROBJECT@meta.data 
ggplot(meta,aes(x=seurat_clusters,y=percent.mt,fill=seurat_clusters)) + geom_boxplot() 


## subset to only cells expressing a certain gene
EpcamExpressing<-subset(x = OBJECT, subset = EPCAM > 0)

##create a table of cells
table(OBJECT@active.ident,OBJECT@meta.data$orig.ident)

##subset object by cluster
Cluster5<-subset(OBJECT,idents=5)

##Correlation plot for 2 genes in a cluster
cluster11<-subset(OBJECT,idents=11)
FeatureScatter(cluster11, feature1="Nox1",feature2="Reg4", group.by="group")

## expression of a particular gene for each cell
EpcamExpressing<-subset(x = OBJECT, subset = EPCAM > 0)

##Correlation plot for 2 genes in a cluster in a specific group
test <- subset(x = OBJECT, subset = group == "YOURGROUP")
FeatureScatter(test, feature1="Nox1",feature2="Reg4")

###get values for each cell for set of genes for any subset
cluster11<-subset(OBJECT,idents=11)
test <- subset(cluster11, subset = group == "CTRL")
data<-FetchData(test, vars = c("ptpro", "eps8"))

## differntial gene analysis between groups
Genotypemarkers<-FindMarkers(OBJECT, ident.1 = "WT", group.by = 'group')

#### differntial gene analysis between groups in a particular cluster , e.g cluster 2
Genotypemarkers<-FindMarkers(OBJECT, ident.1 = "WT", group.by = 'group',subset.ident = "2")

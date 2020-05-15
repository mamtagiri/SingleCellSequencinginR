
myfile<-readRDS("YOURFILE.rds") #% change the path to your rds integrated model file 

clusterid<-c(0,17) #just add the cluster # you want to recluster in the bracket
myreclusteredfile<-recluster(myfile,clusterid,20,0.6)  #change the values in the brackek starting with myfile
                                                  #the clusterid, PCs to use and then the resolution to use
saveRDS(myrecluster,file="Myreclsuteredfile.rds")


recluster <-function(object,clusterid,pc,resolution) {
  object1 <- subset(object, idents=clusterid,invert=TRUE)
  DefaultAssay(object1)<-"RNA"
  object1 <- NormalizeData(object1)
  all.set<-rownames(object1)
  object1 <- ScaleData(object1, display.progress = F,features=all.set)
  object1 <- FindVariableFeatures(object1,selection.method = "vst", nfeatures = 2000)
  object1 <- RunPCA(object1, npcs = pc, verbose = FALSE)
  object1 <- RunUMAP(object1, reduction = "pca", dims = 1:pc)
  object1 <- FindNeighbors(object1, reduction = "pca", dims = 1:pc)
  object1 <- FindClusters(object1, resolution = resolution)
  return(object1)
}

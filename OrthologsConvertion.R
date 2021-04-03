fish<-readRDS("")
DefaultAssay(mice)<-"RNA"
human<-readRDS("")

#read the orthology file
mart_export <- read.delim("orthology FILE", row.names=NULL)
mart1<-mart_export[,c(6,8)]
mart1$id<-paste0(mart1$Human.gene.name,"-",mart1$Gene.name)
mart2<-subset(mart1,!duplicated(mart1$id))


#filter the duplicates
mart2<-filter(mart2,!duplicated(Gene.name))

#subset data to the genes with orthologs
humangenes<-as.vector(unlist(mart2$Gene.name))
fishgenes<-as.vector(unlist(mart2$Human.gene.name))
fish1<-subset(fish, features=micegenes)
human1<-subset(human, features=humangenes)

#convert gene names for fish data
for (i in fishgenes){
  p=which(grepl(pn,mart2$Gene.name))
  p1<-p[1]
  f<-mart2$Human.gene.name[p1]
  pn2<-i
  pn3<-paste0("^",pn2,"$")
  k=grep(pn3,fish1@assays$RNA@counts@Dimnames[[1]])
  k1<-k[1]
  fish1@assays$RNA@counts@Dimnames[[1]][k1]<-as.character(f)}

for (i in fishgenes){
  pn <-paste0("^",i,"$")
  p=which(grepl(pn,mart2$Gene.name))
  p1<-p[1]
  f<-mart2$Human.gene.name[p1]
  pn2<-i
  pn3<-paste0("^",pn2,"$")
  k=grep(pn3,fish1@assays$RNA@data@Dimnames[[1]])
  k1<-k[1]
  fish1@assays$RNA@data@Dimnames[[1]][k1]<-as.character(f)}

#save the rds with gene names changed
saveRDS(fish1,file="FishOrthologs.rds")
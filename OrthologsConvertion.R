##########################################
# Ortholog gene name conversion in seurat
##########################################
library(Seurat)
library(dplyr)

# Add the path to the appropriate rds file
# read the fish rds file after seurat clustering and set asssay to RNA
fish<-readRDS("")
DefaultAssay(mice)<-"RNA"

# read the fish rds file after seurat clustering and set asssay to RNA
human<-readRDS("")
DefaultAssay(human)<-"RNA"

#read the orthology file, add path to the biomart file downloaded from biomart or Zfin

mart_export <- read.delim("ORTHOLOGY FILE", row.names=NULL)

# the mart file has column 6= human gene name, col 8= fish gene names
mart1<-mart_export[,c(6,8)]

# create a column with Humna genename- Fish gene name to remove any duplicates (keep unique pairs)
mart1$id<-paste0(mart1$Human.gene.name,"-",mart1$Gene.name)
mart2<-subset(mart1,!duplicated(mart1$id))


#filter the duplicates for same gene names 
mart2<-filter(mart2,!duplicated(Gene.name))

#subset data to the genes with orthologs
# create a list of all the gene names in the above step (ortholog pairs)
fishgenes<-as.vector(unlist(mart2$Gene.name))
humangenes<-as.vector(unlist(mart2$Human.gene.name))

# subsetting the rds file read in line 9, 13 to only genes that have orthologous pairs, you can skip this step to keep all genes

fish1<-subset(fish, features=micegenes)
human1<-subset(human, features=humangenes)

#convert gene names for fish data
# iterate over each fish gene name, find its position and corresponding ortholog name from the mart file. Search the fish name in the rds file and change to the ortholog just read.  
for (i in fishgenes){
  pn <-paste0("^",i,"$")  # using regex for exact match
  p=which(grepl(pn,mart2$Gene.name)) # find position of the match
  p1<-p[1]  # keeping the first position from p
  f<-mart2$Human.gene.name[p1] # get the human gene name from the same position
  pn2<-i
  pn3<-paste0("^",pn2,"$") # regex exact match for the fish name
  k=grep(pn3,fish1@assays$RNA@counts@Dimnames[[1]]) # grep the position in the rds fish object for the fish gene name
  k1<-k[1] # keep only the first position from above
  fish1@assays$RNA@counts@Dimnames[[1]][k1]<-as.character(f)} # change the exact position with its correspondong human gene name read above

#convert gene names for fish data- NOW in the DATA slot of rds file
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
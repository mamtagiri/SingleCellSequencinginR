
library(gplots)
library(RColorBrewer)
# Get some nicer colours
mypalette <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypalette)
# Set up colour vector for celltype variable
#col.cell <- c("purple","orange")[test$samples]

#load your file with rows as genes and columns as samples

testdata <- read.csv("~/PATH/YOURFILE.csv", row.names=1)

# Plot the heatmap
testnum<-as.matrix(testdata)

heatmap.2(t(testnum1), 
          col=rev(morecols(50)),
          dendrogram = "none",scale="row",margins = c(16,12), keysize = 0.8,cexRow=1.9,sepwidth=c(0.02,0.02),
          sepcolor="Gray",cexCol=1.9,Rowv=FALSE,Colv=FALSE,colsep=1:ncol(testnum1),
          rowsep=1:nrow(testnum1))
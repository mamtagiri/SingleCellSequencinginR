
#######################
tracks=NULL

##Load file with all the columns needed, first columns as chromosome names (or genes names)
n1 <- read.csv("~/PATH/n1.csv")
n1<-n1  %>% arrange(desc("YOUR DIFFERENT COLUMNS")) %>% group_by(avg_logFC)


#make circumferences
genomeChr = as.list(n1$chromosome)
lengthChr =n1$length
names(lengthChr) <- genomeChr

#YOUR first ring
genomeChr1 = as.list(n1$gene.x)
lengthChr1 = n1$length2T
names(lengthChr1) <- genomeChr1


#Second RING

genomeChr2 = as.list(n1$gene.y)
lengthChr2 = n1$length2L
names(lengthChr2) <- genomeChr2

#Third ring
genomeChr3 = as.list(n1$gene.y.y)
lengthChr3 = n1$length2C
names(lengthChr3) <- genomeChr3

#Fourth ring
genomeChr4 = as.list(n1$gene.x.x)
lengthChr4 = n1$length2I
names(lengthChr4) <- genomeChr4



library(dplyr)

# Define boxes positions

boxChromosomes = rep(genomeChr,lengthChr)
boxChromosomes1 = rep(genomeChr1,lengthChr1)
boxChromosomes2 = rep(genomeChr2,lengthChr2)
boxChromosomes3 = rep(genomeChr3,lengthChr3)
boxChromosomes4 = rep(genomeChr4,lengthChr4)

# Define values for two heatmap tracks
boxVal1 = rep(n1$avg_logFC,lengthChr)
#First ring
val2 = rep(n1$lengtht,lengthChr1)
#Second ring
val3 = rep(n1$lengthL,lengthChr2)
#third ring
val4 = rep(n1$lengthC,lengthChr3)
#fourth ring
val5 = rep(n1$lengthI,lengthChr4)


#define tracks

tracks = BioCircosHeatmapTrack("heatmap1", boxChromosomes,1,500,
                               boxVal1, minRadius = 1.0, maxRadius = 1.15,color = c("green", "red"))
tracks = tracks+BioCircosHeatmapTrack("heatmap2", boxChromosomes1,1,500,values = val2,range=c(0,1),
                                      minRadius = 0.9, maxRadius = 0.8,color = c("white","#14C6CC")) #14c6cc
tracks = tracks+BioCircosHeatmapTrack("heatmap3", boxChromosomes2,1,500,values = val3,range=c(0,1),
                                      minRadius = 0.75, maxRadius = 0.65,color = c("white","#D8BFD8"))

tracks = tracks+BioCircosHeatmapTrack("heatmap4", boxChromosomes3,1,500,values = val5,range=c(0,1),
                                      minRadius = 0.45, maxRadius = 0.35,color = c("white","#F58F73"))#f58f73

tracks = tracks+BioCircosHeatmapTrack("heatmap5",boxChromosomes4,1,500,values = val4,range=c(0,1),
                                      minRadius = 0.6, maxRadius = 0.5,color = c("white","#8FBC8F")) ##f58f73

#Put them all together
BioCircos(tracks, genome = as.list(lengthChr), chrPad = 0.02,displayGenomeBorder = TRUE,genomeBorderSize = 0.5, genomeTicksLen = 0, genomeTicksTextSize = 0, genomeTicksScale = 1e+8,genomeLabelTextSize = "15pt", genomeLabelDy = 50,genomeLabelOrientation = 90)

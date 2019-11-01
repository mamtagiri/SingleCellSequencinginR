######## CHANGE the LOCATION#######

load("/Path to.RData")

genestouse <-c("Genes of interest")

DoHeatmap(
  object = seuratobject, 
  genes.use =genestouse, use.scaled = TRUE,
  slim.col.label = TRUE,disp.min = -5, disp.max = 8,col.low = "#FF00FF",
  col.mid = "#000000", col.high = "#FFFF00"
)


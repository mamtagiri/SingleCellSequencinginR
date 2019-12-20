#monocle


library(monocle)


######load the seurat object first ######

seuratobject<-readRDS(YOURPATH/SEURATRDSfile)
data <- as(as.matrix(seuratobject@assays$RNA@data), 'sparseMatrix')

pd <- new('AnnotatedDataFrame', data = seuratobject@meta.data)

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)

FibroMono <- newCellDataSet(data,
                          phenoData  = pd,
                          featureData = fd,expressionFamily = negbinomial.size())


####ESTIMATE DISPERSIONS AND SIZE FACTORS############
FibroMono <- estimateSizeFactors(FibroMono)
FibroMono <- estimateDispersions(FibroMono)


#######FILTER LOW QUALITY GENES/CELLS##########
FibroMono <- detectGenes(FibroMono, min_expr = 0.1)
print(head(fData(FibroMono)))

expressed_genes <- row.names(subset(fData(FibroMono),
                                    num_cells_expressed >= 10))
######DEG TEST##########


diff_test_res <- differentialGeneTest(FibroMono[expressed_genes,],
                                      fullModelFormulaStr = "~timegroup")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

FibroMono<- setOrderingFilter(FibroMono, ordering_genes)
FibroMono <- reduceDimension(FibroMono, max_components = 2,method = 'DDRTree')
FibroMono <- orderCells(FibroMono)
png(filename="Trajectory.png")
plot_ordering_genes(FibroMono)
dev.off()
png(filename="TrajectorybyGroup.png")
plot_cell_trajectory(FibroMono, color_by = "group")
dev.off()


###SIGNIFICANT GENES CHANGING WITH PSEUDOTIME##########
GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$Hours)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
diff_test_res <- differentialGeneTest(FibroMono[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res, qval < 0.1))
FibroMono <- orderCells(FibroMono, root_state = GM_state(FibroMono))

###PLOT####
png(filename="TrajectorybyPSEUDOTIME.png")
plot_cell_trajectory(FibroMono, color_by = "Pseudotime")
dev.off()
png(filename="PSEUDOTIMEHEATMAP.png")
plot_pseudotime_heatmap(FibroMono[sig_gene_names,],
                        num_clusters = 3,
                        cores = 1,
                        show_rownames = T)
dev.off()


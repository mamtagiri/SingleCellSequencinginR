if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')


#make valcano plot from sheet with p value and fold change
EnhancedVolcano(Yourfile,
                lab = Yourfile$GeneSymbol,
                x = 'log2FoldChange',
                y = 'pval',
                xlim = c(-5, 8),
                FCcutoff = 1.0,
                pCutoff = 0.05,
                ylim=c(0,5),
                subtitle = 'Your sub title',
                title = 'your title')

library(calibrate)
#with(subset(Yourfile, pval<.05 ), points(log2FoldChange, -log10(pval), pch=20, col="red"))


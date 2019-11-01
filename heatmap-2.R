library(ggplot2)
library(reshape2)

###change your file location!!!#####
test_heatmap <- read.csv("~/Downloads/test_heatmap.csv")
############################


test = setNames(data.frame(t(test_heatmap[,-1])), test_heatmap[,1])
test$samples<-rownames(test)
data <- ddply(melt(test), .(variable), transform)
ggplot(data, aes(samples, variable )) +
  geom_tile(aes(fill = value), , color = "white") +
  scale_fill_gradient(low = "light blue", high = "orange") +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 12),
        plot.title = element_text(size=16),
        axis.title=element_text(size=14,face="bold"),
        axis.text.x = element_text(angle = 90, hjust = 1)) + ylab("genes")+xlab("samples")+labs(fill = "Scale")

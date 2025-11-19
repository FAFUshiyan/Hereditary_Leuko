library(pheatmap)
he <- read.table("/Users/shiyan/Downloads/median.tsv",header = T,row.names = 1,sep = "\t")
my_colors <- c("#67AAC7", "white", "#BF6863")
my_palette <- colorRampPalette(my_colors)(100)
pheatmap_breaks <- seq(-1, 1, length.out = 101)
#h <- pheatmap(he,scale = "row",cluster_rows = T,cluster_cols = F,color = my_palette,legend_breaks = c(-1,0,1), breaks = pheatmap_breaks)
h <- pheatmap(he,scale = "row",cluster_rows = F,cluster_cols = F,color = my_palette,legend_breaks = c(-1,0,1))
h
ggsave("~/Downloads/heatmap.pdf",plot = h,width = 9,height = 4)

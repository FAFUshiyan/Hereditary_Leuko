library(pheatmap)

# 1. 读入矩阵
he <- read.table(
  "/Users/shiyan/Downloads/human.diff.sort.lipid.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

# 2. 颜色与 break
my_colors <- c("#67AAC7", "white", "#BF6863")
my_palette <- colorRampPalette(my_colors)(100)
pheatmap_breaks <- seq(-2, 2, length.out = 101)

# 3. 按行 scale（模拟 pheatmap 的 scale="row"）
#   注意：如果某些行是常数行，会出现 NA，可以在后面用 na.omit 或者用一个很小的 sd 替代
he_scaled <- t(scale(t(he)))

# 4. 定义原始 gaps_row
gaps_row_vec <- c(12, 31, 65, 90, 96, 104)

nr <- nrow(he_scaled)

# 5. 根据 gaps_row 计算每一块的行索引
#    例如：1:12, 13:31, 32:65, ...
block_bounds <- c(0, gaps_row_vec, nr)
block_indices <- lapply(seq_len(length(block_bounds) - 1), function(i) {
  seq(block_bounds[i] + 1, block_bounds[i + 1])
})

# 6. 在每一块内部做聚类并得到重排后的索引
reordered_idx_list <- lapply(block_indices, function(idx) {
  submat <- he_scaled[idx, , drop = FALSE]
  
  # 若某个块太小或只有一行，直接返回
  if (length(idx) == 1) {
    return(idx)
  }
  
  # 如果出现 NA（比如全 0 行被 scale 出 NA），简单处理一下
  # 把含 NA 的行去掉或用原顺序
  if (any(is.na(submat))) {
    # 这里采用保守做法：不聚类，直接原顺序
    return(idx)
  }
  
  d  <- dist(submat, method = "euclidean")
  hc <- hclust(d, method = "complete")  # 也可以换成 "ward.D2"
  
  # hc$order 是块内的相对顺序，这里转换回原始行号
  idx[hc$order]
})

# 7. 拼成新的整体行顺序
row_order_new <- unlist(reordered_idx_list)

# 8. 用新的顺序重排原始矩阵（不要用 he_scaled，画图时交给 pheatmap 再 scale）
he_reordered <- he[row_order_new, , drop = FALSE]

# 9. 重新计算 gaps_row（按每一块大小的累积和来算）
block_sizes <- sapply(block_indices, length)
gaps_row_new <- cumsum(block_sizes)[-length(block_sizes)]

# 10. 画最终的 heatmap（pheatmap 内部再做 scale="row"）
h <- pheatmap(
  he_reordered,
  scale         = "row",
  cluster_rows  = FALSE,       # 不再让它重新聚类
  cluster_cols  = FALSE,       # 如果想对列聚类可以改 TRUE
  color         = my_palette,
  legend_breaks = c(-2, 0, 2),
  breaks        = pheatmap_breaks,
  gaps_row      = gaps_row_new,
  gaps_col      = c(26),
  border_color  = NA,          # 建议去掉格子边框，更现代
  fontsize_row  = 6,           # 行多时减小字号或直接隐藏
  fontsize_col  = 8
)

# 11. 导出 PDF（推荐用 pdf()，也可以继续 ggsave）
pdf("~/Downloads/lipid.all.heatmap.block_clustered.pdf", width = 9, height = 6)
print(h)
dev.off()
